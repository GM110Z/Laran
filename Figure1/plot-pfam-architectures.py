#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


# Edit these values, choose --palette both/warm/blue/reference, or supply colours.
PALETTES = {
    "both": ("#F25532", "#4D5793", "#758BAE", "#111A45",
             "#554747", "#A7A593", "#FFF4E5", "#EAE0CF"),
    "warm": ("#F25532", "#554747", "#A7A593", "#FFF4E5"),
    "blue": ("#4D5793", "#111A45", "#758BAE", "#EAE0CF"),
    "reference": ("#F5EBDD", "#F2765E", "#315B8C", "#413333"),
}
NO_DOMAINS = "No Pfam domain detected"
PROTEIN_FIELDS = ("protein", "protein_accession", "accession")
HATCHES = ("", "..", "xx", "++", "oo", "\\\\", "--", "**", "||")


def table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = [(key or "").strip().lstrip("#") for key in reader.fieldnames or []]
        rows = [{(key or "").strip().lstrip("#"): (value or "").strip()
                 for key, value in row.items() if key is not None} for row in reader]
    return fields, rows


def pick(row: dict[str, str], names: tuple[str, ...]) -> str:
    return next((row[key] for key in names if row.get(key)), "")


def fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name, pieces = "", []

    def store() -> None:
        if not name:
            return
        sequence = "".join(pieces).upper().rstrip("*")
        if not sequence:
            raise ValueError(f"Empty FASTA sequence: {name}")
        if name in sequences and sequences[name] != sequence:
            raise ValueError(f"Conflicting sequences for FASTA identifier {name}")
        sequences[name] = sequence

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                name, pieces = line[1:].split()[0], []
            elif name:
                pieces.append("".join(line.split()))
            else:
                raise ValueError(f"Sequence before first FASTA header in {path}")
    store()
    if not sequences:
        raise ValueError(f"No FASTA records in {path}")
    return sequences


def require_stage(run: Path, stage: str) -> None:
    path = run / stage / "status.tsv"
    if not path.is_file():
        raise ValueError(f"Missing {path}; cannot verify this stage completed. "
                         "For exported/older results, supply explicit input files without --run.")
    status = dict(line.rstrip("\n").split("\t", 1) for line in path.read_text().splitlines()
                  if "\t" in line and not line.startswith("#"))
    if status.get("status") != "ok":
        raise ValueError(f"{stage} stage is {status.get('status', 'unknown')!r}, not 'ok': {path}")


def coordinates(start: str, end: str, protein: str, sequences: dict[str, str]) -> tuple[int, int]:
    left, right = int(start), int(end)
    if not 1 <= left <= right <= len(sequences[protein]):
        raise ValueError(f"Invalid 1-based coordinates {left}-{right} for {protein} "
                         f"(length {len(sequences[protein])})")
    return left, right


def read_domains(path: Path, sequences: dict[str, str]) -> dict[str, list[dict]]:
    fields, rows = table(path)
    if not set(fields).intersection(PROTEIN_FIELDS) or not {"domain", "start", "end"} <= set(fields):
        raise ValueError(f"Unrecognised domain TSV columns in {path}: {fields}")
    if not {"pfam", "pfam_accession", "database"}.intersection(fields):
        raise ValueError("Domain TSV needs a pfam/pfam_accession or database column to distinguish Pfam hits.")
    result: dict[str, list[dict]] = defaultdict(list)
    seen = set()
    for row in rows:
        protein = pick(row, PROTEIN_FIELDS)
        if protein not in sequences:
            continue  
        accession = pick(row, ("pfam", "pfam_accession"))
        match = re.fullmatch(r"(PF\d+)(?:\.\d+)?", accession, flags=re.I)
        if not match and "pfam" not in row.get("database", "").lower():
            continue  
        name = row["domain"]
        if not name or name == "-":
            raise ValueError(f"Missing domain name for {protein}")
        identity = match.group(1).upper() if match else name
        start, end = coordinates(row["start"], row["end"], protein, sequences)
        key = (protein, identity, start, end)
        if key not in seen:
            seen.add(key)
            result[protein].append(dict(identity=identity, domain=name, start=start, end=end))
    for hits in result.values():
        hits.sort(key=lambda hit: (hit["start"], hit["end"], hit["identity"]))
    return result


def read_tms(path: Path, sequences: dict[str, str], summary: bool) -> dict[str, list[tuple[int, int]]]:
    fields, rows = table(path)
    needed = {"tm_helix_coordinates_aa"} if summary else {"kind", "start", "end"}
    if not set(fields).intersection(PROTEIN_FIELDS) or not needed <= set(fields):
        raise ValueError(f"Unrecognised TM TSV columns in {path}: {fields}")
    result: dict[str, set[tuple[int, int]]] = defaultdict(set)
    for row in rows:
        protein = pick(row, PROTEIN_FIELDS)
        if protein not in sequences:
            continue
        if summary:
            text = row["tm_helix_coordinates_aa"]
            if text in ("", "-", "NA"):
                continue
            pairs = []
            for item in text.split(";"):
                match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", item)
                if not match:
                    raise ValueError(f"Invalid TM coordinates for {protein}: {item!r}")
                pairs.append(match.groups())
        elif row["kind"].lower() in ("tm", "tmhelix", "transmembrane"):
            pairs = [(row["start"], row["end"])]
        else:
            continue
        for start, end in pairs:
            result[protein].add(coordinates(start, end, protein, sequences))
    return {protein: sorted(pairs) for protein, pairs in result.items()}


def group_proteins(sequences: dict[str, str], domains: dict[str, list[dict]],
                   preferred: list[str]) -> list[dict]:
    groups = defaultdict(list)
    for protein in sorted(sequences):
        key = tuple(hit["identity"] for hit in domains.get(protein, []))
        groups[key].append(protein)
    for protein in preferred:
        if protein not in sequences:
            raise ValueError(f"Requested representative is not in the query FASTA: {protein}")
    output = []
    for key, members in groups.items():
        median = statistics.median(len(sequences[protein]) for protein in members)
        choices = [protein for protein in preferred if protein in members]
        if len(set(choices)) > 1:
            raise ValueError(f"More than one forced representative for architecture {key}")
        representative = choices[0] if choices else min(
            members, key=lambda protein: (abs(len(sequences[protein]) - median), protein))
        label = " + ".join(hit["domain"] for hit in domains.get(representative, [])) or NO_DOMAINS
        output.append(dict(key=key, members=members, representative=representative,
                           median_length=median, architecture=label))
    return sorted(output, key=lambda group: (-len(group["members"]), group["key"]))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def domain_styles(domains, palette, highlighted, custom_colours, previous=None):
    """One persistent style per Pfam; exhaust colours before adding hatching."""
    names = {}
    for hits in domains.values():
        for hit in hits:
            names.setdefault(hit["identity"], hit["domain"])
    order = sorted(names)
    if highlighted in order:
        order.remove(highlighted)
        order.insert(0, highlighted)
    styles = {}
    used = {(style["colour"].lower(), style["hatch"]) for style in (previous or {}).values()}
    for identity in order:
        if previous and identity in previous:
            style = dict(previous[identity])
            style["domain"] = names[identity]
            styles[identity] = style
            used.add((style["colour"].lower(), style["hatch"]))
    available = ((colour, hatch) for hatch in HATCHES for colour in palette)
    for identity in order:
        if identity in styles:
            continue
        for colour, hatch in available:
            if (colour.lower(), hatch) not in used:
                break
        else:
            raise ValueError("Too many domain styles for this palette; provide more colours.")
        colour = custom_colours.get(identity, custom_colours.get(names[identity], colour))
        # Explicit colour overrides remain distinguishable by pattern if reused.
        while (colour.lower(), hatch) in used:
            hatch += "."
        used.add((colour.lower(), hatch))
        styles[identity] = dict(domain=names[identity], colour=colour, hatch=hatch)
    for identity, style in styles.items():
        style["colour"] = custom_colours.get(identity, custom_colours.get(names[identity], style["colour"]))
    return styles


def plot(args, sequences, domains, tms, groups, styles) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.ticker import MaxNLocator

    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                                "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    body, outline = args.backbone_colour, "#302B2B"
    has_tm = any(tms.get(group["representative"]) for group in groups)
    legend_rows = (len(styles) + int(has_tm) + 3) // 4
    main_height = max(1.4, 0.5 * len(groups))
    legend_height = max(.30, legend_rows * .30)
    fig = plt.figure(figsize=(args.width, main_height + legend_height + 1.1), layout="constrained")
    grid = fig.add_gridspec(2, 1, height_ratios=(main_height, legend_height), hspace=.06)
    ax = fig.add_subplot(grid[0])
    legend_ax = fig.add_subplot(grid[1])
    legend_ax.set_axis_off()
    longest = max(len(sequences[group["representative"]]) for group in groups)
    used = {}
    for y, group in enumerate(groups):
        protein = group["representative"]
        length = len(sequences[protein])
        ax.add_patch(Rectangle((0, y - .16), length, .32, facecolor=body,
                               edgecolor=outline, linewidth=.8))
        for hit in domains.get(protein, []):
            style = styles[hit["identity"]]
            ax.add_patch(Rectangle((hit["start"] - 1, y - .16), hit["end"] - hit["start"] + 1,
                                   .32, facecolor=style["colour"], edgecolor=outline,
                                   hatch=style["hatch"], linewidth=.45))
            used[hit["identity"]] = style
        for start, end in tms.get(protein, []):
            ax.add_patch(Rectangle((start - 1, y - .16), end - start + 1, .32,
                                   facecolor=args.tm_colour, edgecolor=outline,
                                   hatch="////", linewidth=.65, zorder=3))
        suffix = f"; {protein}" if args.show_accessions else ""
        ax.text(length + longest * .025, y, f"n={len(group['members'])}{suffix}", va="center", fontsize=10)
    ax.set_yticks(range(len(groups)), [group["architecture"] for group in groups], fontsize=10)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xlim(0, longest * (1.48 if args.show_accessions else 1.20))
    ax.set_ylim(len(groups) - .45, -.65)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=7, integer=True))
    ax.set_xlabel("Amino-acid position", fontsize=11)
    ax.spines[["top", "right", "left"]].set_visible(False)
    if args.title:
        ax.set_title(args.title, loc="left", pad=16, fontsize=12, weight="bold")
    handles = [Patch(facecolor=style["colour"], edgecolor=outline,
                     hatch=style["hatch"], label=style["domain"])
               for identity, style in sorted(used.items())]
    if has_tm:
        handles.append(Patch(facecolor=args.tm_colour, edgecolor=outline, hatch="////", label="TM helix"))
    if handles:
        legend_ax.legend(handles=handles, loc="center left", ncol=min(4, len(handles)),
                         frameon=False, fontsize=9)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(args.out / f"architectures.{extension}", dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, help="Existing FlaGs3 run; automatically selects query proteins only")
    parser.add_argument("--fasta", type=Path, help="Explicit protein FASTA; without --run, every unique FASTA ID is a query")
    parser.add_argument("--domains", type=Path, help="FlaGs Pfam TSV; accepts current protein or older #protein headers")
    tm = parser.add_mutually_exclusive_group()
    tm.add_argument("--features", type=Path, help="FlaGs features.tsv for TM overlays")
    tm.add_argument("--tm-summary", type=Path, help="Older summarize_deeptmhmm.py output TSV")
    tm.add_argument("--no-tm", action="store_true", help="Explicitly omit TM overlays")
    parser.add_argument("--out", type=Path, default=Path("architecture_panel"))
    parser.add_argument("--palette", default="both", help="both (all eight supplied colours), warm, blue, reference, or any comma-separated list of domain colours")
    parser.add_argument("--highlight-domain", help="Give this domain the first palette colour; all other domains still get individual styles")
    parser.add_argument("--domain-colour", action="append", default=[], metavar="DOMAIN=HEX", help="Override a domain colour; repeat as needed")
    parser.add_argument("--colour-map", type=Path, help="Reuse domain_colours.tsv from a previous plot to keep Pfam styles fixed across datasets")
    parser.add_argument("--backbone-colour", default="#F7F7F7", help="Protein backbone fill, separate from the domain palette")
    parser.add_argument("--tm-colour", default="none", help="TM overlay fill; default transparent with diagonal hatching")
    parser.add_argument("--representative", action="append", default=[], metavar="ACCESSION", help="Force one representative for its architecture; repeat for different architectures")
    parser.add_argument("--title", default="Protein domain architectures")
    parser.add_argument("--show-accessions", action="store_true")
    parser.add_argument("--width", type=float, default=10.0)
    args = parser.parse_args(argv)
    try:
        from matplotlib.colors import is_color_like
        inputs = {}
        if args.run:
            args.run = args.run.expanduser().resolve()
            if args.out.expanduser().resolve().is_relative_to(args.run):
                raise ValueError("Choose an --out directory outside the FlaGs run to preserve its contents.")
            require_stage(args.run, "extract")
            if not args.domains:
                require_stage(args.run, "domains")
            args.fasta = args.fasta or args.run / "extract" / "proteins.faa"
            args.domains = args.domains or args.run / "domains" / "domains.tsv"
            if not (args.features or args.tm_summary or args.no_tm):
                if (args.run / "features" / "features.tsv").is_file():
                    require_stage(args.run, "features")
                    args.features = args.run / "features" / "features.tsv"
                else:
                    raise ValueError("No FlaGs TM feature table; provide --features/--tm-summary or explicitly use --no-tm.")
        if not args.fasta or not args.domains:
            raise ValueError("Provide --run, or both --fasta and --domains.")
        if args.width <= 0:
            raise ValueError("--width must be positive")
        sequences = fasta(args.fasta)
        if args.run:
            fields, rows = table(args.run / "extract" / "rows.tsv")
            if "accession" not in fields:
                raise ValueError("extract/rows.tsv has no accession column")
            queries = {row["accession"] for row in rows if row.get("accession") not in (None, "", "-")}
            missing = queries - sequences.keys()
            if missing:
                raise ValueError(f"Missing query sequences: {', '.join(sorted(missing)[:10])}")
            sequences = {name: sequences[name] for name in sorted(queries)}
            inputs["query_rows"] = args.run / "extract" / "rows.tsv"
        if not sequences:
            raise ValueError("No query proteins selected")
        domains = read_domains(args.domains, sequences)
        tms = read_tms(args.features or args.tm_summary, sequences, bool(args.tm_summary)) if (args.features or args.tm_summary) else {}
        prevalence = Counter(identity for protein in sequences
                             for identity in {hit["identity"] for hit in domains.get(protein, [])})
        names = {hit["domain"]: hit["identity"] for hits in domains.values() for hit in hits}
        highlighted = names.get(args.highlight_domain, args.highlight_domain)
        if args.highlight_domain and highlighted not in prevalence:
            raise ValueError(f"Highlighted domain not found in query Pfam calls: {args.highlight_domain}")
        if not prevalence:
            print("Warning: no query Pfam hits matched. Check identifiers before interpreting the no-domain group.", file=sys.stderr)
        palette = list(PALETTES[args.palette]) if args.palette in PALETTES else [value.strip() for value in args.palette.split(",")]
        if not palette or not all(is_color_like(value) for value in palette):
            raise ValueError("--palette must be both/warm/blue/reference or valid comma-separated colours (include # for hex colours).")
        if not is_color_like(args.backbone_colour) or not is_color_like(args.tm_colour):
            raise ValueError("Invalid --backbone-colour or --tm-colour")
        custom = {}
        for setting in args.domain_colour:
            name, sep, colour = setting.partition("=")
            if not sep or not name or not is_color_like(colour):
                raise ValueError(f"Invalid --domain-colour: {setting}")
            custom[name] = colour
        previous = {}
        if args.colour_map:
            fields, rows = table(args.colour_map)
            if not {"pfam", "domain", "colour", "hatch"} <= set(fields):
                raise ValueError("--colour-map needs pfam, domain, colour, hatch columns")
            for row in rows:
                if not row["pfam"] or not is_color_like(row["colour"]):
                    raise ValueError("Invalid identity or colour in --colour-map")
                if any(char not in '/\\|-.+xoO*' for char in row["hatch"]):
                    raise ValueError("Unsupported hatch pattern in --colour-map")
                previous[row["pfam"]] = dict(domain=row["domain"], colour=row["colour"], hatch=row["hatch"])
        groups = group_proteins(sequences, domains, args.representative)
        styles = domain_styles(domains, palette, highlighted, custom, previous)
        args.out = args.out.expanduser()
        args.out.mkdir(parents=True, exist_ok=True)
        plot(args, sequences, domains, tms, groups, styles)
        saved_styles = {**previous, **styles}
        write_tsv(args.out / "domain_colours.tsv", ["pfam", "domain", "colour", "hatch"],
                  [dict(pfam=identity, **style) for identity, style in sorted(saved_styles.items())])
        counts, membership = [], []
        for group in groups:
            rep = group["representative"]
            counts.append(dict(architecture=group["architecture"], pfam_order=" + ".join(group["key"]),
                               n=len(group["members"]), representative=rep, representative_length=len(sequences[rep]),
                               median_length=group["median_length"], percent=100 * len(group["members"]) / len(sequences)))
            for protein in group["members"]:
                membership.append(dict(protein=protein, length=len(sequences[protein]), architecture=group["architecture"],
                    representative="yes" if protein == rep else "no",
                    domain_coordinates="; ".join(f"{h['domain']}:{h['start']}-{h['end']}" for h in domains.get(protein, [])),
                    tm_coordinates="; ".join(f"{start}-{end}" for start, end in tms.get(protein, []))))
        write_tsv(args.out / "architecture_counts.tsv", list(counts[0]), counts)
        write_tsv(args.out / "architecture_membership.tsv", list(membership[0]), membership)
        inputs.update(fasta=args.fasta, domains=args.domains)
        if args.features or args.tm_summary:
            inputs["tm"] = args.features or args.tm_summary
        hashes = {}
        for name, path in inputs.items():
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            hashes[name] = {"path": str(path.resolve()), "sha256": digest.hexdigest()}
        metadata = dict(query_proteins=len(sequences), architectures=len(groups),
                        counting_unit="unique query protein accession, not genome occurrences or unique sequences",
                        architecture_rule="ordered Pfam accessions; repeats retained; exact duplicate hits removed",
                        representative_rule="closest to architecture median length; accession breaks ties; explicit overrides allowed",
                        first_colour_domain=highlighted, palette=palette, domain_colour_overrides=custom,
                        domain_styles=styles, backbone_colour=args.backbone_colour, tm_colour=args.tm_colour,
                        colour_rule="one style per Pfam; use every palette colour before reusing colours with different hatching",
                        tm_overlay=bool(args.features or args.tm_summary), coordinate_system="1-based inclusive input",
                        inputs=hashes)
        (args.out / "plot_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"Done: {len(sequences)} unique query proteins, {len(groups)} architectures. Outputs: {args.out.resolve()}")
    except (ValueError, OSError, ImportError) as error:
        parser.exit(2, f"Error: {error}\n")


if __name__ == "__main__":
    main()
