#!/usr/bin/env python3
"""FlankGuard: cached neighbourhood analysis with harmonised system reporting.

Existing contexts and successful scanner runs are reused automatically.
Reports are always rebuilt from raw DefenseFinder/PADLOC results.

--replot-only prohibits downloads and scanner execution.
--retry-failed retries previously failed IPG/context retrievals.

Only Python's standard library is required for this script.
DefenseFinder and PADLOC are external scanner dependencies.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import html
from http.client import IncompleteRead
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UNKNOWN = "Unknown"
ACCESSION_RE = re.compile(
    r"(?:[A-Z]{1,4}_)?[A-Z]{0,3}\d{5,}(?:\.\d+)?", re.I
)
EXCLUDED_SYSTEMS = {
    "vspr": "Very short patch repair; not an antiviral defence system",
    "dms-other": "Unclassified DNA-modification call; excluded from defence summaries",
}
CONTEXT_FILES = (
    "context.json",
    "genes.tsv",
    "neighborhood.fna",
    "neighborhood.faa",
    "neighborhood.gff",
)


def clean(value):
    text = "" if value is None else str(value).strip()
    return text or UNKNOWN


def key(value):
    return re.sub(r"[^a-z0-9]+", "-", str(value).lower()).strip("-")


def local_name(tag):
    return tag.rsplit("}", 1)[-1].lower()


def attr(element, *names):
    values = {key(k): clean(v) for k, v in element.attrib.items()}
    for name in names:
        value = values.get(key(name), UNKNOWN)
        if value != UNKNOWN:
            return value
    return UNKNOWN


def element_text(element, suffix):
    for child in element.iter():
        if local_name(child.tag).endswith(suffix.lower()):
            value = clean(child.text)
            if value != UNKNOWN:
                return value
    return UNKNOWN


def integer(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


def joined(values):
    return "; ".join(dict.fromkeys(
        str(v).strip() for v in values
        if v is not None and str(v).strip()
    ))


def atomic_write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def read_tsv(path, optional=False):
    path = Path(path)
    if optional and not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(
            column for row in rows for column in row
        ))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t",
            extrasaction="ignore", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_json(path, default=None):
    path = Path(path)
    if not path.is_file():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_proteins(path, protein_column=None):
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    if not text.strip():
        raise ValueError(f"Input is empty: {path}")
    try:
        delimiter = csv.Sniffer().sniff(
            text[:65536], delimiters=",\t;|"
        ).delimiter
    except csv.Error:
        delimiter = "\t" if "\t" in text.partition("\n")[0] else ","

    rows = [
        row for row in csv.reader(text.splitlines(), delimiter=delimiter)
        if any(cell.strip() for cell in row)
    ]
    if not rows or len(rows[0]) < 2:
        raise ValueError("Supply a table with at least two columns")

    first = [cell.strip() for cell in rows[0]]
    headers = {key(name): index for index, name in enumerate(first)}
    recognised = [
        "protein-accession", "proteinaccession", "sseqid",
        "subject-accession", "subjectaccession", "subjectid",
        "subject-id", "accession",
    ]
    if protein_column:
        wanted = key(protein_column)
        if wanted not in headers:
            raise ValueError(f"Protein column not found: {protein_column}")
        column = headers[wanted]
        has_header = True
    else:
        matched = next((name for name in recognised if name in headers), None)
        column = headers[matched] if matched else 1
        has_header = matched is not None

    data = rows[1:] if has_header else rows
    proteins = []
    for row in data:
        if len(row) <= column or not row[column].strip():
            continue
        matches = ACCESSION_RE.findall(row[column])
        protein = (
            matches[-1] if matches else row[column].split()[0]
        ).upper()
        if protein not in proteins:
            proteins.append(protein)
    if not proteins:
        raise ValueError("No protein accessions found")
    return proteins, len(data), first[column] if has_header else "column_2"


class NcbiClient:
    def __init__(self, email, api_key, delay, offline=False):
        self.email = email
        self.api_key = api_key
        self.delay = delay
        self.offline = offline
        self.last_request = 0.0

    def fetch(self, cache_file, endpoint, params):
        cache_file = Path(cache_file)
        if cache_file.is_file() and cache_file.stat().st_size:
            return cache_file.read_text(encoding="utf-8")
        if self.offline:
            raise FileNotFoundError(f"Missing cached response: {cache_file}")

        params = {
            **params, "email": self.email, "tool": "flankguard"
        }
        if self.api_key:
            params["api_key"] = self.api_key
        request = Request(
            BASE + endpoint + "?" + urlencode(params),
            headers={"User-Agent": f"FlankGuard/2.0 ({self.email})"},
        )
        for attempt in range(8):
            pause = self.delay - (time.monotonic() - self.last_request)
            if pause > 0:
                time.sleep(pause)
            try:
                with urlopen(request, timeout=120) as response:
                    text = response.read().decode("utf-8")
                self.last_request = time.monotonic()
                atomic_write(cache_file, text)
                return text
            except HTTPError as exc:
                self.last_request = time.monotonic()
                if exc.code not in {429, 500, 502, 503, 504}:
                    raise RuntimeError(f"NCBI HTTP {exc.code}") from exc
                reason = str(exc)
            except (
                URLError, TimeoutError, IncompleteRead, ConnectionError
            ) as exc:
                self.last_request = time.monotonic()
                reason = str(exc)
            if attempt == 7:
                raise RuntimeError(f"NCBI request failed: {reason}")
            wait = 2 ** attempt
            print(f"NCBI request failed: {reason}; retrying in {wait}s")
            time.sleep(wait)


def parse_ipg(xml_text, protein):
    root = ET.fromstring(xml_text)
    mappings = []

    def walk(element, current=None):
        tag = local_name(element.tag)
        if tag in {"protein", "proteinrecord"}:
            current = {
                "source": attr(element, "source", "database", "db"),
                "accession": attr(element, "accver", "accession", "protein_acc"),
                "name": attr(element, "name", "product", "protein_name"),
            }
        if tag in {"cds", "cdsrecord", "nucleotide"}:
            current = current or {}
            source = current.get("source", UNKNOWN)
            if source == UNKNOWN:
                source = attr(element, "source", "database", "db")
            nucleotide = attr(
                element, "accver", "nucleotide_acc", "nucleotide_accession",
                "nucaccver", "nuc_accession", "accession",
            )
            if source.lower() == "refseq" and nucleotide != UNKNOWN:
                mappings.append({
                    "protein_accession": protein,
                    "ipg_protein_accession": current.get("accession", protein),
                    "protein_name": current.get("name", UNKNOWN),
                    "nuccore_accession": nucleotide,
                    "cds_start": attr(element, "start", "from"),
                    "cds_stop": attr(element, "stop", "to", "end"),
                    "strand": attr(element, "strand"),
                    "organism": attr(element, "organism", "org"),
                    "strain": attr(element, "strain"),
                    "assembly": attr(
                        element, "assembly", "assembly_acc", "assembly_accession"
                    ),
                })
        for child in element:
            walk(child, current)

    walk(root)
    unique = {}
    for row in mappings:
        unique[(
            row["nuccore_accession"], row["cds_start"], row["cds_stop"]
        )] = row
    return list(unique.values())


def choose_representative(rows):
    if not rows:
        raise ValueError("No RefSeq IPG mappings")

    def rank(row):
        accession = row["nuccore_accession"].upper()
        accession_rank = (
            0 if accession.startswith("NC_") else
            1 if accession.startswith("NZ_CP") else
            2 if accession.startswith(("NZ_AP", "NZ_BA")) else
            3 if accession.startswith("NZ_") else 4
        )
        coordinates_valid = (
            integer(row["cds_start"]) is not None
            and integer(row["cds_stop"]) is not None
        )
        return (
            not coordinates_valid, accession_rank,
            row["assembly"] == UNKNOWN,
            row["organism"] == UNKNOWN,
            row["strain"] == UNKNOWN,
            accession,
        )

    selected = min(rows, key=rank).copy()
    if integer(selected["cds_start"]) is None:
        raise ValueError("Representative has no valid CDS coordinates")
    if integer(selected["cds_stop"]) is None:
        raise ValueError("Representative has no valid CDS coordinates")
    selected["ipg_refseq_candidates"] = str(len(rows))
    selected["representative_rule"] = (
        "Valid coordinates; NC_, NZ_CP, NZ_AP/NZ_BA, other NZ_; "
        "metadata completeness; lexical tie-break"
    )
    return selected


def requested_intervals(start, stop, flank, length, topology):
    desired_start, desired_stop = start - flank, stop + flank
    if not length:
        return [(max(1, desired_start), desired_stop)], "length unknown"
    if (
        topology.lower() == "circular"
        and desired_stop - desired_start + 1 < length
    ):
        if desired_start < 1 or desired_stop > length:
            return [
                (((desired_start - 1) % length) + 1, length),
                (1, ((desired_stop - 1) % length) + 1),
            ], "complete circular wrap"
        return [(desired_start, desired_stop)], "complete"
    clipped = (max(1, desired_start), min(length, desired_stop))
    status = (
        "complete" if clipped == (desired_start, desired_stop)
        else "truncated at linear replicon edge"
    )
    return [clipped], status


def qualifiers(feature):
    result = defaultdict(list)
    for item in feature.iter():
        if local_name(item.tag) not in {"gbqualifier", "insdqualifier"}:
            continue
        name = element_text(item, "_name")
        value = element_text(item, "_value")
        if name != UNKNOWN and value != UNKNOWN:
            result[name.lower()].append(value)
    return result


def parse_region(xml_text, representative, request_start):
    root = ET.fromstring(xml_text)
    record = next(
        (e for e in root.iter() if local_name(e.tag) in {"gbseq", "insdseq"}),
        None,
    )
    if record is None:
        raise ValueError("No GenBank record returned")
    sequence = element_text(record, "_sequence")
    if sequence == UNKNOWN:
        raise ValueError("No nucleotide sequence returned")
    sequence = re.sub(r"\s+", "", sequence).upper()
    target_ids = {
        representative["protein_accession"].split(".")[0],
        representative["ipg_protein_accession"].split(".")[0],
    }
    target_start, target_stop = sorted([
        int(representative["cds_start"]),
        int(representative["cds_stop"]),
    ])
    genes = []
    for feature in record.iter():
        if local_name(feature.tag) not in {"gbfeature", "insdfeature"}:
            continue
        if element_text(feature, "_key").lower() != "cds":
            continue
        spans, strand = [], "+"
        for interval in feature.iter():
            if local_name(interval.tag) not in {"gbinterval", "insdinterval"}:
                continue
            start = integer(element_text(interval, "_from"))
            stop = integer(element_text(interval, "_to"))
            point = integer(element_text(interval, "_point"))
            if start is None or stop is None:
                start = stop = point
            if start is None:
                continue
            spans.append(tuple(sorted((start, stop))))
            for child in interval:
                if local_name(child.tag).endswith("_iscomp"):
                    value = clean(child.text)
                    if value == UNKNOWN:
                        value = attr(child, "value")
                    if value.lower() in {"true", "1", "yes"}:
                        strand = "-"
        if not spans:
            continue
        q = qualifiers(feature)
        translation = re.sub(r"\s+", "", q.get("translation", [""])[0])
        if not translation:
            continue
        raw_start = min(x[0] for x in spans)
        raw_stop = max(x[1] for x in spans)
        local_start, local_stop = max(1, raw_start), min(len(sequence), raw_stop)
        if local_stop < local_start:
            continue
        genomic_start = request_start + raw_start - 1
        genomic_stop = request_start + raw_stop - 1
        protein_id = q.get("protein_id", q.get("locus_tag", [UNKNOWN]))[0]
        is_target = (
            protein_id.split(".")[0] in target_ids
            or (
                genomic_start <= target_stop
                and genomic_stop >= target_start
            )
        )
        genes.append({
            "local_start": local_start, "local_stop": local_stop,
            "genomic_start": genomic_start, "genomic_stop": genomic_stop,
            "strand": strand, "original_protein_id": protein_id,
            "locus_tag": q.get("locus_tag", [UNKNOWN])[0],
            "product": q.get("product", [UNKNOWN])[0],
            "translation": translation,
            "is_target": "yes" if is_target else "no",
        })
    if not genes:
        raise ValueError("No translated CDS features returned")
    genes.sort(key=lambda g: (g["local_start"], g["local_stop"]))
    return sequence, genes, {
        "returned_accession": element_text(record, "_accession-version"),
        "definition": element_text(record, "_definition"),
        "organism_from_record": element_text(record, "_organism"),
    }


def fetch_context(client, representative, flank, cache):
    accession = representative["nuccore_accession"]
    text = client.fetch(
        cache / "nuccore_summary" / f"{safe_name(accession)}.xml",
        "esummary.fcgi",
        {"db": "nuccore", "id": accession, "retmode": "xml"},
    )
    values = {}
    for element in ET.fromstring(text).iter():
        if local_name(element.tag) == "item":
            name = element.attrib.get("Name", element.attrib.get("name", ""))
            values[name.lower()] = clean(element.text)
    length = integer(values.get("length"))
    topology = values.get("topology", UNKNOWN)
    start, stop = sorted([
        int(representative["cds_start"]), int(representative["cds_stop"])
    ])
    intervals, edge = requested_intervals(start, stop, flank, length, topology)
    sequences, genes, pieces = [], [], []
    offset = 0
    for begin, end in intervals:
        xml_text = client.fetch(
            cache / "regions" / f"{safe_name(accession)}_{begin}_{end}.xml",
            "efetch.fcgi",
            {
                "db": "nuccore", "id": accession,
                "rettype": "gbwithparts", "retmode": "xml",
                "seq_start": begin, "seq_stop": end,
            },
        )
        sequence, features, metadata = parse_region(
            xml_text, representative, begin
        )
        for feature in features:
            feature["local_start"] += offset
            feature["local_stop"] += offset
            genes.append(feature)
        sequences.append(sequence)
        pieces.append(metadata)
        offset += len(sequence)
    for index, gene in enumerate(genes, 1):
        gene["feature_id"] = f"FG_{index:05d}"
    context = {
        **representative, **pieces[0],
        "request_start": str(intervals[0][0]),
        "requested_stop": str(intervals[-1][1]),
        "request_intervals": ";".join(f"{a}-{b}" for a, b in intervals),
        "returned_length_bp": str(offset),
        "replicon_length": str(length) if length else UNKNOWN,
        "replicon_topology": topology,
        "edge_status": edge,
        "coordinate_mode": "slice-relative-1-based",
        "target_cds_found": (
            "yes" if any(g["is_target"] == "yes" for g in genes) else "no"
        ),
        "cds_count": str(len(genes)),
        "flank_bp_each_side": flank,
    }
    return context, "".join(sequences), genes


def wrap(sequence):
    return "\n".join(sequence[i:i + 70] for i in range(0, len(sequence), 70))


def write_context(directory, context, sequence, genes):
    contig = safe_name(
        f"{context['protein_accession']}__{context['nuccore_accession']}"
    )
    atomic_write(directory / "neighborhood.fna", f">{contig}\n{wrap(sequence)}\n")
    atomic_write(
        directory / "neighborhood.faa",
        "".join(
            f">{g['feature_id']} {g['product']}\n{wrap(g['translation'])}\n"
            for g in genes
        ),
    )
    gff = ["##gff-version 3"]
    for gene in genes:
        attributes = (
            f"ID={gene['feature_id']};protein_id={gene['feature_id']};"
            f"original_protein_id={safe_name(gene['original_protein_id'])};"
            f"product={safe_name(gene['product'])}"
        )
        gff.append("\t".join(map(str, [
            contig, "FlankGuard", "CDS",
            gene["local_start"], gene["local_stop"],
            ".", gene["strand"], "0", attributes,
        ])))
    atomic_write(directory / "neighborhood.gff", "\n".join(gff) + "\n")
    write_tsv(directory / "genes.tsv", genes, [
        "feature_id", "local_start", "local_stop",
        "genomic_start", "genomic_stop", "strand",
        "original_protein_id", "locus_tag", "product", "is_target",
    ])
    context["context_dir"] = str(directory)
    atomic_write(directory / "context.json", json.dumps(context, indent=2) + "\n")


def cached_context(directory, protein, flank):
    if not all((directory / name).is_file() for name in CONTEXT_FILES):
        return None
    context = read_json(directory / "context.json")
    if context.get("protein_accession") != protein:
        return None
    start = integer(context.get("cds_start"))
    stop = integer(context.get("cds_stop"))
    if start is None or stop is None:
        return None
    intervals, _ = requested_intervals(
        min(start, stop), max(start, stop), flank,
        integer(context.get("replicon_length")),
        context.get("replicon_topology", UNKNOWN),
    )
    expected = ";".join(f"{a}-{b}" for a, b in intervals)
    recorded_flank = integer(context.get("flank_bp_each_side"))
    if recorded_flank is not None and recorded_flank != flank:
        return None
    if context.get("request_intervals") != expected:
        return None
    context["context_dir"] = str(directory)
    return context


def prepare_contexts(args, proteins, client):
    prior_summary = read_json(args.outdir / "summary.json", {}) or {}
    if prior_summary:
        previous_flank = integer(prior_summary.get("flank_bp_each_side"))
        if previous_flank is not None and previous_flank != args.flank_bp:
            raise RuntimeError(
                "This output folder contains a different flank size. "
                "Use a new output folder for different extraction parameters."
            )

    prior_failures = read_tsv(args.outdir / "failures.tsv", optional=True)
    failed = {
        row["protein_accession"]: row for row in prior_failures
        if row.get("stage") in {"ipg", "context"}
    }
    representatives = {
        row["protein_accession"]: row
        for row in read_tsv(args.outdir / "representatives.tsv", optional=True)
    }
    contexts, failures = [], []
    reused = 0

    for index, protein in enumerate(proteins, 1):
        directory = args.outdir / "neighborhoods" / safe_name(protein)
        context = cached_context(directory, protein, args.flank_bp)
        if context is not None:
            contexts.append(context)
            reused += 1
            continue

        if (
            protein in failed and prior_summary
            and not args.retry_failed
        ):
            failures.append(failed[protein])
            print(f"[Context] Retaining previous retrieval failure: {protein}")
            continue

        if args.replot_only:
            raise RuntimeError(
                f"No compatible completed context for {protein}. "
                "Remove --replot-only to prepare missing inputs."
            )

        if any(
            (directory / tool / ".flankguard_complete").exists()
            for tool in ("defensefinder", "padloc")
        ):
            raise RuntimeError(
                f"Context inputs changed or are incomplete for {protein}, "
                "but scanner results already exist. Use a fresh output folder "
                "to avoid reusing results for different sequences."
            )

        print(f"[Context {index}/{len(proteins)}] Preparing {protein}", flush=True)
        stage = "ipg"
        try:
            representative = representatives.get(protein)
            if representative is None:
                xml_text = client.fetch(
                    args.outdir / "cache" / "ipg" / f"{safe_name(protein)}.xml",
                    "efetch.fcgi",
                    {
                        "db": "protein", "id": protein,
                        "rettype": "ipg", "retmode": "xml",
                    },
                )
                representative = choose_representative(parse_ipg(xml_text, protein))
                representatives[protein] = representative
            stage = "context"
            context, sequence, genes = fetch_context(
                client, representative, args.flank_bp, args.outdir / "cache"
            )
            write_context(directory, context, sequence, genes)
            contexts.append(context)
        except Exception as exc:
            print(f"  FAILED: {exc}", file=sys.stderr)
            failures.append({
                "stage": stage, "protein_accession": protein, "error": str(exc)
            })

    selected_representatives = [
        representatives[p] for p in proteins if p in representatives
    ]
    write_tsv(
        args.outdir / "representatives.tsv", selected_representatives,
        [
            "protein_accession", "ipg_protein_accession", "protein_name",
            "nuccore_accession", "cds_start", "cds_stop", "strand",
            "organism", "strain", "assembly",
            "ipg_refseq_candidates", "representative_rule",
        ],
    )
    write_tsv(args.outdir / "neighborhood_manifest.tsv", contexts)
    print(
        f"Contexts: {reused} reused, {len(contexts) - reused} prepared, "
        f"{len(failures)} unavailable"
    )
    return contexts, selected_representatives, failures, reused


def conda_prefix(conda, environment):
    expanded = os.path.expandvars(os.path.expanduser(environment.strip()))
    if not expanded:
        raise RuntimeError("Conda environment must not be empty")
    path = Path(expanded)
    is_path = (
        path.is_absolute() or "/" in expanded or "\\" in expanded
        or (path / "conda-meta").is_dir()
    )
    if is_path:
        path = path.resolve()
        if not (path / "conda-meta").is_dir():
            raise RuntimeError(
                f"Not a Conda environment root: {path}. "
                "Supply the environment root, not its bin directory."
            )
        selector = ["-p", str(path)]
    else:
        selector = ["-n", expanded]
    return [
        os.path.expandvars(os.path.expanduser(conda)),
        "run", "--no-capture-output", *selector,
        "/bin/bash", "-c",
        'export PATH="$CONDA_PREFIX/bin:$PATH"; exec "$@"',
        "flankguard-environment",
    ]


def scanner_signature(context, tool):
    directory = Path(context["context_dir"])
    parts = [file_hash(directory / "neighborhood.faa")]
    if tool == "padloc":
        parts.append(file_hash(directory / "neighborhood.gff"))
    return hashlib.sha256("|".join(parts).encode()).hexdigest()


def scanner_complete(context, tool):
    directory = Path(context["context_dir"]) / tool
    marker = directory / ".flankguard_complete"
    if not marker.is_file():
        return False
    text = marker.read_text(encoding="utf-8").strip()
    if text.startswith("{"):
        saved = json.loads(text)
        if saved.get("input_sha256") != scanner_signature(context, tool):
            raise RuntimeError(
                f"Scanner inputs changed for {context['protein_accession']} "
                f"({tool}); use a fresh output folder."
            )
    elif text != "complete":
        return False

    if tool == "defensefinder":
        return any(directory.rglob("*systems.tsv"))

    if any(directory.rglob("*_padloc.csv")):
        return True

    # PADLOC can complete successfully without creating a CSV when no
    # systems are found. That is not a failed or incomplete scan.
    log = directory / "flankguard.log"
    return (
        log.is_file()
        and "nothing found" in log.read_text(
            encoding="utf-8", errors="replace"
        ).lower()
    )


def scanner_command(args, context, tool, destination):
    directory = Path(context["context_dir"])
    if tool == "defensefinder":
        return conda_prefix(args.conda, args.defensefinder_env) + [
            "defense-finder", "run", str(directory / "neighborhood.faa"),
            "--out-dir", str(destination),
            "--workers", str(args.threads),
            "--antidefensefinder",
            "--index-dir", str(destination / "index"),
        ]
    return conda_prefix(args.conda, args.padloc_env) + [
        "padloc",
        "--faa", str(directory / "neighborhood.faa"),
        "--gff", str(directory / "neighborhood.gff"),
        "--outdir", str(destination), "--cpu", str(args.threads),
    ]


def preflight(args, tools):
    if args.skip_tool_checks or args.dry_run:
        return
    for tool in sorted(tools):
        if tool == "defensefinder":
            command = conda_prefix(args.conda, args.defensefinder_env) + [
                "defense-finder", "run", "--help"
            ]
        else:
            command = conda_prefix(args.conda, args.padloc_env) + [
                "padloc", "--check-deps"
            ]
        result = subprocess.run(command, capture_output=True, text=True)
        output = (result.stdout or "") + (result.stderr or "")
        if result.returncode or (
            tool == "padloc"
            and re.search(r"warning:.*not installed", output, re.I)
        ):
            raise RuntimeError(f"{tool} preflight failed:\n{output}")


def run_scanner(args, context, tool):
    parent = Path(context["context_dir"])
    destination = parent / tool
    staging = parent / f".{tool}.running-{time.time_ns()}"
    staging.mkdir(parents=True)
    command = scanner_command(args, context, tool, staging)
    log = staging / "flankguard.log"
    print("  $ " + shlex.join(command), flush=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\n" + shlex.join(command) + "\n\nOUTPUT\n")
        handle.flush()
        result = subprocess.run(
            command, stdout=handle, stderr=subprocess.STDOUT, text=True
        )
    if result.returncode:
        raise RuntimeError(
            f"{tool} exited with code {result.returncode}; see {log}"
        )
    atomic_write(
        staging / ".flankguard_complete",
        json.dumps({
            "input_sha256": scanner_signature(context, tool),
            "tool": tool,
        }) + "\n",
    )
    # Replace this scanner's old generated outputs only after success.
    if destination.exists():
        shutil.rmtree(destination)
    staging.rename(destination)


def first_value(row, *names):
    values = {key(k): v for k, v in row.items()}
    for name in names:
        value = values.get(key(name), "")
        if value is not None and str(value).strip():
            return str(value).strip()
    return UNKNOWN


def parse_scanner(context, tool):
    directory = Path(context["context_dir"]) / tool
    calls = []
    common = {
        "protein_accession": context["protein_accession"],
        "nuccore_accession": context["nuccore_accession"],
    }
    if tool == "defensefinder":
        files = sorted(set(directory.rglob("*systems.tsv")))
        for path in files:
            for row in read_tsv(path):
                activity = first_value(row, "activity").lower()
                label = (
                    "AntiDefenseFinder"
                    if "antidefense" in activity or "anti-defense" in activity
                    else "DefenseFinder"
                )
                calls.append({
                    **common, "tool": label,
                    "system": first_value(
                        row, "subtype", "type", "system", "system_name",
                        "model_fqn", "model",
                    ),
                    "system_id": first_value(row, "sys_id", "system_id", "id"),
                    "result_file": str(path),
                    "raw_result": json.dumps(row, sort_keys=True),
                })
    else:
        for path in sorted(directory.rglob("*_padloc.csv")):
            with path.open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            groups = defaultdict(list)
            for row in rows:
                groups[(
                    first_value(row, "system", "system_name"),
                    first_value(row, "system.number", "system_number", "system.id"),
                )].append(row)
            for (system, system_id), members in groups.items():
                calls.append({
                    **common, "tool": "PADLOC", "system": system,
                    "system_id": system_id, "result_file": str(path),
                    "raw_result": json.dumps(members, sort_keys=True),
                })
    return calls


def scan_contexts(args, contexts):
    tools = []
    if not args.skip_defensefinder:
        tools.append("defensefinder")
    if not args.skip_padloc:
        tools.append("padloc")

    completed = {
        (c["protein_accession"], tool): scanner_complete(c, tool)
        for c in contexts for tool in tools
    }
    pending_tools = {
        tool for (_, tool), done in completed.items() if not done
    }
    if args.replot_only and pending_tools:
        raise RuntimeError(
            "Some required scanner results are incomplete. "
            "Remove --replot-only to run missing scans."
        )
    preflight(args, pending_tools)

    def work(context):
        protein = context["protein_accession"]
        calls, failures, status = [], [], {}
        reused, executed = 0, 0
        for tool in tools:
            try:
                if completed[(protein, tool)]:
                    reused += 1
                elif args.dry_run:
                    command = scanner_command(
                        args, context, tool, Path(context["context_dir"]) / tool
                    )
                    print("DRY RUN: " + shlex.join(command), flush=True)
                    status[tool] = "not_run"
                    continue
                else:
                    run_scanner(args, context, tool)
                    executed += 1
                calls.extend(parse_scanner(context, tool))
                status[tool] = "complete"
            except Exception as exc:
                status[tool] = "failed"
                failures.append({
                    "stage": tool, "protein_accession": protein, "error": str(exc)
                })
                print(f"{protein}: {exc}", file=sys.stderr)
        return protein, calls, failures, status, reused, executed

    all_calls, failures, statuses = [], [], {}
    reused_total = executed_total = 0
    if args.jobs == 1:
        results = map(work, contexts)
        for protein, calls, failed, status, reused, executed in results:
            all_calls.extend(calls)
            failures.extend(failed)
            statuses[protein] = status
            reused_total += reused
            executed_total += executed
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [executor.submit(work, c) for c in contexts]
            for future in as_completed(futures):
                protein, calls, failed, status, reused, executed = future.result()
                all_calls.extend(calls)
                failures.extend(failed)
                statuses[protein] = status
                reused_total += reused
                executed_total += executed

    all_calls.sort(key=lambda r: (
        r["protein_accession"], r["tool"], r["system"], r["system_id"]
    ))
    print(f"Scanner runs: {reused_total} reused, {executed_total} executed")
    return all_calls, failures, statuses, reused_total, executed_total


def canonical_system(name):
    value = key(name)

    # The methyltransferase/nuclease prefixes are model detail, not
    # separate top-level RM types. Match IIG before II or I.
    match = re.search(
        r"(?:^|-)rm-(?:type-)?(iig|iii|iv|ii|i)(?=-|$)", value
    )
    if match:
        return "RM type " + match.group(1).upper()

    aliases = {
        "septu": "Septu type I",
        "septu-i": "Septu type I",
        "septu-type-i": "Septu type I",
        "cas-class1-subtype-i-c": "CRISPR-Cas type I-C",
        "cas-type-i-c": "CRISPR-Cas type I-C",
        "cas-i-c": "CRISPR-Cas type I-C",
        "gao-hhe": "Hhe",
        "hhe": "Hhe",
        "gao-ppl": "Ppl",
        "ppl": "Ppl",
        "pago-s1b": "Argonaute type I",
        "argonaute-type-i": "Argonaute type I",
        "argonaute-i": "Argonaute type I",
    }
    if value in aliases:
        return aliases[value]
    value = re.sub(
        r"-type-(iis|iv|iii|ii|i)(?=-|$)",
        lambda m: "-" + ("ii" if m.group(1) == "iis" else m.group(1)),
        value,
    )
    labels = {
        "cbass-i": "CBASS type I", "cbass-ii": "CBASS type II",
        "wadjet-i": "Wadjet type I", "gabija": "Gabija",
    }
    return labels.get(value, name.replace("_", " "))


def call_details(raw_calls, contexts):
    context_map = {c["protein_accession"]: c for c in contexts}
    gene_cache = {}
    details = []

    for call in raw_calls:
        protein = call["protein_accession"]
        context = context_map[protein]
        if protein not in gene_cache:
            rows = read_tsv(Path(context["context_dir"]) / "genes.tsv")
            gene_cache[protein] = (
                {r["feature_id"]: r for r in rows},
                next((r for r in rows if r.get("is_target") == "yes"), None),
            )
        genes, target = gene_cache[protein]
        raw = json.loads(call["raw_result"])

        if call["tool"] == "PADLOC":
            members = raw if isinstance(raw, list) else [raw]
            feature_ids = list(dict.fromkeys(
                r.get("target.name", "") for r in members
                if r.get("target.name")
            ))
            profiles = joined(r.get("hmm.name", "") for r in members)
            accessions = joined(r.get("hmm.accession", "") for r in members)
        else:
            feature_ids = list(dict.fromkeys(re.findall(
                r"FG_\d+", str(raw.get("protein_in_syst", ""))
            )))
            profiles = raw.get("name_of_profiles_in_sys", "")
            accessions = ""

        mapped = [genes[f] for f in feature_ids if f in genes]
        missing = [f for f in feature_ids if f not in genes]
        start = min((int(g["local_start"]) for g in mapped), default=None)
        stop = max((int(g["local_stop"]) for g in mapped), default=None)
        gap, position = "", "unavailable"
        if target is not None and start is not None:
            target_start, target_stop = (
                int(target["local_start"]), int(target["local_stop"])
            )
            if stop < target_start:
                gap = target_start - stop - 1
                position = "upstream" if target["strand"] == "+" else "downstream"
            elif start > target_stop:
                gap = start - target_stop - 1
                position = "downstream" if target["strand"] == "+" else "upstream"
            else:
                gap, position = 0, "overlaps target CDS"

        details.append({
            "protein_accession": protein,
            "nuccore_accession": call["nuccore_accession"],
            "organism": context.get("organism", UNKNOWN),
            "strain": context.get("strain", UNKNOWN),
            "assembly": context.get("assembly", UNKNOWN),
            "program": call["tool"],
            "call_class": (
                "anti-defence" if call["tool"] == "AntiDefenseFinder" else "defence"
            ),
            "source_system": call["system"],
            "canonical_system": canonical_system(call["system"]),
            "system_id": call["system_id"],
            "model_profiles": profiles,
            "model_accessions": accessions,
            "local_feature_ids": joined(feature_ids),
            "mapped_gene_count": len(mapped),
            "coordinate_mapping_status": (
                "complete" if mapped and not missing else
                "partial" if mapped else "unmapped"
            ),
            "system_local_start": start if start is not None else "",
            "system_local_stop": stop if stop is not None else "",
            "position_relative_to_target": position,
            "nearest_distance_to_target_bp": gap,
            "RefSeq_system_protein_accessions": joined(
                g["original_protein_id"] for g in mapped
            ),
            "locus_tags": joined(g["locus_tag"] for g in mapped),
            "gene_products": joined(g["product"] for g in mapped),
            "gene_coordinates_in_neighborhood": joined(
                f"{g['feature_id']}:{g['local_start']}-{g['local_stop']}({g['strand']})"
                for g in mapped
            ),
            "gene_coordinates_on_replicon": joined(
                f"{g['feature_id']}:{g['genomic_start']}-{g['genomic_stop']}({g['strand']})"
                for g in mapped
            ),
            "retrieved_intervals": context.get("request_intervals", ""),
            "edge_status": context.get("edge_status", ""),
            "result_file": call["result_file"],
        })
    return details


def merge_loci(calls):
    groups = defaultdict(list)
    for row in calls:
        groups[(
            row["protein_accession"], row["call_class"],
            key(row["canonical_system"]),
        )].append(row)

    loci = []
    for _, rows in sorted(groups.items()):
        parents = list(range(len(rows)))

        def find(index):
            while parents[index] != index:
                parents[index] = parents[parents[index]]
                index = parents[index]
            return index

        feature_sets = [
            set(r["local_feature_ids"].split("; ")) - {""} for r in rows
        ]
        for left in range(len(rows)):
            for right in range(left + 1, len(rows)):
                # Same canonical type alone is insufficient.
                # Require an actual shared gene, not just overlapping spans.
                if not feature_sets[left] & feature_sets[right]:
                    continue
                # Separate calls from one scanner remain separate unless
                # they contain precisely the same reported gene set.
                if (
                    rows[left]["program"] == rows[right]["program"]
                    and feature_sets[left] != feature_sets[right]
                ):
                    continue
                a, b = find(left), find(right)
                if a != b:
                    parents[b] = a

        clusters = defaultdict(list)
        for index, row in enumerate(rows):
            clusters[find(index)].append(row)

        for cluster in clusters.values():
            first = cluster[0]
            starts = [
                int(r["system_local_start"]) for r in cluster
                if r["system_local_start"] != ""
            ]
            stops = [
                int(r["system_local_stop"]) for r in cluster
                if r["system_local_stop"] != ""
            ]
            distances = [
                int(r["nearest_distance_to_target_bp"]) for r in cluster
                if r["nearest_distance_to_target_bp"] != ""
            ]
            programs = sorted({r["program"] for r in cluster})
            loci.append({
                "protein_accession": first["protein_accession"],
                "nuccore_accession": first["nuccore_accession"],
                "organism": first["organism"],
                "assembly": first["assembly"],
                "call_class": first["call_class"],
                "consensus_system": first["canonical_system"],
                "supporting_programs": joined(programs),
                "source_systems": joined(r["source_system"] for r in cluster),
                "tool_specific_call_count": len(cluster),
                "local_feature_ids": joined(
                    f for r in cluster for f in r["local_feature_ids"].split("; ")
                ),
                "consensus_local_start": min(starts) if starts else "",
                "consensus_local_stop": max(stops) if stops else "",
                "nearest_distance_to_target_bp": (
                    min(distances) if distances else ""
                ),
                "position_relative_to_target": joined(
                    r["position_relative_to_target"] for r in cluster
                ),
                "RefSeq_system_protein_accessions": joined(
                    f for r in cluster
                    for f in r["RefSeq_system_protein_accessions"].split("; ")
                ),
                "gene_coordinates_on_replicon": joined(
                    f for r in cluster
                    for f in r["gene_coordinates_on_replicon"].split("; ")
                ),
                "coordinate_mapping_status": joined(
                    r["coordinate_mapping_status"] for r in cluster
                ),
                "merge_rule": "same canonical system and shared reported gene",
            })

    loci.sort(key=lambda r: (
        r["protein_accession"],
        integer(r["consensus_local_start"])
        if integer(r["consensus_local_start"]) is not None else math.inf,
        r["call_class"], r["consensus_system"],
    ))
    counters = Counter()
    for row in loci:
        counters[row["protein_accession"]] += 1
        row["consensus_locus_id"] = (
            f"{row['protein_accession']}_L{counters[row['protein_accession']]:03d}"
        )
    return loci


def hex_colour(value):
    if not re.fullmatch(r"#[0-9a-fA-F]{6}", value):
        raise ValueError("Supply a colour such as #264979")
    return value


def tint(colour, fraction):
    channels = [int(colour[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(
        f"{round(c + (255 - c) * fraction):02x}" for c in channels
    )


def svg_begin(width, height, title, subtitle):
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;fill:#111}'
        '.title{font-size:24px;font-weight:bold}'
        '.note{font-size:13px;fill:#444}'
        '.label{font-size:14px}.small{font-size:11px}'
        '.value{font-size:13px;font-weight:bold}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="30" y="35" class="title">{html.escape(title)}</text>',
        f'<text x="30" y="60" class="note">{html.escape(subtitle)}</text>',
    ]


def frequency_svg(rows, title, denominator, colour, top, strip=False):
    shown = rows[:top]
    left, top_y, step = 470, 100, 32
    height = top_y + max(1, len(shown)) * step + 45
    parts = svg_begin(
        1100, height, title,
        f"Unique analysed neighbourhoods; denominator n={denominator}.",
    )
    if not shown:
        parts.append('<text x="30" y="100">No retained calls.</text>')
    for index, row in enumerate(shown):
        y = top_y + index * step
        label = row["plot_label"]
        percent = float(row["percent"])
        count = row["neighborhood_count"]
        parts.append(
            f'<text x="{left - 15}" y="{y + 15}" text-anchor="end" '
            f'class="label">{html.escape(label)}</text>'
        )
        if strip:
            fill = tint(colour, 1 - min(1.0, percent / 100))
            parts.append(
                f'<rect x="{left}" y="{y}" width="24" height="24" '
                f'fill="{fill}" stroke="#222"/>'
            )
            value_x = left + 40
        else:
            width = 500 * percent / 100
            parts.append(
                f'<rect x="{left}" y="{y}" width="{width:.2f}" '
                f'height="22" fill="{colour}"/>'
            )
            value_x = left + width + 12
        parts.append(
            f'<text x="{value_x:.2f}" y="{y + 16}" class="value">'
            f'{count} ({percent:.2f}%)</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def presence_heatmap(summary, presence, labels, title, colour, max_rows):
    selected = [
        r for r in summary
        if r["analysis_status"] == "analysed"
        and presence.get(r["protein_accession"], set()) & set(labels)
    ]
    selected.sort(key=lambda r: (
        -len(presence.get(r["protein_accession"], set()) & set(labels)),
        r["organism"], r["protein_accession"],
    ))
    shown = selected[:max_rows]
    cell, left, top = 20, 440, 290
    width = max(1000, left + len(labels) * cell + 30)
    height = top + max(1, len(shown)) * cell + 60
    parts = svg_begin(
        width, height, title,
        f"Showing {len(shown)} of {len(selected)} neighbourhoods "
        "with a displayed system; white means not detected.",
    )
    for column, label in enumerate(labels):
        x = left + column * cell + 10
        parts.append(
            f'<text x="{x}" y="{top - 10}" class="small" '
            f'transform="rotate(-58 {x} {top - 10})">'
            f'{html.escape(label)}</text>'
        )
    for index, row in enumerate(shown):
        y = top + index * cell
        protein = row["protein_accession"]
        text = f"{row['organism']} | {protein}"
        parts.append(
            f'<text x="{left - 10}" y="{y + 14}" text-anchor="end" '
            f'class="small">{html.escape(text)}</text>'
        )
        for column, label in enumerate(labels):
            fill = colour if label in presence.get(protein, set()) else "white"
            parts.append(
                f'<rect x="{left + column * cell}" y="{y}" '
                f'width="{cell}" height="{cell}" fill="{fill}" '
                'stroke="#555" stroke-width=".5"/>'
            )
    if not shown:
        parts.append('<text x="30" y="110">No retained calls.</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def genus_heatmap(summary, presence, labels, target, colour):
    groups = defaultdict(list)
    for row in summary:
        if row["analysis_status"] != "analysed":
            continue
        organism = row["organism"].strip()
        genus = organism.split()[0] if organism else UNKNOWN
        groups[genus].append(row["protein_accession"])
    ordered = sorted(groups, key=lambda g: (-len(groups[g]), g))
    displayed = [(g, groups[g]) for g in ordered[:12]]
    if len(ordered) > 12:
        displayed.append((
            "Other genera",
            [p for g in ordered[12:] for p in groups[g]],
        ))
    labels = labels[:15]
    left, top, cw, ch = 265, 290, 50, 29
    parts = svg_begin(
        1123, 794, f"Systems near {target} across genera",
        "Within-genus prevalence (%); n is the number of analysed neighbourhoods.",
    )
    for column, label in enumerate(labels):
        x = left + column * cw + 20
        parts.append(
            f'<text x="{x}" y="{top - 12}" class="small" '
            f'transform="rotate(-58 {x} {top - 12})">'
            f'{html.escape(label)}</text>'
        )
    for index, (genus, proteins) in enumerate(displayed):
        y = top + index * ch
        parts.append(
            f'<text x="{left - 12}" y="{y + 19}" text-anchor="end" '
            f'class="label">{html.escape(genus)} (n={len(proteins)})</text>'
        )
        for column, label in enumerate(labels):
            percent = 100 * sum(
                label in presence.get(p, set()) for p in proteins
            ) / len(proteins)
            x = left + column * cw
            fill = tint(colour, 1 - percent / 100)
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" '
                f'fill="{fill}" stroke="#555" stroke-width=".5"/>'
            )
            if percent:
                text_colour = "white" if percent >= 55 else "#111"
                parts.append(
                    f'<text x="{x + cw / 2}" y="{y + 19}" '
                    f'text-anchor="middle" '
                    f'style="fill:{text_colour};font-size:11px">{percent:.0f}</text>'
                )
    parts.append("</svg>")
    return "\n".join(parts)


def environment_svg(positive, total, target, colour):
    percent = 100 * positive / total if total else 0
    parts = svg_begin(
        1000, 245, f"{target}: predicted defence/anti-defence environments",
        f"At least one retained call in {positive} of {total} analysed "
        f"neighbourhoods ({percent:.2f}%).",
    )
    parts.extend([
        '<rect x="40" y="105" width="900" height="50" fill="#eeeeee"/>',
        f'<rect x="40" y="105" width="{9 * percent:.2f}" '
        f'height="50" fill="{colour}"/>',
        '<text x="40" y="195" class="note">'
        'VSPR and DMS-other are excluded. Failed/incomplete scans '
        'are not treated as negative results.</text>',
        "</svg>",
    ])
    return "\n".join(parts)


def breadth_svg(summary, loci, target):
    analysed = [r for r in summary if r["analysis_status"] == "analysed"]
    positive = [r for r in analysed if r["defence_environment"] == "yes"]
    counts = Counter(r["protein_accession"] for r in loci)
    multiple = sum(counts[r["protein_accession"]] >= 2 for r in positive)
    genera = {
        r["organism"].split()[0] for r in positive if r["organism"].split()
    }
    parts = svg_begin(
        1000, 270, f"Breadth of system association near {target}",
        "Observed counts; no predefined claim of a majority.",
    )
    for index, text in enumerate([
        f"{len(positive)} of {len(analysed)} analysed neighbourhoods have retained calls.",
        f"{len(genera)} genera are represented among positive neighbourhoods.",
        f"{multiple} positive neighbourhoods contain at least two consensus loci.",
        "Shared-gene evidence is required to merge tool calls into a locus.",
    ]):
        parts.append(
            f'<text x="35" y="{110 + index * 34}" class="label">'
            f'{html.escape(text)}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def rebuild_reports(args, contexts, raw_calls, statuses):
    # Preserve every original scanner call before exclusions.
    write_tsv(args.outdir / "combined_systems.tsv", raw_calls, [
        "protein_accession", "nuccore_accession", "tool", "system",
        "system_id", "result_file", "raw_result",
    ])
    all_details = call_details(raw_calls, contexts)
    excluded, retained = [], []
    for row in all_details:
        reason = EXCLUDED_SYSTEMS.get(key(row["source_system"]))
        if reason:
            excluded.append({**row, "exclusion_reason": reason})
        else:
            retained.append(row)

    write_tsv(args.outdir / "excluded_system_calls.tsv", excluded, [
        "protein_accession", "nuccore_accession", "program",
        "source_system", "canonical_system", "exclusion_reason",
        "local_feature_ids", "result_file",
    ])
    write_tsv(args.outdir / "supplementary_system_calls.tsv", retained)

    by_protein = defaultdict(list)
    for row in retained:
        by_protein[row["protein_accession"]].append(row)

    summary = []
    for context in contexts:
        protein = context["protein_accession"]
        state = statuses.get(protein, {})
        if args.dry_run:
            analysis_status = "dry_run"
        elif not state:
            analysis_status = "no_scanners_selected"
        elif any(v != "complete" for v in state.values()):
            analysis_status = "scanner_incomplete"
        elif context.get("target_cds_found") != "yes":
            analysis_status = "target_not_found"
        else:
            analysis_status = "analysed"

        calls = by_protein[protein]
        labels = {
            (r["call_class"], r["canonical_system"]) for r in calls
        }
        summary.append({
            "protein_accession": protein,
            "nuccore_accession": context["nuccore_accession"],
            "organism": context.get("organism", UNKNOWN),
            "assembly": context.get("assembly", UNKNOWN),
            "analysis_status": analysis_status,
            "scanner_status": joined(f"{t}:{v}" for t, v in sorted(state.items())),
            "defence_environment": (
                "yes" if labels else
                "no" if analysis_status == "analysed" else "unknown"
            ),
            "total_systems": len(labels),
            "defence_system_count": sum(c == "defence" for c, _ in labels),
            "anti_defence_system_count": sum(c == "anti-defence" for c, _ in labels),
            "systems": joined(f"{c}: {name}" for c, name in sorted(labels)),
            "target_cds_found": context.get("target_cds_found", "no"),
            "edge_status": context.get("edge_status", ""),
        })

    analysed_ids = {
        r["protein_accession"] for r in summary
        if r["analysis_status"] == "analysed"
    }
    analysed_calls = [
        r for r in retained if r["protein_accession"] in analysed_ids
    ]
    loci = merge_loci(retained)
    analysed_loci = [
        r for r in loci if r["protein_accession"] in analysed_ids
    ]
    locus_counts = Counter(r["protein_accession"] for r in loci)
    for row in summary:
        row["nonredundant_consensus_locus_count"] = locus_counts[
            row["protein_accession"]
        ]

    write_tsv(args.outdir / "neighborhood_summary.tsv", summary)
    write_tsv(args.outdir / "supplementary_neighborhoods.tsv", summary)
    write_tsv(args.outdir / "supplementary_consensus_loci.tsv", loci, [
        "consensus_locus_id", "protein_accession", "nuccore_accession",
        "organism", "assembly", "call_class", "consensus_system",
        "supporting_programs", "source_systems", "tool_specific_call_count",
        "local_feature_ids", "consensus_local_start", "consensus_local_stop",
        "nearest_distance_to_target_bp", "position_relative_to_target",
        "RefSeq_system_protein_accessions", "gene_coordinates_on_replicon",
        "coordinate_mapping_status", "merge_rule",
    ])

    total = len(analysed_ids)
    sets = defaultdict(set)
    tool_sets = defaultdict(set)
    for row in analysed_calls:
        sets[(row["call_class"], row["canonical_system"])].add(
            row["protein_accession"]
        )
        tool_sets[(row["program"], row["source_system"])].add(
            row["protein_accession"]
        )

    frequencies = [
        {
            "call_class": call_class,
            "consensus_system": system,
            "neighborhood_count": len(proteins),
            "percent_of_analysed_neighborhoods": (
                f"{100 * len(proteins) / total:.2f}" if total else "0.00"
            ),
        }
        for (call_class, system), proteins in sets.items()
    ]
    frequencies.sort(key=lambda r: (
        -r["neighborhood_count"], r["call_class"], r["consensus_system"]
    ))
    write_tsv(args.outdir / "consensus_system_frequencies.tsv", frequencies, [
        "call_class", "consensus_system", "neighborhood_count",
        "percent_of_analysed_neighborhoods",
    ])

    tool_frequencies = [
        {
            "tool": tool, "system": system,
            "neighborhood_count": len(proteins),
            "percent_of_neighborhoods": (
                f"{100 * len(proteins) / total:.2f}" if total else "0.00"
            ),
        }
        for (tool, system), proteins in tool_sets.items()
    ]
    tool_frequencies.sort(key=lambda r: (
        -r["neighborhood_count"], r["tool"], r["system"]
    ))
    write_tsv(args.outdir / "system_frequencies.tsv", tool_frequencies, [
        "tool", "system", "neighborhood_count", "percent_of_neighborhoods",
    ])

    presence, tool_presence = defaultdict(set), defaultdict(set)
    for (call_class, system), proteins in sets.items():
        for protein in proteins:
            presence[protein].add(f"{call_class}: {system}")
    for (tool, system), proteins in tool_sets.items():
        for protein in proteins:
            tool_presence[protein].add(f"{tool}: {system}")

    all_labels = [
        f"{r['call_class']}: {r['consensus_system']}" for r in frequencies
    ]
    matrix = [
        {
            "protein_accession": r["protein_accession"],
            "nuccore_accession": r["nuccore_accession"],
            "organism": r["organism"],
            **{
                label: int(label in presence[r["protein_accession"]])
                for label in all_labels
            },
        }
        for r in summary if r["analysis_status"] == "analysed"
    ]
    write_tsv(args.outdir / "system_presence_matrix.tsv", matrix, [
        "protein_accession", "nuccore_accession", "organism", *all_labels
    ])

    plot_rows = [
        {
            "plot_label": f"{r['call_class']}: {r['consensus_system']}",
            "neighborhood_count": r["neighborhood_count"],
            "percent": r["percent_of_analysed_neighborhoods"],
        }
        for r in frequencies
    ]
    tool_plot_rows = [
        {
            "plot_label": f"{r['tool']}: {r['system']}",
            "neighborhood_count": r["neighborhood_count"],
            "percent": r["percent_of_neighborhoods"],
        }
        for r in tool_frequencies
    ]
    target, colour = args.target_label, args.heatmap_color

    atomic_write(
        args.outdir / "system_frequencies.svg",
        frequency_svg(
            plot_rows, f"Systems near {target}", total, colour,
            args.heatmap_top_systems,
        ),
    )
    for call_class, filename in [
        ("defence", "defence_system_frequencies.svg"),
        ("anti-defence", "antidefence_system_frequencies.svg"),
    ]:
        selected = [
            row for row in plot_rows
            if row["plot_label"].startswith(call_class + ": ")
        ]
        atomic_write(
            args.outdir / filename,
            frequency_svg(
                selected, f"{call_class.capitalize()} systems near {target}",
                total, colour, args.heatmap_top_systems,
            ),
        )
    atomic_write(
        args.outdir / "tool_specific_system_frequencies.svg",
        frequency_svg(
            tool_plot_rows, f"Tool-specific calls near {target}",
            total, colour, args.heatmap_top_systems,
        ),
    )
    atomic_write(
        args.outdir / "tool_specific_system_frequency_strip.svg",
        frequency_svg(
            tool_plot_rows, f"Tool-specific calls near {target}",
            total, colour, args.heatmap_top_systems, strip=True,
        ),
    )
    atomic_write(
        args.outdir / "system_heatmap.svg",
        presence_heatmap(
            summary, presence, all_labels[:args.heatmap_top_systems],
            f"System presence near {target}", colour,
            args.heatmap_max_neighborhoods,
        ),
    )
    tool_labels = [
        f"{r['tool']}: {r['system']}"
        for r in tool_frequencies[:args.heatmap_top_systems]
    ]
    atomic_write(
        args.outdir / "tool_specific_system_heatmap.svg",
        presence_heatmap(
            summary, tool_presence, tool_labels,
            f"Tool-specific system presence near {target}", colour,
            args.heatmap_max_neighborhoods,
        ),
    )
    atomic_write(
        args.outdir / "system_heatmap_A4.svg",
        genus_heatmap(summary, presence, all_labels, target, colour),
    )

    positive = sum(bool(presence[p]) for p in analysed_ids)
    environment_rows = [
        {
            "category": "At least one retained defence/anti-defence system",
            "count": positive,
            "percent": f"{100 * positive / total:.2f}" if total else "0.00",
        },
        {
            "category": "No retained system detected",
            "count": total - positive,
            "percent": (
                f"{100 * (total - positive) / total:.2f}" if total else "0.00"
            ),
        },
    ]
    write_tsv(
        args.outdir / "query_environment_counts.tsv",
        environment_rows, ["category", "count", "percent"],
    )
    atomic_write(
        args.outdir / "defence_environment_percentage.svg",
        environment_svg(positive, total, target, colour),
    )
    atomic_write(
        args.outdir / "association_breadth.svg",
        breadth_svg(summary, analysed_loci, target),
    )
    return {
        "neighborhoods_analysed": total,
        "neighborhoods_excluded_from_denominator": len(contexts) - total,
        "neighborhoods_with_detected_systems": positive,
        "neighborhoods_without_detected_systems": total - positive,
        "raw_system_calls": len(raw_calls),
        "retained_system_calls": len(retained),
        "excluded_system_calls": len(excluded),
        "consensus_loci": len(analysed_loci),
        "counting_unit": "unique analysed representative protein neighbourhood",
        "excluded_systems": EXCLUDED_SYSTEMS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hit_table", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--email", required=True)
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    parser.add_argument("--protein-column")
    parser.add_argument("--flank-bp", type=int, default=30000)
    parser.add_argument(
        "--defensefinder-env", default="defensefinder",
        help="Conda environment name or environment root path",
    )
    parser.add_argument(
        "--padloc-env", default="padloc",
        help="Conda environment name or environment root path",
    )
    parser.add_argument(
        "--conda",
        default=shutil.which("conda") or shutil.which("mamba") or "conda",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--max-proteins", type=int)
    parser.add_argument("--delay", type=float)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-tool-checks", action="store_true")
    parser.add_argument("--skip-defensefinder", action="store_true")
    parser.add_argument("--skip-padloc", action="store_true")
    parser.add_argument(
        "--replot-only", action="store_true",
        help="Require cached contexts/scans; never download or run scanners",
    )
    parser.add_argument(
        "--retry-failed", action="store_true",
        help="Retry previously failed IPG/context retrievals",
    )
    parser.add_argument("--target-label", default="PAPI28")
    parser.add_argument("--heatmap-color", default="#264979")
    parser.add_argument("--heatmap-top-systems", type=int, default=25)
    parser.add_argument("--heatmap-max-neighborhoods", type=int, default=250)
    args = parser.parse_args()

    args.hit_table = args.hit_table.expanduser().resolve()
    args.outdir = args.outdir.expanduser().resolve()
    if not args.hit_table.is_file():
        parser.error(f"Input not found: {args.hit_table}")
    if args.flank_bp < 0 or min(args.threads, args.jobs) < 1:
        parser.error("Invalid flank size, thread count or job count")
    if args.max_proteins is not None and args.max_proteins < 1:
        parser.error("--max-proteins must be positive")
    if min(args.heatmap_top_systems, args.heatmap_max_neighborhoods) < 1:
        parser.error("Heatmap limits must be positive")
    if args.replot_only and (args.retry_failed or args.dry_run):
        parser.error("--replot-only cannot be combined with --retry-failed or --dry-run")
    if args.delay is not None and args.delay < 0:
        parser.error("--delay must be non-negative")

    try:
        hex_colour(args.heatmap_color)
        proteins, input_rows, protein_column = read_proteins(
            args.hit_table, args.protein_column
        )
        if args.max_proteins:
            proteins = proteins[:args.max_proteins]
        args.outdir.mkdir(parents=True, exist_ok=True)
        delay = (
            args.delay if args.delay is not None
            else 0.11 if args.api_key else 0.34
        )
        client = NcbiClient(
            args.email, args.api_key, delay,
            offline=args.offline or args.replot_only,
        )
        contexts, representatives, failures, reused_contexts = prepare_contexts(
            args, proteins, client
        )
        if not contexts:
            raise RuntimeError("No analysable cached or newly extracted contexts")
        raw_calls, scan_failures, statuses, reused_scans, executed_scans = scan_contexts(
            args, contexts
        )
        failures.extend(scan_failures)
        print("Rebuilding harmonised loci, frequencies and plots...", flush=True)
        report = rebuild_reports(args, contexts, raw_calls, statuses)
        write_tsv(args.outdir / "failures.tsv", failures, [
            "stage", "protein_accession", "error"
        ])
        summary = {
            "input": str(args.hit_table),
            "input_sha256": file_hash(args.hit_table),
            "input_rows": input_rows,
            "protein_column": protein_column,
            "unique_proteins_selected": len(proteins),
            "representatives_selected": len(representatives),
            "neighborhoods_extracted": len(contexts),
            "flank_bp_each_side": args.flank_bp,
            "target_label": args.target_label,
            "defensefinder_environment": args.defensefinder_env,
            "padloc_environment": args.padloc_env,
            "contexts_reused": reused_contexts,
            "scanner_runs_reused": reused_scans,
            "scanner_runs_executed": executed_scans,
            "replot_only": args.replot_only,
            "dry_run": args.dry_run,
            "failures": len(failures),
            **report,
        }
        atomic_write(
            args.outdir / "summary.json", json.dumps(summary, indent=2) + "\n"
        )
        print(
            f"Done: {report['neighborhoods_with_detected_systems']} positive / "
            f"{report['neighborhoods_analysed']} analysed neighbourhoods."
        )
        print(f"Results: {args.outdir}")
        return 0
    except (ValueError, RuntimeError, OSError, ET.ParseError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
