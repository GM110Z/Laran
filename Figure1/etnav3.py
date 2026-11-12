#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
from http.client import IncompleteRead
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UNKNOWN = "Unknown"
UNKNOWN_VALUES = {
    "", "missing", "unknown", "not applicable", "not_applicable", "n/a",
    "na", "none", "not provided", "not collected", "not determined", "nd",
    "null", "nan", "-", "--",
}
BLAST13 = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
]
PROTEIN_RE = re.compile(r"(?:[A-Z]{1,4}_)?[A-Z]{0,3}\d{5,}(?:\.\d+)?", re.I)


def clean(value) -> str:
    if value is None:
        return UNKNOWN
    text = str(value).strip()
    return UNKNOWN if text.lower().rstrip(".") in UNKNOWN_VALUES else text


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1].lower()


def normalise_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def pick(mapping: dict[str, str], *names: str) -> str:
    norm = {normalise_key(k): clean(v) for k, v in mapping.items()}
    for name in names:
        value = norm.get(normalise_key(name), UNKNOWN)
        if value != UNKNOWN:
            return value
    return UNKNOWN


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(text, encoding="utf-8")
    temp.replace(path)


class NcbiClient:
    def __init__(self, email: str, api_key: str | None, delay: float, retries: int = 8):
        self.email = email
        self.api_key = api_key
        self.delay = delay
        self.retries = retries
        self.last_request = 0.0

    def get(self, endpoint: str, params: dict[str, object]) -> str:
        query = {**params, "email": self.email, "tool": "blast_ipg_metadata"}
        if self.api_key:
            query["api_key"] = self.api_key
        url = BASE + endpoint + "?" + urlencode(query)
        for attempt in range(self.retries):
            wait = self.delay - (time.monotonic() - self.last_request)
            if wait > 0:
                time.sleep(wait)
            request = Request(url, headers={"User-Agent": f"blast-ipg-metadata/1.0 ({self.email})"})
            try:
                with urlopen(request, timeout=90) as response:
                    text = response.read().decode("utf-8")
                self.last_request = time.monotonic()
                return text
            except HTTPError as exc:
                self.last_request = time.monotonic()
                if exc.code not in {429, 500, 502, 503, 504}:
                    raise RuntimeError(f"NCBI HTTP {exc.code} for {endpoint}") from exc
                reason = f"HTTP {exc.code}"
            except (URLError, TimeoutError, IncompleteRead) as exc:
                self.last_request = time.monotonic()
                reason = str(exc)
            retry_wait = max(1.0, self.delay) * (2**attempt)
            print(f"  NCBI request failed ({reason}); retrying in {retry_wait:.1f}s", file=sys.stderr)
            time.sleep(retry_wait)
        raise RuntimeError(f"NCBI request failed after {self.retries} attempts: {endpoint}")


def sniff_table(path: Path, header_mode: str) -> tuple[list[str], list[dict[str, str]]]:
    sample = path.read_text(encoding="utf-8-sig", errors="replace")[:65536]
    if not sample.strip():
        raise ValueError(f"Input is empty: {path}")
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        delimiter = dialect.delimiter
    except csv.Error:
        delimiter = "\t" if "\t" in sample.partition("\n")[0] else ","
    raw = list(csv.reader(sample.splitlines(), delimiter=delimiter))
    if not raw or len(raw[0]) < 2:
        raise ValueError("The hit list must contain at least two columns")
    first = [x.strip().lower() for x in raw[0]]
    looks_header = (
        first[1] in {"sseqid", "subject", "subject_id", "accession", "protein_accession"}
        or first[0] in {"qseqid", "query", "query_id"}
    )
    has_header = looks_header if header_mode == "auto" else header_mode == "yes"
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
        rows = list(csv.reader(handle, delimiter=delimiter))
    width = max(len(row) for row in rows)
    if has_header:
        headers = [cell.strip() or f"column_{i + 1}" for i, cell in enumerate(rows[0])]
        headers.extend(f"column_{i + 1}" for i in range(len(headers), width))
        data_rows = rows[1:]
    else:
        headers = BLAST13[:width]
        headers.extend(f"extra_{i + 1}" for i in range(len(headers), width))
        data_rows = rows
    records = []
    for line_number, row in enumerate(data_rows, start=2 if has_header else 1):
        if not row or not any(x.strip() for x in row):
            continue
        padded = row + [""] * (width - len(row))
        record = {headers[i]: padded[i].strip() for i in range(width)}
        record["input_line"] = str(line_number)
        records.append(record)
    return headers, records


def extract_protein(value: str) -> str:
    value = clean(value)
    if value == UNKNOWN:
        return UNKNOWN
    candidates = PROTEIN_RE.findall(value)
    return candidates[-1].upper() if candidates else value.split()[0]


def attr(element: ET.Element, *names: str) -> str:
    attrs = {normalise_key(k): v for k, v in element.attrib.items()}
    for name in names:
        value = clean(attrs.get(normalise_key(name)))
        if value != UNKNOWN:
            return value
    return UNKNOWN


def parse_ipg_xml(xml_text: str, requested_protein: str) -> list[dict[str, str]]:
    root = ET.fromstring(xml_text)
    mappings: list[dict[str, str]] = []
    current_protein: dict[str, str] | None = None

    def walk(element: ET.Element, protein: dict[str, str] | None) -> None:
        nonlocal mappings
        tag = local_name(element.tag)
        if tag in {"protein", "proteinrecord"}:
            protein = {
                "ipg_source": attr(element, "source", "database", "db"),
                "ipg_protein_accession": attr(element, "accver", "accession", "protein_acc"),
                "protein_name": attr(element, "name", "product", "protein_name"),
            }
        if tag in {"cds", "cdsrecord", "nucleotide"}:
            source = (protein or {}).get("ipg_source", UNKNOWN)
            source = attr(element, "source", "database", "db") if source == UNKNOWN else source
            nucleotide = attr(
                element, "accver", "nucleotide_acc", "nucleotide_accession",
                "nucaccver", "nuc_accession", "accession",
            )
            if source.lower() == "refseq" and nucleotide != UNKNOWN:
                mappings.append({
                    "blast_protein_accession": requested_protein,
                    "ipg_source": source,
                    "ipg_protein_accession": (protein or {}).get("ipg_protein_accession", requested_protein),
                    "protein_name": (protein or {}).get("protein_name", UNKNOWN),
                    "nuccore_accession": nucleotide,
                    "cds_start": attr(element, "start", "from"),
                    "cds_stop": attr(element, "stop", "to", "end"),
                    "strand": attr(element, "strand"),
                    "organism": attr(element, "organism", "org"),
                    "strain": attr(element, "strain"),
                    "assembly": attr(element, "assembly", "assembly_acc", "assembly_accession"),
                })
        for child in element:
            walk(child, protein)

    walk(root, current_protein)
    if not mappings:
        for element in root.iter():
            nucleotide = attr(element, "nucleotide_acc", "nucleotide_accession", "nucaccver")
            source = attr(element, "source", "database", "db")
            if source.lower() == "refseq" and nucleotide != UNKNOWN:
                mappings.append({
                    "blast_protein_accession": requested_protein,
                    "ipg_source": source,
                    "ipg_protein_accession": attr(element, "protein_acc", "protein_accession"),
                    "protein_name": attr(element, "name", "product"),
                    "nuccore_accession": nucleotide,
                    "cds_start": attr(element, "start", "from"),
                    "cds_stop": attr(element, "stop", "to", "end"),
                    "strand": attr(element, "strand"),
                    "organism": attr(element, "organism", "org"),
                    "strain": attr(element, "strain"),
                    "assembly": attr(element, "assembly", "assembly_acc", "assembly_accession"),
                })
    unique = {}
    for row in mappings:
        key = (row["blast_protein_accession"], row["nuccore_accession"], row["cds_start"], row["cds_stop"])
        unique[key] = row
    return list(unique.values())


def parse_nuccore_xml(xml_text: str) -> dict[str, str]:
    root = ET.fromstring(xml_text)
    record = next((x for x in root.iter() if local_name(x.tag) in {"gbseq", "insdseq"}), root)
    prefix = "gbseq" if local_name(record.tag) == "gbseq" else "insdseq"

    def child_text(suffix: str) -> str:
        target = f"{prefix}_{suffix}".lower()
        for element in record.iter():
            if local_name(element.tag) == target:
                return clean(element.text)
        return UNKNOWN

    qualifiers: dict[str, list[str]] = defaultdict(list)
    in_source = False
    for feature in record.iter():
        if local_name(feature.tag) not in {f"{prefix}_feature", "gbfeature", "insdfeature"}:
            continue
        key = UNKNOWN
        for element in feature:
            if local_name(element.tag).endswith("_key"):
                key = clean(element.text)
                break
        if key.lower() != "source":
            continue
        in_source = True
        name = None
        for element in feature.iter():
            tag = local_name(element.tag)
            if tag.endswith("_name"):
                name = clean(element.text)
            elif tag.endswith("_value") and name and name != UNKNOWN:
                value = clean(element.text)
                if value != UNKNOWN:
                    qualifiers[name.lower()].append(value)
                name = None
    flat = {key: "; ".join(dict.fromkeys(values)) for key, values in qualifiers.items()}
    biosample = UNKNOWN
    for value in flat.get("db_xref", "").split(";"):
        if value.strip().lower().startswith("biosample:"):
            biosample = value.split(":", 1)[1].strip()
            break
    # RefSeq records commonly expose BioSample as a record-level GBXref rather
    # than as a /db_xref qualifier on the source feature.
    if biosample == UNKNOWN:
        for xref in record.iter():
            if local_name(xref.tag) not in {"gbxref", "insdxref"}:
                continue
            database = UNKNOWN
            identifier = UNKNOWN
            for child in xref:
                tag = local_name(child.tag)
                if tag.endswith("_dbname"):
                    database = clean(child.text)
                elif tag.endswith("_id"):
                    identifier = clean(child.text)
            if database.lower() == "biosample" and identifier != UNKNOWN:
                biosample = identifier
                break
    return {
        "nuccore_accession": child_text("accession-version"),
        "nuccore_definition": child_text("definition"),
        "organism": child_text("organism"),
        "taxonomy": child_text("taxonomy"),
        "biosample_accession": biosample,
        "isolation_source_nuccore": pick(flat, "isolation_source", "sample source", "source"),
        "host_nuccore": pick(flat, "host", "host scientific name"),
        "disease_nuccore": pick(flat, "host_disease", "disease"),
        "geo_location_nuccore": pick(flat, "geo_loc_name", "country"),
        "collection_date_nuccore": pick(flat, "collection_date"),
        "strain_nuccore": pick(flat, "strain", "isolate"),
        "all_nuccore_source_qualifiers": "; ".join(f"{k}={v}" for k, v in sorted(flat.items())),
        "nuccore_source_feature_found": "yes" if in_source else "no",
    }


def parse_elink_biosample_uid(xml_text: str) -> str:
    root = ET.fromstring(xml_text)
    for linksetdb in root.iter():
        if local_name(linksetdb.tag) != "linksetdb":
            continue
        db_to = ""
        for child in linksetdb:
            if local_name(child.tag) == "dbto":
                db_to = clean(child.text).lower()
        if db_to and db_to != "biosample":
            continue
        for element in linksetdb.iter():
            if local_name(element.tag) == "id" and clean(element.text) != UNKNOWN:
                return clean(element.text)
    return UNKNOWN


def parse_esearch_uid(xml_text: str) -> str:
    root = ET.fromstring(xml_text)
    for id_list in root.iter():
        if local_name(id_list.tag) == "idlist":
            for element in id_list:
                if local_name(element.tag) == "id":
                    return clean(element.text)
    return UNKNOWN


def parse_biosample_xml(xml_text: str) -> dict[str, str]:
    root = ET.fromstring(xml_text)
    sample = next((x for x in root.iter() if local_name(x.tag) == "biosample"), None)
    if sample is None:
        raise ValueError("No BioSample element in NCBI response")
    attributes: dict[str, str] = {}
    identifiers: dict[str, str] = {}
    for element in sample.iter():
        tag = local_name(element.tag)
        if tag == "attribute":
            value = clean(element.text)
            for key in ("attribute_name", "harmonized_name", "display_name", "name"):
                label = element.attrib.get(key)
                if label and value != UNKNOWN:
                    attributes[label.lower().strip()] = value
        elif tag == "id":
            value = clean(element.text)
            for key in ("db", "db_label", "namespace"):
                label = element.attrib.get(key)
                if label and value != UNKNOWN:
                    identifiers[label.lower().strip()] = value
    title = UNKNOWN
    organism = UNKNOWN
    for element in sample.iter():
        tag = local_name(element.tag)
        if tag == "title" and title == UNKNOWN:
            title = clean(element.text)
        elif tag == "organism" and organism == UNKNOWN:
            organism = attr(element, "taxonomy_name", "organism_name")
    return {
        "biosample_accession": clean(sample.attrib.get("accession")),
        "biosample_uid": clean(sample.attrib.get("id")),
        "biosample_title": title,
        "biosample_organism": organism,
        "geo_location": pick(attributes, "geo_loc_name", "geographic location", "country"),
        "isolation_source": pick(attributes, "isolation_source", "sample source", "source"),
        "host": pick(attributes, "host", "host scientific name"),
        "disease": pick(attributes, "disease", "host_disease", "host health state", "disease state"),
        "collection_date": pick(attributes, "collection_date", "collection date"),
        "strain": pick(attributes, "strain", "isolate"),
        "sample_name": pick(attributes, "sample name", "sample_name"),
        "submitter_id": clean(identifiers.get("submitter id") or identifiers.get("external id")),
        "all_biosample_attributes": "; ".join(f"{k}={v}" for k, v in sorted(attributes.items())),
    }


def genus_from_organism(organism: str) -> str:
    organism = clean(organism).strip("[]")
    if organism == UNKNOWN or not organism:
        return UNKNOWN
    words = organism.split()
    if words[0].lower() in {"uncultured", "unclassified", "environmental", "synthetic"}:
        return UNKNOWN
    if words[0].lower() == "candidatus" and len(words) > 1:
        return "Candidatus " + words[1]
    return words[0]


CLINICAL_PATTERNS = [
    (r"\bclinical\b", "explicit clinical label"),
    (r"\bpatient\b|\bhospital\b|\bnosocomial\b", "patient/hospital wording"),
    (r"\bblood(?:stream)?\b|\bsepsis\b|\bbacteremia\b", "bloodstream specimen/disease"),
    (r"\bburn\b|\bwound\b|\beschar\b|\babscess\b", "burn/wound specimen"),
    (r"\bcystic fibrosis\b", "cystic fibrosis"),
    (r"\bsputum\b|\bbronchoalveolar\b|\brespiratory tract\b|\btracheal aspirate\b", "respiratory specimen"),
    (r"\burine\b|\bcerebrospinal fluid\b|\bcatheter\b|\bbiopsy\b", "clinical specimen"),
]
NONCLINICAL_PATTERNS = [
    (r"\benvironmental\b|\benvironment\b", "explicit environmental label"),
    (r"\bsoil\b|\bsediment\b|\bcompost\b", "soil/sediment source"),
    (r"\bsea water\b|\bseawater\b|\bfreshwater\b|\briver\b|\blake\b", "natural-water source"),
    (r"\bplant\b|\brhizosphere\b|\bleaf\b|\broot\b", "plant-associated source"),
    (r"\bfood\b|\bcheese\b|\bmilk\b|\bmeat\b", "food source"),
    (r"\bhealthy (?:human|person|volunteer|subject)\b", "healthy donor"),
]


def classify_source(row: dict[str, str]) -> tuple[str, str, str, str]:
    fields = [
        row.get("isolation_source", UNKNOWN), row.get("isolation_source_nuccore", UNKNOWN),
        row.get("host", UNKNOWN), row.get("host_nuccore", UNKNOWN),
        row.get("disease", UNKNOWN), row.get("disease_nuccore", UNKNOWN),
        row.get("biosample_title", UNKNOWN), row.get("sample_name", UNKNOWN),
        row.get("nuccore_definition", UNKNOWN), row.get("all_biosample_attributes", ""),
    ]
    evidence = " | ".join(x for x in fields if clean(x) != UNKNOWN)
    lowered = evidence.lower()
    clinical_matches = [(pattern, label) for pattern, label in CLINICAL_PATTERNS if re.search(pattern, lowered)]
    nonclinical_matches = [(pattern, label) for pattern, label in NONCLINICAL_PATTERNS if re.search(pattern, lowered)]
    if clinical_matches:
        clinical, clinical_rule = "Clinical", clinical_matches[0][1]
    elif nonclinical_matches:
        clinical, clinical_rule = "Non-clinical", nonclinical_matches[0][1]
    else:
        clinical, clinical_rule = UNKNOWN, "no decisive keyword"

    if re.search(r"\bburn\b|\bwound\b|\beschar\b", lowered):
        category, source_rule = "Burn/wound", "burn, wound, or eschar keyword"
    elif re.search(r"\bcystic fibrosis\b", lowered) or (
        re.search(r"\bcf\b", lowered) and re.search(r"\blung\b|\bsputum\b|\brespiratory\b", lowered)
    ):
        category, source_rule = "CF respiratory", "cystic fibrosis respiratory wording"
    elif re.search(r"\blung\b|\bsputum\b|\brespiratory\b|\bbronchoalveolar\b|\btracheal\b", lowered):
        category, source_rule = "Other respiratory", "respiratory wording without explicit CF"
    elif clinical == "Clinical":
        category, source_rule = "Other clinical", clinical_rule
    elif clinical == "Non-clinical":
        category, source_rule = "Non-clinical", clinical_rule
    else:
        category, source_rule = UNKNOWN, "no decisive source keyword"
    return clinical, clinical_rule, category, source_rule


def read_overrides(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "nuccore_accession" not in reader.fieldnames:
            raise ValueError("Override TSV requires a nuccore_accession column")
        return {row["nuccore_accession"].strip(): row for row in reader if row.get("nuccore_accession", "").strip()}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def count_rows(values: Iterable[str], total: int) -> list[dict[str, object]]:
    counts = Counter(clean(value) for value in values)
    return [
        {"category": category, "count": count, "percent": f"{100 * count / total:.2f}" if total else "0.00"}
        for category, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]


PALETTE = ["#19647E", "#28AFB0", "#F4D35E", "#EE964B", "#C44536", "#6C5B7B", "#5B8E7D"]


def svg_bar_panel(title: str, rows: list[dict[str, object]], x: int, y: int, width: int, max_rows: int) -> tuple[str, int]:
    shown = rows[:max_rows]
    height = 52 + 31 * max(1, len(shown))
    maximum = max([int(row["count"]) for row in shown] or [1])
    label_width = min(220, max([len(str(row["category"])) * 7 for row in shown] or [80]))
    bar_x = x + label_width + 20
    bar_width = width - label_width - 85
    parts = [f'<text x="{x}" y="{y + 20}" class="panel-title">{html.escape(title)}</text>']
    for index, row in enumerate(shown):
        yy = y + 43 + index * 31
        value = int(row["count"])
        length = max(2, round(bar_width * value / maximum))
        label = html.escape(str(row["category"]))
        parts.extend([
            f'<text x="{bar_x - 10}" y="{yy + 15}" text-anchor="end" class="label">{label}</text>',
            f'<rect x="{bar_x}" y="{yy}" width="{length}" height="20" rx="4" fill="{PALETTE[index % len(PALETTE)]}"/>',
            f'<text x="{bar_x + length + 8}" y="{yy + 15}" class="value">{value}</text>',
        ])
    return "\n".join(parts), height


def make_svg(genus: list[dict[str, object]], clinical: list[dict[str, object]], sources: list[dict[str, object]], total: int, top_genera: int) -> str:
    width = 1080
    genus_h = 52 + 31 * min(top_genera, max(1, len(genus)))
    lower_h = max(52 + 31 * max(1, len(clinical)), 52 + 31 * max(1, len(sources)))
    height = 145 + genus_h + lower_h + 70
    panel1, _ = svg_bar_panel(f"Top genera — unique RefSeq nucleotide records (n={total})", genus, 55, 105, 970, top_genera)
    panel2, _ = svg_bar_panel("Clinical status", clinical, 55, 130 + genus_h, 455, 8)
    panel3, _ = svg_bar_panel("Isolation-source category", sources, 565, 130 + genus_h, 460, 8)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
  .title {{ font: 700 28px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; fill:#17324D }}
  .subtitle {{ font: 14px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; fill:#60758A }}
  .panel-title {{ font: 700 18px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; fill:#17324D }}
  .label {{ font: 13px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; fill:#334E68 }}
  .value {{ font: 700 13px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; fill:#334E68 }}
</style>
<rect width="100%" height="100%" rx="18" fill="#F7FAFC"/>
<text x="55" y="48" class="title">RefSeq IPG hit landscape</text>
<text x="55" y="75" class="subtitle">Counts are deduplicated by nucleotide accession; Unknown is retained rather than guessed.</text>
{panel1}
{panel2}
{panel3}
<text x="55" y="{height - 28}" class="subtitle">Clinical/source groups are keyword-assisted screening labels. Review the evidence columns before biological interpretation.</text>
</svg>'''


def html_table(title: str, rows: list[dict[str, object]]) -> str:
    body = "".join(
        f"<tr><td>{html.escape(str(r['category']))}</td><td>{r['count']}</td><td>{r['percent']}%</td></tr>" for r in rows
    )
    return f"<section><h2>{html.escape(title)}</h2><table><thead><tr><th>Category</th><th>Count</th><th>Percent</th></tr></thead><tbody>{body}</tbody></table></section>"


def make_html(total: int, proteins: int, mappings: int, genus, clinical, sources) -> str:
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>RefSeq IPG hit report</title><style>
body{{margin:0;background:#edf3f7;color:#17324d;font:15px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}}main{{max-width:1120px;margin:35px auto;padding:0 24px}}
.hero,section{{background:white;border-radius:16px;padding:24px;margin:18px 0;box-shadow:0 8px 26px #17324d12}}h1,h2{{margin-top:0}}.metrics{{display:flex;gap:16px;flex-wrap:wrap}}.metric{{background:#f1f7f8;padding:14px 20px;border-radius:12px}}.metric b{{font-size:24px;display:block}}
img{{width:100%;height:auto;border-radius:16px}}table{{border-collapse:collapse;width:100%}}th,td{{text-align:left;padding:9px 12px;border-bottom:1px solid #dce6ed}}th{{color:#60758a}}.note{{color:#60758a}}
</style></head><body><main><div class="hero"><h1>RefSeq IPG hit report</h1><div class="metrics"><div class="metric"><b>{proteins}</b>unique BLAST proteins</div><div class="metric"><b>{mappings}</b>RefSeq IPG mappings</div><div class="metric"><b>{total}</b>unique nucleotide records</div></div><p class="note">Counts below use unique nucleotide accessions. Automated clinical/source labels are screening aids and retain uncertain records as Unknown.</p></div>
<img src="hit_summary.svg" alt="Bar plots summarising genus, clinical status, and isolation source">
{html_table('Genus counts', genus)}{html_table('Clinical-status counts', clinical)}{html_table('Isolation-source counts', sources)}
</main></body></html>'''


def cache_fetch(client: NcbiClient, path: Path, endpoint: str, params: dict[str, object], offline: bool) -> str:
    if path.exists() and path.stat().st_size:
        return path.read_text(encoding="utf-8")
    if offline:
        raise FileNotFoundError(f"Missing cache file in --offline mode: {path}")
    text = client.get(endpoint, params)
    atomic_write(path, text)
    return text


def fetch_metadata(client: NcbiClient, accession: str, cache_dir: Path, offline: bool) -> dict[str, str]:
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", accession)
    xml = cache_fetch(client, cache_dir / "nuccore" / f"{safe}.xml", "efetch.fcgi", {
        # One base is sufficient: the record metadata and source feature are
        # retained, while large CDS feature tables and sequences are avoided.
        "db": "nuccore", "id": accession, "rettype": "gb", "retmode": "xml",
        "seq_start": 1, "seq_stop": 1,
    }, offline)
    row = parse_nuccore_xml(xml)
    if row["nuccore_accession"] == UNKNOWN:
        row["nuccore_accession"] = accession
    biosample_acc = row["biosample_accession"]
    biosample_uid = UNKNOWN
    if biosample_acc != UNKNOWN:
        search_xml = cache_fetch(client, cache_dir / "biosample_search" / f"{safe}.xml", "esearch.fcgi", {
            "db": "biosample", "term": f"{biosample_acc}[Accession]", "retmax": 1, "retmode": "xml",
        }, offline)
        biosample_uid = parse_esearch_uid(search_xml)
    if biosample_uid == UNKNOWN:
        link_xml = cache_fetch(client, cache_dir / "biosample_link" / f"{safe}.xml", "elink.fcgi", {
            "dbfrom": "nuccore", "db": "biosample", "id": accession, "retmode": "xml",
        }, offline)
        biosample_uid = parse_elink_biosample_uid(link_xml)
    biosample = {}
    if biosample_uid != UNKNOWN:
        bs_xml = cache_fetch(client, cache_dir / "biosample" / f"{biosample_uid}.xml", "efetch.fcgi", {
            "db": "biosample", "id": biosample_uid, "retmode": "xml",
        }, offline)
        biosample = parse_biosample_xml(bs_xml)
    else:
        biosample = {
            "biosample_accession": biosample_acc, "biosample_uid": UNKNOWN,
            "biosample_title": UNKNOWN, "biosample_organism": UNKNOWN,
            "geo_location": UNKNOWN, "isolation_source": UNKNOWN, "host": UNKNOWN,
            "disease": UNKNOWN, "collection_date": UNKNOWN, "strain": UNKNOWN,
            "sample_name": UNKNOWN, "submitter_id": UNKNOWN, "all_biosample_attributes": "",
        }
    row.update(biosample)
    row["organism"] = clean(row.get("biosample_organism")) if clean(row.get("biosample_organism")) != UNKNOWN else row["organism"]
    row["genus"] = genus_from_organism(row["organism"])
    clinical, clinical_rule, category, source_rule = classify_source(row)
    row.update({
        "clinical_status": clinical, "clinical_rule": clinical_rule,
        "source_category": category, "source_rule": source_rule,
        "classification_evidence": " | ".join(clean(row.get(key)) for key in (
            "isolation_source", "isolation_source_nuccore", "host", "host_nuccore",
            "disease", "disease_nuccore", "biosample_title", "sample_name",
        ) if clean(row.get(key)) != UNKNOWN),
        "error": "",
    })
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hitlist", type=Path, help="BLAST table; protein accessions must be in column 2")
    parser.add_argument("outdir", type=Path, help="Output folder (also holds the resumable XML cache)")
    parser.add_argument("--email", required=True, help="Contact email required by NCBI E-utilities")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"), help="Optional NCBI API key (or NCBI_API_KEY)")
    parser.add_argument("--header", choices=("auto", "yes", "no"), default="auto", help="Whether the input has a header")
    parser.add_argument("--delay", type=float, default=None, help="Seconds between NCBI requests (default: 0.34, or 0.11 with API key)")
    parser.add_argument("--top-genera", type=int, default=20, help="Number of genera displayed in the plot")
    parser.add_argument("--max-proteins", type=int, help="Process only the first N unique proteins (useful for a test run)")
    parser.add_argument("--max-nuccore", type=int, help="Fetch metadata for only the first N unique nucleotide records (smoke tests only)")
    parser.add_argument("--offline", action="store_true", help="Read existing XML cache only; make no NCBI requests")
    parser.add_argument("--overrides", type=Path, help="Optional TSV with nuccore_accession and curated classification columns")
    args = parser.parse_args()
    if not args.hitlist.is_file():
        parser.error(f"Hit list not found: {args.hitlist}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.outdir / "cache"
    delay = args.delay if args.delay is not None else (0.11 if args.api_key else 0.34)
    client = NcbiClient(args.email, args.api_key, delay)

    headers, blast_rows = sniff_table(args.hitlist, args.header)
    subject_column = headers[1]
    for row in blast_rows:
        row["protein_accession"] = extract_protein(row.get(subject_column, ""))
    proteins = list(dict.fromkeys(row["protein_accession"] for row in blast_rows if row["protein_accession"] != UNKNOWN))
    if args.max_proteins:
        proteins = proteins[: args.max_proteins]
    selected = set(proteins)
    normalized_hits = [row for row in blast_rows if row["protein_accession"] in selected]
    write_tsv(args.outdir / "blast_hits.tsv", normalized_hits, ["input_line", "protein_accession", *headers])
    print(f"Input rows: {len(blast_rows)}; unique proteins selected: {len(proteins)}", flush=True)

    mappings: list[dict[str, str]] = []
    failures = []
    for index, protein in enumerate(proteins, 1):
        print(f"[IPG {index}/{len(proteins)}] {protein}", flush=True)
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", protein)
        try:
            xml = cache_fetch(client, cache_dir / "ipg" / f"{safe}.xml", "efetch.fcgi", {
                "db": "protein", "id": protein, "rettype": "ipg", "retmode": "xml",
            }, args.offline)
            found = parse_ipg_xml(xml, protein)
            mappings.extend(found)
            if not found:
                failures.append({"stage": "ipg", "accession": protein, "error": "No RefSeq CDS mappings in IPG report"})
        except Exception as exc:
            print(f"  FAILED: {exc}", file=sys.stderr)
            failures.append({"stage": "ipg", "accession": protein, "error": str(exc)})
    mapping_fields = [
        "blast_protein_accession", "ipg_source", "ipg_protein_accession", "protein_name",
        "nuccore_accession", "cds_start", "cds_stop", "strand", "organism", "strain", "assembly",
    ]
    write_tsv(args.outdir / "ipg_refseq_hits.tsv", mappings, mapping_fields)

    nucleotide_accessions = list(dict.fromkeys(row["nuccore_accession"] for row in mappings))
    all_nucleotide_count = len(nucleotide_accessions)
    if args.max_nuccore:
        nucleotide_accessions = nucleotide_accessions[: args.max_nuccore]
    metadata = []
    for index, accession in enumerate(nucleotide_accessions, 1):
        print(f"[Metadata {index}/{len(nucleotide_accessions)}] {accession}", flush=True)
        try:
            metadata.append(fetch_metadata(client, accession, cache_dir, args.offline))
        except Exception as exc:
            print(f"  FAILED: {exc}", file=sys.stderr)
            metadata.append({
                "nuccore_accession": accession, "organism": UNKNOWN, "genus": UNKNOWN,
                "clinical_status": UNKNOWN, "clinical_rule": "metadata retrieval failed",
                "source_category": UNKNOWN, "source_rule": "metadata retrieval failed",
                "classification_evidence": "", "error": str(exc),
            })
            failures.append({"stage": "metadata", "accession": accession, "error": str(exc)})

    overrides = read_overrides(args.overrides)
    for row in metadata:
        override = overrides.get(row["nuccore_accession"])
        if not override:
            continue
        for column in ("clinical_status", "source_category", "clinical_rule", "source_rule"):
            if clean(override.get(column)) != UNKNOWN:
                row[column] = override[column].strip()
        row["classification_override"] = "yes"
        row["override_note"] = clean(override.get("note"))
    metadata_fields = sorted({key for row in metadata for key in row})
    preferred = [
        "nuccore_accession", "organism", "genus", "assembly", "biosample_accession", "biosample_uid",
        "clinical_status", "clinical_rule", "source_category", "source_rule", "classification_evidence",
        "isolation_source", "host", "disease", "geo_location", "collection_date", "strain",
        "biosample_title", "sample_name", "nuccore_definition", "error",
    ]
    metadata_fields = [x for x in preferred if x in metadata_fields] + [x for x in metadata_fields if x not in preferred]
    write_tsv(args.outdir / "nuccore_metadata.tsv", metadata, metadata_fields)

    meta_by_accession = {row["nuccore_accession"]: row for row in metadata}
    enriched = []
    for mapping in mappings:
        row = dict(mapping)
        row.update({f"metadata_{k}" if k in row else k: v for k, v in meta_by_accession.get(mapping["nuccore_accession"], {}).items()})
        enriched.append(row)
    write_tsv(args.outdir / "refseq_hits_enriched.tsv", enriched)
    write_tsv(args.outdir / "failures.tsv", failures, ["stage", "accession", "error"])

    total = len(metadata)
    genus = count_rows((row.get("genus", UNKNOWN) for row in metadata), total)
    clinical = count_rows((row.get("clinical_status", UNKNOWN) for row in metadata), total)
    sources = count_rows((row.get("source_category", UNKNOWN) for row in metadata), total)
    write_tsv(args.outdir / "genus_counts.tsv", genus, ["category", "count", "percent"])
    write_tsv(args.outdir / "clinical_status_counts.tsv", clinical, ["category", "count", "percent"])
    write_tsv(args.outdir / "source_category_counts.tsv", sources, ["category", "count", "percent"])
    atomic_write(args.outdir / "hit_summary.svg", make_svg(genus, clinical, sources, total, args.top_genera))
    atomic_write(args.outdir / "report.html", make_html(total, len(proteins), len(mappings), genus, clinical, sources))
    summary = {
        "input": str(args.hitlist.resolve()), "input_rows": len(blast_rows), "unique_proteins": len(proteins),
        "refseq_ipg_mappings": len(mappings), "unique_nuccore_records_in_ipg": all_nucleotide_count,
        "unique_nuccore_records_with_metadata": total,
        "metadata_failures": sum(bool(row.get("error")) for row in metadata),
        "counting_unit": "unique nuccore_accession", "clinical_status_counts": dict(Counter(row.get("clinical_status", UNKNOWN) for row in metadata)),
        "source_category_counts": dict(Counter(row.get("source_category", UNKNOWN) for row in metadata)),
    }
    atomic_write(args.outdir / "summary.json", json.dumps(summary, indent=2) + "\n")
    print(f"Done. Report: {args.outdir / 'report.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
