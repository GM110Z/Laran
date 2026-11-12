1. **plot-pfam-architectures.py**

Plots Pfam architectures, showing one representative protein and the count for each architecture.
Chooses the protein closest to each group’s median length and exports PDF, SVG, PNG and supporting tables.

## Options

| Option | Input / example | Description |
| :--- | :--- | :--- |
| `--run` | `results_20261110_195331` |PFAM results folder. |
| `--out` | `papi28_architectures_8colours` | Output folder|
| `--title` | `"PAPI28 domain architectures"` | Figure title in quotes. |
| `--palette` | `'#F25532,#4D5793,#758BAE'` | Comma-separated colours in quotes. Supply enough colours for every distinct Pfam type to avoid patterns. |
| `--no-tm` | No value | Omit transmembrane-helix overlays. |
| `--representative` | `PROTEIN_ACCESSION` | Force this protein as its architecture’s representative. Repeat for different architectures. |
| `--domain-colour` | `'FunZ_C=#F25532'` | Set a specific domain’s colour. Repeat for other domains. |
| `--show-accessions` | No value | Print representative accessions beside the counts. |

2. **flankguard.py**

Identifies defence and anti-defence systems around representative protein hits using DefenseFinder, AntiDefenseFinder and PADLOC. Reuses completed extractions and scans, harmonizes system names, merges calls sharing genes and rebuilds tables and SVG plots. Excludes VSPR and DMS-other from defence summaries while preserving raw calls.

## Options

| Option | Input / example | Description |
| :--- | :--- | :--- |
| `hit_table` | `PAPI28.csv` | Required positional input. BLAST table containing protein accessions. |
| `outdir` | `flankguard_results` | Required positional output folder. Reuse the same folder to resume completed work. |
| `--email` | `"name@example.com"` | Required contact email for NCBI requests. |
| `--api-key` | `"YOUR_NCBI_API_KEY"` | Optional NCBI API key. Also reads the `NCBI_API_KEY` environment variable. |
| `--protein-column` | `sseqid` | Header of the protein-accession column. Otherwise detected automatically, falling back to column 2. |
| `--flank-bp` | `30000` | Bases extracted on each side of the target gene. Default: 30,000. Use a new output folder when changing this value. |
| `--defensefinder-env` | `/path/to/defensefinder` | DefenseFinder Conda environment name or full environment root path. Default: `defensefinder`. |
| `--padloc-env` | `/path/to/padloc` | PADLOC Conda environment name or full environment root path. Default: `padloc`. |
| `--conda` | `/path/to/conda` | Conda or Mamba executable. Detected automatically where possible. |
| `--threads` | `4` | CPUs allocated to each scanner. Default: 1. |
| `--jobs` | `2` | Neighbourhoods processed concurrently. Default: 1. |
| `--max-proteins` | `3` | Process only the first N unique proteins for a trial run. |
| `--replot-only` | No value | Rebuild merged calls, tables and figures from completed cached results without downloads or scanner execution. |
| `--retry-failed` | No value | Retry previously failed protein mappings or neighbourhood retrievals. Cannot be combined with `--replot-only`. |
| `--offline` | No value | Use cached NCBI responses only. Missing scanner runs may still execute locally. |
| `--dry-run` | No value | Prepare inputs and print missing scanner commands without executing them. May download uncached inputs. |
| `--skip-tool-checks` | No value | Skip scanner dependency checks. |
| `--skip-defensefinder` | No value | Omit DefenseFinder and AntiDefenseFinder, including their cached calls. |
| `--skip-padloc` | No value | Omit PADLOC, including its cached calls. |
| `--target-label` | `"PAPI28"` | Target name displayed in figure titles. Default: `PAPI28`. |
| `--heatmap-color` | `'#264979'` | Six-digit hexadecimal plot colour. Default: `#264979`. |
| `--heatmap-top-systems` | `25` | Maximum systems displayed in frequency plots and neighbourhood heatmaps. Default: 25. |
| `--heatmap-max-neighborhoods` | `250` | Maximum rows displayed in neighbourhood heatmaps. Default: 250. |
| `--delay` | `0.34` | Seconds between NCBI requests. Default: 0.34, or 0.11 with an API key. |

3. **etnav3.py**

Retrieves all available RefSeq genomic mappings and source metadata for proteins in a BLAST table. Summarizes genus, clinical status and isolation source across unique nucleotide records. Exports TSV tables, an SVG figure, an HTML report and a JSON summary. Clinical and source classifications use metadata keywords and support curated corrections.

## Options

| Option | Input / example | Description |
| :--- | :--- | :--- |
| `hitlist` | `PAPI28.csv` | Required positional input. BLAST table with protein accessions in column 2. |
| `outdir` | `etna_results` | Required positional output folder. Also stores cached NCBI responses for reuse. |
| `--email` | `"name@example.com"` | Required contact email for NCBI requests. |
| `--api-key` | `"YOUR_NCBI_API_KEY"` | Optional NCBI API key. Also reads the `NCBI_API_KEY` environment variable. |
| `--header` | `auto` | Whether the input contains a header: `auto`, `yes` or `no`. Default: `auto`. |
| `--top-genera` | `20` | Number of genera displayed in the summary figure. Default: 20. |
| `--max-proteins` | `3` | Process only the first N unique proteins for a trial run. |
| `--max-nuccore` | `10` | Retrieve metadata for only the first N unique nucleotide records. For trial runs; produces an incomplete metadata summary. |
| `--offline` | No value | Read cached NCBI responses only, without network requests. |
| `--overrides` | `curated_sources.tsv` | Optional TSV containing `nuccore_accession` and corrections to `clinical_status`, `source_category`, `clinical_rule` or `source_rule`. Also accepts a `note` column. |
| `--delay` | `0.34` | Seconds between NCBI requests. Default: 0.34, or 0.11 with an API key. |
