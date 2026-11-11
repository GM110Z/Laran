1. **plot-pfam-architectures.py**

Plots Pfam architectures, showing one representative protein and the count for each architecture.
Chooses the protein closest to each group’s median length and exports PDF, SVG, PNG and supporting tables.

## Options

| Option | Input / example | Description |
| :--- | :--- | :--- |
| `--run` | `results_20261110_195331` | Existing FlaGs3 results folder. |
| `--out` | `papi28_architectures_8colours` | Output folder, outside the FlaGs3 results folder. |
| `--title` | `"PAPI28 domain architectures"` | Figure title in quotes. |
| `--palette` | `'#F25532,#4D5793,#758BAE'` | Comma-separated colours in quotes. Supply enough colours for every distinct Pfam type to avoid patterns. |
| `--no-tm` | No value | Omit transmembrane-helix overlays. |
| `--representative` | `PROTEIN_ACCESSION` | Force this protein as its architecture’s representative. Repeat for different architectures. |
| `--domain-colour` | `'FunZ_C=#F25532'` | Set a specific domain’s colour. Repeat for other domains. |
| `--show-accessions` | No value | Print representative accessions beside the counts. |
