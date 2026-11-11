1. **plot-pfam-architectures.py**

Plots Pfam architectures, showing one representative protein and the count for each architecture.
Chooses the protein closest to each group’s median length and exports PDF, SVG, PNG and supporting tables.


Parameter	What to enter
--run	folder with pfam table.
--out	Where to save figures and tables. 
--title	Figure title in quotes.
--palette	Quoted, comma-separated hex colours. Supply enough colours for all distinct Pfam types to avoid dots.
--no-tm	Omit TM helices; takes no value (if protein has none)
--representative	Optional accession to force as a group’s representative.
--domain-colour	Optional specific colour, e.g. 'FunZ_C=#F25532'.
--show-accessions	Print representative accessions beside counts.
