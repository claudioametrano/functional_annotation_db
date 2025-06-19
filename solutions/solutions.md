#### TASK1
- No, the reference structure shows the entire protein, while the ~200 aa are only one subunit
- 8 identical subunits
- try copy/paste 8 time the monomeric subunit and see if Alphafold manage to predict something similar to the original pore (the 8 subunits polymerize due to H bonds)
#### TASK2
- Checks dependancies
- tRNA prediction (aragorn)
- rRNA prediction (Barrnap)
- CRISPR repeats
- Predicting CDS (Prodigal)
- blastp for initial match wit IS (insertion sequence: bacterial mobile element) ans AMR (anti microbial resistance) databases
- blastp for initial match with Swissprot
- Hmmer3 for the protein which did not received annotation (the majority) with blastp
- Creation of files compatible for NCBI submission (tbl2asn and sed)
- File output


#### TASK3

```bash
$ grep -E "tRNA-[A-Z][a-z][a-z]\(" results/Synechocystis_GCF_000009725.1/Synechocystis.tsv
```
add | wc -l to count them

Total codon  = 4^3 = 64
Sense codon = 64 -3 stop codons = 61
The **third base of the codon** (often called the **"wobble position"**) can pair with multiple bases in the tRNA anticodon, so less that 61 tRNA are often enough!

#### TASK4
A.

```bash
$ grep -c "GO:0009399" results/eggnog_annotation/*.annotations
```
results/eggnog_annotation/nostoc.emapper.annotations:31

results/eggnog_annotation/prochlorococcus.emapper.annotations:2

results/eggnog_annotation/synechocystis.emapper.annotations:3

B. *Prochlorococcus* lacks Phycobilisomes, it has indeed other adaptations to low light (not this particular strain though, which is adapted to high light) 

```bash
$ grep -c   "Phycobilisome" results/eggnog_annotation/*.annotations
```
results/eggnog_annotation/nostoc.emapper.annotations:14

results/eggnog_annotation/prochlorococcus.emapper.annotations:1

results/eggnog_annotation/synechocystis.emapper.annotations:13

```bash
$ grep -c   "GO:0030089" results/eggnog_annotation/*.annotations
```
results/eggnog_annotation/nostoc.emapper.annotations:25

results/eggnog_annotation/prochlorococcus.emapper.annotations:2

results/eggnog_annotation/synechocystis.emapper.annotations:17

