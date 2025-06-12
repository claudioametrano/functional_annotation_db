# Functional annotation databases

### Genomic workflow: a recap
![workflow](images/genomic_workflow.png)
from [Jung et al. 2020](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008325)modif.
In brief:
- Sample -> DNA extraction
- Library preparation
- Sequencing/base-calling
- Raw reads QC and filtering
- De novo assembly (or mapping on reference)
- Polishing and scaffolding
- Structural annotation -> non-coding and coding region: gene prediction (*ab initio* + evidence-based)
- **functional annotation**
- Visualization and submission
### What do we mean by genome annotation?
The term **"genome annotation"** is typically used to refer to two discrete bioinformatics tasks. The first of these is structural annotation, namely the prediction of the genomic intervals comprised of functional genome features: genes, transcripts associated with those genes, the exons comprising those transcripts, and for protein-coding transcripts their coding sequences (CDS) and untranslated regions (UTRs). The second  is functional annotation, namely assigning gene symbols (names) and putative functions to the structurally identified gene models, in order to associate biological information to genes and corresponding proteins.
We will focus here on functional annotation, displaying relevant database and software resources to assign amino-acid sequences to a putative function.

![func. annot.](images/functional_annot.png)
from [Del Angel et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5850084/)

Some of the method are not that different from the one you applied in day 1 and 2 (BLAST),
what makes the difference here is the completeness (in term of phylogenetic breadth) and accuracy of the annotations stored in databases, also in term of delimited orthogroup.

### Database and tools for functional annotation
Many tool and database with annotated genes/genes families are available to perform functional annotation. Also many of these databases have their own integrated tool to perform annotation (e.g. EggNOG and EggNOG mapper) and web-servers for quick search, but: 
- They do not necessarily "speak the same language": different gene identifiers, organization and hierarchy of the database and targeted feature (e.g. full coding sequence, domains)
- They do not cover necessarily the same set of genes
- Lack of annotation for a relevant fraction of genes outside core metabolism (especially in non-model organisms)
**However there is a good interchange  among databases, which also report the accession and the classificatio (e.g. GO terms) of other relevant databases for the same protein.**

### Relevant resources
| Database               | Number entries                                                                                                                                                                  | classification methos                                                            | tools/software                                                                | URL                                                              | Notes                                                                                                               |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- | ---------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **UniProtKB** 2025 _02 | 573 230 **Swiss-Prot** (reviewed); 252 188 522 **TrEMBL** proteins                                                                                                              |                                                                                  | No single pipeline, but UniRule, BLAST+, DIAMOND and REST                     | [https://www.uniprot.org](https://www.uniprot.org)               | Sequence records split into manually-curated and automatically annotated sets; cross-links to nearly every other DB |
| **InterPro** rel. 93   | ≈ 48 k protein signatures                                                                                                                                                       | Hidden-Markov profiles/domains mapped onto proteins; GO terms present            | **InterProScan 5** (CLI/Web) does the search & GO transfer                    | [https://www.ebi.ac.uk/interpro](https://www.ebi.ac.uk/interpro) | One-stop domain scan; covers > 81 % of UniProtKB                                                                    |
| **Pfam** v37.3         | 24 424 families                                                                                                                                                                 | HMM families grouped into structural “clans”                                     |                                                                               | now hosted by Interpro                                           | Lean subset of InterPro—handy for rapid scans & benchmarking                                                        |
| **Gene Ontology (GO)** | 185 k controlled terms (BP/MF/CC) (June 2025 snapshot)([geneontology.org](https://geneontology.org/stats.html?utm_source=chatgpt.com "Release statistics - Gene Ontology"))     | Structured ontology; annotations stored separately                               | **AmiGO** (browse), **OBO-Edit** (curate); hundreds of enrichment tools       | [https://geneontology.org](https://geneontology.org)             | Standard vocabulary for functional enrichment; CC-BY 4.0                                                            |
| **KEGG** rel. 115      | 27 663 K numbers; 580 pathway maps([genome.jp](https://www.genome.jp/kegg/docs/statistics.html?utm_source=chatgpt.com "KEGG - Current Statistics"))                             | Genes grouped into KO orthologs; KO nodes populate manually drawn pathway graphs | **KAAS / BlastKOALA** assign K numbers                                        | [https://www.genome.jp/kegg](https://www.genome.jp/kegg)         | Gold-standard metabolic & signaling network backbone; licence required for bulk                                     |
| **Reactome** v90       | 2 742 curated human pathways (15 492 reactions)([reactome.org](https://reactome.org/about/news/261-v90-news?utm_source=chatgpt.com "V90 Released - Reactome Pathway Database")) | Expert pathway events → reactions → complexes                                    | Web analysis service; **ReactomeFIViz** Cytoscape app                         | [https://reactome.org](https://reactome.org)                     |                                                                                                                     |
| **MetaCyc** v28.5      | 3 153 pathways; 19 020 reactions; 19 372 metabolites([metacyc.org](https://metacyc.org/ "MetaCyc: Metabolic Pathways From all Domains of Life"))                                | Experimentally validated metabolic pathways across all life                      | **Pathway Tools** predicts networks & builds PGDBs                            | [https://metacyc.org](https://metacyc.org)                       | Ideal reference for gap-filling genome-scale models                                                                 |
| **eggNOG** v6.0        | 17 032 907 orthologous groups across 12 535 genomes                                                                                                                             | Hierarchical orthology clusters with functional tags                             | **eggNOG-mapper v2** (DIAMOND/MMseqs2 backend)                                | [https://eggnog6.embl.de](https://eggnog6.embl.de)               | Fast GO/KEGG/COG transfer; powers many metagenome pipelines                                                         |
| **OrthoDB** v12        | 162 M genes from 5 827 euks + 18 158 proks + 7 962 viruses([orthodb.org](https://www.orthodb.org/?utm_source=chatgpt.com "OrthoDB"))                                            | Ortholog groups at multiple phylogenetic depths                                  | Supplies lineage datasets for **BUSCO**, OrthoDB API                          | [https://www.orthodb.org](https://www.orthodb.org)               | Deepest taxonomic sampling; excellent for completeness checks                                                       |
| **Rfam** v15.0         | 4 178 non-coding RNA families([rfam.org](https://rfam.org/?utm_source=chatgpt.com "Rfam: The RNA families database"))                                                           | Covariance-model families (Infernal)                                             | **cmscan** / **Infernal** package                                             | [https://rfam.org](https://rfam.org)                             | Adds structured RNA & ribozyme annotation to genomes                                                                |
| **AlphaFold DB** 2024  | > 214 M predicted 3-D structures([alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/?utm_source=chatgpt.com "AlphaFold Protein Structure Database"))                             | Per-residue pLDDT confidence; complexes in v3                                    | **AlphaFold Server** / API; compatible with **FoldSeek** for structure search | [https://alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk)       | Structure-based function clues for “hypothetical” proteins                                                          |
**UniProt**: Swiss-Prot (manually curated: ~600000 entries) + TrEMBL (Unreviewed: ~250 milion entries).

InterPro

Gene Ontology (GO): Describes the function of gene products, using a dictionary of standardized words, it is a series of relations between controlled
vocabulary terms, they are available both in human readable format and as codes ()
https://geneontology.org/



### Wandering through databases 
Let's try some of the annotation resources with this not yet annotated protein that could be one of the output of our structural annotation workflow

> protein1
MKAISRVLIAMVAAIAALFTSTGTSHAGLDNELSLVDGQDRTLTVQQWDTFLNGVFPLDRNRLTREWFHSGRAKYIVAGPGADEFEGTLELGYQIGFPWSLGVGINFSYTTPNILIDDGDITAPPFGLNSVITPNLFPGVSISADLGNGPGIQEVATFSVDVSGAEGGVAVSNAHGTVTGAAGGVLLRPFARLIASTGDSVTTYGEPWNMN