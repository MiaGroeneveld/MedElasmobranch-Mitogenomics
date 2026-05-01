# MedElasmoBranch-Mitogenomics

Bioinformatics pipeline for mitochondrial genome assembly, annotation and characterisation of Mediterranean elasmobranchs.

---

## Associated publication

> **Novel mitogenomic resources for Mediterranean elasmobranchs (Spiga et al.)**
> *(in preparation)*

---

## Repository contents

| File | Description |
|------|-------------|
| `pipeline.md` | Full annotated pipeline with all scripts and usage instructions |
| `draw_trna.py` | Python script for tRNA secondary structure visualisation using RNAplot |

---

## Pipeline overview

1. **Mitogenome assembly** — GetOrganelle v1.7.7
2. **Annotation** — MitoAnnotator / MITOS2 / Geneious Prime
3. **Assembly statistics** — R (ggplot2): base composition, nucleotide skew
4. **Control Region characterisation** — Tandem Repeats Finder
5. **tRNA structure prediction** — MiTFi + RNAplot (ViennaRNA)
6. **Codon usage** — DAMBE v7.3.32 (RSCU)
7. **Ka/Ks analysis** — MEGA 11 + R (ggplot2)
8. **DNA extraction protocol comparison** — BWA + SAMtools
9. **Genome size estimation** — Jellyfish + GSET (MATLAB)

---

## Dependencies

### Command line
- [GetOrganelle v1.7.7](https://github.com/kinggerm/getorganelle)
- [MiTFi](https://github.com/RemiAllio/MitoFinder/tree/master/mitfi)
- [ViennaRNA / RNAplot](https://www.tbi.univie.ac.at/RNA/)
- [BWA v0.7.13](https://github.com/lh3/bwa)
- [SAMtools v1.19.2](https://www.htslib.org/)
- [Jellyfish v2.3.1](https://github.com/gmarcais/Jellyfish)
- [fastp](https://github.com/OpenGene/fastp)
- [GSET](https://github.com/Xingyu-Liao/GSET)

### GUI / web
- [Bandage](https://rrwick.github.io/Bandage/)
- [MitoFish / MitoAnnotator](https://mitofish.aori.u-tokyo.ac.jp/annotation/input/)
- [MITOS2 on Galaxy Europe](https://usegalaxy.eu)
- [Geneious Prime v2024.0.5](https://www.geneious.com/)
- [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html)
- [MEGA 11](https://www.megasoftware.net/)
- [DAMBE v7.3.32](http://dambe.bio.uottawa.ca/)

### R packages
```r
install.packages(c("readxl", "tidyverse", "ggplot2", "viridis",
                   "dplyr", "tidyr", "writexl", "forcats", "gridExtra"))
```

### Python packages
```bash
pip install viennarna
```

---

### R scripts
All R code is embedded in `pipeline.md` with inline explanations

---

## Citation

If you use this pipeline, please cite the associated publication (above) and the relevant tools listed under Dependencies.

---
