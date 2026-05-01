# Mitogenome assembly, annotation and characterisation pipeline for Mediterranean elasmobranchs
### Article title: Novel mitogenomic resources for Mediterranean elasmobranchs

---

## 1. Paired-end mitochondrial genome assembly using GetOrganelle v1.7.7 (Jin et al. 2020)
[https://github.com/kinggerm/getorganelle]

Some samples had to be individually optimised. This can be run locally on an Ubuntu system (A) or on an HPC (B). Use the `animal_mt` database.

### (A) Single sample
Usage: `bash getorganelle.sh <sample_name> <R1.fastq.gz> <R2.fastq.gz> <seed.fasta>`

```bash
SAMPLE=$1
R1=$2
R2=$3
SEED=$4

mkdir -p ${SAMPLE}
cd ${SAMPLE}

echo "=== Running GetOrganelle for $SAMPLE ==="
get_organelle_from_reads.py \
  -1 ../${R1} \
  -2 ../${R2} \
  -s ../${SEED} \
  -o go_${SAMPLE} \
  -F animal_mt \
  -t 8 \
  -R 20 \
  -k 21,45,65,85,105

# rename main output
cd go_${SAMPLE}
mv *.path_sequence.fasta ${SAMPLE}.complete.path_sequence.fasta
mv *.fastg ${SAMPLE}.assembly_graph.fastg
mv *.selected_graph.gfa ${SAMPLE}.selected_graph.gfa
mv get_org.log.txt ${SAMPLE}.get_org.log.txt

echo "=== All done for $SAMPLE! ==="
```

### (A) Batch job
Usage: `bash getorganelle.sh <input_directory> <output_directory> <seed_directory>`

```bash
RAW_DIR=$1
OUT_DIR=$2
SEED_DIR=$3
THREADS=8
ROUNDS=20
KMERS="21,45,65,85,105"

mkdir -p "$OUT_DIR"
FAILED="$OUT_DIR/failed_samples.txt"

for R1 in "${RAW_DIR}"/*_R1_001.fastq*; do
    R2=${R1/_R1_/_R2_}
    [[ ! -f "$R2" ]] && echo "Missing R2 for $R1, skipping..." && continue

    sample=$(basename "$R1" | sed -E 's/_R1.*$//')
    species=$(basename "$R1" | awk -F'[-_]' '{print $1 "-" $2}')
    seed="${SEED_DIR}/${species}_seed.fasta"

    [[ ! -f "$seed" ]] && echo "No seed for $species, skipping..." && continue

    sample_out="${OUT_DIR}/${sample}_go"
    mkdir -p "$sample_out"

    echo "=== Running GetOrganelle for $sample ==="
    if ! get_organelle_from_reads.py \
        -1 "$R1" \
        -2 "$R2" \
        -s "$seed" \
        -o "$sample_out" \
        -F animal_mt \
        -t $THREADS \
        -R $ROUNDS \
        -k $KMERS \
        --overwrite \
        > "${sample_out}/${sample}_getorganelle.log" 2>&1; then

        echo "GetOrganelle crashed for $sample - check log!"
        echo "$sample" >> "$FAILED"
        continue
    fi

    # rename main outputs
    cd "$sample_out"
    if ls *.path_sequence.fasta 1> /dev/null 2>&1; then
        mv *.path_sequence.fasta "${sample}.complete.path_sequence.fasta"
        mv *.fastg "${sample}.assembly_graph.fastg" 2>/dev/null || true
        mv *.selected_graph.gfa "${sample}.selected_graph.gfa" 2>/dev/null || true
        echo "Finished $sample successfully"
    else
        echo "No FASTA output for $sample"
        echo "$sample" >> "$FAILED"
    fi
    cd - > /dev/null
done

echo "=== DONE ==="
```

### (B) HPC

```bash
get_organelle_from_reads.py \
  -1 /path/to/raw_data/Sample_R1.fastq.gz \
  -2 /path/to/raw_data/Sample_R2.fastq.gz \
  -s /path/to/seeds/Sample_seed.fasta \
  -o /path/to/output_folder/Sample_go \
  -F animal_mt \
  -t 12 \
  -R 20 \
  -k 21,45,65,85,105
```

Visualise `.gfa` graph in Bandage (Wick et al. 2015).

---

## 2. Annotation
Annotation using MitoAnnotator (Iwasaki et al. 2013) on the MitoFish webserver [https://mitofish.aori.u-tokyo.ac.jp/annotation/input/] and MITOS2 (Donath et al. 2019) on the GalaxyEuro platform [https://usegalaxy.eu]. Use MITOS2 to guide manual ATP8 annotation in Geneious Prime v2024.0.5 (Kearse et al. 2012), if missed, and perform manual curation/checks.

---

## 3. Assembly statistics and visualisation
Calculate assembly statistics based on one representative individual per species in R. Obtain values from Geneious and calculate skews: AT skew = (A−T)/(A+T), GC skew = (G−C)/(G+C).

```r
library(readxl)
library(tidyverse)
library(viridis)
library(ggplot2)

col_w <- 80 / 25.4

df <- read_excel("Assembly-statistics.xlsx", sheet = 1)
df$Species <- factor(df$Species, levels = unique(df$Species))

df_long <- df %>%
  pivot_longer(
    cols = c(`A (%)`, `T (%)`, `G (%)`, `C (%)`),
    names_to = "Base",
    values_to = "Percent"
  ) %>%
  mutate(
    `Genome size (bp)` = as.numeric(`Genome size (bp)`),
    Percent = as.numeric(Percent),
    Base_bp = `Genome size (bp)` * (Percent / 100)
  )

df_long$Base    <- factor(df_long$Base, levels = c("A (%)", "T (%)", "G (%)", "C (%)"))
df_long$Species <- factor(df_long$Species, levels = unique(df_long$Species))

# Figure 1 — base composition (legend at bottom)
p1 <- ggplot(df_long, aes(x = Species, y = Base_bp, fill = Base)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(
    values = c(
      "A (%)" = "#66c2a5",
      "T (%)" = "#fc8d62",
      "G (%)" = "#8da0cb",
      "C (%)" = "#e78ac3"
    )
  ) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme_bw(base_size = 7) +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 7),
    legend.title     = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 5),
    legend.key.size  = unit(0.3, "cm"),
    legend.margin    = margin(t = -8),
    panel.grid.minor = element_blank(),
    plot.title       = element_blank()
  ) +
  labs(y = "Genome size (bp)")

ggsave("compositionsize.png", plot = p1, width = col_w, height = 100 / 25.4,
       dpi = 600, units = "in")
ggsave("compositionsize.pdf", plot = p1, width = col_w, height = 100 / 25.4,
       units = "in")

# Figure 1 — base composition (legend inside top left)
p1 <- ggplot(df_long, aes(x = Species, y = Base_bp, fill = Base)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(
    values = c(
      "A (%)" = "#66c2a5",
      "T (%)" = "#fc8d62",
      "G (%)" = "#8da0cb",
      "C (%)" = "#e78ac3"
    )
  ) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme_bw(base_size = 7) +
  theme(
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 7),
    legend.title      = element_blank(),
    legend.position   = c(0.10, 0.90),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.text       = element_text(size = 6),
    legend.key.size   = unit(0.3, "cm"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_blank()
  ) +
  labs(y = "Genome size (bp)")

ggsave("compositionsize_T.png", plot = p1, width = col_w, height = 100 / 25.4,
       dpi = 600, units = "in")
ggsave("compositionsize_T.pdf", plot = p1, width = col_w, height = 100 / 25.4,
       units = "in")

# Figure 2 — nucleotide skew
ray_species <- c(
  "A. bovinus", "B. lata", "D. marmorata", "G. cemiculus",
  "L. melitensis", "R. asterias", "R. montagui", "R. rhinobatos",
  "R. alba", "T. nobiliana", "T. torpedo"
)

df <- df %>%
  mutate(
    Clade = ifelse(Species %in% ray_species, "Batoidea (rays)", "Selachii (sharks)"),
    Clade = factor(Clade, levels = c("Selachii (sharks)", "Batoidea (rays)"))
  )

df$Species <- factor(df$Species, levels = unique(df$Species))

shark_species <- df %>%
  filter(Clade == "Selachii (sharks)") %>%
  pull(Species) %>% unique()

ray_species_in_df <- df %>%
  filter(Clade == "Batoidea (rays)") %>%
  pull(Species) %>% unique()

shark_cols <- colorRampPalette(c("#1b9e77", "#66c2a5", "#3288bd"))(length(shark_species))
ray_cols   <- colorRampPalette(c("#d7301f", "#fc8d59", "#fee08b"))(length(ray_species_in_df))

species_cols <- c(
  setNames(shark_cols, shark_species),
  setNames(ray_cols,   ray_species_in_df)
)

p <- ggplot(df, aes(x = `AT skew`, y = `GC skew`, color = Species)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = species_cols) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10, face = "italic"),
    plot.title       = element_blank(),
    plot.subtitle    = element_text(size = 11)
  ) +
  labs(x = "AT skew", y = "GC skew")

ggsave("nucleotideskew.png", plot = p, width = 180/25.4, height = 8,
       dpi = 600, units = "in")
ggsave("nucleotideskew.pdf", plot = p, width = 180/25.4, height = 8,
       units = "in")
```

---

## 4. Control Region characterisation
Using Tandem Repeats Finder (Benson 1999) [https://tandem.bu.edu/trf/trf.html]. Input individual Control Region sequences as `.fasta` files.

---

## 5. tRNA structure prediction and visualisation

tRNA genes predicted using MiTFi (Jühling et al. 2012) distributed as part of MitoFinder 
(Allio et al. 2020) [https://github.com/RemiAllio/MitoFinder/tree/master/mitfi]. Ensure all dependencies are installed locally before running. MiTFi was run locally on WSL - the script below is a batch wrapper around the original MiTFi command to process multiple `.fasta` files automatically.

Usage: `bash mitfi_run.sh`

```bash
mkdir -p mitfi_out

# loop over all .fasta files in input directory
for f in *.fasta; do
    [ -e "$f" ] || continue

    out="mitfi_out/${f%.fasta}.struct.fa"
    echo "Running MiTFi on $f  ->  $out"
    java -jar mitfi.jar -fasta -structure "$f" > "$out"
done

echo "All done. Results saved in mitfi_out/"
```

tRNA secondary structures visualised using RNAplot from the ViennaRNA package (Lorenz et al. 2011) via a custom Python script. The script reads MiTFi `.struct.fa` output files, runs RNAplot on each tRNA, post-processes the SVG output, and arranges all tRNA structures into a grid figure per species.

Usage: `python3 draw_trna.py <species.struct.fa>`

```python
import os
import sys
import re
import subprocess
import tempfile
import math

def parse_struct_fa(filepath):
    trnas = []
    with open(filepath) as f:
        lines = [l.rstrip() for l in f if l.strip()]
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header = lines[i]
            seq = lines[i+1] if i+1 < len(lines) else ''
            db  = lines[i+2] if i+2 < len(lines) else ''
            db  = re.sub(r'[^.()\[\]]', '', db)
            minlen = min(len(seq), len(db))
            trnas.append((header, seq[:minlen], db[:minlen]))
            i += 3
        else:
            i += 1
    return trnas

def get_tRNA_name(header):
    parts = header.split('|')
    aa_map = {
        'A':'Alanine','C':'Cysteine','D':'Aspartic acid','E':'Glutamic acid',
        'F':'Phenylalanine','G':'Glycine','H':'Histidine','I':'Isoleucine',
        'K':'Lysine','L1':'Leucine-1','L2':'Leucine-2','M':'Methionine',
        'N':'Asparagine','P':'Proline','Q':'Glutamine','R':'Arginine',
        'S1':'Serine-1','S2':'Serine-2','T':'Threonine','V':'Valine',
        'W':'Tryptophan','Y':'Tyrosine'
    }
    if len(parts) >= 7:
        return aa_map.get(parts[6], parts[6])
    return header

def get_species_name(filepath):
    base = os.path.splitext(os.path.basename(filepath))[0]
    name = re.sub(r'[-_](EL|RR|PAL)\w+.*$', '', base)
    return name.replace('-', ' ')

def rnaplot_svg(seq, db, name, workdir, idx=0):
    safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', name) + f'_{idx}'
    infile = os.path.join(workdir, f'{safe_name}.fa')
    with open(infile, 'w') as f:
        f.write(f'>{safe_name}\n{seq}\n{db}\n')
    subprocess.run(
        ['RNAplot', '-f', 'svg', '-t', '1', '-i', infile],
        cwd=workdir,
        capture_output=True
    )
    outfile = os.path.join(workdir, f'{safe_name}_ss.svg')
    if os.path.exists(outfile):
        return outfile
    return None

def post_process_svg(svgfile):
    with open(svgfile, 'r', encoding='utf-8') as f:
        content = f.read()
    content = re.sub(r'<\?xml[^>]+\?>', '', content).strip()
    content = re.sub(r'<script[^>]*>.*?</script>', '', content, flags=re.DOTALL)
    content = re.sub(r'(stroke:\s*)grey', r'\1#d0d0d0', content)
    content = re.sub(r'(stroke:\s*)red', r'\1#7aade0', content)
    inner = re.sub(r'^<svg[^>]*>', '', content).strip()
    inner = re.sub(r'</svg>\s*$', '', inner).strip()
    return inner

def make_grid(svg_files, labels, output_path, ncols=8, species_name=''):
    src_size  = 452
    cell_size = 200
    scale     = cell_size / src_size
    label_h   = 20
    title_h   = 50
    pad       = 15

    n = len(svg_files)
    nrows = math.ceil(n / ncols)
    total_w = pad + ncols * (cell_size + pad)
    total_h = title_h + nrows * (cell_size + label_h + pad) + pad

    lines = []
    lines.append('<?xml version="1.0" encoding="UTF-8"?>')
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
                 f'width="{total_w}" height="{total_h}" version="1.1">')
    lines.append(f'  <rect width="{total_w}" height="{total_h}" fill="white"/>')
    lines.append(f'  <text x="{total_w//2}" y="35" text-anchor="middle" '
                 f'font-size="20" font-style="italic" font-family="Arial" '
                 f'fill="black">{species_name}</text>')

    for idx, (svgf, label) in enumerate(zip(svg_files, labels)):
        row, col = divmod(idx, ncols)
        x = pad + col * (cell_size + pad)
        y = title_h + row * (cell_size + label_h + pad)

        if svgf and os.path.exists(svgf):
            inner = post_process_svg(svgf)
            lines.append(f'  <g transform="translate({x},{y}) scale({scale:.6f})">')
            lines.append(inner)
            lines.append('  </g>')

        lines.append(f'  <text x="{x + cell_size//2}" y="{y + cell_size + 15}" '
                     f'text-anchor="middle" font-size="11" font-family="Arial" '
                     f'fill="black">{label}</text>')

    lines.append('</svg>')

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'Saved: {output_path}')

def process_species(filepath):
    trnas = parse_struct_fa(filepath)
    if not trnas:
        print(f'No tRNAs found in {filepath}')
        return

    species = get_species_name(filepath)
    base = os.path.splitext(filepath)[0]
    output_svg = base + '_tRNA_grid.svg'

    with tempfile.TemporaryDirectory() as workdir:
        svg_files = []
        labels = []
        for header, seq, db in trnas:
            name = get_tRNA_name(header)
            print(f'  Plotting {name}...')
            svgf = rnaplot_svg(seq, db, name, workdir, idx=len(svg_files))
            svg_files.append(svgf)
            labels.append(name)
        make_grid(svg_files, labels, output_svg, ncols=8, species_name=species)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python3 draw_trna.py <species.struct.fa>')
        sys.exit(1)
    process_species(sys.argv[1])
```

---

## 6. Relative synonymous codon usage (RSCU)
Calculate RSCU of the PCGs in DAMBE v7.3.32 (Xia 2018). Export PCGs from Geneious to use as input.

---

## 7. Ka/Ks analysis

Calculate Ks (synonymous substitutions) and Ka (non-synonymous substitutions) pairwise values between all PCGs in MEGA 11 (Tamura et al. 2021) using the Nei-Gojobori model. Copy results into Excel for input.

```r
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
library(forcats)
library(gridExtra)

ks_file <- "Synonymous_substitutions_Ks.xlsx"
ka_file <- "Nonsynonymous_substitutions_Ka.xlsx"

make_ka_ks_table_for_gene <- function(gene, ks_file, ka_file) {
  ks <- read_excel(ks_file, sheet = gene, na = c("", "?"))
  ka <- read_excel(ka_file, sheet = gene, na = c("", "?"))
  names(ks)[1] <- "speciesA"
  names(ka)[1] <- "speciesA"
  ks_long <- ks %>%
    pivot_longer(cols = -speciesA, names_to = "speciesB", values_to = "Ks") %>%
    filter(!is.na(Ks))
  ka_long <- ka %>%
    pivot_longer(cols = -speciesA, names_to = "speciesB", values_to = "Ka") %>%
    filter(!is.na(Ka))
  out <- ks_long %>%
    full_join(ka_long, by = c("speciesA", "speciesB")) %>%
    filter(speciesA != speciesB) %>%
    mutate(`Ka/Ks` = Ka / Ks) %>%
    relocate(Ks, Ka, `Ka/Ks`, .after = speciesB) %>%
    mutate(Gene = gene, .before = 1)
  return(out)
}

# use function for each gene
ATP6_table  <- make_ka_ks_table_for_gene("ATP6",  ks_file, ka_file)
ATP8_table  <- make_ka_ks_table_for_gene("ATP8",  ks_file, ka_file)
COI_table   <- make_ka_ks_table_for_gene("COI",   ks_file, ka_file)
COII_table  <- make_ka_ks_table_for_gene("COII",  ks_file, ka_file)
COIII_table <- make_ka_ks_table_for_gene("COIII", ks_file, ka_file)
CytB_table  <- make_ka_ks_table_for_gene("CytB",  ks_file, ka_file)
ND1_table   <- make_ka_ks_table_for_gene("ND1",   ks_file, ka_file)
ND2_table   <- make_ka_ks_table_for_gene("ND2",   ks_file, ka_file)
ND3_table   <- make_ka_ks_table_for_gene("ND3",   ks_file, ka_file)
ND4_table   <- make_ka_ks_table_for_gene("ND4",   ks_file, ka_file)
ND4L_table  <- make_ka_ks_table_for_gene("ND4L",  ks_file, ka_file)
ND5_table   <- make_ka_ks_table_for_gene("ND5",   ks_file, ka_file)
ND6_table   <- make_ka_ks_table_for_gene("ND6",   ks_file, ka_file)

# manually add results into one file ("KaKs_all_vs_all.xlsx") with a separate sheet for each gene

# read combined file — only KaKs_all_vs_all.xlsx needed from here onwards
file <- "KaKs_all_vs_all.xlsx"
genes <- excel_sheets(file)

all_kaks <- lapply(genes, function(g) {
  df <- read_excel(file, sheet = g)
  if ("Ka/Ks" %in% names(df)) df <- df %>% rename(Ka_Ks = `Ka/Ks`)
  df %>% mutate(Gene = g)
}) %>%
  bind_rows() %>%
  filter(!is.na(Ka_Ks), is.finite(Ka_Ks))

# rename genes to standard nomenclature
all_kaks <- all_kaks %>%
  mutate(Gene = recode(Gene,
    "COI"   = "cox1", "COII"  = "cox2", "COIII" = "cox3",
    "CytB"  = "cob",  "ATP6"  = "atp6", "ATP8"  = "atp8",
    "ND1"   = "nd1",  "ND2"   = "nd2",  "ND3"   = "nd3",
    "ND4"   = "nd4",  "ND4L"  = "nd4l", "ND5"   = "nd5",
    "ND6"   = "nd6"
  ))

full_w <- 180 / 25.4

kaks_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "italic", size = 10),
    axis.text.y     = element_text(size = 10),
    axis.title.y    = element_text(size = 11),
    plot.title      = element_text(size = 12, face = "bold")
  )

# combined — no title (main text figure)
all <- ggplot(all_kaks, aes(x = Gene, y = Ka_Ks, fill = Gene)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.9) +
  kaks_theme +
  theme(plot.title = element_blank()) +
  ylab("Ka/Ks") +
  xlab(NULL)

ggsave("KaKs_all.png", plot = all, width = full_w, height = 150/25.4, dpi = 600, units = "in")
ggsave("KaKs_all.pdf", plot = all, width = full_w, height = 150/25.4, units = "in")

# filter rays and sharks
ray_species <- c(
  "A-bovinus-EL1536", "B-lata-RR002810", "D-marmorata-RR001836_T",
  "G-cemiculus-PAL07", "L-melitensis-RR002790", "R-asterias-EL3420",
  "R-montagui-RR002299", "R-rhinobatos-EL1571", "R-alba-021",
  "T-nobiliana-EL0254_minimal", "T-torpedo-RR003058"
)

all_kaks_rays   <- all_kaks %>% filter(speciesA %in% ray_species, speciesB %in% ray_species)
all_kaks_sharks <- all_kaks %>% filter(!(speciesA %in% ray_species), !(speciesB %in% ray_species))

# rays (supplementary)
raja <- ggplot(all_kaks_rays, aes(x = Gene, y = Ka_Ks, fill = Gene)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.9) +
  kaks_theme +
  ylab("Ka/Ks") +
  xlab(NULL) +
  ggtitle("Figure S45: Pairwise Ka/Ks of 13 PCGs in 11 species of rays")

ggsave("KaKs_raja.png", plot = raja, width = full_w, height = 150/25.4, dpi = 600, units = "in")
ggsave("KaKs_raja.pdf", plot = raja, width = full_w, height = 150/25.4, units = "in")

# sharks (supplementary)
squali <- ggplot(all_kaks_sharks, aes(x = Gene, y = Ka_Ks, fill = Gene)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.9) +
  kaks_theme +
  ylab("Ka/Ks") +
  xlab(NULL) +
  ggtitle("Figure S46: Pairwise Ka/Ks of 13 PCGs in 11 species of sharks")

ggsave("KaKs_squali.png", plot = squali, width = full_w, height = 150/25.4, dpi = 600, units = "in")
ggsave("KaKs_squali.pdf", plot = squali, width = full_w, height = 150/25.4, units = "in")

# combined supplementary PDF (rays + sharks stacked)
pdf("KaKs_supp_combined.pdf", width = full_w, height = 2 * 150/25.4)
grid.arrange(raja, squali, nrow = 2)
dev.off()
```

---

## 8. Comparison of DNA extraction protocols
Comparison of (a) standard salting-out protocol and (b) exonuclease-treated salting-out protocol. Map raw reads to the respective assembled mitogenomes using BWA v0.7.13 (Li & Durbin 2009) and SAMtools v1.19.2 (Li et al. 2009). Calculate the percentage of mapped reads for each sample.

Note: ensure FASTA names match the raw FASTQ prefix (`_T` suffix for treated samples).

```bash
THREADS=12
FASTA_DIR="/path/to/assembled.fasta"

echo -e "sample\ttotal_reads\tmapped_reads\tpct_mapped" > mt_summary.tsv

for R1 in *_R1.fastq.gz; do
    [ -e "$R1" ] || continue

    PREFIX="${R1%_R1.fastq.gz}"
    R2="${PREFIX}_R2.fastq.gz"
    REF="${FASTA_DIR}/${PREFIX}.fasta"
    BAM="${PREFIX}.mtmap.bam"

    if [ ! -f "$R2" ]; then
        echo "WARNING: Missing R2 for ${PREFIX}, skipping." >&2
        continue
    fi

    if [ ! -f "$REF" ]; then
        echo "WARNING: Missing fasta for ${PREFIX}, skipping." >&2
        continue
    fi

    # index the mitogenome if needed
    if [ ! -f "${REF}.bwt" ]; then
        bwa index "$REF"
    fi

    # map raw reads to the assembled mitogenome
    bwa mem -t $THREADS "$REF" "$R1" "$R2" | \
        samtools view -@ $THREADS -b -o "$BAM" -

    # count total and mapped reads
    total=$(samtools view -@ $THREADS -c "$BAM")
    mapped=$(samtools view -@ $THREADS -c -F 4 "$BAM")
    pct=$(awk -v m="$mapped" -v t="$total" 'BEGIN {printf "%.4f", (m/t)*100}')

    echo -e "${PREFIX}\t${total}\t${mapped}\t${pct}" >> mt_summary.tsv
done

echo "Done. Summary written to mt_summary.tsv"
```

---

## 9. Genome size estimation
Estimate genome size using Jellyfish v2.3.1 (Marçais & Kingsford 2011) and the Genome Size Estimation Tool (GSET; Liao et al. 2024) in MATLAB.

### Trim reads using fastp

```bash
fastp \
  -i /path/to/input_folder/sample_R1.fastq.gz \
  -I /path/to/input_folder/sample_R2.fastq.gz \
  -o /path/to/trimmed_folder/sample_R1.trimmed.fastq.gz \
  -O /path/to/trimmed_folder/sample_R2.trimmed.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --trim_poly_g \
  --cut_tail \
  -w 12 \
  -h /path/to/output_folder/sample_fastp.html \
  -j /path/to/output_folder/sample_fastp.json
```

The `.html` file contains the L-value (mean read length after trimming).

### Run Jellyfish to obtain k-mer histogram

```bash
R1=/path/to/trimmed_folder/sample_R1.trimmed.fastq.gz
R2=/path/to/trimmed_folder/sample_R2.trimmed.fastq.gz

# count k-mers
zcat $R1 $R2 | jellyfish count \
  -C \
  -m 21 \
  -s 5G \
  -t 12 \
  -o sample_21mer.jf \
  /dev/fd/0

# generate histogram
jellyfish histo \
  -t 12 \
  sample_21mer.jf \
  > sample_21mer.histo
```

### Run GSET in MATLAB
[https://github.com/Xingyu-Liao/GSET] — ensure dependencies are downloaded.

```matlab
% usage: GSE('sample_kmer.histo', L, K)
% L = mean read length from fastp, K = k-mer size (21)
GSE('sample_21mer.histo', L, 21);
```
