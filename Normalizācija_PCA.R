
# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager")
} else {
  message("BiocManager is already installed.")
}

# Check if DESeq2 is installed
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  message("Installing DESeq2...")
  BiocManager::install("DESeq2", dependecies = TRUE)
} else {
  message("DESeq2 is already installed.")
}

# Load DESeq2
library(DESeq2)

message("DESeq2 loaded successfully!")

meta <- read_excel("Metatabula_baka.xlsx")

#Ielādē gene count matrix
counts <- read.delim(
  "Gene.counts.strand2.tsv",
  check.names = FALSE
)

dim(meta)
dim(counts)

head(meta)
counts[1:6, 1:6]

#paskatāmies pirmās rindas
head(counts$Gene, 10)

gene_counts <- counts[grepl("^ENSG", counts$Gene), ]

dim(gene_counts)
head(gene_counts$Gene)

rownames(gene_counts) <- gene_counts$Gene
gene_counts$Gene <- NULL

#Cik kopā readu ir katrā paraugā
library_sizes <- colSums(gene_counts)
#kopsavilkums par library sizes
summary(library_sizes)
#Paraugi sakārtoti no mazākā līdz lielākajam readu skaitam
sort(library_sizes)

#FILTRĒŠANA

# Izveidojam pilno DESeq2 objektu no sākotnējā gene_counts
dds_full <- DESeqDataSetFromMatrix(
  countData = gene_counts,
  colData = meta,
  design = ~ 1
)

# Pārbaudām sākotnējo izmēru
dim(dds_full)

# Saglabājam darba kopiju
dds <- dds_full

# Pārbaudām, cik gēniem visos 72 paraugos count = 0
sum(rowSums(counts(dds)) == 0)

# Maigā filtrēšana:
# izmetam tikai tos gēnus, kuriem visos 72 sample count = 0
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]

# Pārbaudām, cik gēnu palika pēc filtrēšanas
dim(dds)

# Pārbaudām, vai tiešām vairs nav gēnu ar 0 visos sample
sum(rowSums(counts(dds)) == 0)




# =========================
# PCA PIRMS NORMALIZĀCIJAS
# =========================
library(ggplot2)
library(patchwork)

# ==========================================
# PCA PIRMS NORMALIZĀCIJAS
# DIETA + FĀZE
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)

# 1. RAW counts log-transformācija
raw_log <- log2(counts(dds) + 1)

# 2. Aprēķina variāciju katram gēnam
rv_raw <- rowVars(raw_log)

# 3. Izvēlas 500 visvariablākos gēnus
select_raw <- order(
  rv_raw,
  decreasing = TRUE
)[1:500]

# 4. PCA tikai ar top 500 gēniem
pca_raw <- prcomp(
  t(raw_log[select_raw, ]),
  scale. = FALSE
)

# 5. PCA koordinātes + metadata
pca_df <- data.frame(
  PC1 = pca_raw$x[,1],
  PC2 = pca_raw$x[,2],
  PC3 = pca_raw$x[,3],
  PC4 = pca_raw$x[,4],
  diet = meta$diet,
  phase = meta$phase
)

# 6. Apvienojam diet + phase vienā grupā
pca_df$diet_phase <- interaction(
  pca_df$diet,
  pca_df$phase,
  sep = " "
)

# 7. PC1 vs PC2
p12 <- ggplot(
  pca_df,
  aes(PC1, PC2, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC2",
    color = "Diet + phase"
  )

# 8. PC1 vs PC3
p13 <- ggplot(
  pca_df,
  aes(PC1, PC3, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC3",
    color = "Diet + phase"
  )

# 9. PC1 vs PC4
p14 <- ggplot(
  pca_df,
  aes(PC1, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC4",
    color = "Diet + phase"
  )

# 10. PC2 vs PC3
p23 <- ggplot(
  pca_df,
  aes(PC2, PC3, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC3",
    color = "Diet + phase"
  )

# 11. PC2 vs PC4
p24 <- ggplot(
  pca_df,
  aes(PC2, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC4",
    color = "Diet + phase"
  )

# 12. PC3 vs PC4
p34 <- ggplot(
  pca_df,
  aes(PC3, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC3 vs PC4",
    color = "Diet + phase"
  )

# 13. Saliek visus grafikus kopā
(p12 | p13 | p14) /
  (p23 | p24 | p34) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")



# ==========================================
# PCA PIRMS NORMALIZĀCIJAS
# DZIMUMS
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)

# 1. RAW counts log-transformācija
raw_log <- log2(counts(dds) + 1)

# 2. Aprēķina variāciju katram gēnam
rv_raw <- rowVars(raw_log)

# 3. Izvēlas 500 visvariablākos gēnus
select_raw <- order(
  rv_raw,
  decreasing = TRUE
)[1:500]

# 4. PCA tikai ar top 500 gēniem
pca_raw <- prcomp(
  t(raw_log[select_raw, ]),
  scale. = FALSE
)

# 5. PCA koordinātes + dzimums
pca_df_sex <- data.frame(
  PC1 = pca_raw$x[,1],
  PC2 = pca_raw$x[,2],
  PC3 = pca_raw$x[,3],
  PC4 = pca_raw$x[,4],
  sex = meta$sex
)

# PC1 vs PC2
p12_sex <- ggplot(
  pca_df_sex,
  aes(PC1, PC2, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC2",
    color = "Sex"
  )

# PC1 vs PC3
p13_sex <- ggplot(
  pca_df_sex,
  aes(PC1, PC3, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC3",
    color = "Sex"
  )

# PC1 vs PC4
p14_sex <- ggplot(
  pca_df_sex,
  aes(PC1, PC4, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC4",
    color = "Sex"
  )

# PC2 vs PC3
p23_sex <- ggplot(
  pca_df_sex,
  aes(PC2, PC3, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC3",
    color = "Sex"
  )

# PC2 vs PC4
p24_sex <- ggplot(
  pca_df_sex,
  aes(PC2, PC4, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC4",
    color = "Sex"
  )

# PC3 vs PC4
p34_sex <- ggplot(
  pca_df_sex,
  aes(PC3, PC4, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC3 vs PC4",
    color = "Sex"
  )

# Visi kopā
(p12_sex | p13_sex | p14_sex) /
  (p23_sex | p24_sex | p34_sex) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ==========================================
# PCA PIRMS NORMALIZĀCIJAS
# TIKAI SIEVIETES
# DIETA + FĀZE
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)

# 1. Atlasām tikai sieviešu sample
female_samples <- meta$sex == "F"

# 2. RAW counts tikai sievietēm
female_counts <- counts(dds)[, female_samples]

# 3. Log transformācija
female_log <- log2(female_counts + 1)

# 4. Atrodam top 500 visvariablākos gēnus
rv_female <- rowVars(female_log)

select_female <- order(
  rv_female,
  decreasing = TRUE
)[1:500]

# 5. PCA tikai sievietēm
pca_female <- prcomp(
  t(female_log[select_female, ]),
  scale. = FALSE
)

# 6. Metadata tikai sievietēm
meta_female <- meta[female_samples, ]

# 7. PCA koordinātes + diet + phase
pca_female_df <- data.frame(
  PC1 = pca_female$x[,1],
  PC2 = pca_female$x[,2],
  PC3 = pca_female$x[,3],
  PC4 = pca_female$x[,4],
  diet = meta_female$diet,
  phase = meta_female$phase
)

pca_female_df$diet_phase <- interaction(
  pca_female_df$diet,
  pca_female_df$phase,
  sep = " "
)

# PC1 vs PC2
p12_f <- ggplot(pca_female_df, aes(PC1, PC2, color = diet_phase)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC1 vs PC2", color = "Diet + phase")

# PC1 vs PC3
p13_f <- ggplot(pca_female_df, aes(PC1, PC3, color = diet_phase)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC1 vs PC3", color = "Diet + phase")

# PC1 vs PC4
p14_f <- ggplot(pca_female_df, aes(PC1, PC4, color = diet_phase)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC1 vs PC4", color = "Diet + phase")

# PC2 vs PC3
p23_f <- ggplot(pca_female_df, aes(PC2, PC3, color = diet_phase)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC2 vs PC3", color = "Diet + phase")

# PC2 vs PC4
p24_f <- ggplot(pca_female_df, aes(PC2, PC4, color = diet_phase)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC2 vs PC4", color = "Diet + phase")

# PC3 vs PC4
p34_f <- ggplot(pca_female_df, aes(PC3, PC4, color = diet_phase)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC3 vs PC4", color = "Diet + phase")

# Visi kopā
(p12_f | p13_f | p14_f) /
  (p23_f | p24_f | p34_f) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# ==========================================
# PCA PIRMS NORMALIZĀCIJAS
# VECUMS
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)

# 1. RAW counts log-transformācija
raw_log <- log2(counts(dds) + 1)

# 2. Variācija katram gēnam
rv_raw <- rowVars(raw_log)

# 3. Top 500 visvariablākie gēni
select_raw <- order(
  rv_raw,
  decreasing = TRUE
)[1:500]

# 4. PCA
pca_raw <- prcomp(
  t(raw_log[select_raw, ]),
  scale. = FALSE
)

# 5. PCA koordinātes + vecums
pca_df_age <- data.frame(
  PC1 = pca_raw$x[,1],
  PC2 = pca_raw$x[,2],
  PC3 = pca_raw$x[,3],
  PC4 = pca_raw$x[,4],
  age = meta$age
)

# PC1 vs PC2
p12_age <- ggplot(pca_df_age, aes(PC1, PC2, color = age)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC1 vs PC2", color = "Age")

# PC1 vs PC3
p13_age <- ggplot(pca_df_age, aes(PC1, PC3, color = age)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC1 vs PC3", color = "Age")

# PC1 vs PC4
p14_age <- ggplot(pca_df_age, aes(PC1, PC4, color = age)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC1 vs PC4", color = "Age")

# PC2 vs PC3
p23_age <- ggplot(pca_df_age, aes(PC2, PC3, color = age)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC2 vs PC3", color = "Age")

# PC2 vs PC4
p24_age <- ggplot(pca_df_age, aes(PC2, PC4, color = age)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC2 vs PC4", color = "Age")

# PC3 vs PC4
p34_age <- ggplot(pca_df_age, aes(PC3, PC4, color = age)) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(title = "PC3 vs PC4", color = "Age")

# Visi kopā
(p12_age | p13_age | p14_age) /
  (p23_age | p24_age | p34_age) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# ==========================================
# PCA PIRMS NORMALIZĀCIJAS
# GP0071 SAMPLE NOSAUKUMI
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)

# RAW counts log-transformācija
raw_log <- log2(counts(dds) + 1)

# Top 500 visvariablākie gēni
rv_raw <- rowVars(raw_log)
select_raw <- order(rv_raw, decreasing = TRUE)[1:500]

# PCA
pca_raw <- prcomp(
  t(raw_log[select_raw, ]),
  scale. = FALSE
)

# PCA koordinātes + sample + patient
pca_check <- data.frame(
  PC1 = pca_raw$x[,1],
  PC2 = pca_raw$x[,2],
  PC3 = pca_raw$x[,3],
  PC4 = pca_raw$x[,4],
  sample = meta$sample,
  patient = meta$patient
)

# Iezīmējam tikai GP0071
pca_check$highlight <- ifelse(
  pca_check$patient == "GP0071",
  "GP0071",
  "Other"
)

#Pārbaudu konkrētu pacientu
# PC1 vs PC2
ggplot(pca_check, aes(PC1, PC2, color = highlight)) +
  geom_point(size = 2.5) +
  geom_text(
    data = subset(pca_check, patient == "GP0071"),
    aes(label = sample),
    vjust = -0.7,
    size = 3,
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    title = "GP0071 samples — PC1 vs PC2",
    color = ""
  )

#Pārbaudu konkrētu pacientu
ggplot(pca_check, aes(PC1, PC3, color = highlight)) +
  geom_point(size = 2.5) +
  geom_text(
    data = subset(pca_check, patient == "GP0071"),
    aes(label = sample),
    vjust = -0.7,
    size = 3,
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    title = "GP0071 samples — PC1 vs PC3",
    color = ""
  )

#Pārbaudu konkrētu pacientu
ggplot(pca_check, aes(PC1, PC4, color = highlight)) +
  geom_point(size = 2.5) +
  geom_text(
    data = subset(pca_check, patient == "GP0071"),
    aes(label = sample),
    vjust = -0.7,
    size = 3,
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    title = "GP0071 samples — PC1 vs PC4",
    color = ""
  )

# =========================
#NORMALIZĀCIJA
# =========================

library(DESeq2)

# Pirms normalizācijas
raw_total <- colSums(counts(dds))

# Pēc normalizācijas
norm_total <- colSums(counts(dds, normalized = TRUE))

# Pārbaude
head(raw_total)
head(norm_total)
# Izveidojam DESeq2 objektu no gene counts + metadata
dds <- DESeqDataSetFromMatrix(
  countData = gene_counts,
  colData = meta,
  design = ~ 1
)

# Maigā filtrēšana:
# izmetam tikai gēnus, kuriem visos 72 sample count = 0
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]

# Aprēķinām sequencing depth / library size korekcijas faktorus
dds <- estimateSizeFactors(dds)

# Apskatām korekcijas faktorus
sizeFactors(dds)

# Iegūstam normalizētos gene counts
norm_counts <- counts(dds, normalized = TRUE)

# Pārbaudām izmēru
dim(norm_counts)

#Comparing normalized counts vs not normalized
# Kopējais gene count katram sample PIRMS normalizācijas
raw_total <- colSums(counts(dds))

# Kopējais gene count katram sample PĒC normalizācijas
norm_total <- colSums(counts(dds, normalized = TRUE))

# Saliekam tabulā
norm_check <- data.frame(
  sample = colnames(dds),
  raw = raw_total,
  normalized = norm_total
)

# Pārveidojam ērtāk grafikam
library(tidyr)
library(ggplot2)

norm_check_long <- pivot_longer(
  norm_check,
  cols = c(raw, normalized),
  names_to = "data",
  values_to = "total_counts"
)

# Grafiks
ggplot(
  norm_check_long,
  aes(x = reorder(sample, total_counts),
      y = total_counts,
      color = data,
      group = data)
) +
  geom_point(size = 2) +
  geom_line() +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Sample",
    y = "Total gene counts",
    color = "",
    title = "Library size before and after normalization"
  )


#Checking how normalisation worked
# Size factors
summary(sizeFactors(dds))

# Raw vs normalized library sizes
raw_lib <- colSums(counts(dds))
norm_lib <- colSums(counts(dds, normalized = TRUE))

summary(raw_lib)
summary(norm_lib)

# Relatīvā izkliede
sd(raw_lib) / mean(raw_lib)
sd(norm_lib) / mean(norm_lib)




#PCA PĒC NORMALIZĀCIJAS

# ==========================================
# PCA PĒC NORMALIZĀCIJAS
# DIĒTA + FĀZE
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)
library(matrixStats)

# 1. Normalizēto counts log-transformācija
norm_log <- log2(norm_counts + 1)

# 2. Aprēķina variāciju katram gēnam
rv_norm <- rowVars(norm_log)

# 3. Izvēlas 500 visvariablākos gēnus
select_norm <- order(
  rv_norm,
  decreasing = TRUE
)[1:500]

# 4. PCA tikai ar top 500 gēniem
pca_norm <- prcomp(
  t(norm_log[select_norm, ]),
  scale. = FALSE
)

# 5. PCA koordinātes + metadata
pca_df_norm <- data.frame(
  PC1 = pca_norm$x[,1],
  PC2 = pca_norm$x[,2],
  PC3 = pca_norm$x[,3],
  PC4 = pca_norm$x[,4],
  diet = meta$diet,
  phase = meta$phase
)

# 6. Apvienojam diet + phase vienā grupā
pca_df_norm$diet_phase <- interaction(
  pca_df_norm$diet,
  pca_df_norm$phase,
  sep = " "
)

# 7. PC1 vs PC2
p12_norm <- ggplot(
  pca_df_norm,
  aes(PC1, PC2, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC2",
    color = "Diet + phase"
  )

# 8. PC1 vs PC3
p13_norm <- ggplot(
  pca_df_norm,
  aes(PC1, PC3, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC3",
    color = "Diet + phase"
  )

# 9. PC1 vs PC4
p14_norm <- ggplot(
  pca_df_norm,
  aes(PC1, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC4",
    color = "Diet + phase"
  )

# 10. PC2 vs PC3
p23_norm <- ggplot(
  pca_df_norm,
  aes(PC2, PC3, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC3",
    color = "Diet + phase"
  )

# 11. PC2 vs PC4
p24_norm <- ggplot(
  pca_df_norm,
  aes(PC2, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC4",
    color = "Diet + phase"
  )

# 12. PC3 vs PC4
p34_norm <- ggplot(
  pca_df_norm,
  aes(PC3, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC3 vs PC4",
    color = "Diet + phase"
  )

# 13. Saliek visus grafikus kopā
(p12_norm | p13_norm | p14_norm) /
  (p23_norm | p24_norm | p34_norm) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# ==========================================
# PCA PĒC NORMALIZĀCIJAS
# TIKAI SIEVIETES
# DIĒTA + FĀZE
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)
library(matrixStats)

# 1. Atlasām tikai sieviešu sample
female_samples <- meta$sex == "F"

# 2. Normalizētie counts tikai sievietēm
female_counts_norm <- norm_counts[, female_samples]

# 3. Log transformācija
female_log_norm <- log2(female_counts_norm + 1)

# 4. Atrodam top 500 visvariablākos gēnus
rv_female_norm <- rowVars(female_log_norm)

select_female_norm <- order(
  rv_female_norm,
  decreasing = TRUE
)[1:500]

# 5. PCA tikai sievietēm
pca_female_norm <- prcomp(
  t(female_log_norm[select_female_norm, ]),
  scale. = FALSE
)

# 6. Metadata tikai sievietēm
meta_female <- meta[female_samples, ]

# 7. PCA koordinātes + diet + phase
pca_female_df_norm <- data.frame(
  PC1 = pca_female_norm$x[,1],
  PC2 = pca_female_norm$x[,2],
  PC3 = pca_female_norm$x[,3],
  PC4 = pca_female_norm$x[,4],
  diet = meta_female$diet,
  phase = meta_female$phase
)

# 8. Apvienojam diet + phase vienā grupā
pca_female_df_norm$diet_phase <- interaction(
  pca_female_df_norm$diet,
  pca_female_df_norm$phase,
  sep = " "
)

# 9. PC1 vs PC2
p12_f_norm <- ggplot(
  pca_female_df_norm,
  aes(PC1, PC2, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC2",
    color = "Diet + phase"
  )

# 10. PC1 vs PC3
p13_f_norm <- ggplot(
  pca_female_df_norm,
  aes(PC1, PC3, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC3",
    color = "Diet + phase"
  )

# 11. PC1 vs PC4
p14_f_norm <- ggplot(
  pca_female_df_norm,
  aes(PC1, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC4",
    color = "Diet + phase"
  )

# 12. PC2 vs PC3
p23_f_norm <- ggplot(
  pca_female_df_norm,
  aes(PC2, PC3, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC3",
    color = "Diet + phase"
  )

# 13. PC2 vs PC4
p24_f_norm <- ggplot(
  pca_female_df_norm,
  aes(PC2, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC4",
    color = "Diet + phase"
  )

# 14. PC3 vs PC4
p34_f_norm <- ggplot(
  pca_female_df_norm,
  aes(PC3, PC4, color = diet_phase)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC3 vs PC4",
    color = "Diet + phase"
  )

# 15. Visi grafiki kopā
(p12_f_norm | p13_f_norm | p14_f_norm) /
  (p23_f_norm | p24_f_norm | p34_f_norm) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# ==========================================
# PCA PĒC NORMALIZĀCIJAS
# VECUMS (VISI DZIMUMI)
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)
library(matrixStats)

# 1. Normalizēto counts log-transformācija
norm_log <- log2(norm_counts + 1)

# 2. Variācija katram gēnam
rv_norm <- rowVars(norm_log)

# 3. Top 500 visvariablākie gēni
select_norm <- order(
  rv_norm,
  decreasing = TRUE
)[1:500]

# 4. PCA
pca_norm <- prcomp(
  t(norm_log[select_norm, ]),
  scale. = FALSE
)

# 5. PCA koordinātes + vecums
pca_df_age_norm <- data.frame(
  PC1 = pca_norm$x[,1],
  PC2 = pca_norm$x[,2],
  PC3 = pca_norm$x[,3],
  PC4 = pca_norm$x[,4],
  age = as.numeric(meta$age)
)

# 6. PC1 vs PC2
p12_age_norm <- ggplot(
  pca_df_age_norm,
  aes(PC1, PC2, color = age)
) +
  geom_point(size = 2.5) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(
    title = "PC1 vs PC2",
    color = "Age"
  )

# 7. PC1 vs PC3
p13_age_norm <- ggplot(
  pca_df_age_norm,
  aes(PC1, PC3, color = age)
) +
  geom_point(size = 2.5) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(
    title = "PC1 vs PC3",
    color = "Age"
  )

# 8. PC1 vs PC4
p14_age_norm <- ggplot(
  pca_df_age_norm,
  aes(PC1, PC4, color = age)
) +
  geom_point(size = 2.5) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(
    title = "PC1 vs PC4",
    color = "Age"
  )

# 9. PC2 vs PC3
p23_age_norm <- ggplot(
  pca_df_age_norm,
  aes(PC2, PC3, color = age)
) +
  geom_point(size = 2.5) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(
    title = "PC2 vs PC3",
    color = "Age"
  )

# 10. PC2 vs PC4
p24_age_norm <- ggplot(
  pca_df_age_norm,
  aes(PC2, PC4, color = age)
) +
  geom_point(size = 2.5) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(
    title = "PC2 vs PC4",
    color = "Age"
  )

# 11. PC3 vs PC4
p34_age_norm <- ggplot(
  pca_df_age_norm,
  aes(PC3, PC4, color = age)
) +
  geom_point(size = 2.5) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(
    title = "PC3 vs PC4",
    color = "Age"
  )

# 12. Visi grafiki kopā
(p12_age_norm | p13_age_norm | p14_age_norm) /
  (p23_age_norm | p24_age_norm | p34_age_norm) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ==========================================
# PCA PĒC NORMALIZĀCIJAS
# KATRAM PACIENTAM SAVA KRĀSA
# KATRAM PARAUGAM SAVS SAMPLE ID
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)

# 1. PCA koordinātes + patient + sample
pca_patient_df <- data.frame(
  PC1 = pca_norm$x[,1],
  PC2 = pca_norm$x[,2],
  PC3 = pca_norm$x[,3],
  PC4 = pca_norm$x[,4],
  sample = meta$sample,
  patient = factor(meta$patient)
)

# 2. Funkcija PCA grafiku veidošanai
plot_patient_pca <- function(x, y) {
  
  ggplot(
    pca_patient_df,
    aes(
      x = .data[[x]],
      y = .data[[y]],
      color = patient
    )
  ) +
    geom_point(size = 2.5) +
    geom_text(
      aes(label = sample),
      vjust = -0.7,
      size = 2.5,
      show.legend = FALSE
    ) +
    theme_classic() +
    labs(
      title = paste(x, "vs", y),
      color = "Patient"
    )
}

# 3. Seši PCA grafiki
p12_patient <- plot_patient_pca("PC1", "PC2")
p13_patient <- plot_patient_pca("PC1", "PC3")
p14_patient <- plot_patient_pca("PC1", "PC4")

p23_patient <- plot_patient_pca("PC2", "PC3")
p24_patient <- plot_patient_pca("PC2", "PC4")
p34_patient <- plot_patient_pca("PC3", "PC4")

# 4. Visi grafiki kopā
(p12_patient | p13_patient | p14_patient) /
  (p23_patient | p24_patient | p34_patient) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# ==========================================
# PCA PĒC NORMALIZĀCIJAS
# DZIMUMS (VISI PARAUGI)
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(ggplot2)
library(patchwork)

# 1. PCA koordinātes + dzimums
pca_df_sex_norm <- data.frame(
  PC1 = pca_norm$x[,1],
  PC2 = pca_norm$x[,2],
  PC3 = pca_norm$x[,3],
  PC4 = pca_norm$x[,4],
  sex = factor(meta$sex)
)

# 2. PC1 vs PC2
p12_sex_norm <- ggplot(
  pca_df_sex_norm,
  aes(PC1, PC2, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC2",
    color = "Sex"
  )

# 3. PC1 vs PC3
p13_sex_norm <- ggplot(
  pca_df_sex_norm,
  aes(PC1, PC3, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC3",
    color = "Sex"
  )

# 4. PC1 vs PC4
p14_sex_norm <- ggplot(
  pca_df_sex_norm,
  aes(PC1, PC4, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC1 vs PC4",
    color = "Sex"
  )

# 5. PC2 vs PC3
p23_sex_norm <- ggplot(
  pca_df_sex_norm,
  aes(PC2, PC3, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC3",
    color = "Sex"
  )

# 6. PC2 vs PC4
p24_sex_norm <- ggplot(
  pca_df_sex_norm,
  aes(PC2, PC4, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC2 vs PC4",
    color = "Sex"
  )

# 7. PC3 vs PC4
p34_sex_norm <- ggplot(
  pca_df_sex_norm,
  aes(PC3, PC4, color = sex)
) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    title = "PC3 vs PC4",
    color = "Sex"
  )

# 8. Visi grafiki kopā
(p12_sex_norm | p13_sex_norm | p14_sex_norm) /
  (p23_sex_norm | p24_sex_norm | p34_sex_norm) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ==========================================
# PCA PĒC NORMALIZĀCIJAS — BATCH EFFECT
# TOP 500 VISVARIABLĀKIE GĒNI
# ==========================================

library(readxl)
library(ggplot2)
library(patchwork)

# 1. Ielādējam jauno metadatu tabulu
meta <- read_excel("Metatabula_baka.xlsx")

# 2. Sakārtojam metadata pēc count matrix paraugu secības
meta <- meta[match(colnames(norm_counts), meta$sample), ]

stopifnot(identical(colnames(norm_counts), meta$sample))

# 3. PCA koordinātes + batch
pca_df_batch <- data.frame(
  PC1 = pca_norm$x[,1],
  PC2 = pca_norm$x[,2],
  PC3 = pca_norm$x[,3],
  PC4 = pca_norm$x[,4],
  batch = factor(meta$batch)
)

# 4. Funkcija PCA grafikiem
plot_batch_pca <- function(x, y) {
  
  ggplot(
    pca_df_batch,
    aes(
      x = .data[[x]],
      y = .data[[y]],
      color = batch
    )
  ) +
    geom_point(size = 2.5) +
    theme_classic() +
    labs(
      title = paste(x, "vs", y),
      color = "Batch"
    )
}

# 5. Seši PCA grafiki
p12_batch <- plot_batch_pca("PC1", "PC2")
p13_batch <- plot_batch_pca("PC1", "PC3")
p14_batch <- plot_batch_pca("PC1", "PC4")

p23_batch <- plot_batch_pca("PC2", "PC3")
p24_batch <- plot_batch_pca("PC2", "PC4")
p34_batch <- plot_batch_pca("PC3", "PC4")

# 6. Visi grafiki kopā
(p12_batch | p13_batch | p14_batch) /
  (p23_batch | p24_batch | p34_batch) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")