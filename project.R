#set the working directory
#setwd("~/Documents/MBBC/AOm/project") 


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

install.packages("readxl")

# Load DESeq2
library(DESeq2)
library(readxl)

message("DESeq2 loaded successfully!")

#Load meta data 
meta <- read_excel("Metatabula_baka.xlsx")

#Load gene count matrix
counts <- read.delim(
  "Gene.counts.strand2.tsv",
  check.names = FALSE
)

#Check the dimensions of each matrix
dim(meta)
dim(counts)


# Inspect metadata
head(meta)
colnames(meta)
str(meta)

# Inspect count matrix
head(counts)
colnames(counts)
str(counts)


head(meta[[1]])
head(colnames(counts)[-1])

# remove ".R.L" from the gene count column names
colnames(counts) <- sub("\\.R\\.L$", "", colnames(counts))

# creates an array with the names of the samples (already w the matching names with the metadata)
sample_ids <- colnames(counts)[-1]

# re order the sample ids in metada to be the same as in the data matrix (same order, ex:GP0017_1, GP0017_2... instead of GP0017_4, GP0017_2...)
# no need because it was already in the same order
#meta <- meta[match(sample_ids, meta[[1]]), , drop = FALSE]


#check if it is ordered, need to return TRUE
identical(
  colnames(counts)[-1],
  as.character(meta[[1]])
)


#drop all the samples in which 'diet order' column have LC-IF (do not start with the Intermittent fasting diet)

#show us a table with the diet order and diet values
table(meta[["diet order"]], useNA = "ifany")

table(meta[["diet"]], useNA = "ifany")

#create a condition w the values to keep
#we only want individuals that went through the intermittent fasting first and only the intermittent diet values
keep_rows <- meta[["diet order"]] != "LC-IF" & meta[["diet"]] != "LC"
keep_rows

#create a new metadata file containing only the samples we want (filtering with the keep_rows condition)
filtermeta <- meta[keep_rows, , drop = FALSE]

#create a new data file filtering it with the metadata 
#as.character converts the metadata ids into text so we can compare wit the counts columns
filteredcounts <- counts[, c("Gene", as.character(filtermeta[[1]])), drop = FALSE]

#now we have all the samples we're going to work with in filteredcounts


#now to remove all the genes that have 0 expression before and after (both) the diet we create another filter

#first create an array with the sum of all genes in every sample (the sum of every row)
gene_totals <- rowSums(filteredcounts[, -1])
#check how many genes were checked, -4 that were the first 4 lines because they had information about the N_unmapped, N_multimapped etc
length(gene_totals)-4

#check the dimension of filtered counts before
dim(filteredcounts)

#now we apply the condition that the gene_totals should have a value > 0
filteredcounts <- filteredcounts[gene_totals > 0, ]

#check the dimension of filtered counts after
dim(filteredcounts)

#it went from 58306 to 30862

write.table(
  filteredcounts,
  file = "filtered_gene_counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

install.packages("writexl")
library(writexl)

write_xlsx(
  filtermeta,
  "filtered_metadata.xlsx"
)


