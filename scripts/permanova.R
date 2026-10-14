#!/usr/bin/env Rscript

library(vegan)
library(permute)

# Check args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: Rscript permanova.R <distance-matrix.tsv>")
}

distance_file <- args[1]

# Metadata
metadata <- read.delim(
    "data/metadata/sample-metadata.tsv",
    row.names = 1,
    sep = "\t",
    comment.char = "#"
)

# Distance matrix
distance_matrix <- read.delim(
    distance_file,
    row.names = 1,
    sep = "\t"
)

distance_matrix <- as.dist(distance_matrix)

# Restricted permutations within PatientID
permutation_control <- how(
    within = Within(type = "free"),
    plots = Plots(
        strata = metadata$PatientID,
        type = "none"
    ),
    nperm = 999
)

# PERMANOVA
result <- adonis2(
    distance_matrix ~ Tissue,
    data = metadata,
    permutations = permutation_control
)

model <- result[1, ]

cat("R2\tF\tp\n")
cat(
    model$R2, "\t",
    model$F, "\t",
    model$`Pr(>F)`,
    "\n",
    sep = ""
)
