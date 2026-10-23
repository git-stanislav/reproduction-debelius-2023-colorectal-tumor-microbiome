#!/usr/bin/env Rscript

library(vegan)

# 1. Check command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: Rscript permanova.R <distance-matrix.tsv>")
}

distance_file <- args[1]

# 2. Load metadata
metadata <- read.delim(
    "data/metadata/sample-metadata.tsv",
    row.names = 1,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
)

# 3. Load the distance matrix
distance_matrix <- read.delim(
    distance_file,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
)

distance_matrix <- as.matrix(distance_matrix)

# 4. Check for the required columns
required_columns <- c("host_subject_id", "tissue_type", "long_survival")

if (!all(required_columns %in% colnames(metadata))) {
    stop("Metadata must contain host_subject_id, tissue_type and long_survival columns.")
}

# 5. Check sample correspondence
if (!setequal(rownames(distance_matrix), colnames(distance_matrix))) {
    stop("Row and column sample IDs in the distance matrix do not match.")
}


print(rownames(distance_matrix))
print(rownames(metadata))

if (!setequal(rownames(distance_matrix), rownames(metadata))) {
    stop("Sample IDs in the distance matrix and metadata do not match.")
}

# 6. Align metadata with the order of the distance matrix
metadata <- metadata[
    rownames(distance_matrix),
    ,
    drop = FALSE
]

# 7. Sort by patient and tissue type
ord <- order(metadata$host_subject_id, metadata$tissue_type)

metadata <- metadata[ord, , drop = FALSE]

# The rows and columns of the matrix must retain the same order
sample_order <- rownames(metadata)

distance_matrix <- distance_matrix[
    sample_order,
    sample_order,
    drop = FALSE
]

# 8. Check pair integrity
if (anyNA(metadata$host_subject_id) ||
    anyNA(metadata$tissue_type) ||
    anyNA(metadata$long_survival)) {
    stop("host_subject_id, tissue_type and long_survival must not contain missing values.")
}

if (any(table(metadata$host_subject_id) != 2)) {
    stop("Each patient must have exactly two samples.")
}

tissues <- split(
    as.character(metadata$tissue_type),
    metadata$host_subject_id
)

valid_tissues <- vapply(
    tissues,
    function(x) setequal(x, c("normal tissue", "tumour tissue")),
    logical(1)
)

if (!all(valid_tissues)) {
    stop("Each patient must have one normal and one tumor sample.")
}

survivals <- split(
    as.character(metadata$long_survival),
    metadata$host_subject_id
)

valid_survival <- vapply(
    survivals,
    function(x) length(unique(x)) == 1,
    logical(1)
)

if (!all(valid_survival)) {
    stop("long_survival must be identical for both samples of each patient.")
}

# 9. Identify the indices of each patient's samples
patients <- unique(metadata$host_subject_id)
n_patients <- length(patients)

pair_indices <- lapply(
    patients,
    function(patient) which(metadata$host_subject_id == patient)
)

# 10. Set the random number generator
set.seed(42)

nperm <- 999

# 11. Generate permutations of entire patient pairs
permutations <- t(replicate(nperm, {
    patient_order <- sample(seq_len(n_patients))

    unlist(
        pair_indices[patient_order],
        use.names = FALSE
    )
}))

# 12. Convert the distance matrix to a dist object
distance_matrix <- as.dist(distance_matrix)

# 13. Perform PERMANOVA
result <- adonis2(
    distance_matrix ~ long_survival,
    data = metadata,
    permutations = permutations
)

# 14. Print the results
print(result)

