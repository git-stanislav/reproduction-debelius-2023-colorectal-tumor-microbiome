from pathlib import Path
import sys

import pandas as pd


# Input files
fastq_path = Path(sys.argv[1])
metadata_path = Path(sys.argv[2])
output_path = Path(sys.argv[3])

# Read all columns as strings to preserve identifiers and values.
fastq = pd.read_csv(fastq_path, sep="\t", dtype=str)
metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)

# Check the required columns.
required_fastq = {"sample_alias"}
required_metadata = {"host_subject_id", "tissue_type"}

if not required_fastq.issubset(fastq.columns):
    raise ValueError(f"Missing columns in FASTQ table: {required_fastq - set(fastq.columns)}")

if not required_metadata.issubset(metadata.columns):
    raise ValueError(
        f"Missing columns in metadata table: "
        f"{required_metadata - set(metadata.columns)}"
    )

# Parse sample aliases such as:
# short_indv145_tumor
# short_indv145_normal
parsed = fastq["sample_alias"].str.extract(
    r"^[^_]+_indv(?P<host_subject_id>\d+)_(?P<tissue_alias>tumou?r|normal)$",
    flags=2,  # case-insensitive
)

# Do not proceed if any alias has an unexpected format.
invalid = parsed["host_subject_id"].isna() | parsed["tissue_alias"].isna()

if invalid.any():
    print("Could not parse these sample_alias values:")
    print(fastq.loc[invalid, "sample_alias"].to_string(index=False))
    raise SystemExit("Fix the aliases or update the regular expression.")

# Convert the alias tissue suffix to the vocabulary used by the metadata.
tissue_map = {
    "tumor": "tumour tissue",
    "tumour": "tumour tissue",
    "normal": "normal tissue",
}

fastq["host_subject_id"] = parsed["host_subject_id"]
fastq["tissue_type"] = parsed["tissue_alias"].str.lower().map(tissue_map)

# Normalize the metadata join keys.
metadata["host_subject_id"] = metadata["host_subject_id"].str.strip()
metadata["tissue_type"] = metadata["tissue_type"].str.strip().str.lower()

# Each patient/tissue combination should identify at most one metadata row.
keys = ["host_subject_id", "tissue_type"]
duplicates = metadata.duplicated(subset=keys, keep=False)

if duplicates.any():
    print("Duplicate keys in paired_metadata.tsv:")
    print(metadata.loc[duplicates, keys].sort_values(keys).to_string(index=False))
    raise SystemExit("Cannot safely merge: metadata keys are not unique.")

# Keep every row in the FASTQ table, even if no metadata match is found.
merged = fastq.merge(
    metadata,
    on=keys,
    how="left",
    validate="many_to_one",
    indicator=True,
    suffixes=("", "_paired"),
)

# Report the merge quality.
print(f"FASTQ rows: {len(fastq)}")
print("Merge status:")
print(merged["_merge"].value_counts(dropna=False).to_string())

unmatched = merged.loc[
    merged["_merge"] == "left_only",
    ["sample_alias", "host_subject_id", "tissue_type"],
]

if not unmatched.empty:
    print("\nFASTQ samples without a metadata match:")
    print(unmatched.to_string(index=False))

# Sample ID column must be the first.
merged.insert(0, "sample-id", merged.pop("sample-id"))

# Keep the match status as an audit column.
merged = merged.rename(columns={"_merge": "merge_status"})

# In QIIME 2 'sample_name' conflicts with a name reserved for the ID column header.
merged = merged.rename(columns={"sample_name": "sample_label"})

output_path.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(output_path, sep="\t", index=False)

print(f"\nSaved merged metadata to: {output_path}")


# Create QIIME 2 manifest.
manifest_path = Path(sys.argv[4])
raw_dir = Path("data/raw/fastq").resolve()

# Detect the FASTQ run accession column.
accession_col = next(
    (col for col in ("run_accession", "accession") if col in merged.columns),
    None,
)

if accession_col is None:
    raise ValueError(
        "Cannot find the run accession column. "
        f"Available columns: {list(merged.columns)}"
    )

accessions = merged[accession_col].str.strip()
sample_ids = merged["sample-id"].str.strip()

if accessions.isna().any() or accessions.eq("").any():
    raise ValueError("Some samples have no run accession.")

if sample_ids.duplicated().any():
    duplicates = sample_ids[sample_ids.duplicated(keep=False)]
    raise ValueError(f"Duplicate sample IDs:\n{duplicates.to_string(index=False)}")

manifest = pd.DataFrame({
    "sample-id": sample_ids,
    "forward-absolute-filepath": accessions.map(
        lambda acc: str(raw_dir / f"{acc}_1.fastq.gz")
    ),
    "reverse-absolute-filepath": accessions.map(
        lambda acc: str(raw_dir / f"{acc}_2.fastq.gz")
    ),
})

# Check that every FASTQ pair exists.
fastq_paths = (
    manifest["forward-absolute-filepath"].tolist()
    + manifest["reverse-absolute-filepath"].tolist()
)
missing = [path for path in fastq_paths if not Path(path).is_file()]

if missing:
    raise FileNotFoundError(
        "Some FASTQ files are missing:\n" + "\n".join(missing)
    )

manifest.to_csv(manifest_path, sep="\t", index=False)

print(f"Saved QIIME 2 manifest to: {manifest_path}")

