from pathlib import Path

import pandas as pd


TMPDIR = Path.home() / "qiime2-tmp"
TMPDIR.mkdir(parents=True, exist_ok=True)

shell.prefix(f"export TMPDIR='{TMPDIR}'; ")


def get_all_fastq_files(wildcards):
    checkpoint_output = checkpoints.download_fastq_urls.get().output[0]
    df = pd.read_csv(checkpoint_output, sep="\t")

    accessions = df["run_accession"].dropna().astype(str).drop_duplicates().tolist()

    return [
        f"data/raw/fastq/{accession}_{read}.fastq.gz"
        for accession in accessions
        for read in ("1", "2")
    ]


def get_fastq_url(wildcards):
    checkpoint_output = checkpoints.download_fastq_urls.get().output[0]
    df = pd.read_csv(checkpoint_output, sep="\t")

    rows = df.loc[
        df["run_accession"].astype(str) == wildcards.accession,
        "fastq_ftp",
    ]

    if rows.empty:
        raise ValueError(
            f"No ENA FASTQ entry found for run accession {wildcards.accession}"
        )

    urls = str(rows.iloc[0]).split(";")
    read_index = int(wildcards.read) - 1

    if read_index >= len(urls) or not urls[read_index].strip():
        raise ValueError(
            f"No FASTQ URL for read {wildcards.read} of run {wildcards.accession}"
        )

    url = urls[read_index].strip()

    if not url.startswith(("ftp://", "http://", "https://")):
        url = "https://" + url

    return url


BETA_SCORES = [
    "braycurtis",
    "jaccard",
    "aitchison",
    "weighted_unifrac",
    "unweighted_unifrac",
]


rule all:
    input:
        "results/visualizations/jaccard-tissue-significance.qzv",
        "results/visualizations/braycurtis-tissue-significance.qzv",
        "results/visualizations/aitchison-tissue-significance.qzv",
        "results/visualizations/unweighted_unifrac-tissue-significance.qzv",
        "results/visualizations/weighted_unifrac-tissue-significance.qzv",
        "results/visualizations/jaccard-emperor.qzv",
        "results/visualizations/braycurtis-emperor.qzv",
        "results/visualizations/aitchison-emperor.qzv",
        "results/visualizations/unweighted_unifrac-emperor.qzv",
        "results/visualizations/weighted_unifrac-emperor.qzv",
        "report/latex/tables/permanova-results.tex",


rule permanova_to_latex:
    input:
        expand(
            "results/statistics/permanova/{score}-matrix-permanova.txt",
            score=BETA_SCORES,
        ),
    output:
        "report/latex/tables/permanova-results.tex",
    shell:
        """
        python3 scripts/permanova-to-latex.py {input} {output}
        """


rule permanova_r:
    input:
        "results/dissimilarity-matrices/{score}-matrix.tsv",
    output:
        "results/statistics/permanova/{score}-matrix-permanova.txt",
    shell:
        """
        mkdir -p results/statistics/permanova

        Rscript scripts/permanova.R {input} > {output}
        """


rule tissue_significance:
    input:
        score="results/diversity/{file}-diversity.qza",
        sample_metadata="data/metadata/sample-metadata.tsv",
    output:
        "results/visualizations/{file}-tissue-significance.qzv",
    shell:
        """
        qiime diversity beta-group-significance \
        --i-distance-matrix {input.score} \
        --m-metadata-file {input.sample_metadata} \
        --m-metadata-column Tissue \
        --p-method permanova \
        --o-visualization {output}
        """


rule emperor:
    input:
        pcoa="results/diversity/{file}-pcoa.qza",
        sample_metadata="data/metadata/sample-metadata.tsv",
    output:
        "results/visualizations/{file}-emperor.qzv",
    shell:
        """
        mkdir -p results/visualizations

        qiime emperor plot \
            --i-pcoa {input.pcoa} \
            --m-metadata-file {input.sample_metadata} \
            --o-visualization {output}
        """


rule pcoa:
    input:
        "results/diversity/{file}-diversity.qza",
    output:
        "results/diversity/{file}-pcoa.qza",
    shell:
        """
        qiime diversity pcoa \
        --i-distance-matrix {input} \
        --o-pcoa {output}
        """


rule export_diversity_scores:
    input:
        "results/diversity/{score}-diversity.qza",
    output:
        "results/dissimilarity-matrices/{score}-matrix.tsv",
    shell:
        """
        mkdir -p results/dissimilarity-matrices

        tmpdir=$(mktemp -d)

        qiime tools export \
            --input-path {input} \
            --output-path "$tmpdir"

        mv "$tmpdir/distance-matrix.tsv" {output}

        rmdir "$tmpdir"
        """


rule evaluate_diversity_phylogeny:
    input:
        feature_table_2500="work/qiime2/feature-table-2500.qza",
        rooted_tree="work/qiime2/sepp-tree.qza",
    output:
        unweighted_unifrac="results/diversity/unweighted_unifrac-diversity.qza",
        weighted_unifrac="results/diversity/weighted_unifrac-diversity.qza",
    shell:
        """
        qiime diversity beta-phylogenetic \
        --i-table {input.feature_table_2500} \
        --i-phylogeny {input.rooted_tree} \
        --p-metric unweighted_unifrac \
        --o-distance-matrix {output.unweighted_unifrac}
    
        qiime diversity beta-phylogenetic \
        --i-table {input.feature_table_2500} \
        --i-phylogeny {input.rooted_tree} \
        --p-metric weighted_unifrac \
        --o-distance-matrix {output.weighted_unifrac}
        """


rule evaluate_diversity:
    input:
        feature_table="work/qiime2/sepp-filtered-table.qza",
        feature_table_2500="work/qiime2/feature-table-2500.qza",
    output:
        bray_curtis="results/diversity/braycurtis-diversity.qza",
        aitchison="results/diversity/aitchison-diversity.qza",
        jaccard="results/diversity/jaccard-diversity.qza",
    shell:
        """
        mkdir -p results/diversity

        qiime diversity beta \
            --i-table {input.feature_table_2500} \
            --p-metric jaccard \
            --o-distance-matrix {output.jaccard}

        qiime diversity beta \
            --i-table {input.feature_table_2500} \
            --p-metric braycurtis \
            --o-distance-matrix {output.bray_curtis}

        qiime diversity beta \
            --i-table {input.feature_table} \
            --p-metric aitchison \
            --o-distance-matrix {output.aitchison}
        """


rule rarefy:
    input:
        "work/qiime2/sepp-filtered-table.qza",
    output:
        "work/qiime2/feature-table-2500.qza",
    shell:
        """
        qiime feature-table rarefy \
            --i-table {input} \
            --p-sampling-depth 2500 \
            --o-rarefied-table {output}
        """


rule filter_sepp_features:
    input:
        table="work/qiime2/taxonomy-filtered-table.qza",
        tree="work/qiime2/sepp-tree.qza",
    output:
        filtered_table="work/qiime2/sepp-filtered-table.qza",
        removed_table="work/qiime2/sepp-removed-table.qza",
    shell:
        """
        qiime fragment-insertion filter-features \
            --i-table {input.table} \
            --i-tree {input.tree} \
            --o-filtered-table {output.filtered_table} \
            --o-removed-table {output.removed_table}
        """


rule phylogeny:
    input:
        sequences="work/qiime2/asv-sequences.qza",
        sepp_reference="resources/sepp/sepp-refs-silva-128.qza",
    output:
        tree="work/qiime2/sepp-tree.qza",
        placements="work/qiime2/sepp-placements.qza",
    log:
        "logs/phylogeny.log",
    shell:
        """
        qiime fragment-insertion sepp \
            --i-representative-sequences {input.sequences} \
            --i-reference-database {input.sepp_reference} \
            --o-tree {output.tree} \
            --o-placements {output.placements} \
            --p-threads 8 \
            > {log} 2>&1
        """


rule get_sepp_silva:
    output:
        "resources/sepp/sepp-refs-silva-128.qza",
    shell:
        """
        mkdir -p resources/sepp

        wget \
        -O {output} \
        "https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-silva-128.qza"
        """


rule feature_table_summarize:
    input:
        feature_table="work/qiime2/taxonomy-filtered-table.qza",
        sample_metadata="data/metadata/sample-metadata.tsv",
    output:
        dir=directory("feature-table-summarize-output-dir"),
        summary="work/qiime2/visualizations/feature-table-summary.qzv",
    shell:
        """
        qiime feature-table summarize \
        --i-table {input.feature_table} \
        --m-metadata-file {input.sample_metadata} \
        --output-dir {output.dir} \
        --o-summary {output.summary}
        """


rule taxonomy_visualization:
    input:
        taxonomy="work/qiime2/taxonomy.qza",
        feature_table="work/qiime2/taxonomy-filtered-table.qza",
        sample_metadata="data/metadata/sample-metadata.tsv",
    output:
        taxonomy="work/qiime2/visualizations/taxonomy.qzv",
        taxa_barplot="work/qiime2/visualizations/taxa-barplot.qzv",
    shell:
        """
        qiime metadata tabulate \
            --m-input-file {input.taxonomy} \
            --o-visualization {output.taxonomy}
            
        qiime taxa barplot \
            --i-table {input.feature_table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.sample_metadata} \
            --o-visualization {output.taxa_barplot}
        """


rule filter_taxonomy:
    input:
        table="work/qiime2/feature-table.qza",
        taxonomy="work/qiime2/taxonomy.qza",
    output:
        "work/qiime2/taxonomy-filtered-table.qza",
    shell:
        """
        qiime taxa filter-table \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --p-include "p__" \
            --p-mode contains \
            --o-filtered-table {output}
        """


rule classify_asvs:
    input:
        classifier="resources/reference/silva-128-341-805-classifier.qza",
        sequences="work/qiime2/asv-sequences.qza",
    output:
        "work/qiime2/taxonomy.qza",
    log:
        "logs/classify-asvs.log",
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.sequences} \
            --o-classification {output} \
            > {log} 2>&1
        """


rule train_classifier:
    input:
        sequences="resources/reference/silva-128-341-805-derep.qza",
        taxonomy="resources/reference/silva-128-tax-derep.qza",
    output:
        "resources/reference/silva-128-341-805-classifier.qza",
    log:
        "logs/train-classifier.log",
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads {input.sequences} \
            --i-reference-taxonomy {input.taxonomy} \
            --o-classifier {output} \
            > {log} 2>&1
        """


rule dereplicate_silva:
    input:
        sequences="resources/reference/silva-128-seqs-341-805.qza",
        taxonomy="resources/reference/silva-128-tax-ba.qza",
    output:
        sequences="resources/reference/silva-128-341-805-derep.qza",
        taxonomy="resources/reference/silva-128-tax-derep.qza",
    shell:
        """
        qiime rescript dereplicate \
            --i-sequences {input.sequences} \
            --i-taxa {input.taxonomy} \
            --p-mode uniq \
            --o-dereplicated-sequences {output.sequences} \
            --o-dereplicated-taxa {output.taxonomy}
        """


rule extract_341_805:
    input:
        "resources/reference/silva-128-seqs-dna-ba-cull.qza",
    output:
        "resources/reference/silva-128-seqs-341-805.qza",
    shell:
        """
        qiime feature-classifier extract-reads \
            --i-sequences {input} \
            --p-f-primer CCTACGGGNGGCWGCAG \
            --p-r-primer GGACTACHVGGGTATCTAAT \
            --p-n-jobs 8 \
            --o-reads {output}
        """


rule cull_silva:
    input:
        "resources/reference/silva-128-seqs-dna-ba.qza",
    output:
        "resources/reference/silva-128-seqs-dna-ba-cull.qza",
    shell:
        """
        qiime rescript cull-seqs \
            --i-sequences {input} \
            --p-num-degenerates 5 \
            --p-homopolymer-length 8 \
            --o-clean-sequences {output}
        """


rule filter_silva_sequences:
    input:
        sequences="resources/reference/silva-128-seqs-dna.qza",
        taxonomy="resources/reference/silva-128-tax-ba.qza",
    output:
        "resources/reference/silva-128-seqs-dna-ba.qza",
    shell:
        """
        qiime feature-table filter-seqs \
            --i-data {input.sequences} \
            --m-metadata-file {input.taxonomy} \
            --o-filtered-data {output}
        """


rule silva_reverse_transcribe:
    input:
        "resources/reference/silva-128-seqs.qza",
    output:
        "resources/reference/silva-128-seqs-dna.qza",
    shell:
        """
        qiime rescript reverse-transcribe \
            --i-rna-sequences {input} \
            --o-dna-sequences {output}
        """


rule filter_silva_taxonomy:
    input:
        "resources/reference/silva-128-tax.qza",
    output:
        "resources/reference/silva-128-tax-ba.qza",
    shell:
        """
        qiime rescript filter-taxa \
            --i-taxonomy {input} \
            --p-include Bacteria Archaea \
            --p-exclude Chloroplast Mitochondria \
            --o-filtered-taxonomy {output}
        """


rule get_silva:
    output:
        seqs="resources/reference/silva-128-seqs.qza",
        tax="resources/reference/silva-128-tax.qza",
    shell:
        """
  mkdir -p resources/reference
  
  qiime rescript get-silva-data \
    --p-version 128 \
    --p-target SSURef_NR99 \
    --o-silva-sequences {output.seqs} \
    --o-silva-taxonomy {output.tax}
  """


rule denoising_visualization:
    input:
        denoising_stats="work/qiime2/dada2-stats.qza",
        feature_table="work/qiime2/feature-table.qza",
        sample_metadata="data/metadata/sample-metadata.tsv",
        asv_sequences="work/qiime2/asv-sequences.qza",
    output:
        dada2_stats_summ="work/qiime2/visualizations/dada2-stats-summ.qzv",
        feature_table_summ="work/qiime2/visualizations/feature-table-summ.qzv",
        feature_freqs="work/qiime2/feature-frequencies.qza",
        sample_freqs="work/qiime2/sample-frequencies.qza",
        asv_sequences_summ="work/qiime2/visualizations/asv-sequences-summ.qzv",
    shell:
        """
  qiime metadata tabulate --m-input-file {input.denoising_stats} --o-visualization {output.dada2_stats_summ}
  
  qiime feature-table summarize --i-table {input.feature_table} --m-metadata-file {input.sample_metadata} --o-feature-frequencies {output.feature_freqs} --o-sample-frequencies {output.sample_freqs} --o-summary {output.feature_table_summ}
  
  qiime feature-table tabulate-seqs --i-data {input.asv_sequences} --o-visualization {output.asv_sequences_summ}
  """


rule denoising:
    input:
        "work/qiime2/demux-paired-end-trimmed.qza",
    output:
        asv_sequences="work/qiime2/asv-sequences.qza",
        feature_table="work/qiime2/feature-table.qza",
        denoising_stats="work/qiime2/dada2-stats.qza",
        base_transition_stats="work/qiime2/base-transition-stats.qza",
    log:
        "logs/denoising.log",
    shell:
        """
  qiime dada2 denoise-paired --i-demultiplexed-seqs {input} --p-trunc-len-f 265 --p-trunc-len-r 225 --p-min-overlap 30 --p-n-threads 8 --o-representative-sequences {output.asv_sequences} --o-table {output.feature_table} --o-denoising-stats {output.denoising_stats} --o-base-transition-stats {output.base_transition_stats} > {log} 2>&1
  """


rule remove_primers:
    input:
        "work/qiime2/demux-paired-end.qza",
    output:
        "work/qiime2/demux-paired-end-trimmed.qza",
    log:
        "logs/remove-primers.log",
    shell:
        """
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences {input} \
    --p-cores 8 \
    --p-front-f CCTACGGGNGGCWGCAG \
    --p-front-r GGACTACHVGGGTATCTAAT \
    --o-trimmed-sequences {output} \
    > {log} 2>&1
  """


rule demux:
    input:
        "work/metadata/manifest.tsv",
    output:
        demux="work/qiime2/demux-paired-end.qza",
        demux_summ="work/qiime2/visualizations/demux-paired-end-summary.qzv",
    log:
        import_log="logs/demux-import.log",
        summarize_log="logs/demux-summarize.log",
    shell:
        """
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path {input} \
    --output-path {output.demux} \
    --input-format PairedEndFastqManifestPhred33V2 \
    > {log.import_log} 2>&1
     
  mkdir -p work/qiime2/visualizations
  
  qiime demux summarize \
    --i-data {output.demux} \
    --o-visualization {output.demux_summ} \
    > {log.summarize_log} 2>&1
  """


rule create_manifest:
    input:
        metadata="data/metadata/sample-metadata.tsv",
        fastq=get_all_fastq_files,
    output:
        "work/metadata/manifest.tsv",
    shell:
        """
  mkdir -p work/metadata

  RAW="$(realpath data/raw/fastq)"

  {{
    printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n'

    tail -n +3 {input.metadata} |
    while IFS=$'\t' read -r sample accession patient tissue survival; do
        r1="$RAW/${{accession}}_1.fastq.gz"
        r2="$RAW/${{accession}}_2.fastq.gz"

        if [[ ! -f "$r1" || ! -f "$r2" ]]; then
            echo "ERROR: No pair of FASTQ files for $sample ($accession)" >&2
            exit 1
        fi

        printf '%s\t%s\t%s\n' "$sample" "$r1" "$r2"
    done
  }} > {output}
  """


rule extract_metadata:
    input:
        "work/metadata/fastq-urls.tsv",
    output:
        "data/metadata/sample-metadata.tsv",
    shell:
        """
         mkdir -p data/metadata

         awk -F'\t' '
         BEGIN{{
           OFS="\t"
           print "SampleID\tAccession\tPatientID\tTissue\tSurvival"
           print "#q2:types\tcategorical\tcategorical\tcategorical\tcategorical"
         }} FNR>1 && !seen[$3]++ {{
         n=split($3, a, "_")
         if (n < 3) {{
           print "ERROR: Cannot parse sample_alias: " $3 > "/dev/stderr"
           exit 1
         }}
         print $3, $1, a[2], a[3], a[1]
         }}' {input} > {output}
         """


rule download_fastq:
    output:
        "data/raw/fastq/{accession}_{read}.fastq.gz",
    wildcard_constraints:
        read="[12]",
    params:
        url=get_fastq_url,
    log:
        "logs/download-fastq/{accession}_{read}.log",
    shell:
        """
        mkdir -p data/raw/fastq logs/download-fastq

        curl -L --fail --retry 10 --retry-delay 5 --retry-all-errors --continue-at - \
            -o {output} "{params.url}" > {log} 2>&1
        
        gzip -t {output} >> {log} 2>&1
        """


checkpoint download_fastq_urls:
    output:
        "work/metadata/fastq-urls.tsv",
    log:
        "logs/download-fastq-urls.log",
    shell:
        """
        mkdir -p work/metadata logs

        url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB57580&result=read_run&fields=run_accession,sample_accession,sample_alias,sample_title,experiment_accession,experiment_alias,experiment_title,run_alias,library_name,instrument_model,fastq_ftp,fastq_md5,fastq_bytes'

        curl -L --fail --retry 10 --retry-delay 5 --retry-all-errors \
            "$url" -o {output} > {log} 2>&1
        """
