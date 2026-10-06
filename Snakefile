from pathlib import Path

import pandas as pd


TMPDIR = Path.home() / "qiime2-tmp"
TMPDIR.mkdir(parents=True, exist_ok=True)

shell.prefix(f"export TMPDIR='{TMPDIR}'; ")


def get_all_fastq_files(wildcards):
    checkpoint_output = checkpoints.download_fastq_urls.get().output[0]
    df = pd.read_csv(checkpoint_output, sep="\t")

    accessions = (
        df["run_accession"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    return [
        f"raw-data/{accession}_{read}.fastq.gz"
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
    "bray-curtis",
    "jaccard",
    "aitchison",
    "weighted-unifrac",
    "unweighted-unifrac"
]


rule all:
    input:
        "jaccard-pcoa.qza",
        "bray-curtis-pcoa.qza",
        "aitchison-pcoa.qza",
        "unweighted-unifrac-pcoa.qza",
        "weighted-unifrac-pcoa.qza",
        "jaccard-tissue-significance.qzv",
        "bray-curtis-tissue-significance.qzv",
        "aitchison-tissue-significance.qzv",
        "unweighted-unifrac-tissue-significance.qzv",
        "weighted-unifrac-tissue-significance.qzv",
        "jaccard-emperor.qzv",
        "bray-curtis-emperor.qzv",
        "aitchison-emperor.qzv",
        "unweighted-unifrac-emperor.qzv",
        "weighted-unifrac-emperor.qzv",
        "report/tex-tables/permanova-results.tex"

# expand("raw-data/{file}.fastq.gz", file=FILES)

rule permanova_to_latex:
    input:
        expand("results/{score}-matrix-permanova.txt",
            score=BETA_SCORES)
    output:
        "report/tex-tables/permanova-results.tex"
    shell:
        """
        python3 scripts/permanova-to-latex.py {input} {output}
        """

rule permanova_r:
    input:
        "export-diversity-scores/{score}-matrix.tsv"
    output:
        "results/{score}-matrix-permanova.txt"
    shell:
        """
        Rscript scripts/permanova.R {input} > {output}
        """

rule tissue_significance:
    input:
        score="{file}-diversity.qza",
        sample_metadata="sample-metadata.tsv"
    output:
        "{file}-tissue-significance.qzv"
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
        pcoa="{file}-pcoa.qza",
        sample_metadata="sample-metadata.tsv"
    output:
        "{file}-emperor.qzv"
    shell:
        """
        qiime emperor plot \
            --i-pcoa {input.pcoa} \
            --m-metadata-file {input.sample_metadata} \
            --o-visualization {output}
        """

rule pcoa:
    input:
        "{file}-diversity.qza"
    output:
        "{file}-pcoa.qza"
    shell:
        """
        qiime diversity pcoa \
        --i-distance-matrix {input} \
        --o-pcoa {output}
        """

rule export_diversity_scores:
    input:
        "{score}-diversity.qza"
    output:
        "export-diversity-scores/{score}-matrix.tsv"
    shell:
        """
        tmpdir=$(mktemp -d)

        qiime tools export \
            --input-path {input} \
            --output-path "$tmpdir"

        mv "$tmpdir/distance-matrix.tsv" {output}

        rmdir "$tmpdir"
        """

rule evaluate_diversity_phylogeny:
    input:
        feature_table_2500="feature-table-2500.qza",
        rooted_tree="rooted-tree.qza"
    output:
        unweighted_unifrac="unweighted-unifrac-diversity.qza",
        weighted_unifrac="weighted-unifrac-diversity.qza"
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

rule phylogeny:
    input:
        "asv-sequences.qza"
    output:
        aligned_rep_seqs="aligned-rep-seqs.qza",
        masked_aligned_rep_seqs="masked-aligned-rep-seqs.qza",
        unrooted_tree="unrooted-tree.qza",
        rooted_tree="rooted-tree.qza"
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences {input} \
        --o-alignment {output.aligned_rep_seqs} \
        --o-masked-alignment {output.masked_aligned_rep_seqs} \
        --o-tree {output.unrooted_tree} \
        --o-rooted-tree {output.rooted_tree}
        """

rule evaluate_diversity:
    input:
        feature_table="feature-table.qza",
        feature_table_2500="feature-table-2500.qza"
    output:
        bray_curtis="bray-curtis-diversity.qza",
        aitchison="aitchison-diversity.qza",
        jaccard="jaccard-diversity.qza"
    shell:
        """
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
        "feature-table.qza"
    output:
        "feature-table-2500.qza"
    shell:
        """
        qiime feature-table rarefy \
            --i-table {input} \
            --p-sampling-depth 2500 \
            --o-rarefied-table {output}
        """

rule feature_table_summarize:
    input:
        feature_table="feature-table.qza",
        sample_metadata="sample-metadata.tsv"
    output:
        dir=directory("feature-table-summarize-output-dir"),
        summary="feature-table-summary.qzv"
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
        taxonomy="taxonomy.qza",
        feature_table="feature-table.qza",
        sample_metadata="sample-metadata.tsv"

    output:
        taxonomy="taxonomy.qzv",
        taxa_barplot="taxa-barplot.qzv"
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

rule classify_asvs:
    input:
        classifier="reference/silva-128-341-805-classifier.qza",
        sequences="asv-sequences.qza"
    output:
        "taxonomy.qza"
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.sequences} \
            --o-classification {output}
        """

rule train_classifier:
    input:
        sequences="reference/silva-128-341-805-derep.qza",
        taxonomy="reference/silva-128-tax-derep.qza"
    output:
        "reference/silva-128-341-805-classifier.qza"
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads {input.sequences} \
            --i-reference-taxonomy {input.taxonomy} \
            --o-classifier {output}
        """

rule dereplicate_silva:
    input:
        sequences="reference/silva-128-seqs-341-805.qza",
        taxonomy="reference/silva-128-tax-ba.qza"
    output:
        sequences="reference/silva-128-341-805-derep.qza",
        taxonomy="reference/silva-128-tax-derep.qza"
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
        "reference/silva-128-seqs-dna-ba-cull.qza"
    output:
        "reference/silva-128-seqs-341-805.qza"
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
        "reference/silva-128-seqs-dna-ba.qza"
    output:
        "reference/silva-128-seqs-dna-ba-cull.qza"
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
        sequences="reference/silva-128-seqs-dna.qza",
        taxonomy="reference/silva-128-tax-ba.qza"
    output:
        "reference/silva-128-seqs-dna-ba.qza"
    shell:
        """
        qiime feature-table filter-seqs \
            --i-data {input.sequences} \
            --m-metadata-file {input.taxonomy} \
            --o-filtered-data {output}
        """

rule silva_reverse_transcribe:
    input:
        "reference/silva-128-seqs.qza"
    output:
        "reference/silva-128-seqs-dna.qza"
    shell:
        """
        qiime rescript reverse-transcribe \
            --i-rna-sequences {input} \
            --o-dna-sequences {output}
        """

rule filter_silva_taxonomy:
    input: "reference/silva-128-tax.qza"
    output: "reference/silva-128-tax-ba.qza"
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
        seqs="reference/silva-128-seqs.qza",
        tax="reference/silva-128-tax.qza"
    shell: """
  mkdir -p reference
  
  qiime rescript get-silva-data \
    --p-version 128 \
    --p-target SSURef_NR99 \
    --o-silva-sequences {output.seqs} \
    --o-silva-taxonomy {output.tax}
  """

rule denoising_visualization:
    input:
        denoising_stats="dada2-stats.qza",
        feature_table="feature-table.qza",
        sample_metadata="sample-metadata.tsv",
        asv_sequences="asv-sequences.qza"
    output:
        dada2_stats_summ="dada2-stats-summ.qzv",
        feature_table_summ="feature-table-summ.qzv",
        feature_freqs="feature-frequencies.qza",
        sample_freqs="sample-frequencies.qza",
        asv_sequences_summ="asv-sequences-summ.qzv"
    shell: """
  qiime metadata tabulate --m-input-file {input.denoising_stats} --o-visualization {output.dada2_stats_summ}
  
  qiime feature-table summarize --i-table {input.feature_table} --m-metadata-file {input.sample_metadata} --o-feature-frequencies {output.feature_freqs} --o-sample-frequencies {output.sample_freqs} --o-summary {output.feature_table_summ}
  
  qiime feature-table tabulate-seqs --i-data {input.asv_sequences} --o-visualization {output.asv_sequences_summ}
  """

rule denoising:
    input: "demux-paired-end-trimmed.qza"
    output:
        asv_sequences="asv-sequences.qza",
        feature_table="feature-table.qza",
        denoising_stats="dada2-stats.qza",
        base_transition_stats="base-transition-stats.qza",
    shell: """
  qiime dada2 denoise-paired --i-demultiplexed-seqs {input} --p-trunc-len-f 265 --p-trunc-len-r 225 --p-min-overlap 30 --p-n-threads 8 --o-representative-sequences {output.asv_sequences} --o-table {output.feature_table} --o-denoising-stats {output.denoising_stats} --o-base-transition-stats {output.base_transition_stats}
  """

rule remove_primers:
    input: "demux-paired-end.qza"
    output: "demux-paired-end-trimmed.qza"
    shell: """
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences {input} \
    --p-cores 8 \
    --p-front-f CCTACGGGNGGCWGCAG \
    --p-front-r GGACTACHVGGGTATCTAAT \
    --o-trimmed-sequences {output}
  """

rule demux:
    input: "manifest.tsv"
    output:
        demux="demux-paired-end.qza",
        demux_summ="demux-paired-end-summary.qzv"
    log: "logs/demux_summarize.log"
    shell: """
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path {input} \
    --output-path {output.demux} \
    --input-format PairedEndFastqManifestPhred33V2
    
  mkdir -p logs
  
  qiime demux summarize \
    --i-data {output.demux} \
    --o-visualization {output.demux_summ} \
    > {log} 2>&1
  """

rule create_manifest:
    input:
        metadata="sample-metadata.tsv",
        fastq=get_all_fastq_files
    output:
        "manifest.tsv"
    shell: """
  RAW="$(realpath raw-data)"

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
        "fastq-urls.tsv"
    output:
        "sample-metadata.tsv"
    shell: """
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
        "raw-data/{accession}_{read}.fastq.gz"
    wildcard_constraints:
        read="[12]"
    params:
        url=get_fastq_url
    shell:
        """
        mkdir -p raw-data
        curl -L --fail --retry 10 --retry-delay 5 --retry-all-errors --continue-at - \
            -o {output} "{params.url}"
        gzip -t {output}
        """

checkpoint download_fastq_urls:
    output:
        "fastq-urls.tsv"
    shell:
        """
        url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB57580&result=read_run&fields=run_accession,sample_accession,sample_alias,sample_title,experiment_accession,experiment_alias,experiment_title,run_alias,library_name,instrument_model,fastq_ftp,fastq_md5,fastq_bytes'

        curl -L --fail --retry 10 --retry-delay 5 --retry-all-errors \
            "$url" -o {output}
        """
