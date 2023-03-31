rule all:
    input:
        auspice = "auspice/ToBRFV_20220412.json"

input_fasta = "data/ToBRFV_20220412.fa",
input_metadata = "data/meta_ToBRFV_20220412.tsv",
reference = "config/KT383474_NCBI-edit.gb",
lat_longs = "config/lat_longs.tsv",
description = "config/description.md",
auspice_config = "config/auspice_config.json"

#rule filter:
#    message:
#        """
#        Filtering to
#          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
#          - from {params.min_date} onwards
#          - excluding strains in {input.exclude}
#        """
#    input:
#        sequences = input_fasta,
#        metadata = input_metadata,
#        exclude = dropped_strains
#    output:
#        sequences = "results/filtered.fasta"
#    params:
#        group_by = "country year month",
#        sequences_per_group = 20,
#        min_date = 2012
#    shell:
#        """
#        augur filter \
#            --sequences {input.sequences} \
#            --metadata {input.metadata} \
#            --exclude {input.exclude} \
#            --output {output.sequences} \
#            --group-by {params.group_by} \
#            --sequences-per-group {params.sequences_per_group} \
#            --min-date {params.min_date}
#        """

rule align:
    message:
        """
        Aligning sequences to {input}
          - filling gaps with N
        """
    input:
        sequences = input_fasta
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --output {output.alignment} \
            --method mafft \
            --nthreads auto
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --method raxml \
            --output {output.tree}\
            --substitution-model GTR \
            --nthreads auto
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = input_metadata
    output:
        tree = "results/refined_tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference}
        """

rule traits:
    message: "Infer ancestral traits from an existing phylogenetic tree and the metadata annotating each tip of the tree"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata
    output:
        traits = "results/traits.json"
    params:
        columns = "continent country state municipality host"

    shell:
        """
        augur traits \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --output {output.traits} \
        --columns {params.columns}\
        --confidence
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment
    output:
        node_data = "results/nt_muts.json",
        sequences = "results/nucleotide-mutations.fasta"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --output-sequences {output.sequences}\
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule export:
    message: "Exporting data files fo auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.traits,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        description = description,
        lat_longs = lat_longs,
        auspice_config = auspice_config
    output:
        auspice = rules.all.input.auspice,
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --output {output.auspice} \
        """
