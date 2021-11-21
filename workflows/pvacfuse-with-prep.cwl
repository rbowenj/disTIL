class: Workflow
cwlVersion: v1.0
id: pvacfuse_with_prep
doc: >-
  # About this workflow

  This workflow runs all pVACfuse input file preparation steps followed by
  pVACfuse prediction of neoepitopes.


  ## Steps

  1. `agfusion` - for more details, this tool can be found in the repo.

  2. `pvacfuse` v2.0.4
label: pvacfuse-with-prep
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: epitope_lengths_class_i
    type: 'int[]?'
    'sbg:exposed': true
  - id: epitope_lengths_class_ii
    type: 'int[]?'
    'sbg:exposed': true
  - id: binding_threshold
    type: int?
    'sbg:exposed': true
  - id: percentile_threshold
    type: int?
    'sbg:exposed': true
  - id: iedb_retries
    type: int?
    'sbg:exposed': true
  - id: keep_tmp_files
    type: boolean?
    'sbg:exposed': true
  - id: net_chop_method
    type:
      - 'null'
      - type: enum
        symbols:
          - cterm
          - 20s
        name: net_chop_method
    'sbg:exposed': true
  - id: netmhc_stab
    type: boolean?
    'sbg:exposed': true
  - id: top_score_metric
    type:
      - 'null'
      - type: enum
        symbols:
          - lowest
          - median
        name: top_score_metric
    'sbg:exposed': true
  - id: net_chop_threshold
    type: float?
    'sbg:exposed': true
  - id: run_reference_proteome_similarity
    type: boolean?
    'sbg:exposed': true
  - id: additional_report_columns
    type:
      - 'null'
      - type: enum
        symbols:
          - sample_name
        name: additional_report_columns
    'sbg:exposed': true
  - id: fasta_size
    type: int?
    'sbg:exposed': true
  - id: downstream_sequence_length
    type: string?
    'sbg:exposed': true
  - id: exclude_nas
    type: boolean?
    'sbg:exposed': true
  - id: n_threads
    type: int?
    'sbg:exposed': true
  - id: alleles
    'sbg:fileTypes': TXT
    type: File
    label: HLA Alleles
    doc: The text file containing HLA alleles for the tumour sample.
    'sbg:x': 185.5625
    'sbg:y': 93
  - id: sample_name
    type: string
    label: Tumour Sample Name
    doc: The tumour sample name.
    'sbg:x': 0
    'sbg:y': 0
  - id: ref_genome
    type:
      type: enum
      symbols:
        - hg38
        - hg19
      name: ref_genome
    label: Reference Genome Build
    doc: >-
      The referece genome build used to generate the fusion variant calls. This
      is used to determine which agfusion database to use (release 87 or 75).
    'sbg:x': 0
    'sbg:y': 107
  - id: fusion_tsv
    'sbg:fileTypes': TSV
    type: File
    label: Fusion TSV
    doc: The fusion TSV to be annotated and used for neoepitope prediction.
    'sbg:x': 0
    'sbg:y': 214
  - id: fusion_caller
    type:
      type: enum
      symbols:
        - starfusion
        - jaffa
        - bellerophontes
        - breakfusion
        - chimerascan
        - chimerscope
        - defuse
        - ericscript
        - fusioncatcher
        - fusionhunter
        - fusionmap
        - fusioninspector
        - infusion
        - mapsplice
        - tophatfusion
      name: fusion_caller
    'sbg:x': 0
    'sbg:y': 321
  - id: prediction_algorithms
    type: 'string[]'
    'sbg:exposed': true
outputs:
  - id: pvacfuse_predictions
    outputSource:
      - pvacfuse/pvacfuse_predictions
    type: Directory?
    label: pVACfuse Predictions
    doc: The directory containing all pVACfuse output files.
    'sbg:x': 912.2371826171875
    'sbg:y': 0
  - id: mhc_ii_aggregated_report
    outputSource:
      - pvacfuse/mhc_ii_aggregated_report
    type: File?
    label: MHC-II Aggregated Report
    doc: The aggregated report for all MHC-II-binding predicted neoepitopes.
    'sbg:x': 912.2371826171875
    'sbg:y': 107
  - id: mhc_i_aggregated_report
    outputSource:
      - pvacfuse/mhc_i_aggregated_report
    type: File?
    label: MHC I Aggregated Report
    doc: The aggregated report for all MHC-I-binding predicted neoepitopes.
    'sbg:x': 912.2371826171875
    'sbg:y': 214
  - id: combined_aggregated_report
    outputSource:
      - pvacfuse/combined_aggregated_report
    type: File?
    label: Combined Aggregated Report
    doc: >-
      The report containing all (MHC-I and/or MHC-II binding) predicted
      neoepitopes from pVACfuse.
    'sbg:x': 912.2371826171875
    'sbg:y': 321
  - id: agfusion_output
    outputSource:
      - agfusion/output_dir
    type: Directory?
    label: AGFusion Output
    doc: AGFusion output directory
    'sbg:x': 481.8861083984375
    'sbg:y': 214
  - id: mhc_ii_filtered_epitopes
    outputSource:
      - pvacfuse/mhc_ii_filtered_epitopes
    type: File?
    label: MHC-II Filtered Neoepitope Report
    doc: Filtered MHC-II-binding predicted neoepitopes.
    'sbg:x': 850.0859375
    'sbg:y': -317.598388671875
  - id: mhc_i_filtered_epitopes
    outputSource:
      - pvacfuse/mhc_i_filtered_epitopes
    type: File?
    label: MHC-I Filtered Neoepitope Report
    doc: Filtered MHC-I-binding predicted neoepitopes.
    'sbg:x': 1041.0859375
    'sbg:y': -217.59840393066406
  - id: combined_filtered_epitopes
    outputSource:
      - pvacfuse/combined_filtered_epitopes
    type: File?
    label: Combined Filtered Neoepitope Report
    doc: Filtered MHC-I- and/or MHC-II-binding predicted neoepitopes.
    'sbg:x': 1239.28271484375
    'sbg:y': -86.59840393066406
steps:
  - id: agfusion
    in:
      - id: fusion_tsv
        source: fusion_tsv
      - id: fusion_caller
        source: fusion_caller
      - id: ref_genome
        source: ref_genome
    out:
      - id: output_dir
    run: ../tools/agfusion.cwl
    label: agfusion
    doc: Run agfusion annotation of gene fusions.
    'sbg:x': 185.5625
    'sbg:y': 214
  - id: pvacfuse
    in:
      - id: input_file
        source: agfusion/output_dir
      - id: sample_name
        source: sample_name
      - id: alleles
        source: alleles
      - id: prediction_algorithms
        default:
          - all
        source:
          - prediction_algorithms
      - id: epitope_lengths_class_i
        source:
          - epitope_lengths_class_i
      - id: epitope_lengths_class_ii
        source:
          - epitope_lengths_class_ii
      - id: binding_threshold
        source: binding_threshold
      - id: percentile_threshold
        source: percentile_threshold
      - id: iedb_retries
        source: iedb_retries
      - id: keep_tmp_files
        source: keep_tmp_files
      - id: net_chop_method
        source: net_chop_method
      - id: netmhc_stab
        source: netmhc_stab
      - id: top_score_metric
        source: top_score_metric
      - id: net_chop_threshold
        source: net_chop_threshold
      - id: run_reference_proteome_similarity
        source: run_reference_proteome_similarity
      - id: additional_report_columns
        source: additional_report_columns
      - id: fasta_size
        source: fasta_size
      - id: downstream_sequence_length
        source: downstream_sequence_length
      - id: exclude_nas
        source: exclude_nas
      - id: n_threads
        source: n_threads
    out:
      - id: mhc_i_all_epitopes
      - id: mhc_i_aggregated_report
      - id: mhc_i_filtered_epitopes
      - id: mhc_ii_all_epitopes
      - id: mhc_ii_aggregated_report
      - id: mhc_ii_filtered_epitopes
      - id: combined_all_epitopes
      - id: combined_aggregated_report
      - id: combined_filtered_epitopes
      - id: mhc_i
      - id: mhc_ii
      - id: pvacfuse_predictions
    run: ../tools/pvacfuse.cwl
    label: pvacfuse
    'sbg:x': 481.8861083984375
    'sbg:y': 30
requirements: []
