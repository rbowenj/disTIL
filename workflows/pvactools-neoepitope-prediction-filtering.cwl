class: Workflow
cwlVersion: v1.0
id: pvactools_neoepitope_prediction_filtering
doc: >-
  # About this workflow

  This workflow runs pVACseq and pVACfuse (including all input file
  pre-processing and annotation) disTIL pVACfuse STAR annotation.
label: pvactools-neoepitope-prediction-filtering
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: alleles
    'sbg:fileTypes': TXT
    type: File
    label: HLA Alleles
    doc: The text file containing HLA alleles for the tumour sample.
    'sbg:x': 0
    'sbg:y': 1390.1875
  - id: fusion_tsv
    'sbg:fileTypes': TSV
    type: File
    label: Input Fusion TSV
    doc: The input TSV of gene fusions called from the tumour sample.
    'sbg:x': 0
    'sbg:y': 1283.25
  - id: gene_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Gene Expression File
    doc: Gene level expression file.
    'sbg:x': 0
    'sbg:y': 1176.3125
  - id: input_bam_dna
    'sbg:fileTypes': BAM
    type: File
    label: DNA BAM
    doc: The DNA BAM for the sample of interest.
    'sbg:x': 0
    'sbg:y': 1069.375
  - id: input_bam_rna
    'sbg:fileTypes': BAM
    type: File
    label: RNA BAM
    doc: The RNA BAM for the sample of interest.
    'sbg:x': 0
    'sbg:y': 962.4375
  - id: input_vcf
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    label: Input VCF
    'sbg:x': 0
    'sbg:y': 855.5
  - id: ref_genome_dna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: DNA Reference Genome
    doc: Reference sequence used to align the DNA BAM.
    'sbg:x': 0
    'sbg:y': 534.6875
  - id: ref_genome_rna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: RNA Reference Genome
    doc: Reference sequence used to align the RNA BAM.
    'sbg:x': 0
    'sbg:y': 427.75
  - id: transcript_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Transcript Expression File
    doc: Transcript level expression file.
    'sbg:x': 0
    'sbg:y': 213.875
  - id: vep_cache
    'sbg:fileTypes': TAR.GZ
    type: File
    label: VEP Cache
    doc: A VEP cache TAR file.
    'sbg:x': 0
    'sbg:y': 106.9375
  - id: vep_plugin_files
    'sbg:fileTypes': PM
    type: 'File[]'
    label: VEP Plugin Files
    doc: VEP plugin files to use in VCF annotation.
    'sbg:x': 0
    'sbg:y': 0
  - id: phased_proximal_variants_vcf
    type: File?
    label: Phased Proximal Variants VCF
    doc: VCF of phased proximal variants (see pVACseq documentation).
    'sbg:x': 0
    'sbg:y': 641.625
  - id: pvacfuse_binding_threshold
    type: int?
    'sbg:exposed': true
  - id: pvacfuse_percentile_threshold
    type: int?
    'sbg:exposed': true
  - id: pvacfuse_iedb_retries
    type: int?
    'sbg:exposed': true
  - id: pvacfuse_keep_tmp_files
    type: boolean?
    'sbg:exposed': true
  - id: pvacfuse_top_score_metric
    type:
      - 'null'
      - type: enum
        symbols:
          - lowest
          - median
        name: top_score_metric
    'sbg:exposed': true
  - id: pvacfuse_net_chop_threshold
    type: float?
    'sbg:exposed': true
  - id: pvacfuse_run_reference_proteome_similarity
    type: boolean?
    'sbg:exposed': true
  - id: pvacfuse_additional_report_columns
    type:
      - 'null'
      - type: enum
        symbols:
          - sample_name
        name: additional_report_columns
    'sbg:exposed': true
  - id: pvacfuse_fasta_size
    type: int?
    'sbg:exposed': true
  - id: pvacfuse_downstream_sequence_length
    type: string?
    'sbg:exposed': true
  - id: pvacfuse_exclude_nas
    type: boolean?
    'sbg:exposed': true
  - id: n_threads
    type: int?
    'sbg:exposed': true
  - id: ref_genome
    type:
      type: enum
      symbols:
        - hg38
        - hg19
      name: ref_genome
    'sbg:exposed': true
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
    'sbg:exposed': true
  - id: pvacfuse_epitope_lengths_class_i
    type: 'int[]?'
    'sbg:exposed': true
  - id: pvacfuse_epitope_lengths_class_ii
    type: 'int[]?'
    'sbg:exposed': true
  - id: pvacseq_epitope_lengths_class_i
    type: 'int[]?'
    'sbg:exposed': true
  - id: pvacseq_epitope_lengths_class_ii
    type: 'int[]?'
    'sbg:exposed': true
  - id: pvacfuse_net_chop_method
    type:
      - 'null'
      - type: enum
        symbols:
          - cterm
          - 20s
        name: net_chop_method
    'sbg:exposed': true
  - id: pvacfuse_netmhc_stab
    type: boolean?
    'sbg:exposed': true
  - id: pvacseq_run_reference_proteome_similarity
    type: boolean?
    'sbg:exposed': true
  - id: gene_id_column
    type: string?
    'sbg:exposed': true
  - id: gene_expression_column
    type: string?
    'sbg:exposed': true
  - id: gene_quant_algo
    type:
      type: enum
      symbols:
        - kallisto
        - stringtie
        - cufflinks
        - custom
      name: gene_quant_algo
    'sbg:exposed': true
  - id: transcript_quant_algo
    type:
      type: enum
      symbols:
        - kallisto
        - stringtie
        - cufflinks
        - custom
      name: transcript_quant_algo
    'sbg:exposed': true
  - id: pvacseq_min_base_qual
    type: int?
    'sbg:exposed': true
  - id: pvacseq_min_mapping_qual
    type: int?
    'sbg:exposed': true
  - id: cache_version
    type: int
    'sbg:exposed': true
  - id: pvacseq_net_chop_method
    type:
      - 'null'
      - type: enum
        symbols:
          - cterm
          - 20s
        name: net_chop_method
    'sbg:exposed': true
  - id: pvacseq_netmhc_stab
    type: boolean?
    'sbg:exposed': true
  - id: pvacseq_additional_report_columns
    type:
      - 'null'
      - type: enum
        symbols:
          - sample_name
        name: additional_report_columns
    'sbg:exposed': true
  - id: pvacseq_allele_specific_binding_thresholds
    type: boolean?
    'sbg:exposed': true
  - id: pvacseq_binding_threshold
    type: int?
    'sbg:exposed': true
  - id: pvacseq_downstream_sequence_length
    type: string?
    'sbg:exposed': true
  - id: pvacseq_exclude_nas
    type: boolean?
    'sbg:exposed': true
  - id: pvacseq_expn_val
    type: float?
    'sbg:exposed': true
  - id: pvacseq_fasta_size
    type: int?
    'sbg:exposed': true
  - id: pvacseq_iedb_retries
    type: int?
    'sbg:exposed': true
  - id: pvacseq_keep_tmp_files
    type: boolean?
    'sbg:exposed': true
  - id: pvacseq_maximum_transcript_support_level
    type:
      - 'null'
      - type: enum
        symbols:
          - '1'
          - '2'
          - '3'
          - '4'
          - '5'
        name: maximum_transcript_support_level
    'sbg:exposed': true
  - id: pvacseq_minimum_fold_change
    type: float?
    'sbg:exposed': true
  - id: pvacseq_net_chop_threshold
    type: float?
    'sbg:exposed': true
  - id: normal_cov
    type: int?
    'sbg:exposed': true
  - id: normal_vaf
    type: float?
    'sbg:exposed': true
  - id: pvacseq_percentile_threshold
    type: int?
    'sbg:exposed': true
  - id: pvacseq_prediction_algorithms
    type: 'string[]'
    'sbg:exposed': true
  - id: tdna_cov
    type: int?
    'sbg:exposed': true
  - id: tdna_vaf
    type: float?
    'sbg:exposed': true
  - id: pvacseq_top_score_metric
    type:
      - 'null'
      - type: enum
        symbols:
          - lowest
          - median
        name: top_score_metric
    'sbg:exposed': true
  - id: trna_cov
    type: int?
    'sbg:exposed': true
  - id: trna_vaf
    type: float?
    'sbg:exposed': true
  - id: pvacfuse_prediction_algorithms
    type: 'string[]'
    'sbg:exposed': true
  - id: transcript_id_column
    type: string?
    'sbg:exposed': true
  - id: transcript_expression_column
    type: string?
    'sbg:exposed': true
  - id: normal_sample_name
    type: string?
    label: Normal Sample Name
    doc: The name of the matched normal sample in the input VCF.
    'sbg:x': 0
    'sbg:y': 748.5625
  - id: sample_name
    type: string
    label: Tumour Sample Name
    doc: >-
      The name of the tumour sample in the input VCF to use for neoepitope
      prediction.
    'sbg:x': 0
    'sbg:y': 320.8125
outputs:
  - id: pvacfuse_star_filtered_
    outputSource:
      - pvacfuse_star_filter_mhc_i/pvacfuse_star_filtered
    'sbg:fileTypes': TSV
    type: File
    label: pVACfuse MHC-I Star-Annotated Filtered Report
    'sbg:x': 1472.73681640625
    'sbg:y': 748.5625
  - id: pvacfuse_star_filtered_mhc_ii
    outputSource:
      - pvacfuse_star_filter_mhc_ii/pvacfuse_star_filtered
    'sbg:fileTypes': TSV
    type: File
    label: pVACfuse MHC-II Star-Annotated Filtered Report
    'sbg:x': 1472.73681640625
    'sbg:y': 641.625
  - id: agfusion_output
    outputSource:
      - pvactools_neoepitope_prediction/agfusion_output
    type: Directory?
    label: AGFusion Output
    'sbg:x': 990.9447631835938
    'sbg:y': 1083.375
  - id: pvacseq_predictions
    outputSource:
      - pvactools_neoepitope_prediction/pvacseq_predictions
    type: Directory
    label: pVACseq Predictions
    'sbg:x': 990.9447631835938
    'sbg:y': 413.75
  - id: pvacfuse_predictions
    outputSource:
      - pvactools_neoepitope_prediction/pvacfuse_predictions
    type: Directory?
    label: pVACfuse Predictions
    'sbg:x': 990.9447631835938
    'sbg:y': 976.4375
  - id: pvacseq_mhc_ii_filtered_epitopes
    outputSource:
      - pvactools_neoepitope_prediction/mhc_ii_filtered_epitopes
    type: File?
    label: pVACseq MHC-II Filtered Report
    'sbg:x': 990.9447631835938
    'sbg:y': 520.6875
  - id: pvacseq_mhc_i_filtered_epitopes
    outputSource:
      - pvactools_neoepitope_prediction/mhc_i_filtered_epitopes
    type: File?
    label: pVACseq MHC-I Filtered Report
    'sbg:x': 990.9447631835938
    'sbg:y': 627.625
  - id: rc_exp_annotated_vcf
    outputSource:
      - pvactools_neoepitope_prediction/rc_exp_annotated_vcf
    'sbg:fileTypes': VCF.GZ
    type: File
    label: Fully Annotated VCF
    'sbg:x': 990.9447631835938
    'sbg:y': 306.8125
steps:
  - id: pvactools_neoepitope_prediction
    in:
      - id: pvacfuse_binding_threshold
        source: pvacfuse_binding_threshold
      - id: pvacfuse_percentile_threshold
        source: pvacfuse_percentile_threshold
      - id: pvacfuse_iedb_retries
        source: pvacfuse_iedb_retries
      - id: pvacfuse_keep_tmp_files
        source: pvacfuse_keep_tmp_files
      - id: pvacfuse_top_score_metric
        source: pvacfuse_top_score_metric
      - id: pvacfuse_net_chop_threshold
        source: pvacfuse_net_chop_threshold
      - id: pvacfuse_run_reference_proteome_similarity
        source: pvacfuse_run_reference_proteome_similarity
      - id: pvacfuse_additional_report_columns
        source: pvacfuse_additional_report_columns
      - id: pvacfuse_fasta_size
        source: pvacfuse_fasta_size
      - id: pvacfuse_downstream_sequence_length
        source: pvacfuse_downstream_sequence_length
      - id: pvacfuse_exclude_nas
        source: pvacfuse_exclude_nas
      - id: n_threads
        source: n_threads
      - id: ref_genome
        source: ref_genome
      - id: fusion_caller
        source: fusion_caller
      - id: pvacfuse_epitope_lengths_class_i
        source:
          - pvacfuse_epitope_lengths_class_i
      - id: pvacfuse_epitope_lengths_class_ii
        source:
          - pvacfuse_epitope_lengths_class_ii
      - id: pvacseq_epitope_lengths_class_i
        source:
          - pvacseq_epitope_lengths_class_i
      - id: pvacseq_epitope_lengths_class_ii
        source:
          - pvacseq_epitope_lengths_class_ii
      - id: pvacfuse_net_chop_method
        source: pvacfuse_net_chop_method
      - id: pvacfuse_netmhc_stab
        source: pvacfuse_netmhc_stab
      - id: pvacseq_run_reference_proteome_similarity
        source: pvacseq_run_reference_proteome_similarity
      - id: gene_id_column
        source: gene_id_column
      - id: gene_expression_column
        source: gene_expression_column
      - id: gene_quant_algo
        source: gene_quant_algo
      - id: transcript_quant_algo
        source: transcript_quant_algo
      - id: pvacseq_min_base_qual
        source: pvacseq_min_base_qual
      - id: pvacseq_min_mapping_qual
        source: pvacseq_min_mapping_qual
      - id: cache_version
        source: cache_version
      - id: pvacseq_net_chop_method
        source: pvacseq_net_chop_method
      - id: pvacseq_netmhc_stab
        source: pvacseq_netmhc_stab
      - id: fusion_tsv
        source: fusion_tsv
      - id: alleles
        source: alleles
      - id: vep_plugin_files
        source:
          - vep_plugin_files
      - id: vep_cache
        source: vep_cache
      - id: transcript_expression_file
        source: transcript_expression_file
      - id: ref_genome_rna
        source: ref_genome_rna
      - id: ref_genome_dna
        source: ref_genome_dna
      - id: input_vcf
        source: input_vcf
      - id: input_bam_rna
        source: input_bam_rna
      - id: input_bam_dna
        source: input_bam_dna
      - id: gene_expression_file
        source: gene_expression_file
      - id: sample_name
        source: sample_name
      - id: normal_sample_name
        source: normal_sample_name
      - id: phased_proximal_variants_vcf
        source: phased_proximal_variants_vcf
      - id: pvacseq_additional_report_columns
        source: pvacseq_additional_report_columns
      - id: pvacseq_allele_specific_binding_thresholds
        source: pvacseq_allele_specific_binding_thresholds
      - id: pvacseq_binding_threshold
        source: pvacseq_binding_threshold
      - id: pvacseq_downstream_sequence_length
        source: pvacseq_downstream_sequence_length
      - id: pvacseq_exclude_nas
        source: pvacseq_exclude_nas
      - id: pvacseq_expn_val
        source: pvacseq_expn_val
      - id: pvacseq_fasta_size
        source: pvacseq_fasta_size
      - id: pvacseq_iedb_retries
        source: pvacseq_iedb_retries
      - id: pvacseq_keep_tmp_files
        source: pvacseq_keep_tmp_files
      - id: pvacseq_maximum_transcript_support_level
        source: pvacseq_maximum_transcript_support_level
      - id: pvacseq_minimum_fold_change
        source: pvacseq_minimum_fold_change
      - id: pvacseq_net_chop_threshold
        source: pvacseq_net_chop_threshold
      - id: normal_cov
        source: normal_cov
      - id: normal_vaf
        source: normal_vaf
      - id: pvacseq_percentile_threshold
        source: pvacseq_percentile_threshold
      - id: pvacseq_prediction_algorithms
        source:
          - pvacseq_prediction_algorithms
      - id: tdna_cov
        source: tdna_cov
      - id: tdna_vaf
        source: tdna_vaf
      - id: pvacseq_top_score_metric
        source: pvacseq_top_score_metric
      - id: trna_cov
        source: trna_cov
      - id: trna_vaf
        source: trna_vaf
      - id: pvacfuse_prediction_algorithms
        source:
          - pvacfuse_prediction_algorithms
      - id: transcript_id_column
        source: transcript_id_column
      - id: transcript_expression_column
        source: transcript_expression_column
    out:
      - id: rc_exp_annotated_vcf
      - id: pvacseq_predictions
      - id: pvacseq_mhc_ii_aggregated_report
      - id: pvacseq_mhc_i_aggregated_report
      - id: pvacseq_combined_aggregated_report
      - id: pvacfuse_predictions
      - id: pvacfuse_mhc_ii_aggregated_report
      - id: pvacfuse_mhc_i_aggregated_report
      - id: pvacfuse_combined_aggregated_report
      - id: mhc_ii_filtered_epitopes
      - id: mhc_i_filtered_epitopes
      - id: combined_filtered_epitopes
      - id: agfusion_output
      - id: mhc_i_filtered_epitopes_1
      - id: mhc_ii_filtered_epitopes_1
      - id: combined_filtered_epitopes_1
    run: ./pvactools-neoepitope-prediction.cwl
    label: pvactools-neoepitope-prediction
    'sbg:x': 294.484375
    'sbg:y': 590.09375
  - id: pvacfuse_star_filter_mhc_ii
    in:
      - id: pvacfuse_report
        source: pvactools_neoepitope_prediction/mhc_ii_filtered_epitopes_1
      - id: star_fusions
        source: fusion_tsv
    out:
      - id: pvacfuse_star_filtered
    run: ../tools/pvacfuse-star-filter.cwl
    label: pvacfuse-star-filter-mhc-ii
    'sbg:x': 990.9447631835938
    'sbg:y': 741.5625
  - id: pvacfuse_star_filter_mhc_i
    in:
      - id: pvacfuse_report
        source: pvactools_neoepitope_prediction/mhc_i_filtered_epitopes_1
      - id: star_fusions
        source: fusion_tsv
    out:
      - id: pvacfuse_star_filtered
    run: ../tools/pvacfuse-star-filter.cwl
    label: pvacfuse-star-filter-mhc-i
    'sbg:x': 990.9447631835938
    'sbg:y': 862.5
requirements:
  - class: SubworkflowFeatureRequirement
