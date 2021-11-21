class: Workflow
cwlVersion: v1.0
id: pvactools_neoepitope_prediction
doc: >-
  # About this workflow

  This workflow runs neoepitope prediction using pVACseq and pVACfuse, including
  all necessary input file preparation steps.  

  pVACtools version 2.0.3 is used in this workflow.


  **Note that if you do not have all the required inputs for this workflow, you
  can run the individual tools and subworkflows which are also in the project
  repo.**


  ## Steps

  ### pVACseq Input File Preparation

  1. Annotates the input VCF using **vep**.

  2. Decomposes the input VCF using **vt-decompose**.

  3. Generates SNV and indel readcounts from a BAM using **bam-readcount** via
  an adapted version of the bam-readcount-helper from mgibio.

  4. Annotates the VCF with SNV readcounts using **vatools
  vcf-readcount-annotator**.

  5. Annotates the VCF with indel readcounts using **vatools
  vcf-readcount-annotator**.

  6. Annotates the VCF with gene level expression values using **vatools
  vcf-expression-annotator**.

  7. Annotates the VCF with transcript level expression values using **vatools
  vcf-expression-annotator**.


  ### pVACfuse Input File Preparation

  1. Annotates the input fusion TSV using **agfusion**.


  ### Neoepitope Prediction

  1. pVACseq for prediction of indel- and SNV-derived neoepitopes.

  2. pVACfuse for prediction of fusion-derived neoepitopes.


  ## Documentation

  - [pvactools](https://pvactools.readthedocs.io/en/latest/)

  - [vep](https://www.ensembl.org/info/docs/tools/vep/index.html)

  - [vt-decompose](https://genome.sph.umich.edu/wiki/Vt)

  - [bam-readcount](https://github.com/genome/bam-readcount)

  -
  [bam-readcount-helper](https://github.com/genome/docker-bam_readcount_helper-cwl)

  - [vatools
  vcf-readcount-annotator](https://vatools.readthedocs.io/en/latest/vcf_readcount_annotator.html)

  - [vatools
  vcf-expression-annotator](https://vatools.readthedocs.io/en/latest/vcf_expression_annotator.html)

  - [agfusion](https://github.com/murphycj/AGFusion)
label: pvactools-neoepitope-prediction
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
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
  - id: fusion_tsv
    'sbg:fileTypes': TSV
    type: File
    label: Input Fusion TSV
    doc: The input TSV of gene fusions called from the tumour sample.
    'sbg:x': 0
    'sbg:y': 1508.9375
  - id: alleles
    'sbg:fileTypes': TXT
    type: File
    label: HLA Alleles
    doc: The text file containing HLA alleles for the tumour sample.
    'sbg:x': 0
    'sbg:y': 1615.71875
  - id: vep_plugin_files
    'sbg:fileTypes': PM
    type: 'File[]'
    label: VEP Plugin Files
    doc: VEP plugin files to use in VCF annotation.
    'sbg:x': 0
    'sbg:y': 227.5625
  - id: vep_cache
    'sbg:fileTypes': TAR.GZ
    type: File
    label: VEP Cache
    doc: A VEP cache TAR file.
    'sbg:x': 0
    'sbg:y': 334.34375
  - id: transcript_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Transcript Expression File
    doc: Transcript level expression file.
    'sbg:x': 0
    'sbg:y': 441.125
  - id: ref_genome_rna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: RNA Reference Genome
    doc: Reference sequence used to align the RNA BAM.
    secondaryFiles:
      - .fai
    'sbg:x': 0
    'sbg:y': 654.6875
  - id: ref_genome_dna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: DNA Reference Genome
    doc: Reference sequence used to align the DNA BAM.
    secondaryFiles:
      - .fai
    'sbg:x': 0
    'sbg:y': 761.46875
  - id: input_vcf
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    label: Input VCF
    doc: The VCF to be annotated and used for neoepitope prediction.
    secondaryFiles:
      - .tbi
    'sbg:x': 0
    'sbg:y': 1081.8125
  - id: input_bam_rna
    'sbg:fileTypes': BAM
    type: File
    label: RNA BAM
    doc: The RNA BAM for the sample of interest.
    secondaryFiles:
      - .bai
    'sbg:x': 0
    'sbg:y': 1188.59375
  - id: input_bam_dna
    'sbg:fileTypes': BAM
    type: File
    label: DNA BAM
    doc: The DNA BAM for the sample of interest.
    secondaryFiles:
      - .bai
    'sbg:x': 0
    'sbg:y': 1295.375
  - id: gene_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Gene Expression File
    doc: Gene level expression file.
    'sbg:x': 0
    'sbg:y': 1402.15625
  - id: sample_name
    type: string
    label: Tumour Sample Name
    doc: >-
      The name of the tumour sample in the input VCF to use for neoepitope
      prediction.
    'sbg:x': 0
    'sbg:y': 547.90625
  - id: normal_sample_name
    type: string?
    label: Normal Sample Name
    doc: The name of the matched normal sample in the input VCF.
    'sbg:x': 0
    'sbg:y': 975.03125
  - id: phased_proximal_variants_vcf
    type: File?
    label: Phased Proximal Variants VCF
    doc: VCF of phased proximal variants (see pVACseq documentation).
    'sbg:x': 0
    'sbg:y': 868.25
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
  - id: n_junction_reads
    type: int?
    'sbg:exposed': true
  - id: n_spanning_frags
    type: int?
    'sbg:exposed': true
outputs:
  - id: rc_exp_annotated_vcf
    outputSource:
      - pvacseq_with_prep/rc_exp_annotated_vcf
    'sbg:fileTypes': VCF.GZ
    type: File
    label: Annotated VCF
    doc: 'Input VCF annotated with VEP, readcounts, and expression data.'
    'sbg:x': 893.3243408203125
    'sbg:y': 0
  - id: pvacseq_predictions
    outputSource:
      - pvacseq_with_prep/pvacseq_predictions
    type: Directory
    label: pVACseq Predictions
    doc: The directory containing all pVACseq output files.
    'sbg:x': 893.3243408203125
    'sbg:y': 106.78125
  - id: pvacseq_mhc_ii_aggregated_report
    outputSource:
      - pvacseq_with_prep/mhc_ii_aggregated_report
    type: File?
    label: pVACseq MHC-II Aggregated Report
    doc: >-
      The aggregated report for all pVACseq MHC-II-binding predicted
      neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 213.5625
  - id: pvacseq_mhc_i_aggregated_report
    outputSource:
      - pvacseq_with_prep/mhc_i_aggregated_report
    type: File?
    label: pVACseq MHC-I Aggregated Report
    doc: The aggregated report for all pVACseq MHC-I-binding predicted neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 320.34375
  - id: pvacseq_combined_aggregated_report
    outputSource:
      - pvacseq_with_prep/combined_aggregated_report
    type: File?
    label: pVACseq Combined Aggregated Report
    doc: >-
      The report containing all (MHC-I and/or MHC-II binding) predicted
      neoepitopes from pVACseq.
    'sbg:x': 893.3243408203125
    'sbg:y': 427.125
  - id: pvacfuse_predictions
    outputSource:
      - pvacfuse_with_prep/pvacfuse_predictions
    type: Directory?
    label: pVACfuse Predictions
    doc: The directory containing all pVACfuse output files.
    'sbg:x': 893.3243408203125
    'sbg:y': 775.46875
  - id: pvacfuse_mhc_ii_aggregated_report
    outputSource:
      - pvacfuse_with_prep/mhc_ii_aggregated_report
    type: File?
    label: pVACfuse MHC-II Aggregated Report
    doc: >-
      The aggregated report for all pVACfuse MHC-II-binding predicted
      neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 882.25
  - id: pvacfuse_mhc_i_aggregated_report
    outputSource:
      - pvacfuse_with_prep/mhc_i_aggregated_report
    type: File?
    label: pVACfuse MHC-I Aggregated Report
    doc: >-
      The aggregated report for all pVACfuse MHC-I-binding predicted
      neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 989.03125
  - id: pvacfuse_combined_aggregated_report
    outputSource:
      - pvacfuse_with_prep/combined_aggregated_report
    type: File?
    label: pVACfuse Combined Aggregated Report
    doc: >-
      The report containing all (MHC-I and/or MHC-II binding) predicted
      neoepitopes from pVACfuse.
    'sbg:x': 893.3243408203125
    'sbg:y': 1095.8125
  - id: mhc_ii_filtered_epitopes
    outputSource:
      - pvacseq_with_prep/mhc_ii_filtered_epitopes
    type: File?
    label: pVACseq MHC-II Filtered Neoepitope Report
    doc: Filtered MHC-II-binding predicted neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 1309.375
  - id: mhc_i_filtered_epitopes
    outputSource:
      - pvacseq_with_prep/mhc_i_filtered_epitopes
    type: File?
    label: pVACseq MHC-I Filtered Neoepitope Report
    doc: Filtered MHC-I-binding predicted neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 1522.9375
  - id: combined_filtered_epitopes
    outputSource:
      - pvacseq_with_prep/combined_filtered_epitopes
    type: File?
    label: pVACseq Combined Filtered Neoepitope Report
    doc: Filtered MHC-I- and/or MHC-II-binding predicted neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 1736.5
  - id: agfusion_output
    outputSource:
      - pvacfuse_with_prep/agfusion_output
    type: Directory?
    label: AGFusion Output
    doc: Output directory from AGFusion.
    'sbg:x': 893.3243408203125
    'sbg:y': 1843.28125
  - id: mhc_i_filtered_epitopes_1
    outputSource:
      - pvacfuse_with_prep/mhc_i_filtered_epitopes
    type: File?
    label: pVACfuse MHC-I Filtered Neoepitope Report
    doc: Filtered MHC-I-binding predicted neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 1416.15625
  - id: mhc_ii_filtered_epitopes_1
    outputSource:
      - pvacfuse_with_prep/mhc_ii_filtered_epitopes
    type: File?
    label: pVACfuse MHC-II Filtered Neoepitope Report
    doc: Filtered MHC-II-binding predicted neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 1202.59375
  - id: combined_filtered_epitopes_1
    outputSource:
      - pvacfuse_with_prep/combined_filtered_epitopes
    type: File?
    label: pVACfuse Combined Filtered Neoepitope Report
    doc: Filtered MHC-I- and/or MHC-II-binding predicted neoepitopes.
    'sbg:x': 893.3243408203125
    'sbg:y': 1629.71875
  - id: pvacfuse_star_filtered_mhc_ii
    outputSource:
      - pvacfuse_star_filter_mhc_ii/pvacfuse_star_filtered
    'sbg:fileTypes': TSV
    type: File
    label: pVACfuse MHC-II Star-Annotated Filtered Report
    'sbg:x': 1375.116455078125
    'sbg:y': 868.25
  - id: pvacfuse_star_filtered_mhc_i
    outputSource:
      - pvacfuse_star_filter_mhc_i/pvacfuse_star_filtered
    'sbg:fileTypes': TSV
    type: File
    label: pVACfuse MHC-I Star-Annotated Filtered Report
    'sbg:x': 1375.116455078125
    'sbg:y': 975.03125
steps:
  - id: pvacfuse_with_prep
    in:
      - id: epitope_lengths_class_i
        source:
          - pvacfuse_epitope_lengths_class_i
      - id: epitope_lengths_class_ii
        source:
          - pvacfuse_epitope_lengths_class_ii
      - id: binding_threshold
        source: pvacfuse_binding_threshold
      - id: percentile_threshold
        source: pvacfuse_percentile_threshold
      - id: iedb_retries
        source: pvacfuse_iedb_retries
      - id: keep_tmp_files
        source: pvacfuse_keep_tmp_files
      - id: net_chop_method
        source: pvacfuse_net_chop_method
      - id: netmhc_stab
        default: false
        source: pvacfuse_netmhc_stab
      - id: top_score_metric
        source: pvacfuse_top_score_metric
      - id: net_chop_threshold
        source: pvacfuse_net_chop_threshold
      - id: run_reference_proteome_similarity
        source: pvacfuse_run_reference_proteome_similarity
      - id: additional_report_columns
        source: pvacfuse_additional_report_columns
      - id: fasta_size
        source: pvacfuse_fasta_size
      - id: downstream_sequence_length
        source: pvacfuse_downstream_sequence_length
      - id: exclude_nas
        source: pvacfuse_exclude_nas
      - id: n_threads
        source: n_threads
      - id: alleles
        source: alleles
      - id: sample_name
        source: sample_name
      - id: ref_genome
        source: ref_genome
      - id: fusion_tsv
        source: fusion_tsv
      - id: fusion_caller
        default: starfusion
        source: fusion_caller
      - id: prediction_algorithms
        source:
          - pvacfuse_prediction_algorithms
    out:
      - id: pvacfuse_predictions
      - id: mhc_ii_aggregated_report
      - id: mhc_i_aggregated_report
      - id: combined_aggregated_report
      - id: agfusion_output
      - id: mhc_ii_filtered_epitopes
      - id: mhc_i_filtered_epitopes
      - id: combined_filtered_epitopes
    run: ./pvacfuse-with-prep.cwl
    label: pvacfuse-with-prep
    doc: Predict fusion-derived neoepitopes with pVACfuse.
    'sbg:x': 294.53125
    'sbg:y': 975.03125
  - id: pvacseq_with_prep
    in:
      - id: epitope_lengths_class_i
        source:
          - pvacseq_epitope_lengths_class_i
      - id: epitope_lengths_class_ii
        source:
          - pvacseq_epitope_lengths_class_ii
      - id: net_chop_method
        source: pvacseq_net_chop_method
      - id: sample_name
        source: sample_name
      - id: input_vcf
        source: input_vcf
      - id: ref_genome_dna
        source: ref_genome_dna
      - id: ref_genome_rna
        source: ref_genome_rna
      - id: transcript_expression_file
        source: transcript_expression_file
      - id: input_bam_rna
        source: input_bam_rna
      - id: input_bam_dna
        source: input_bam_dna
      - id: gene_expression_file
        source: gene_expression_file
      - id: id_column
        source: gene_id_column
      - id: gene_quant_algo
        source: gene_quant_algo
      - id: expression_column
        source: gene_expression_column
      - id: transcript_quant_algo
        source: transcript_quant_algo
      - id: min_base_qual
        source: pvacseq_min_base_qual
      - id: min_mapping_qual
        source: pvacseq_min_mapping_qual
      - id: vep_plugin_files
        source:
          - vep_plugin_files
      - id: vep_cache
        source: vep_cache
      - id: cache_version
        source: cache_version
      - id: alleles
        source: alleles
      - id: additional_report_columns
        source: pvacseq_additional_report_columns
      - id: allele_specific_binding_thresholds
        source: pvacseq_allele_specific_binding_thresholds
      - id: binding_threshold
        source: pvacseq_binding_threshold
      - id: downstream_sequence_length
        source: pvacseq_downstream_sequence_length
      - id: exclude_nas
        source: pvacseq_exclude_nas
      - id: expn_val
        source: pvacseq_expn_val
      - id: fasta_size
        source: pvacseq_fasta_size
      - id: iedb_retries
        source: pvacseq_iedb_retries
      - id: keep_tmp_files
        source: pvacseq_keep_tmp_files
      - id: maximum_transcript_support_level
        source: pvacseq_maximum_transcript_support_level
      - id: minimum_fold_change
        source: pvacseq_minimum_fold_change
      - id: n_threads
        source: n_threads
      - id: net_chop_threshold
        source: pvacseq_net_chop_threshold
      - id: netmhc_stab
        source: pvacseq_netmhc_stab
      - id: normal_cov
        source: normal_cov
      - id: normal_sample_name
        source: normal_sample_name
      - id: phased_proximal_variants_vcf
        source: phased_proximal_variants_vcf
      - id: normal_vaf
        source: normal_vaf
      - id: percentile_threshold
        source: pvacseq_percentile_threshold
      - id: prediction_algorithms
        source:
          - pvacseq_prediction_algorithms
      - id: run_reference_proteome_similarity
        source: pvacseq_run_reference_proteome_similarity
      - id: tdna_cov
        source: tdna_cov
      - id: tdna_vaf
        source: tdna_vaf
      - id: top_score_metric
        source: pvacseq_top_score_metric
      - id: trna_cov
        source: trna_cov
      - id: trna_vaf
        source: trna_vaf
      - id: transcript_id_column
        source: transcript_id_column
      - id: transcript_expression_column
        source: transcript_expression_column
    out:
      - id: pvacseq_predictions
      - id: combined_aggregated_report
      - id: mhc_i_aggregated_report
      - id: mhc_ii_aggregated_report
      - id: rc_exp_annotated_vcf
      - id: mhc_ii_filtered_epitopes
      - id: mhc_i_filtered_epitopes
      - id: combined_filtered_epitopes
    run: ./pvacseq-with-prep.cwl
    label: pvacseq-with-prep
    doc: Predict indel- and SNV-derived neoepitopes with pVACseq.
    'sbg:x': 294.53125
    'sbg:y': 735.25
  - id: pvacfuse_star_filter_mhc_ii
    in:
      - id: pvacfuse_report
        source: pvacfuse_with_prep/mhc_ii_filtered_epitopes
      - id: star_fusions
        source: fusion_tsv
      - id: n_junction_reads
        source: n_junction_reads
      - id: n_spanning_frags
        source: n_spanning_frags
    out:
      - id: pvacfuse_star_filtered
    run: ../tools/pvacfuse-star-filter.cwl
    label: pvacfuse-star-filter-mhc-ii
    'sbg:x': 893.3243408203125
    'sbg:y': 540.90625
  - id: pvacfuse_star_filter_mhc_i
    in:
      - id: pvacfuse_report
        source: pvacfuse_with_prep/mhc_i_filtered_epitopes
      - id: star_fusions
        source: fusion_tsv
      - id: n_junction_reads
        source: n_junction_reads
      - id: n_spanning_frags
        source: n_spanning_frags
    out:
      - id: pvacfuse_star_filtered
    run: ../tools/pvacfuse-star-filter.cwl
    label: pvacfuse-star-filter-mhc-i
    'sbg:x': 893.3243408203125
    'sbg:y': 661.6875
requirements:
  - class: SubworkflowFeatureRequirement
'sbg:toolAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
