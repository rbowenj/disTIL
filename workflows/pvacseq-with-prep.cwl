class: Workflow
cwlVersion: v1.0
id: pvacseq_with_prep
doc: >-
  # About this workflow

  This workflow runs all pVACseq input file preparation steps followed by
  pVACseq prediction of neoepitopes.


  ## Steps

  1. `pvacseq-vcf-prep` - for more details, this subworkflow can be found in the
  repo (includes filtering for pass variants).

  2. `pvacseq` v2.0.3
label: pvacseq-with-prep
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: epitope_lengths_class_i
    type: 'int[]?'
    'sbg:exposed': true
  - id: epitope_lengths_class_ii
    type: 'int[]?'
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
  - id: sample_name
    type: string
    label: Tumour Sample Name
    doc: >-
      The name of the tumour sample in the input VCF to use for neoepitope
      prediction.
    'sbg:x': 0
    'sbg:y': 321.421875
  - id: input_vcf
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    label: Input VCF
    doc: The VCF to be annotated and used for neoepitope prediction.
    secondaryFiles:
      - .tbi
    'sbg:x': 0
    'sbg:y': 642.84375
  - id: ref_genome_dna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: DNA Reference Genome
    doc: Reference sequence used to align the DNA BAM.
    secondaryFiles:
      - .fai
    'sbg:x': 0
    'sbg:y': 535.703125
  - id: ref_genome_rna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: RNA Reference Genome
    doc: Reference sequence used to align the RNA BAM.
    secondaryFiles:
      - .fai
    'sbg:x': 0
    'sbg:y': 428.5625
  - id: transcript_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Transcript Expression File
    doc: Transcript level expression file.
    'sbg:x': 0
    'sbg:y': 214.28125
  - id: input_bam_rna
    'sbg:fileTypes': BAM
    type: File
    label: RNA BAM
    doc: The RNA BAM for the sample of interest.
    secondaryFiles:
      - .bai
    'sbg:x': 0
    'sbg:y': 749.984375
  - id: input_bam_dna
    'sbg:fileTypes': BAM
    type: File
    label: DNA BAM
    doc: The DNA BAM for the sample of interest.
    secondaryFiles:
      - .bai
    'sbg:x': 0
    'sbg:y': 857.125
  - id: gene_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Gene Expression File
    doc: Gene level expression file.
    'sbg:x': 0
    'sbg:y': 964.265625
  - id: id_column
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
  - id: expression_column
    type: string?
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
  - id: min_base_qual
    type: int?
    'sbg:exposed': true
  - id: min_mapping_qual
    type: int?
    'sbg:exposed': true
  - id: vep_plugin_files
    'sbg:fileTypes': PM
    type: 'File[]'
    label: VEP Plugin Files
    doc: VEP plugin files to use in VCF annotation.
    'sbg:x': 0
    'sbg:y': 0
  - id: vep_cache
    'sbg:fileTypes': TAR.GZ
    type: File
    label: VEP Cache
    doc: A VEP cache TAR file.
    'sbg:x': 0
    'sbg:y': 107.140625
  - id: cache_version
    type: int
    label: VEP Cache Version
    doc: The version of the VEP cache file.
    'sbg:x': 0
    'sbg:y': 1071.40625
  - id: alleles
    'sbg:fileTypes': TXT
    type: File
    label: HLA Alleles
    doc: The text file containing HLA alleles for the tumour sample.
    'sbg:x': 254.53125
    'sbg:y': 696.4140625
  - id: additional_report_columns
    type:
      - 'null'
      - type: enum
        symbols:
          - sample_name
        name: additional_report_columns
    'sbg:exposed': true
  - id: allele_specific_binding_thresholds
    type: boolean?
    'sbg:exposed': true
  - id: binding_threshold
    type: int?
    'sbg:exposed': true
  - id: downstream_sequence_length
    type: string?
    'sbg:exposed': true
  - id: exclude_nas
    type: boolean?
    'sbg:exposed': true
  - id: expn_val
    type: float?
    'sbg:exposed': true
  - id: fasta_size
    type: int?
    'sbg:exposed': true
  - id: iedb_retries
    type: int?
    'sbg:exposed': true
  - id: keep_tmp_files
    type: boolean?
    'sbg:exposed': true
  - id: maximum_transcript_support_level
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
  - id: minimum_fold_change
    type: float?
    'sbg:exposed': true
  - id: n_threads
    type: int?
    'sbg:exposed': true
  - id: net_chop_threshold
    type: float?
    'sbg:exposed': true
  - id: netmhc_stab
    type: boolean?
    'sbg:exposed': true
  - id: normal_cov
    type: int?
    'sbg:exposed': true
  - id: normal_sample_name
    type: string?
    label: Normal Sample Name
    doc: Name of the matched normal sample in the input VCF.
    'sbg:x': 254.53125
    'sbg:y': 589.2734375
  - id: phased_proximal_variants_vcf
    type: File?
    label: Phased Proximal Variants VCF
    doc: A VCF of phased proximal variants (see pVACtools documentation).
    'sbg:x': 254.53125
    'sbg:y': 482.1328125
  - id: normal_vaf
    type: float?
    'sbg:exposed': true
  - id: percentile_threshold
    type: int?
    'sbg:exposed': true
  - id: prediction_algorithms
    type: 'string[]'
    'sbg:exposed': true
  - id: run_reference_proteome_similarity
    type: boolean?
    'sbg:exposed': true
  - id: tdna_cov
    type: int?
    'sbg:exposed': true
  - id: tdna_vaf
    type: float?
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
  - id: trna_cov
    type: int?
    'sbg:exposed': true
  - id: trna_vaf
    type: float?
    'sbg:exposed': true
  - id: transcript_id_column
    type: string?
    'sbg:exposed': true
  - id: transcript_expression_column
    type: string?
    'sbg:exposed': true
outputs:
  - id: pvacseq_predictions
    outputSource:
      - pvacseq/pvacseq_predictions
    type: Directory
    label: pVACseq Predictions
    doc: The directory containing all pVACseq output files.
    'sbg:x': 1342.25927734375
    'sbg:y': 214.28125
  - id: combined_aggregated_report
    outputSource:
      - pvacseq/combined_aggregated_report
    type: File?
    label: Combined Aggregated Report
    doc: >-
      The report containing all (MHC-I and/or MHC-II binding) predicted
      neoepitopes from pVACseq.
    'sbg:x': 1342.25927734375
    'sbg:y': 857.125
  - id: mhc_i_aggregated_report
    outputSource:
      - pvacseq/mhc_i_aggregated_report
    type: File?
    label: MHC-I Aggregated Report
    doc: The aggregated report for all MHC-II-binding predicted neoepitopes.
    'sbg:x': 1342.25927734375
    'sbg:y': 642.84375
  - id: mhc_ii_aggregated_report
    outputSource:
      - pvacseq/mhc_ii_aggregated_report
    type: File?
    label: MHC-II Aggregated Report
    doc: The aggregated report for all MHC-II-binding predicted neoepitopes.
    'sbg:x': 1342.25927734375
    'sbg:y': 428.5625
  - id: rc_exp_annotated_vcf
    outputSource:
      - pvacseq_vcf_prep/rc_exp_annotated_vcf
    'sbg:fileTypes': VCF.GZ
    type: File
    label: Fully Annotated VCF
    doc: >-
      The input VCF which has been annotated with VEP, decomposed, then
      annotated with readcount and expression data.
    secondaryFiles:
      - .tbi
    'sbg:x': 841.328125
    'sbg:y': 419.1328125
  - id: mhc_ii_filtered_epitopes
    outputSource:
      - pvacseq/mhc_ii_filtered_epitopes
    type: File?
    label: MHC-II Filtered Neoepitope Report
    doc: Filtered MHC-II-binding predicted neoepitopes.
    'sbg:x': 1342.25927734375
    'sbg:y': 321.421875
  - id: mhc_i_filtered_epitopes
    outputSource:
      - pvacseq/mhc_i_filtered_epitopes
    type: File?
    label: MHC-I Filtered Neoepitope Report
    doc: Filtered MHC-I-binding predicted neoepitopes.
    'sbg:x': 1342.25927734375
    'sbg:y': 535.703125
  - id: combined_filtered_epitopes
    outputSource:
      - pvacseq/combined_filtered_epitopes
    type: File?
    label: Combined Filtered Neoepitope Report
    doc: Filtered MHC-I- and/or MHC-II-binding predicted neoepitopes.
    'sbg:x': 1342.25927734375
    'sbg:y': 749.984375
steps:
  - id: pvacseq
    in:
      - id: additional_report_columns
        source: additional_report_columns
      - id: allele_specific_binding_thresholds
        source: allele_specific_binding_thresholds
      - id: alleles
        source: alleles
      - id: binding_threshold
        source: binding_threshold
      - id: downstream_sequence_length
        source: downstream_sequence_length
      - id: epitope_lengths_class_i
        source:
          - epitope_lengths_class_i
      - id: epitope_lengths_class_ii
        source:
          - epitope_lengths_class_ii
      - id: exclude_nas
        source: exclude_nas
      - id: expn_val
        source: expn_val
      - id: fasta_size
        source: fasta_size
      - id: iedb_retries
        source: iedb_retries
      - id: input_vcf
        source: pvacseq_vcf_prep/rc_exp_annotated_vcf
      - id: keep_tmp_files
        source: keep_tmp_files
      - id: maximum_transcript_support_level
        source: maximum_transcript_support_level
      - id: minimum_fold_change
        source: minimum_fold_change
      - id: n_threads
        source: n_threads
      - id: net_chop_method
        source: net_chop_method
      - id: net_chop_threshold
        source: net_chop_threshold
      - id: netmhc_stab
        source: netmhc_stab
      - id: normal_cov
        source: normal_cov
      - id: normal_sample_name
        source: normal_sample_name
      - id: normal_vaf
        source: normal_vaf
      - id: percentile_threshold
        source: percentile_threshold
      - id: phased_proximal_variants_vcf
        source: phased_proximal_variants_vcf
      - id: prediction_algorithms
        default:
          - all
        source:
          - prediction_algorithms
      - id: run_reference_proteome_similarity
        source: run_reference_proteome_similarity
      - id: sample_name
        source: sample_name
      - id: tdna_cov
        source: tdna_cov
      - id: tdna_vaf
        source: tdna_vaf
      - id: top_score_metric
        source: top_score_metric
      - id: trna_cov
        source: trna_cov
      - id: trna_vaf
        source: trna_vaf
    out:
      - id: combined_aggregated_report
      - id: combined_all_epitopes
      - id: combined_filtered_epitopes
      - id: mhc_i_aggregated_report
      - id: mhc_i_all_epitopes
      - id: mhc_i_filtered_epitopes
      - id: mhc_ii_aggregated_report
      - id: mhc_ii_all_epitopes
      - id: mhc_ii_filtered_epitopes
      - id: pvacseq_predictions
    run: ../tools/pvacseq.cwl
    label: pVACseq
    doc: Run pVACseq to predict SNV/indel-derived neoepitopes.
    'sbg:x': 841.328125
    'sbg:y': 589.2734375
  - id: pvacseq_vcf_prep
    in:
      - id: min_base_qual
        source: min_base_qual
      - id: min_mapping_qual
        source: min_mapping_qual
      - id: sample_name
        source: sample_name
      - id: ref_genome_dna
        source: ref_genome_dna
      - id: ref_genome_rna
        source: ref_genome_rna
      - id: input_vcf
        source: input_vcf
      - id: input_bam_rna
        source: input_bam_rna
      - id: input_bam_dna
        source: input_bam_dna
      - id: id_column
        source: id_column
      - id: expression_column
        source: expression_column
      - id: gene_expression_file
        source: gene_expression_file
      - id: gene_quant_algo
        source: gene_quant_algo
      - id: transcript_expression_file
        source: transcript_expression_file
      - id: transcript_id_column
        source: transcript_id_column
      - id: transcript_expression_column
        source: transcript_expression_column
      - id: transcript_quant_algo
        source: transcript_quant_algo
      - id: vep_cache
        source: vep_cache
      - id: vep_plugin_files
        source:
          - vep_plugin_files
      - id: cache_version
        source: cache_version
    out:
      - id: dna_annot_zipped
      - id: dna_rna_annot_zipped
      - id: snv_readcount_rna
      - id: indel_readcount_rna
      - id: snv_readcount_dna
      - id: indel_readcount_dna
      - id: rc_exp_annotated_vcf
      - id: vep_annot
    run: ./pvacseq-vcf-prep.cwl
    label: pvacseq-vcf-prep
    'sbg:x': 254.53125
    'sbg:y': 304.9921875
requirements:
  - class: SubworkflowFeatureRequirement
