class: Workflow
cwlVersion: v1.0
id: hlahd_pvactools
label: hlahd-pvactools
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: tumour_rna_reads_2
    type: File
    'sbg:x': 0
    'sbg:y': 906.46875
  - id: tumour_rna_reads_1
    type: File
    'sbg:x': 0
    'sbg:y': 1013.0625
  - id: tumour_dna_reads
    'sbg:fileTypes': >-
      FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA, FQ, FASTQ.GZ, FQ.GZ,
      FASTQ.BZ2, FQ.BZ2
    type: 'File[]'
    'sbg:x': 0
    'sbg:y': 1226.25
  - id: normal_dna_reads
    'sbg:fileTypes': >-
      FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA, FQ, FASTQ.GZ, FQ.GZ,
      FASTQ.BZ2, FQ.BZ2
    type: 'File[]'
    'sbg:x': 0
    'sbg:y': 1652.765625
  - id: bowtie_index_archive
    'sbg:fileTypes': TAR
    type: File
    'sbg:x': 0
    'sbg:y': 1865.953125
  - id: minimum_read_length
    type: int
    'sbg:exposed': true
  - id: minimum_read_length_1
    type: int
    'sbg:exposed': true
  - id: minimum_read_length_2
    type: int
    'sbg:exposed': true
  - id: normal_dna_outdir_name
    type: string
    'sbg:x': 0
    'sbg:y': 1759.359375
  - id: tumour_dna_outdir_name
    type: string
    'sbg:x': 0
    'sbg:y': 1332.84375
  - id: tumour_rna_outdir_name
    type: string
    'sbg:x': 0
    'sbg:y': 1119.65625
  - id: sample_id
    type: string
    'sbg:x': 0
    'sbg:y': 1439.4375
  - id: epitope_lengths_class_i_2
    type: 'int[]?'
    'sbg:exposed': true
  - id: epitope_lengths_class_ii_2
    type: 'int[]?'
    'sbg:exposed': true
  - id: binding_threshold_2
    type: int?
    'sbg:exposed': true
  - id: percentile_threshold_2
    type: int?
    'sbg:exposed': true
  - id: iedb_retries_2
    type: int?
    'sbg:exposed': true
  - id: keep_tmp_files_2
    type: boolean?
    'sbg:exposed': true
  - id: net_chop_method_2
    type:
      - 'null'
      - type: enum
        symbols:
          - cterm
          - 20s
        name: net_chop_method
    'sbg:exposed': true
  - id: netmhc_stab_2
    type: boolean?
    'sbg:exposed': true
  - id: top_score_metric_2
    type:
      - 'null'
      - type: enum
        symbols:
          - lowest
          - median
        name: top_score_metric
    'sbg:exposed': true
  - id: net_chop_threshold_2
    type: float?
    'sbg:exposed': true
  - id: additional_report_columns_2
    type:
      - 'null'
      - type: enum
        symbols:
          - sample_name
        name: additional_report_columns
    'sbg:exposed': true
  - id: run_reference_proteome_similarity_2
    type: boolean?
    'sbg:exposed': true
  - id: fasta_size_2
    type: int?
    'sbg:exposed': true
  - id: downstream_sequence_length_2
    type: string?
    'sbg:exposed': true
  - id: exclude_nas_2
    type: boolean?
    'sbg:exposed': true
  - id: n_threads_2
    type: int?
    'sbg:exposed': true
  - id: ref_genome
    type:
      - 'null'
      - type: enum
        symbols:
          - hg38
          - hg19
        name: ref_genome
    'sbg:exposed': true
  - id: sample_name_1
    type: string
    'sbg:exposed': true
  - id: cache_version_1
    type: int
    'sbg:exposed': true
  - id: merged_1
    type: boolean
    'sbg:exposed': true
  - id: symbol_1
    type: boolean?
    'sbg:exposed': true
  - id: biotype_1
    type: boolean?
    'sbg:exposed': true
  - id: numbers_1
    type: boolean?
    'sbg:exposed': true
  - id: canonical_1
    type: boolean?
    'sbg:exposed': true
  - id: total_length_1
    type: boolean?
    'sbg:exposed': true
  - id: sift_1
    type:
      - 'null'
      - type: enum
        symbols:
          - p
          - s
          - b
        name: sift
    'sbg:exposed': true
  - id: polyphen_1
    type:
      - 'null'
      - type: enum
        symbols:
          - p
          - s
          - b
        name: polyphen
    'sbg:exposed': true
  - id: terms_1
    type:
      - 'null'
      - type: enum
        symbols:
          - SO
          - display
          - NCBI
        name: terms
    'sbg:exposed': true
  - id: epitope_lengths_class_i_3
    type: 'int[]?'
    'sbg:exposed': true
  - id: epitope_lengths_class_ii_3
    type: 'int[]?'
    'sbg:exposed': true
  - id: binding_threshold_3
    type: int?
    'sbg:exposed': true
  - id: percentile_threshold_3
    type: int?
    'sbg:exposed': true
  - id: allele_specific_binding_thresholds_1
    type: boolean?
    'sbg:exposed': true
  - id: iedb_retries_3
    type: int?
    'sbg:exposed': true
  - id: keep_tmp_files_3
    type: boolean?
    'sbg:exposed': true
  - id: normal_sample_name_1
    type: string?
    'sbg:exposed': true
  - id: net_chop_method_3
    type:
      - 'null'
      - type: enum
        symbols:
          - cterm
          - 20s
        name: net_chop_method
    'sbg:exposed': true
  - id: netmhc_stab_3
    type: boolean?
    'sbg:exposed': true
  - id: run_reference_proteome_similarity_3
    type: boolean?
    'sbg:exposed': true
  - id: top_score_metric_3
    type:
      - 'null'
      - type: enum
        symbols:
          - lowest
          - median
        name: top_score_metric
    'sbg:exposed': true
  - id: net_chop_threshold_3
    type: float?
    'sbg:exposed': true
  - id: additional_report_columns_3
    type:
      - 'null'
      - type: enum
        symbols:
          - sample_name
        name: additional_report_columns
    'sbg:exposed': true
  - id: fasta_size_3
    type: int?
    'sbg:exposed': true
  - id: downstream_sequence_length_3
    type: string?
    'sbg:exposed': true
  - id: exclude_nas_3
    type: boolean?
    'sbg:exposed': true
  - id: minimum_fold_change_1
    type: float?
    'sbg:exposed': true
  - id: normal_cov_1
    type: int?
    'sbg:exposed': true
  - id: tdna_cov_1
    type: int?
    'sbg:exposed': true
  - id: trna_cov_1
    type: int?
    'sbg:exposed': true
  - id: normal_vaf_1
    type: float?
    'sbg:exposed': true
  - id: tdna_vaf_1
    type: float?
    'sbg:exposed': true
  - id: trna_vaf_1
    type: float?
    'sbg:exposed': true
  - id: expn_val_1
    type: float?
    'sbg:exposed': true
  - id: maximum_transcript_support_level_1
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
  - id: n_threads_3
    type: int?
    'sbg:exposed': true
  - id: fusion_caller_1
    type:
      - 'null'
      - type: enum
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
    'sbg:x': 952.7314453125
    'sbg:y': 1886.8125
  - id: fusion_tsv_1
    'sbg:fileTypes': tsv
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 1780.21875
  - id: input_file_1
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    'sbg:x': 952.7314453125
    'sbg:y': 1098.796875
  - id: phased_proximal_variants_vcf_1
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 779.015625
  - id: prediction_algorithms_2
    type: 'string[]'
    'sbg:x': 952.7314453125
    'sbg:y': 672.3515625
  - id: ref_genome_1
    'sbg:fileTypes': 'FASTA, FA'
    type: File
    'sbg:x': 0
    'sbg:y': 1546.1015625
  - id: vep_cache_1
    'sbg:fileTypes': TAR.GZ
    type: File
    'sbg:x': 0
    'sbg:y': 799.875
  - id: vep_plugin_files_1
    'sbg:fileTypes': PM
    type: 'File[]?'
    'sbg:x': 0
    'sbg:y': 693.2109375
outputs:
  - id: tumour_rna_hla-hd_results
    outputSource:
      - hlatyping_bowtie_samtools/tumour_rna_hla-hd_results
    type: Directory
    'sbg:x': 952.7314453125
    'sbg:y': 245.90625
  - id: hlad_reads1
    outputSource:
      - hlatyping_bowtie_samtools/hlad_reads1
    'sbg:fileTypes': FASTQ
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 1460.4375
  - id: hla_reads2
    outputSource:
      - hlatyping_bowtie_samtools/hla_reads2
    'sbg:fileTypes': FASTQ
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 1673.625
  - id: filtered_BAM
    outputSource:
      - hlatyping_bowtie_samtools/filtered_BAM
    'sbg:fileTypes': 'BAM, SAM, CRAM'
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 2100
  - id: bowtie_sam
    outputSource:
      - hlatyping_bowtie_samtools/bowtie_sam
    'sbg:fileTypes': SAM
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 2313.1875
  - id: aligned_reads_only
    outputSource:
      - hlatyping_bowtie_samtools/aligned_reads_only
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
    type: 'File[]?'
    'sbg:x': 952.7314453125
    'sbg:y': 2526.5859375
  - id: unaligned_reads_only
    outputSource:
      - hlatyping_bowtie_samtools/unaligned_reads_only
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
    type: 'File[]?'
    'sbg:x': 952.7314453125
    'sbg:y': 139.2421875
  - id: tumour_dna_hla-hd_results
    outputSource:
      - hlatyping_bowtie_samtools/tumour_dna_hla-hd_results
    type: Directory
    'sbg:x': 952.7314453125
    'sbg:y': 459.09375
  - id: hlad_reads2
    outputSource:
      - hlatyping_bowtie_samtools/hlad_reads2
    'sbg:fileTypes': FASTQ
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 1353.84375
  - id: hla_reads3
    outputSource:
      - hlatyping_bowtie_samtools/hla_reads3
    'sbg:fileTypes': FASTQ
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 1567.03125
  - id: filtered_BAM_1
    outputSource:
      - hlatyping_bowtie_samtools/filtered_BAM_1
    'sbg:fileTypes': 'BAM, SAM, CRAM'
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 1993.40625
  - id: bowtie_sam_1
    outputSource:
      - hlatyping_bowtie_samtools/bowtie_sam_1
    'sbg:fileTypes': SAM
    type: File?
    'sbg:x': 952.7314453125
    'sbg:y': 2206.59375
  - id: aligned_reads_only_1
    outputSource:
      - hlatyping_bowtie_samtools/aligned_reads_only_1
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
    type: 'File[]?'
    'sbg:x': 952.7314453125
    'sbg:y': 2419.8515625
  - id: unaligned_reads_only_1
    outputSource:
      - hlatyping_bowtie_samtools/unaligned_reads_only_1
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
    type: 'File[]?'
    'sbg:x': 952.7314453125
    'sbg:y': 32.5078125
  - id: normal_dna_hla-hd_results
    outputSource:
      - hlatyping_bowtie_samtools/normal_dna_hla-hd_results
    type: Directory
    'sbg:x': 952.7314453125
    'sbg:y': 885.609375
  - id: tumour_dna_hla-hd_final_result
    outputSource:
      - hlatyping_bowtie_samtools/tumour_dna_hla-hd_final_result
    type: File
    'sbg:x': 952.7314453125
    'sbg:y': 565.6875
  - id: tumour_rna_hla-hd_final
    outputSource:
      - hlatyping_bowtie_samtools/tumour_rna_hla-hd_final
    type: File
    'sbg:x': 952.7314453125
    'sbg:y': 352.5
  - id: normal_dna_hla-hd_final_result
    outputSource:
      - hlatyping_bowtie_samtools/normal_dna_hla-hd_final_result
    type: File
    'sbg:x': 952.7314453125
    'sbg:y': 992.203125
  - id: consensus_txt
    outputSource:
      - hlahd_consensus/consensus_txt
    type: File
    'sbg:x': 1265.7939453125
    'sbg:y': 1279.546875
  - id: consensus_json
    outputSource:
      - hlahd_consensus/consensus_json
    type: File
    'sbg:x': 1265.7939453125
    'sbg:y': 1386.2109375
  - id: output_gz_tbi
    outputSource:
      - pvacseq_pvacfuse/output_gz_tbi
    type: File
    'sbg:x': 2004.48388671875
    'sbg:y': 213.2578125
  - id: vep_vcf
    outputSource:
      - pvacseq_pvacfuse/vep_vcf
    'sbg:fileTypes': VCF
    type: File
    'sbg:x': 2004.48388671875
    'sbg:y': 0
  - id: output_gz
    outputSource:
      - pvacseq_pvacfuse/output_gz
    type: File
    'sbg:x': 2004.48388671875
    'sbg:y': 319.9921875
  - id: mhc_i_all_epitopes
    outputSource:
      - pvacseq_pvacfuse/mhc_i_all_epitopes
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 1492.875
  - id: mhc_i_aggregated_report
    outputSource:
      - pvacseq_pvacfuse/mhc_i_aggregated_report
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 1706.2734375
  - id: mhc_i_filtered_epitopes
    outputSource:
      - pvacseq_pvacfuse/mhc_i_filtered_epitopes
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 1279.6875
  - id: mhc_ii_all_epitopes
    outputSource:
      - pvacseq_pvacfuse/mhc_ii_all_epitopes
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 746.4375
  - id: mhc_ii_aggregated_report
    outputSource:
      - pvacseq_pvacfuse/mhc_ii_aggregated_report
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 959.8359375
  - id: mhc_ii_filtered_epitopes
    outputSource:
      - pvacseq_pvacfuse/mhc_ii_filtered_epitopes
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 533.25
  - id: combined_all_epitopes
    outputSource:
      - pvacseq_pvacfuse/combined_all_epitopes
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 2239.3125
  - id: combined_aggregated_report
    outputSource:
      - pvacseq_pvacfuse/combined_aggregated_report
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 2452.7109375
  - id: combined_filtered_epitopes
    outputSource:
      - pvacseq_pvacfuse/combined_filtered_epitopes
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 2026.125
  - id: pvacseq_predictions
    outputSource:
      - pvacseq_pvacfuse/pvacseq_predictions
    type: Directory
    'sbg:x': 2004.48388671875
    'sbg:y': 106.59375
  - id: mhc_i
    outputSource:
      - pvacseq_pvacfuse/mhc_i
    type: Directory?
    'sbg:x': 2004.48388671875
    'sbg:y': 1812.9375
  - id: mhc_ii
    outputSource:
      - pvacseq_pvacfuse/mhc_ii
    type: Directory?
    'sbg:x': 2004.48388671875
    'sbg:y': 1066.5
  - id: combined
    outputSource:
      - pvacseq_pvacfuse/combined
    type: Directory?
    'sbg:x': 2004.48388671875
    'sbg:y': 2559.3046875
  - id: mhc_ii_filtered_epitopes_1
    outputSource:
      - pvacseq_pvacfuse/mhc_ii_filtered_epitopes_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 426.65625
  - id: mhc_ii_all_epitopes_1
    outputSource:
      - pvacseq_pvacfuse/mhc_ii_all_epitopes_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 639.84375
  - id: mhc_ii_aggregated_report_1
    outputSource:
      - pvacseq_pvacfuse/mhc_ii_aggregated_report_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 853.1015625
  - id: mhc_i_filtered_epitopes_1
    outputSource:
      - pvacseq_pvacfuse/mhc_i_filtered_epitopes_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 1173.09375
  - id: mhc_i_all_epitopes_1
    outputSource:
      - pvacseq_pvacfuse/mhc_i_all_epitopes_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 1386.28125
  - id: mhc_i_aggregated_report_1
    outputSource:
      - pvacseq_pvacfuse/mhc_i_aggregated_report_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 1599.5390625
  - id: combined_filtered_epitopes_1
    outputSource:
      - pvacseq_pvacfuse/combined_filtered_epitopes_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 1919.53125
  - id: combined_all_epitopes_1
    outputSource:
      - pvacseq_pvacfuse/combined_all_epitopes_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 2132.71875
  - id: combined_aggregated_report_1
    outputSource:
      - pvacseq_pvacfuse/combined_aggregated_report_1
    type: File?
    'sbg:x': 2004.48388671875
    'sbg:y': 2345.9765625
steps:
  - id: pvacseq_pvacfuse
    in:
      - id: epitope_lengths_class_i_2
        source:
          - epitope_lengths_class_i_2
      - id: epitope_lengths_class_ii_2
        source:
          - epitope_lengths_class_ii_2
      - id: binding_threshold_2
        source: binding_threshold_2
      - id: percentile_threshold_2
        source: percentile_threshold_2
      - id: iedb_retries_2
        source: iedb_retries_2
      - id: keep_tmp_files_2
        source: keep_tmp_files_2
      - id: net_chop_method_2
        source: net_chop_method_2
      - id: netmhc_stab_2
        source: netmhc_stab_2
      - id: top_score_metric_2
        source: top_score_metric_2
      - id: net_chop_threshold_2
        source: net_chop_threshold_2
      - id: run_reference_proteome_similarity_2
        source: run_reference_proteome_similarity_2
      - id: additional_report_columns_2
        source: additional_report_columns_2
      - id: fasta_size_2
        source: fasta_size_2
      - id: downstream_sequence_length_2
        source: downstream_sequence_length_2
      - id: exclude_nas_2
        source: exclude_nas_2
      - id: n_threads_2
        source: n_threads_2
      - id: ref_genome
        source: ref_genome
      - id: fusion_caller_1
        source: fusion_caller_1
      - id: alleles_1
        source: hlahd_consensus/consensus_string
      - id: fusion_tsv_1
        source: fusion_tsv_1
      - id: prediction_algorithms_2
        source:
          - prediction_algorithms_2
      - id: sample_name_1
        source: sample_name_1
      - id: input_file_1
        source: input_file_1
      - id: ref_genome_1
        source: ref_genome_1
      - id: vep_cache_1
        source: vep_cache_1
      - id: cache_version_1
        source: cache_version_1
      - id: merged_1
        source: merged_1
      - id: symbol_1
        source: symbol_1
      - id: biotype_1
        source: biotype_1
      - id: numbers_1
        source: numbers_1
      - id: canonical_1
        source: canonical_1
      - id: total_length_1
        source: total_length_1
      - id: sift_1
        source: sift_1
      - id: polyphen_1
        source: polyphen_1
      - id: terms_1
        source: terms_1
      - id: epitope_lengths_class_i_3
        source:
          - epitope_lengths_class_i_3
      - id: epitope_lengths_class_ii_3
        source:
          - epitope_lengths_class_ii_3
      - id: binding_threshold_3
        source: binding_threshold_3
      - id: percentile_threshold_3
        source: percentile_threshold_3
      - id: allele_specific_binding_thresholds_1
        source: allele_specific_binding_thresholds_1
      - id: iedb_retries_3
        source: iedb_retries_3
      - id: keep_tmp_files_3
        source: keep_tmp_files_3
      - id: normal_sample_name_1
        source: normal_sample_name_1
      - id: net_chop_method_3
        source: net_chop_method_3
      - id: netmhc_stab_3
        source: netmhc_stab_3
      - id: run_reference_proteome_similarity_3
        source: run_reference_proteome_similarity_3
      - id: top_score_metric_3
        source: top_score_metric_3
      - id: net_chop_threshold_3
        source: net_chop_threshold_3
      - id: additional_report_columns_3
        source: additional_report_columns_3
      - id: fasta_size_3
        source: fasta_size_3
      - id: downstream_sequence_length_3
        source: downstream_sequence_length_3
      - id: exclude_nas_3
        source: exclude_nas_3
      - id: minimum_fold_change_1
        source: minimum_fold_change_1
      - id: normal_cov_1
        source: normal_cov_1
      - id: tdna_cov_1
        source: tdna_cov_1
      - id: trna_cov_1
        source: trna_cov_1
      - id: normal_vaf_1
        source: normal_vaf_1
      - id: tdna_vaf_1
        source: tdna_vaf_1
      - id: trna_vaf_1
        source: trna_vaf_1
      - id: expn_val_1
        source: expn_val_1
      - id: maximum_transcript_support_level_1
        source: maximum_transcript_support_level_1
      - id: n_threads_3
        source: n_threads_3
      - id: vep_plugin_files_1
        source:
          - vep_plugin_files_1
      - id: prediction_algorithms_3
        source:
          - prediction_algorithms_2
      - id: phased_proximal_variants_vcf_1
        source: phased_proximal_variants_vcf_1
    out:
      - id: output_gz_tbi
      - id: vep_vcf
      - id: output_gz
      - id: mhc_i_all_epitopes
      - id: mhc_i_aggregated_report
      - id: mhc_i_filtered_epitopes
      - id: mhc_ii_all_epitopes
      - id: mhc_ii_aggregated_report
      - id: mhc_ii_filtered_epitopes
      - id: combined_all_epitopes
      - id: combined_aggregated_report
      - id: combined_filtered_epitopes
      - id: pvacseq_predictions
      - id: mhc_i_all_epitopes_1
      - id: mhc_i_aggregated_report_1
      - id: mhc_i_filtered_epitopes_1
      - id: mhc_ii_all_epitopes_1
      - id: mhc_ii_aggregated_report_1
      - id: mhc_ii_filtered_epitopes_1
      - id: combined_all_epitopes_1
      - id: combined_aggregated_report_1
      - id: combined_filtered_epitopes_1
      - id: mhc_i
      - id: mhc_ii
      - id: combined
    run:
      class: Workflow
      cwlVersion: v1.0
      id: mwonge/mwtest/pvacseq-pvacfuse/4
      label: pvacseq-pvacfuse
      $namespaces:
        sbg: 'https://sevenbridges.com'
      inputs:
        - id: epitope_lengths_class_i_2
          type: 'int[]?'
          'sbg:exposed': true
        - id: epitope_lengths_class_ii_2
          type: 'int[]?'
          'sbg:exposed': true
        - id: binding_threshold_2
          type: int?
          'sbg:exposed': true
        - id: percentile_threshold_2
          type: int?
          'sbg:exposed': true
        - id: iedb_retries_2
          type: int?
          'sbg:exposed': true
        - id: keep_tmp_files_2
          type: boolean?
          'sbg:exposed': true
        - id: net_chop_method_2
          type:
            - 'null'
            - type: enum
              symbols:
                - cterm
                - 20s
              name: net_chop_method
          'sbg:exposed': true
        - id: netmhc_stab_2
          type: boolean?
          'sbg:exposed': true
        - id: top_score_metric_2
          type:
            - 'null'
            - type: enum
              symbols:
                - lowest
                - median
              name: top_score_metric
          'sbg:exposed': true
        - id: net_chop_threshold_2
          type: float?
          'sbg:exposed': true
        - id: run_reference_proteome_similarity_2
          type: boolean?
          'sbg:exposed': true
        - id: additional_report_columns_2
          type:
            - 'null'
            - type: enum
              symbols:
                - sample_name
              name: additional_report_columns
          'sbg:exposed': true
        - id: fasta_size_2
          type: int?
          'sbg:exposed': true
        - id: downstream_sequence_length_2
          type: string?
          'sbg:exposed': true
        - id: exclude_nas_2
          type: boolean?
          'sbg:exposed': true
        - id: n_threads_2
          type: int?
          'sbg:exposed': true
        - id: ref_genome
          type:
            - 'null'
            - type: enum
              symbols:
                - hg38
                - hg19
              name: ref_genome
          'sbg:exposed': true
        - id: fusion_caller_1
          type:
            - 'null'
            - type: enum
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
          'sbg:y': 1711.25
        - id: alleles_1
          type: string
          'sbg:x': 0
          'sbg:y': 1818.203125
        - id: fusion_tsv_1
          'sbg:fileTypes': tsv
          type: File?
          'sbg:x': 0
          'sbg:y': 1604.296875
        - id: prediction_algorithms_2
          type: 'string[]'
          'sbg:x': 0
          'sbg:y': 1283.4375
        - id: sample_name_1
          type: string
          'sbg:x': 0
          'sbg:y': 962.578125
        - id: input_file_1
          'sbg:fileTypes': 'VCF, VCF.GZ'
          type: File
          'sbg:x': 0
          'sbg:y': 1497.34375
        - id: ref_genome_1
          'sbg:fileTypes': 'FASTA, FA'
          type: File
          'sbg:x': 0
          'sbg:y': 1069.53125
        - id: vep_cache_1
          'sbg:fileTypes': TAR.GZ
          type: File
          'sbg:x': 0
          'sbg:y': 855.625
        - id: cache_version_1
          type: int
          'sbg:exposed': true
        - id: merged_1
          type: boolean
          'sbg:exposed': true
        - id: symbol_1
          type: boolean?
          'sbg:exposed': true
        - id: biotype_1
          type: boolean?
          'sbg:exposed': true
        - id: numbers_1
          type: boolean?
          'sbg:exposed': true
        - id: canonical_1
          type: boolean?
          'sbg:exposed': true
        - id: total_length_1
          type: boolean?
          'sbg:exposed': true
        - id: sift_1
          type:
            - 'null'
            - type: enum
              symbols:
                - p
                - s
                - b
              name: sift
          'sbg:exposed': true
        - id: polyphen_1
          type:
            - 'null'
            - type: enum
              symbols:
                - p
                - s
                - b
              name: polyphen
          'sbg:exposed': true
        - id: terms_1
          type:
            - 'null'
            - type: enum
              symbols:
                - SO
                - display
                - NCBI
              name: terms
          'sbg:exposed': true
        - id: epitope_lengths_class_i_3
          type: 'int[]?'
          'sbg:exposed': true
        - id: epitope_lengths_class_ii_3
          type: 'int[]?'
          'sbg:exposed': true
        - id: binding_threshold_3
          type: int?
          'sbg:exposed': true
        - id: percentile_threshold_3
          type: int?
          'sbg:exposed': true
        - id: allele_specific_binding_thresholds_1
          type: boolean?
          'sbg:exposed': true
        - id: iedb_retries_3
          type: int?
          'sbg:exposed': true
        - id: keep_tmp_files_3
          type: boolean?
          'sbg:exposed': true
        - id: normal_sample_name_1
          type: string?
          'sbg:exposed': true
        - id: net_chop_method_3
          type:
            - 'null'
            - type: enum
              symbols:
                - cterm
                - 20s
              name: net_chop_method
          'sbg:exposed': true
        - id: netmhc_stab_3
          type: boolean?
          'sbg:exposed': true
        - id: run_reference_proteome_similarity_3
          type: boolean?
          'sbg:exposed': true
        - id: top_score_metric_3
          type:
            - 'null'
            - type: enum
              symbols:
                - lowest
                - median
              name: top_score_metric
          'sbg:exposed': true
        - id: net_chop_threshold_3
          type: float?
          'sbg:exposed': true
        - id: additional_report_columns_3
          type:
            - 'null'
            - type: enum
              symbols:
                - sample_name
              name: additional_report_columns
          'sbg:exposed': true
        - id: fasta_size_3
          type: int?
          'sbg:exposed': true
        - id: downstream_sequence_length_3
          type: string?
          'sbg:exposed': true
        - id: exclude_nas_3
          type: boolean?
          'sbg:exposed': true
        - id: minimum_fold_change_1
          type: float?
          'sbg:exposed': true
        - id: normal_cov_1
          type: int?
          'sbg:exposed': true
        - id: tdna_cov_1
          type: int?
          'sbg:exposed': true
        - id: trna_cov_1
          type: int?
          'sbg:exposed': true
        - id: normal_vaf_1
          type: float?
          'sbg:exposed': true
        - id: tdna_vaf_1
          type: float?
          'sbg:exposed': true
        - id: trna_vaf_1
          type: float?
          'sbg:exposed': true
        - id: expn_val_1
          type: float?
          'sbg:exposed': true
        - id: maximum_transcript_support_level_1
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
        - id: n_threads_3
          type: int?
          'sbg:exposed': true
        - id: vep_plugin_files_1
          'sbg:fileTypes': PM
          type: 'File[]?'
          'sbg:x': 0
          'sbg:y': 748.671875
        - id: prediction_algorithms_3
          type: 'string[]'
          'sbg:x': 0
          'sbg:y': 1176.484375
        - id: phased_proximal_variants_vcf_1
          type: File?
          'sbg:x': 0
          'sbg:y': 1390.390625
      outputs:
        - id: output_gz_tbi
          outputSource:
            - vep_pvactools_1/output_gz_tbi
          type: File
          'sbg:x': 864.0404052734375
          'sbg:y': 213.90625
        - id: vep_vcf
          outputSource:
            - vep_pvactools_1/vep_vcf
          'sbg:fileTypes': VCF
          type: File
          'sbg:x': 864.0404052734375
          'sbg:y': 0
        - id: output_gz
          outputSource:
            - vep_pvactools_1/output_gz
          type: File
          'sbg:x': 864.0404052734375
          'sbg:y': 320.859375
        - id: mhc_i_all_epitopes
          outputSource:
            - vep_pvactools_1/mhc_i_all_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 1497.34375
        - id: mhc_i_aggregated_report
          outputSource:
            - vep_pvactools_1/mhc_i_aggregated_report
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 1711.25
        - id: mhc_i_filtered_epitopes
          outputSource:
            - vep_pvactools_1/mhc_i_filtered_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 1283.4375
        - id: mhc_ii_all_epitopes
          outputSource:
            - vep_pvactools_1/mhc_ii_all_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 748.671875
        - id: mhc_ii_aggregated_report
          outputSource:
            - vep_pvactools_1/mhc_ii_aggregated_report
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 962.578125
        - id: mhc_ii_filtered_epitopes
          outputSource:
            - vep_pvactools_1/mhc_ii_filtered_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 534.765625
        - id: combined_all_epitopes
          outputSource:
            - vep_pvactools_1/combined_all_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 2246.015625
        - id: combined_aggregated_report
          outputSource:
            - vep_pvactools_1/combined_aggregated_report
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 2459.921875
        - id: combined_filtered_epitopes
          outputSource:
            - vep_pvactools_1/combined_filtered_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 2032.109375
        - id: pvacseq_predictions
          outputSource:
            - vep_pvactools_1/pvacseq_predictions
          type: Directory
          'sbg:x': 864.0404052734375
          'sbg:y': 106.953125
        - id: mhc_i_all_epitopes_1
          outputSource:
            - agfusion_pvacfuse_1/mhc_i_all_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 1390.390625
        - id: mhc_i_aggregated_report_1
          outputSource:
            - agfusion_pvacfuse_1/mhc_i_aggregated_report
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 1604.296875
        - id: mhc_i_filtered_epitopes_1
          outputSource:
            - agfusion_pvacfuse_1/mhc_i_filtered_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 1176.484375
        - id: mhc_ii_all_epitopes_1
          outputSource:
            - agfusion_pvacfuse_1/mhc_ii_all_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 641.71875
        - id: mhc_ii_aggregated_report_1
          outputSource:
            - agfusion_pvacfuse_1/mhc_ii_aggregated_report
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 855.625
        - id: mhc_ii_filtered_epitopes_1
          outputSource:
            - agfusion_pvacfuse_1/mhc_ii_filtered_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 427.8125
        - id: combined_all_epitopes_1
          outputSource:
            - agfusion_pvacfuse_1/combined_all_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 2139.0625
        - id: combined_aggregated_report_1
          outputSource:
            - agfusion_pvacfuse_1/combined_aggregated_report
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 2352.96875
        - id: combined_filtered_epitopes_1
          outputSource:
            - agfusion_pvacfuse_1/combined_filtered_epitopes
          type: File?
          'sbg:x': 864.0404052734375
          'sbg:y': 1925.15625
        - id: mhc_i
          outputSource:
            - agfusion_pvacfuse_1/mhc_i
          type: Directory?
          'sbg:x': 864.0404052734375
          'sbg:y': 1818.203125
        - id: mhc_ii
          outputSource:
            - agfusion_pvacfuse_1/mhc_ii
          type: Directory?
          'sbg:x': 864.0404052734375
          'sbg:y': 1069.53125
        - id: combined
          outputSource:
            - agfusion_pvacfuse_1/combined
          type: Directory?
          'sbg:x': 864.0404052734375
          'sbg:y': 2566.875
      steps:
        - id: agfusion_pvacfuse_1
          in:
            - id: sample_name
              source: sample_name_1
            - id: alleles
              source: alleles_1
            - id: prediction_algorithms
              source:
                - prediction_algorithms_2
            - id: epitope_lengths_class_i
              source:
                - epitope_lengths_class_i_2
            - id: epitope_lengths_class_ii
              source:
                - epitope_lengths_class_ii_2
            - id: binding_threshold
              source: binding_threshold_2
            - id: percentile_threshold
              source: percentile_threshold_2
            - id: iedb_retries
              source: iedb_retries_2
            - id: keep_tmp_files
              source: keep_tmp_files_2
            - id: net_chop_method
              source: net_chop_method_2
            - id: netmhc_stab
              source: netmhc_stab_2
            - id: top_score_metric
              source: top_score_metric_2
            - id: net_chop_threshold
              source: net_chop_threshold_2
            - id: run_reference_proteome_similarity
              source: run_reference_proteome_similarity_2
            - id: additional_report_columns
              source: additional_report_columns_2
            - id: fasta_size
              source: fasta_size_2
            - id: downstream_sequence_length
              source: downstream_sequence_length_2
            - id: exclude_nas
              source: exclude_nas_2
            - id: n_threads
              source: n_threads_2
            - id: ref_genome
              source: ref_genome
            - id: fusion_tsv
              source: fusion_tsv_1
            - id: fusion_caller
              source: fusion_caller_1
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
            - id: combined
          run:
            class: Workflow
            cwlVersion: v1.0
            id: mwonge/mwtest/agfusion-pvacfuse/11
            label: agfusion-pvacfuse
            $namespaces:
              sbg: 'https://sevenbridges.com'
            inputs:
              - id: sample_name
                type: string
                'sbg:x': 0
                'sbg:y': 481.640625
              - id: alleles
                type: string
                'sbg:x': 161.296875
                'sbg:y': 528.171875
              - id: prediction_algorithms
                type: 'string[]'
                'sbg:exposed': true
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
              - id: ref_genome
                type:
                  - 'null'
                  - type: enum
                    symbols:
                      - hg38
                      - hg19
                    name: ref_genome
                'sbg:exposed': true
              - id: fusion_tsv
                'sbg:fileTypes': tsv
                type: File?
                'sbg:x': 0
                'sbg:y': 588.671875
              - id: fusion_caller
                type:
                  - 'null'
                  - type: enum
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
                'sbg:y': 695.703125
            outputs:
              - id: mhc_i_all_epitopes
                outputSource:
                  - pvacfuse/mhc_i_all_epitopes
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 535.15625
              - id: mhc_i_aggregated_report
                outputSource:
                  - pvacfuse/mhc_i_aggregated_report
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 642.1875
              - id: mhc_i_filtered_epitopes
                outputSource:
                  - pvacfuse/mhc_i_filtered_epitopes
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 428.125
              - id: mhc_ii_all_epitopes
                outputSource:
                  - pvacfuse/mhc_ii_all_epitopes
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 107.03125
              - id: mhc_ii_aggregated_report
                outputSource:
                  - pvacfuse/mhc_ii_aggregated_report
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 214.0625
              - id: mhc_ii_filtered_epitopes
                outputSource:
                  - pvacfuse/mhc_ii_filtered_epitopes
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 0
              - id: combined_all_epitopes
                outputSource:
                  - pvacfuse/combined_all_epitopes
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 963.28125
              - id: combined_aggregated_report
                outputSource:
                  - pvacfuse/combined_aggregated_report
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 1070.3125
              - id: combined_filtered_epitopes
                outputSource:
                  - pvacfuse/combined_filtered_epitopes
                type: File?
                'sbg:x': 817.9067993164062
                'sbg:y': 856.25
              - id: mhc_i
                outputSource:
                  - pvacfuse/mhc_i
                type: Directory?
                'sbg:x': 817.9067993164062
                'sbg:y': 749.21875
              - id: mhc_ii
                outputSource:
                  - pvacfuse/mhc_ii
                type: Directory?
                'sbg:x': 817.9067993164062
                'sbg:y': 321.09375
              - id: combined
                outputSource:
                  - pvacfuse/combined
                type: Directory?
                'sbg:x': 817.9067993164062
                'sbg:y': 1177.34375
            steps:
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
                  - id: combined
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/pvacfuse/10
                  baseCommand:
                    - ln
                    - '-s'
                  inputs:
                    - id: input_file
                      type: Directory
                      inputBinding:
                        position: 1
                        shellQuote: false
                    - id: sample_name
                      type: string
                      inputBinding:
                        position: 2
                        shellQuote: false
                    - id: alleles
                      type: string
                      inputBinding:
                        position: 3
                        prefix: ''
                        separate: false
                        shellQuote: false
                    - id: prediction_algorithms
                      type: 'string[]'
                      inputBinding:
                        position: 4
                        shellQuote: false
                    - id: epitope_lengths_class_i
                      type: 'int[]?'
                      inputBinding:
                        position: 0
                        prefix: '-e1'
                        itemSeparator: ','
                        shellQuote: false
                    - id: epitope_lengths_class_ii
                      type: 'int[]?'
                      inputBinding:
                        position: 0
                        prefix: '-e2'
                        itemSeparator: ','
                        shellQuote: false
                    - id: binding_threshold
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '-b'
                        shellQuote: false
                    - id: percentile_threshold
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '--percentile-threshold'
                        shellQuote: false
                    - id: iedb_retries
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '-r'
                        shellQuote: false
                    - id: keep_tmp_files
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '-k'
                        shellQuote: false
                    - id: net_chop_method
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - cterm
                            - 20s
                          name: net_chop_method
                      inputBinding:
                        position: 0
                        prefix: '--net-chop-method'
                        shellQuote: false
                    - id: netmhc_stab
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '--netmhc-stab'
                        shellQuote: false
                    - id: top_score_metric
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - lowest
                            - median
                          name: top_score_metric
                      inputBinding:
                        position: 0
                        prefix: '-m'
                        shellQuote: false
                    - id: net_chop_threshold
                      type: float?
                      inputBinding:
                        position: 0
                        prefix: '--net-chop-threshold'
                        shellQuote: false
                    - id: run_reference_proteome_similarity
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '--run-reference-proteome-similarity'
                        shellQuote: false
                    - id: additional_report_columns
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - sample_name
                          name: additional_report_columns
                      inputBinding:
                        position: 0
                        prefix: '-a'
                        shellQuote: false
                    - id: fasta_size
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '-s'
                        shellQuote: false
                    - id: downstream_sequence_length
                      type: string?
                      inputBinding:
                        position: 0
                        prefix: '-d'
                        shellQuote: false
                    - id: exclude_nas
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '--exclude-NAs'
                        shellQuote: false
                    - default: 8
                      id: n_threads
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '--n-threads'
                        shellQuote: false
                  outputs:
                    - id: mhc_i_all_epitopes
                      type: File?
                      outputBinding:
                        glob: MHC_Class_I/$(inputs.sample_name).all_epitopes.tsv
                    - id: mhc_i_aggregated_report
                      type: File?
                      outputBinding:
                        glob: >-
                          MHC_Class_I/$(inputs.sample_name).all_epitopes.aggregated.tsv
                    - id: mhc_i_filtered_epitopes
                      type: File?
                      outputBinding:
                        glob: MHC_Class_I/$(inputs.sample_name).filtered.tsv
                    - id: mhc_ii_all_epitopes
                      type: File?
                      outputBinding:
                        glob: MHC_Class_II/$(inputs.sample_name).all_epitopes.tsv
                    - id: mhc_ii_aggregated_report
                      type: File?
                      outputBinding:
                        glob: >-
                          MHC_Class_II/$(inputs.sample_name).all_epitopes.aggregated.tsv
                    - id: mhc_ii_filtered_epitopes
                      type: File?
                      outputBinding:
                        glob: MHC_Class_II/$(inputs.sample_name).filtered.tsv
                    - id: combined_all_epitopes
                      type: File?
                      outputBinding:
                        glob: combined/$(inputs.sample_name).all_epitopes.tsv
                    - id: combined_aggregated_report
                      type: File?
                      outputBinding:
                        glob: >-
                          combined/$(inputs.sample_name).all_epitopes.aggregated.tsv
                    - id: combined_filtered_epitopes
                      type: File?
                      outputBinding:
                        glob: combined/$(inputs.sample_name).filtered.tsv
                    - id: mhc_i
                      type: Directory?
                      outputBinding:
                        glob: MHC_Class_I
                    - id: mhc_ii
                      type: Directory?
                      outputBinding:
                        glob: MHC_Class_II
                    - id: combined
                      type: Directory?
                      outputBinding:
                        glob: combined
                  label: pvacfuse
                  arguments:
                    - position: 0
                      shellQuote: false
                      valueFrom: $TMPDIR
                    - /tmp/pvacseq
                    - position: 0
                      shellQuote: false
                      valueFrom: ' && '
                    - export
                    - TMPDIR=/tmp/pvacseq
                    - position: 0
                      shellQuote: false
                      valueFrom: ' && '
                    - /usr/local/bin/pvacfuse
                    - run
                    - '--iedb-install-directory'
                    - /opt/iedb
                    - position: 5
                      shellQuote: false
                      valueFrom: $(runtime.outdir)
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: 16000
                      coresMin: $(inputs.n_threads)
                    - class: DockerRequirement
                      dockerPull: 'griffithlab/pvactools:2.0.3'
                    - class: InlineJavascriptRequirement
                  'sbg:projectName': zcc-cavatica
                  'sbg:revisionsInfo':
                    - 'sbg:revision': 0
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626149875
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: ./pvacfuse.cwl
                    - 'sbg:revision': 1
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626150485
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 2
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626150538
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 3
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626150592
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 4
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626150643
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 5
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626150663
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 6
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626150822
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 7
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626151410
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 8
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627532510
                      'sbg:revisionNotes': ''
                    - 'sbg:revision': 9
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1628034419
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: ./pvacfuse.cwl
                    - 'sbg:revision': 10
                      'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1628034751
                      'sbg:revisionNotes': Updated pVACtools version to 2.0.3
                  'sbg:image_url': null
                  'sbg:appVersion':
                    - v1.0
                  'sbg:id': mwonge/mwtest/pvacfuse/10
                  'sbg:revision': 10
                  'sbg:revisionNotes': Updated pVACtools version to 2.0.3
                  'sbg:modifiedOn': 1628034751
                  'sbg:modifiedBy': rbowen_james
                  'sbg:createdOn': 1626149875
                  'sbg:createdBy': rbowen_james
                  'sbg:project': mwonge/mwtest
                  'sbg:sbgMaintained': false
                  'sbg:validationErrors': []
                  'sbg:contributors':
                    - rbowen_james
                  'sbg:latestRevision': 10
                  'sbg:publisher': sbg
                  'sbg:content_hash': >-
                    abbc4ec71a2dc198deb720f77604d9d9fd3f81e5a09adffd1cecbd55aa67d28f8
                label: pvacfuse
                'sbg:x': 383.5577697753906
                'sbg:y': 511.671875
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
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/agfusion/14
                  baseCommand: []
                  inputs:
                    - id: fusion_tsv
                      type: File?
                      inputBinding:
                        position: 1
                        prefix: '-f'
                        shellQuote: false
                        valueFrom: infile
                      'sbg:fileTypes': tsv
                    - id: fusion_caller
                      type:
                        - 'null'
                        - type: enum
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
                      inputBinding:
                        position: 2
                        prefix: '-a'
                        shellQuote: false
                    - id: ref_genome
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - hg38
                            - hg19
                          name: ref_genome
                      inputBinding:
                        position: 0
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.ref_genome == "hg38") {
                                  return '&& agfusion download -s homo_sapiens -r 87 -d .'
                              } else if (inputs.ref_genome == "hg19") {
                                  return '&& agfusion download -s homo_sapiens -r 75 -d .'
                              }
                          }
                  outputs:
                    - id: output_dir
                      type: Directory?
                      outputBinding:
                        glob: agfusion-output
                  label: agfusion
                  arguments:
                    - position: 3
                      prefix: '-db'
                      shellQuote: false
                      valueFrom: |-
                        ${
                            if (inputs.ref_genome == "hg38") {
                                return "agfusion.homo_sapiens.87.db"
                            } else if (inputs.ref_genome == "hg19") {
                                return "agfusion.homo_sapiens.75.db"
                            }
                        }
                    - position: 4
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: '--middlestar'
                    - position: 5
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: '--noncanonical'
                    - position: 3
                      prefix: '-o'
                      shellQuote: false
                      valueFrom: ./agfusion-output
                    - position: 0
                      prefix: ''
                      shellQuote: false
                      valueFrom: >-
                        touch infile && if head -n 1 $(inputs.fusion_tsv.path) |
                        grep 'est_J'; then cut -f1-3,6-
                        $(inputs.fusion_tsv.path) > infile ; else cat
                        $(inputs.fusion_tsv.path) > infile ; fi && head infile
                        1>&2
                    - position: 1
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: '&& agfusion batch'
                  requirements:
                    - class: ShellCommandRequirement
                    - class: DockerRequirement
                      dockerPull: 'rachelbj/agfusion:1.1'
                    - class: InlineJavascriptRequirement
                  'sbg:appVersion':
                    - v1.0
                  'sbg:content_hash': >-
                    aacef118efebe45b6a8043d29b795f43d291caee94c0bb3d714f740cb32ef9b8e
                  'sbg:contributors':
                    - rbowen_james
                  'sbg:createdBy': rbowen_james
                  'sbg:createdOn': 1626132773
                  'sbg:id': mwonge/mwtest/agfusion/14
                  'sbg:image_url': null
                  'sbg:latestRevision': 14
                  'sbg:modifiedBy': rbowen_james
                  'sbg:modifiedOn': 1626845542
                  'sbg:project': mwonge/mwtest
                  'sbg:projectName': zcc-cavatica
                  'sbg:publisher': sbg
                  'sbg:revision': 14
                  'sbg:revisionNotes': null
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626132773
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626133908
                      'sbg:revision': 1
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626136375
                      'sbg:revision': 2
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626136866
                      'sbg:revision': 3
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626137191
                      'sbg:revision': 4
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626138935
                      'sbg:revision': 5
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626139315
                      'sbg:revision': 6
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626139355
                      'sbg:revision': 7
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626139386
                      'sbg:revision': 8
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626149500
                      'sbg:revision': 9
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626582610
                      'sbg:revision': 10
                      'sbg:revisionNotes': Added fusion callers
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626583052
                      'sbg:revision': 11
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626583631
                      'sbg:revision': 12
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626845196
                      'sbg:revision': 13
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626845542
                      'sbg:revision': 14
                      'sbg:revisionNotes': null
                  'sbg:sbgMaintained': false
                  'sbg:validationErrors': []
                label: agfusion
                'sbg:x': 161.296875
                'sbg:y': 642.1875
            requirements: []
            'sbg:appVersion':
              - v1.0
            'sbg:id': mwonge/mwtest/agfusion-pvacfuse/11
            'sbg:revision': 11
            'sbg:revisionNotes': Updated pVACtools version to 2.0.3
            'sbg:modifiedOn': 1628035130
            'sbg:modifiedBy': rbowen_james
            'sbg:createdOn': 1626581605
            'sbg:createdBy': rbowen_james
            'sbg:project': mwonge/mwtest
            'sbg:projectName': zcc-cavatica
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
            'sbg:contributors':
              - rbowen_james
            'sbg:latestRevision': 11
            'sbg:revisionsInfo':
              - 'sbg:revision': 0
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626581605
                'sbg:revisionNotes': null
              - 'sbg:revision': 1
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626582635
                'sbg:revisionNotes': Added fusion callers
              - 'sbg:revision': 2
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626582787
                'sbg:revisionNotes': ''
              - 'sbg:revision': 3
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626583088
                'sbg:revisionNotes': null
              - 'sbg:revision': 4
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626583283
                'sbg:revisionNotes': null
              - 'sbg:revision': 5
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626583646
                'sbg:revisionNotes': null
              - 'sbg:revision': 6
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626845228
                'sbg:revisionNotes': null
              - 'sbg:revision': 7
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626845560
                'sbg:revisionNotes': null
              - 'sbg:revision': 8
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627532556
                'sbg:revisionNotes': null
              - 'sbg:revision': 9
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627532744
                'sbg:revisionNotes': null
              - 'sbg:revision': 10
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627533023
                'sbg:revisionNotes': null
              - 'sbg:revision': 11
                'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1628035130
                'sbg:revisionNotes': Updated pVACtools version to 2.0.3
            'sbg:image_url': >-
              https://cavatica.sbgenomics.com/ns/brood/images/mwonge/mwtest/agfusion-pvacfuse/11.png
            'sbg:publisher': sbg
            'sbg:content_hash': afdf8212af094eed4a8e26404077911e765b3a474d02972919b000f02b61ef0d9
          label: agfusion-pvacfuse
          'sbg:x': 313.03125
          'sbg:y': 1336.9140625
        - id: vep_pvactools_1
          in:
            - id: vep_cache
              source: vep_cache_1
            - id: ref_genome
              source: ref_genome_1
            - id: input_file
              source: input_file_1
            - id: vep_plugin_files
              source:
                - vep_plugin_files_1
            - id: cache_version
              source: cache_version_1
            - id: merged
              source: merged_1
            - id: symbol
              source: symbol_1
            - id: biotype
              source: biotype_1
            - id: numbers
              source: numbers_1
            - id: canonical
              source: canonical_1
            - id: total_length
              source: total_length_1
            - id: sift
              source: sift_1
            - id: polyphen
              source: polyphen_1
            - id: terms
              source: terms_1
            - id: epitope_lengths_class_i
              source:
                - epitope_lengths_class_i_3
            - id: epitope_lengths_class_ii
              source:
                - epitope_lengths_class_ii_3
            - id: binding_threshold
              source: binding_threshold_3
            - id: percentile_threshold
              source: percentile_threshold_3
            - id: allele_specific_binding_thresholds
              source: allele_specific_binding_thresholds_1
            - id: iedb_retries
              source: iedb_retries_3
            - id: keep_tmp_files
              source: keep_tmp_files_3
            - id: normal_sample_name
              source: normal_sample_name_1
            - id: net_chop_method
              source: net_chop_method_3
            - id: netmhc_stab
              source: netmhc_stab_3
            - id: run_reference_proteome_similarity
              source: run_reference_proteome_similarity_3
            - id: top_score_metric
              source: top_score_metric_3
            - id: net_chop_threshold
              source: net_chop_threshold_3
            - id: additional_report_columns
              source: additional_report_columns_3
            - id: fasta_size
              source: fasta_size_3
            - id: downstream_sequence_length
              source: downstream_sequence_length_3
            - id: exclude_nas
              source: exclude_nas_3
            - id: minimum_fold_change
              source: minimum_fold_change_1
            - id: normal_cov
              source: normal_cov_1
            - id: tdna_cov
              source: tdna_cov_1
            - id: trna_cov
              source: trna_cov_1
            - id: normal_vaf
              source: normal_vaf_1
            - id: tdna_vaf
              source: tdna_vaf_1
            - id: trna_vaf
              source: trna_vaf_1
            - id: expn_val
              source: expn_val_1
            - id: maximum_transcript_support_level
              source: maximum_transcript_support_level_1
            - id: n_threads
              source: n_threads_3
            - id: alleles
              source: alleles_1
            - id: phased_proximal_variants_vcf
              source: phased_proximal_variants_vcf_1
            - id: prediction_algorithms
              source:
                - prediction_algorithms_3
          out:
            - id: output_gz_tbi
            - id: vep_vcf
            - id: output_gz
            - id: mhc_i_all_epitopes
            - id: mhc_i_aggregated_report
            - id: mhc_i_filtered_epitopes
            - id: mhc_ii_all_epitopes
            - id: mhc_ii_aggregated_report
            - id: mhc_ii_filtered_epitopes
            - id: combined_all_epitopes
            - id: combined_aggregated_report
            - id: combined_filtered_epitopes
            - id: pvacseq_predictions
          run:
            class: Workflow
            cwlVersion: v1.0
            id: mwonge/mwtest/vep-pvactools/4
            label: vep-pvactools
            $namespaces:
              sbg: 'https://sevenbridges.com'
            inputs:
              - id: vep_cache
                'sbg:fileTypes': TAR.GZ
                type: File
                'sbg:x': 0
                'sbg:y': 427.625
              - id: ref_genome
                'sbg:fileTypes': 'FASTA, FA'
                type: File
                'sbg:x': 0
                'sbg:y': 534.53125
              - id: input_file
                'sbg:fileTypes': 'VCF, VCF.GZ'
                type: File
                'sbg:x': 0
                'sbg:y': 641.4375
              - id: vep_plugin_files
                'sbg:fileTypes': PM
                type: 'File[]?'
                'sbg:x': 0
                'sbg:y': 320.71875
              - id: cache_version
                type: int
                'sbg:exposed': true
              - id: merged
                type: boolean
                'sbg:exposed': true
              - id: symbol
                type: boolean?
                'sbg:exposed': true
              - id: biotype
                type: boolean?
                'sbg:exposed': true
              - id: numbers
                type: boolean?
                'sbg:exposed': true
              - id: canonical
                type: boolean?
                'sbg:exposed': true
              - id: total_length
                type: boolean?
                'sbg:exposed': true
              - id: sift
                type:
                  - 'null'
                  - type: enum
                    symbols:
                      - p
                      - s
                      - b
                    name: sift
                'sbg:exposed': true
              - id: polyphen
                type:
                  - 'null'
                  - type: enum
                    symbols:
                      - p
                      - s
                      - b
                    name: polyphen
                'sbg:exposed': true
              - id: terms
                type:
                  - 'null'
                  - type: enum
                    symbols:
                      - SO
                      - display
                      - NCBI
                    name: terms
                'sbg:exposed': true
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
              - id: allele_specific_binding_thresholds
                type: boolean?
                'sbg:exposed': true
              - id: iedb_retries
                type: int?
                'sbg:exposed': true
              - id: keep_tmp_files
                type: boolean?
                'sbg:exposed': true
              - id: normal_sample_name
                type: string?
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
              - id: run_reference_proteome_similarity
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
              - id: minimum_fold_change
                type: float?
                'sbg:exposed': true
              - id: normal_cov
                type: int?
                'sbg:exposed': true
              - id: tdna_cov
                type: int?
                'sbg:exposed': true
              - id: trna_cov
                type: int?
                'sbg:exposed': true
              - id: normal_vaf
                type: float?
                'sbg:exposed': true
              - id: tdna_vaf
                type: float?
                'sbg:exposed': true
              - id: trna_vaf
                type: float?
                'sbg:exposed': true
              - id: expn_val
                type: float?
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
              - id: n_threads
                type: int?
                'sbg:exposed': true
              - id: alleles
                type: string
                'sbg:x': 489.11181640625
                'sbg:y': 701.890625
              - id: phased_proximal_variants_vcf
                type: File?
                'sbg:x': 489.11181640625
                'sbg:y': 474.078125
              - id: prediction_algorithms
                type: 'string[]'
                'sbg:x': 489.11181640625
                'sbg:y': 367.171875
            outputs:
              - id: output_gz_tbi
                outputSource:
                  - bgzip/output_gz_tbi
                type: File
                'sbg:x': 782.31494140625
                'sbg:y': 481.078125
              - id: vep_vcf
                outputSource:
                  - vep_with_plugins/vep_vcf
                'sbg:fileTypes': VCF
                type: File
                'sbg:x': 489.11181640625
                'sbg:y': 260.265625
              - id: output_gz
                outputSource:
                  - bgzip/output_gz
                type: File
                'sbg:x': 782.31494140625
                'sbg:y': 587.984375
              - id: mhc_i_all_epitopes
                outputSource:
                  - pvacseq/mhc_i_all_epitopes
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 534.53125
              - id: mhc_i_aggregated_report
                outputSource:
                  - pvacseq/mhc_i_aggregated_report
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 641.4375
              - id: mhc_i_filtered_epitopes
                outputSource:
                  - pvacseq/mhc_i_filtered_epitopes
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 427.625
              - id: mhc_ii_all_epitopes
                outputSource:
                  - pvacseq/mhc_ii_all_epitopes
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 213.8125
              - id: mhc_ii_aggregated_report
                outputSource:
                  - pvacseq/mhc_ii_aggregated_report
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 320.71875
              - id: mhc_ii_filtered_epitopes
                outputSource:
                  - pvacseq/mhc_ii_filtered_epitopes
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 106.90625
              - id: combined_all_epitopes
                outputSource:
                  - pvacseq/combined_all_epitopes
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 855.25
              - id: combined_aggregated_report
                outputSource:
                  - pvacseq/combined_aggregated_report
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 962.15625
              - id: combined_filtered_epitopes
                outputSource:
                  - pvacseq/combined_filtered_epitopes
                type: File?
                'sbg:x': 1291.4244384765625
                'sbg:y': 748.34375
              - id: pvacseq_predictions
                outputSource:
                  - pvacseq/pvacseq_predictions
                type: Directory
                'sbg:x': 1291.4244384765625
                'sbg:y': 0
            steps:
              - id: vep_with_plugins
                in:
                  - id: input_file
                    source: input_file
                  - id: vep_cache
                    source: vep_cache
                  - id: ref_genome
                    source: ref_genome
                  - id: vep_plugin_files
                    source:
                      - vep_plugin_files
                  - id: cache_version
                    default: 104
                    source: cache_version
                  - id: merged
                    default: true
                    source: merged
                  - id: symbol
                    default: true
                    source: symbol
                  - id: biotype
                    default: true
                    source: biotype
                  - id: numbers
                    default: true
                    source: numbers
                  - id: canonical
                    default: true
                    source: canonical
                  - id: total_length
                    default: true
                    source: total_length
                  - id: sift
                    default: b
                    source: sift
                  - id: polyphen
                    default: b
                    source: polyphen
                  - id: terms
                    default: SO
                    source: terms
                out:
                  - id: vep_vcf
                  - id: vep_stats
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/vep-with-plugins/26
                  baseCommand: []
                  inputs:
                    - id: input_file
                      type: File
                      inputBinding:
                        position: 3
                        prefix: '--input_file'
                        shellQuote: false
                      label: Input File
                      doc: Input VCF to VEP annotate. Can be bgzipped.
                      'sbg:fileTypes': 'VCF, VCF.GZ'
                    - id: vep_cache
                      type: File
                      label: VEP Cache
                      doc: >-
                        VEP cache supplied as a tar.gz. Caches can be downloaded
                        from [VEP release page](http://ftp.ensembl.org/pub/). It
                        is recommended that the cache version used matches the
                        VEP release number (note that this app uses the latest
                        VEP release. See
                        [here](https://hub.docker.com/r/ensemblorg/ensembl-vep/tags?page=1&ordering=last_updated)
                        for more info).
                      'sbg:fileTypes': TAR.GZ
                    - id: ref_genome
                      type: File
                      inputBinding:
                        position: 3
                        prefix: '--fasta'
                        shellQuote: false
                      label: Reference Genome
                      doc: >-
                        The reference genome FASTA file to use to look up
                        reference sequence.
                      'sbg:fileTypes': 'FASTA, FA'
                    - id: vep_plugin_files
                      type: 'File[]?'
                      label: VEP Plugin Files
                      doc: >-
                        Optional VEP plugin files to use when annotating the
                        input VCF.
                      'sbg:fileTypes': PM
                    - id: cache_version
                      type: int
                      inputBinding:
                        position: 3
                        prefix: '--cache_version'
                        shellQuote: false
                      label: Cache Version
                      doc: >-
                        Version number of the cache used. It is recommended that
                        the cache version used matches the VEP release number
                        (note that this app uses the latest VEP release. See
                        [here](https://hub.docker.com/r/ensemblorg/ensembl-vep/tags?page=1&ordering=last_updated)
                        for more info).
                    - id: merged
                      type: boolean
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.merged == true) {
                                  return '--merged'
                              } else {
                                  return ''
                              }
                              
                          }
                      label: Merged
                      doc: >-
                        Use the merged Ensembl and RefSeq cache. Consequences
                        are flagged with the SOURCE of each transcript used.

                        NOTE: This flag MUST be used if the cache is merged.
                    - 'sbg:toolDefaultValue': 'No'
                      id: symbol
                      type: boolean?
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.symbol == true) {
                                  return '--symbol'
                              } else {
                                  return ''
                              }
                          }
                      label: Symbol
                      doc: >-
                        Adds the gene symbol (e.g. HGNC) (where available) to
                        the output.
                    - 'sbg:toolDefaultValue': 'No'
                      id: biotype
                      type: boolean?
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.biotype == true) {
                                  return '--biotype'
                              } else {
                                  return ''
                              }
                              
                          }
                      label: Biotype
                      doc: >-
                        Adds the biotype of the transcript or regulatory
                        feature.
                    - 'sbg:toolDefaultValue': 'No'
                      id: numbers
                      type: boolean?
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.numbers == true) {
                                  return '--numbers'
                              } else {
                                  return ''
                              }
                          }
                      label: Numbers
                      doc: >-
                        Adds affected exon and intron numbering to to output.
                        Format is Number/Total.
                    - 'sbg:toolDefaultValue': 'No'
                      id: canonical
                      type: boolean?
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.canonical == true) {
                                  return '--canonical'
                              } else {
                                  return ''
                              }
                          }
                      label: Canonical
                      doc: >-
                        Adds a flag indicating if the transcript is the
                        canonical transcript for the gene.
                    - 'sbg:toolDefaultValue': 'No'
                      id: total_length
                      type: boolean?
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.total_length == true) {
                                  return '--total_length'
                              } else {
                                  return ''
                              }
                          }
                      label: Total Length
                      doc: 'Give cDNA, CDS and protein positions as Position/Length.'
                    - 'sbg:toolDefaultValue': 'No'
                      id: sift
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - p
                            - s
                            - b
                          name: sift
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.sift == null) {
                                  return ''
                              } else {
                                  return '--sift ' + inputs.sift
                              }
                          }
                      label: Sift
                      doc: >-
                        Species limited SIFT predicts whether an amino acid
                        substitution affects protein function based on sequence
                        homology and the physical properties of amino acids. VEP
                        can output the prediction term, score or both.
                    - 'sbg:toolDefaultValue': 'No'
                      id: polyphen
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - p
                            - s
                            - b
                          name: polyphen
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.polyphen == null) {
                                  return ''
                              } else {
                                  return '--polyphen ' + inputs.polyphen
                              }
                          }
                      label: Polyphen
                      doc: >-
                        Human only PolyPhen is a tool which predicts possible
                        impact of an amino acid substitution on the structure
                        and function of a human protein using straightforward
                        physical and comparative considerations. VEP can output
                        the prediction term, score or both.
                    - 'sbg:toolDefaultValue': SO
                      id: terms
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - SO
                            - display
                            - NCBI
                          name: terms
                      inputBinding:
                        position: 6
                        separate: false
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (inputs.terms == null) {
                                  return ''
                              } else {
                                  return '--terms ' + inputs.terms
                              }
                          }
                      label: Terms
                      doc: >-
                        The type of consequence terms to output. The Ensembl
                        terms are described
                        [here](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences).
                        The Sequence Ontology is a joint effort by genome
                        annotation centres to standardise descriptions of
                        biological sequences.
                  outputs:
                    - id: vep_vcf
                      doc: VEP annotated VCF file.
                      label: Annotated VCF
                      type: File
                      outputBinding:
                        glob: '*.vep.vcf'
                      'sbg:fileTypes': VCF
                    - id: vep_stats
                      doc: Stats file produced by VEP.
                      label: VEP Stats
                      type: File?
                      outputBinding:
                        glob: '*.vep.html'
                      'sbg:fileTypes': HTML
                  doc: >-
                    #About VEP

                    VEP determines the effect of your variants (SNPs,
                    insertions, deletions, CNVs or structural variants) on
                    genes, transcripts, and protein sequence, as well as
                    regulatory regions.


                    Simply input the coordinates of your variants and the
                    nucleotide changes to find out the:

                    Genes and Transcripts affected by the variants

                    Location of the variants (e.g. upstream of a transcript, in
                    coding sequence, in non-coding RNA, in regulatory regions)

                    Consequence of your variants on the protein sequence (e.g.
                    stop gained, missense, stop lost, frameshift), see variant
                    consequences

                    Known variants that match yours, and associated minor allele
                    frequencies from the 1000 Genomes Project

                    SIFT and PolyPhen-2 scores for changes to protein sequence

                    ... And more! See data types, versions.


                    # About this CWL tool

                    This tool is intended for VEP annotation of VCF files prior
                    to neoantigen prediction using pVACseq
                    ([info](https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html)).
                    All VEP options required for pVACseq are exposed as app
                    settings, plus some additional options. 


                    ## Cache

                    The VEP cache must be supplied as a `tar.gz`. Caches can be
                    downloaded from [VEP release
                    page](http://ftp.ensembl.org/pub/). It is recommended that
                    the cache version used matches the VEP release number (note
                    that this app uses the latest VEP release. See
                    [here](https://hub.docker.com/r/ensemblorg/ensembl-vep/tags?page=1&ordering=last_updated)
                    to find out the latest VEP release number). As of July 2021,
                    the latest VEP release was version 104.


                    ### Merged

                    If the cache version used is merged, the `--merged` flag
                    must be used. This is achieved by setting the 'merged' input
                    option to 'Yes'. 


                    ## Plugins

                    pVACseq requires the use of the Frameshift and Wildtype
                    plugins, available for download
                    [here](https://github.com/griffithlab/pVACtools/tree/master/tools/pvacseq/VEP_plugins).
                  label: vep-with-plugins
                  arguments:
                    - position: 4
                      prefix: '--dir_cache'
                      shellQuote: false
                      valueFrom: |-
                        ${
                            return './cache'
                        }
                    - position: 4
                      prefix: '--output_file'
                      shellQuote: false
                      valueFrom: |-
                        ${
                            var in_file = inputs.input_file.basename
                            
                            if (in_file.endsWith(".gz")) {
                                var out_file = in_file.replace('.vcf.gz', '.vep.vcf')
                            } else {
                                var out_file = in_file.replace('.vcf', '.vep.vcf')
                            }
                            
                            return out_file
                        }
                    - position: 5
                      prefix: '--stats_file'
                      shellQuote: false
                      valueFrom: |-
                        ${
                            var in_file = inputs.input_file.basename
                            
                            if (in_file.endsWith(".gz")) {
                                var out_file = in_file.replace('.vcf.gz', '.vep.html')
                            } else {
                                var out_file = in_file.replace('.vcf', '.vep.html')
                            }
                            
                            return out_file
                        }
                    - position: 5
                      prefix: '--species'
                      shellQuote: false
                      valueFrom: homo_sapiens
                    - position: 3
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: >-
                        --cache --format vcf --vcf --offline --fork 8
                        --no_progress --tsl --hgvs --shift_hgvs 1
                    - position: 0
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: |-
                        mkdir ./plugins ${
                          if (inputs.vep_plugin_files == null) {
                            return ""
                          }

                          let mv_cmd = "";
                          for (var i=0; i < inputs.vep_plugin_files.length; i++) {
                            mv_cmd = mv_cmd + " && mv " + inputs.vep_plugin_files[i].path + " ./plugins/";
                          }
                          return mv_cmd;
                        } && ls ./plugins 1>&2 &&
                    - position: 2
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: vep
                    - position: 7
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: |2-
                         ${
                            let plugin_cmd = "";
                            for (var i=0; i < inputs.vep_plugin_files.length; i++) {
                                plugin_split = inputs.vep_plugin_files[i].path.split('/')
                                plugin = plugin_split[plugin_split.length-1].split('.')[0]
                                plugin_cmd = plugin_cmd + "--plugin " + plugin + " ";
                            }
                            return plugin_cmd;
                        }
                    - position: 1
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: |-
                        ${
                          var cache_bundle = inputs.vep_cache.basename
                          return 'mkdir ./cache' + ' && tar -xvf ' + cache_bundle + ' -C ./cache &&'
                        }
                    - position: 10
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: |-
                        ${
                            if (inputs.vep_plugin_files.length != 0) {
                                return '--dir_plugins ./plugins'
                            } else {
                                return ''
                            }
                        }
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: 52000
                      coresMin: 2
                    - class: DockerRequirement
                      dockerPull: 'ensemblorg/ensembl-vep:latest'
                    - class: InitialWorkDirRequirement
                      listing:
                        - $(inputs.vep_plugin_files)
                        - $(inputs.vep_cache)
                    - class: InlineJavascriptRequirement
                  'sbg:appVersion':
                    - v1.0
                  'sbg:categories':
                    - Annotation
                    - VCF Processing
                  'sbg:content_hash': >-
                    a2740e89ad63e3f545ded291672692e7f26815fb735a1201a9bdaa9456e0cf8d5
                  'sbg:contributors':
                    - rbowen_james
                  'sbg:createdBy': rbowen_james
                  'sbg:createdOn': 1625298547
                  'sbg:id': mwonge/mwtest/vep-with-plugins/26
                  'sbg:image_url': null
                  'sbg:latestRevision': 26
                  'sbg:license': Apache 2.0
                  'sbg:modifiedBy': rbowen_james
                  'sbg:modifiedOn': 1625452274
                  'sbg:project': mwonge/mwtest
                  'sbg:projectName': zcc-cavatica
                  'sbg:publisher': sbg
                  'sbg:revision': 26
                  'sbg:revisionNotes': Documentation complete
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625298547
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625299437
                      'sbg:revision': 1
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625300385
                      'sbg:revision': 2
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625300852
                      'sbg:revision': 3
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625301288
                      'sbg:revision': 4
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625301940
                      'sbg:revision': 5
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625302471
                      'sbg:revision': 6
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625303318
                      'sbg:revision': 7
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625307811
                      'sbg:revision': 8
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625361326
                      'sbg:revision': 9
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625362194
                      'sbg:revision': 10
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625364060
                      'sbg:revision': 11
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625365221
                      'sbg:revision': 12
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625376316
                      'sbg:revision': 13
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625376668
                      'sbg:revision': 14
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625376809
                      'sbg:revision': 15
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625377569
                      'sbg:revision': 16
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625378719
                      'sbg:revision': 17
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625379429
                      'sbg:revision': 18
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625380179
                      'sbg:revision': 19
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625387871
                      'sbg:revision': 20
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625387964
                      'sbg:revision': 21
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625398755
                      'sbg:revision': 22
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625399676
                      'sbg:revision': 23
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625399871
                      'sbg:revision': 24
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625401218
                      'sbg:revision': 25
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625452274
                      'sbg:revision': 26
                      'sbg:revisionNotes': Documentation complete
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': Ensembl
                  'sbg:toolkit': VEP
                  'sbg:toolkitVersion': latest
                  'sbg:validationErrors': []
                  'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
                label: vep-with-plugins
                'sbg:x': 177.5
                'sbg:y': 406.625
              - id: bgzip
                in:
                  - id: in_file
                    source: vep_with_plugins/vep_vcf
                  - id: tabix_command
                    default: tabix
                out:
                  - id: output_gz
                  - id: output_gz_tbi
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/bgzip/83
                  baseCommand: []
                  inputs:
                    - id: in_file
                      type: File
                      doc: Any file type to zip
                    - default: ''
                      id: tabix_command
                      type: string?
                      doc: >-
                        Optional, command to tabix after zip (e.g. 'tabix -p
                        vcf' or 'tabix -s1 -b2 -e2')
                    - default: ''
                      id: sort_command
                      type: string?
                      doc: >-
                        Optional, the sort command run before the zip (e.g.
                        bcftools sort). Must be an inplace sort (output is same
                        file name)
                  outputs:
                    - id: output_gz
                      type: File
                      outputBinding:
                        glob: '*.gz'
                    - id: output_gz_tbi
                      type: File
                      outputBinding:
                        glob: '*.gz.[tc][bs]i'
                  doc: >

                    # Samtools, Bcftools, HTSlib


                    These apps all utilise a samtools docker (v1.10) in order to
                    perform a variety of

                    samtools transforms and functions.


                    The docker contains (v1.10) of samtools, htslib and bcftools


                    ## Documentation


                    (Release Page)[http://www.htslib.org/download/]


                    - (BCFtools)[http://www.htslib.org/doc/bcftools.html]

                    - (bgzip)[http://www.htslib.org/doc/bgzip.html]

                    - (HTSfile)[http://www.htslib.org/doc/htsfile.html]

                    - (Samtools)[http://www.htslib.org/doc/samtools.html]

                    - (Tabix)[http://www.htslib.org/doc/tabix.html]


                    ## CWL Tools


                    ### Bamname


                    Uses samtools view to extract the sample name from the BAM
                    header


                    ### Gunzip


                    Unzips a \*.gz file that was zipped by bgzip


                    ### Tabix


                    bgzips and tabix indexes a vcf using htslib tools bgzip and
                    tabix


                    ## Dependencies


                    Access to the docker:

                    `pgc-images.sbgenomics.com/syan/samtools:1.10`


                    The docker contains the following:


                    - wget

                    - htslib (1.10.2)
                      - bgzip
                      - tabix
                    - bcftools (1.10.2)

                    - samtools (1.10)


                    ## Instance Recommendations


                    Notes:


                    - Requirements will depend entirely on your inputs, sizes,
                    and what tools you want to use

                    - The following is just a suggestion


                    - **ramMin** 1000

                    - **coresMin** 1

                    - **tmpdirMin** 10000

                    - **outdirMin** 1000
                  label: bgzip
                  arguments:
                    - position: 0
                      shellQuote: false
                      valueFrom: |
                        ${
                          let cmd = ""
                          if (inputs.sort_command) {
                            cmd += `${inputs.sort_command} ${inputs.in_file.basename}`

                            if (!inputs.in_file.basename.endsWith(".gz") || inputs.tabix_command) {
                              cmd += " && "
                            }
                          }
                          return cmd
                        }
                    - position: 1
                      shellQuote: false
                      valueFrom: |
                        ${
                          let cmd = "";
                          if (!inputs.in_file.basename.endsWith(".gz")) {
                            cmd += `bgzip ${inputs.in_file.basename}`
                            if (inputs.tabix_command) {
                              cmd += " && "
                            }
                          }
                          return cmd
                        }
                    - position: 2
                      shellQuote: false
                      valueFrom: |
                        ${
                          if (inputs.tabix_command) {
                            return `${inputs.tabix_command} *.gz`
                          }
                        }
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: 10000
                      coresMin: 4
                      outdirMin: 500000
                      tmpdirMin: 600000
                    - class: DockerRequirement
                      dockerPull: 'pgc-images.sbgenomics.com/syan/samtools:production'
                    - class: InitialWorkDirRequirement
                      listing:
                        - $(inputs.in_file)
                    - class: InlineJavascriptRequirement
                  'sbg:appVersion':
                    - v1.0
                  'sbg:content_hash': >-
                    a66ceb4796c5999a01ff1d64b9acf2ef0395a5e9100c92bdd4390de017d318c05
                  'sbg:contributors':
                    - syan
                    - lcui
                  'sbg:createdBy': syan
                  'sbg:createdOn': 1604525818
                  'sbg:id': mwonge/mwtest/bgzip/83
                  'sbg:image_url': null
                  'sbg:latestRevision': 83
                  'sbg:license': MIT License (Expat)
                  'sbg:links':
                    - id: 'http://www.htslib.org/doc/'
                      label: documentation
                    - id: 'http://www.htslib.org'
                      label: release
                    - id: 'https://tldrlegal.com/license/mit-license#summary'
                      label: license
                  'sbg:modifiedBy': syan
                  'sbg:modifiedOn': 1618553628
                  'sbg:project': mwonge/mwtest
                  'sbg:projectName': zcc-cavatica
                  'sbg:publisher': sbg
                  'sbg:revision': 83
                  'sbg:revisionNotes': |-
                    Uploaded using sbpack v2020.06.18. 
                    Source: 
                    repo: git@bitbucket.org:cciacb/cwl.git
                    file: tools/samtools/bgzip.cwl
                    commit: 5615bb2
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604525818
                      'sbg:revision': 0
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: 
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604525898
                      'sbg:revision': 1
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: 
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604526216
                      'sbg:revision': 2
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: 
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604526285
                      'sbg:revision': 3
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: 
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604534358
                      'sbg:revision': 4
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604541128
                      'sbg:revision': 5
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604543046
                      'sbg:revision': 6
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604552839
                      'sbg:revision': 7
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e63e5c4
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604619612
                      'sbg:revision': 8
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604620324
                      'sbg:revision': 9
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604893290
                      'sbg:revision': 10
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 02cfa7f
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604893619
                      'sbg:revision': 11
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604893762
                      'sbg:revision': 12
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604894010
                      'sbg:revision': 13
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604894579
                      'sbg:revision': 14
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: (uncommitted file)
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1604964302
                      'sbg:revision': 15
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 824c601
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605049094
                      'sbg:revision': 16
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 7862097
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605082172
                      'sbg:revision': 17
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: bc4fdf8
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605085484
                      'sbg:revision': 18
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 7d7ee70
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605133896
                      'sbg:revision': 19
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 7d7ee70
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605484201
                      'sbg:revision': 20
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 541cebe
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605572336
                      'sbg:revision': 21
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 47dbb50
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605576087
                      'sbg:revision': 22
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 38786f7
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605665282
                      'sbg:revision': 23
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 7718959
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605665515
                      'sbg:revision': 24
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 7718959
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1605678240
                      'sbg:revision': 25
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: b47363e
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1605678517
                      'sbg:revision': 26
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 6da6c0f
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605681517
                      'sbg:revision': 27
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 8b652ac
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605741714
                      'sbg:revision': 28
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 8b652ac
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1605742263
                      'sbg:revision': 29
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 0b4cf01
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1606864100
                      'sbg:revision': 30
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 85328e2
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1606864391
                      'sbg:revision': 31
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: a0d4aff
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1606971684
                      'sbg:revision': 32
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 684d1ed
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607035058
                      'sbg:revision': 33
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 6d018e4
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607040613
                      'sbg:revision': 34
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: a98fa70
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607049954
                      'sbg:revision': 35
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 578b0de
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607316129
                      'sbg:revision': 36
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: ef0c2a5
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607319791
                      'sbg:revision': 37
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 7a8240d
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1607902133
                      'sbg:revision': 38
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: d59f4db
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607906090
                      'sbg:revision': 39
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 5836465
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607906827
                      'sbg:revision': 40
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: ef6ed2b
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1607907893
                      'sbg:revision': 41
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 7f0bdb4
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1607911306
                      'sbg:revision': 42
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: d04ae7d
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1607914409
                      'sbg:revision': 43
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 18dfb49
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1608075465
                      'sbg:revision': 44
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 18dfb49
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1608077672
                      'sbg:revision': 45
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 56beeb1
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1610424224
                      'sbg:revision': 46
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: ecdc528
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1611035954
                      'sbg:revision': 47
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 373907c
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1611098610
                      'sbg:revision': 48
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 9e52f26
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613085598
                      'sbg:revision': 49
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 6e0b99b
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613525098
                      'sbg:revision': 50
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 3e6f0d4
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1613624656
                      'sbg:revision': 51
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 773b440
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1613627542
                      'sbg:revision': 52
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e5bdc2b
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1613687493
                      'sbg:revision': 53
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 32bb0c3
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613704587
                      'sbg:revision': 54
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e670ba9
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613704907
                      'sbg:revision': 55
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e670ba9
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613705148
                      'sbg:revision': 56
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e6a7d65
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613705288
                      'sbg:revision': 57
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e6a7d65
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613705406
                      'sbg:revision': 58
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: befd8b1
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1613706012
                      'sbg:revision': 59
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 48965cf
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1613718333
                      'sbg:revision': 60
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 87d785c
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1613968190
                      'sbg:revision': 61
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 57998bf
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1613969213
                      'sbg:revision': 62
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: faf8958
                    - 'sbg:modifiedBy': lcui
                      'sbg:modifiedOn': 1614034778
                      'sbg:revision': 63
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: abe98a3
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614048474
                      'sbg:revision': 64
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.10.05. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e706c86
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614049803
                      'sbg:revision': 65
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: ddf9af6
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614298538
                      'sbg:revision': 66
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: ae3a543
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614304889
                      'sbg:revision': 67
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: b48195d
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614305581
                      'sbg:revision': 68
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: b1bc92f
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614557653
                      'sbg:revision': 69
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: c6d3eec
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614655430
                      'sbg:revision': 70
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 4edb2c7
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614657413
                      'sbg:revision': 71
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: e8821d7
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1614728728
                      'sbg:revision': 72
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 8aaf9f9
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1615166222
                      'sbg:revision': 73
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: b6911aa
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1615177323
                      'sbg:revision': 74
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: cf2acf6
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1615256201
                      'sbg:revision': 75
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 9637fb6
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1615335986
                      'sbg:revision': 76
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 19905f3
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1616024266
                      'sbg:revision': 77
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 489075b
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1616369886
                      'sbg:revision': 78
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: fcbfa76
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1616470567
                      'sbg:revision': 79
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: a46e093
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1616538066
                      'sbg:revision': 80
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 4ee52e6
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1617254150
                      'sbg:revision': 81
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 4726466
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1618536815
                      'sbg:revision': 82
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: f3f5fd6
                    - 'sbg:modifiedBy': syan
                      'sbg:modifiedOn': 1618553628
                      'sbg:revision': 83
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.06.18. 
                        Source: 
                        repo: git@bitbucket.org:cciacb/cwl.git
                        file: tools/samtools/bgzip.cwl
                        commit: 5615bb2
                  'sbg:sbgMaintained': false
                  'sbg:validationErrors': []
                label: bgzip
                'sbg:x': 489.11181640625
                'sbg:y': 587.984375
              - id: vcf_tumour_name
                in:
                  - id: vcf
                    source: input_file
                out:
                  - id: tumour_name
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/vcf-tumour-name/7
                  baseCommand: []
                  inputs:
                    - id: vcf
                      type: File
                      'sbg:fileTypes': VCF
                  outputs:
                    - id: tumour_name
                      type: string
                      outputBinding:
                        loadContents: true
                        glob: name.txt
                        outputEval: '$(self[0].contents.trim())'
                  doc: >-
                    # vcf-tumour-name

                    Extracts the tumour sample name from a VCF.

                    Uses a script called `get_tumour_name.py` via a Docker
                    container.


                    The tumour sample is identified by:


                    1. Extracting the sample names from the `#CHROM` header
                    line.

                    2. Choosing the sample name that **does not end in
                    "G[0-9]"**. I.e. assume that the sample ending in G is the
                    germline sample, and the other sample must be tumour.



                    *Note:* This method of tumour sample identification was
                    chosen since the order of tumour and normal sample in the
                    VCF may not be uniform (also could have a single sample
                    VCF), and because tumour samples can end in D/P/R etc., but
                    germline is known to always end with G.


                    ## Docker

                    Docker image:
                    `pgc-images.sbgenomics.com/rbowen_james/vcf-tumour-name`


                    ## Inputs

                    `vcf`: The VCF to extract the tumour sample name from. Can
                    be single or multi-sample. **Must contain only one tumour
                    sample.**


                    ## Outputs

                    `tumour_name`: The tumour sample name as a string.
                  label: vcf-tumour-name
                  arguments:
                    - position: 0
                      prefix: ''
                      shellQuote: false
                      valueFrom: >-
                        python3 /app/get_tumour_name.py $(inputs.vcf.path) >
                        name.txt
                  requirements:
                    - class: ShellCommandRequirement
                    - class: DockerRequirement
                      dockerPull: >-
                        pgc-images.sbgenomics.com/rbowen_james/vcf-tumour-name:latest
                    - class: InlineJavascriptRequirement
                  'sbg:appVersion':
                    - v1.0
                  'sbg:content_hash': >-
                    a8c22ba11fb17947632e42fa7f33c80c02728115b7500ca32c2c3cb52bfd9dbf8
                  'sbg:contributors':
                    - rbowen_james
                  'sbg:createdBy': rbowen_james
                  'sbg:createdOn': 1614666562
                  'sbg:id': mwonge/mwtest/vcf-tumour-name/7
                  'sbg:image_url': null
                  'sbg:latestRevision': 7
                  'sbg:modifiedBy': rbowen_james
                  'sbg:modifiedOn': 1626063981
                  'sbg:project': mwonge/mwtest
                  'sbg:projectName': zcc-cavatica
                  'sbg:publisher': sbg
                  'sbg:revision': 7
                  'sbg:revisionNotes': ''
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1614666562
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1614666868
                      'sbg:revision': 1
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1614667522
                      'sbg:revision': 2
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1614722294
                      'sbg:revision': 3
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1615261751
                      'sbg:revision': 4
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1615262123
                      'sbg:revision': 5
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625804786
                      'sbg:revision': 6
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626063981
                      'sbg:revision': 7
                      'sbg:revisionNotes': ''
                  'sbg:sbgMaintained': false
                  'sbg:validationErrors': []
                label: vcf-tumour-name
                'sbg:x': 177.5
                'sbg:y': 534.53125
              - id: pvacseq
                in:
                  - id: input_vcf
                    source: bgzip/output_gz
                  - id: sample_name
                    source: vcf_tumour_name/tumour_name
                  - id: alleles
                    source: alleles
                  - id: prediction_algorithms
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
                  - id: allele_specific_binding_thresholds
                    source: allele_specific_binding_thresholds
                  - id: iedb_retries
                    source: iedb_retries
                  - id: keep_tmp_files
                    source: keep_tmp_files
                  - id: normal_sample_name
                    source: normal_sample_name
                  - id: net_chop_method
                    source: net_chop_method
                  - id: netmhc_stab
                    source: netmhc_stab
                  - id: run_reference_proteome_similarity
                    source: run_reference_proteome_similarity
                  - id: top_score_metric
                    source: top_score_metric
                  - id: net_chop_threshold
                    source: net_chop_threshold
                  - id: additional_report_columns
                    source: additional_report_columns
                  - id: fasta_size
                    source: fasta_size
                  - id: downstream_sequence_length
                    source: downstream_sequence_length
                  - id: exclude_nas
                    source: exclude_nas
                  - id: phased_proximal_variants_vcf
                    source: phased_proximal_variants_vcf
                  - id: minimum_fold_change
                    source: minimum_fold_change
                  - id: normal_cov
                    source: normal_cov
                  - id: tdna_cov
                    source: tdna_cov
                  - id: trna_cov
                    source: trna_cov
                  - id: normal_vaf
                    source: normal_vaf
                  - id: tdna_vaf
                    source: tdna_vaf
                  - id: trna_vaf
                    source: trna_vaf
                  - id: expn_val
                    source: expn_val
                  - id: maximum_transcript_support_level
                    source: maximum_transcript_support_level
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
                  - id: pvacseq_predictions
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/pvacseq/21
                  baseCommand:
                    - ln
                    - '-s'
                  inputs:
                    - id: input_vcf
                      type: File
                      inputBinding:
                        position: 1
                        shellQuote: false
                      secondaryFiles:
                        - .tbi
                    - id: sample_name
                      type: string
                      inputBinding:
                        position: 2
                        shellQuote: false
                    - id: alleles
                      type: string
                      inputBinding:
                        position: 3
                        prefix: ''
                        separate: false
                        itemSeparator: ','
                        shellQuote: false
                    - id: prediction_algorithms
                      type: 'string[]'
                      inputBinding:
                        position: 4
                        shellQuote: false
                    - id: epitope_lengths_class_i
                      type: 'int[]?'
                      inputBinding:
                        position: 0
                        prefix: '-e1'
                        itemSeparator: ','
                        shellQuote: false
                    - id: epitope_lengths_class_ii
                      type: 'int[]?'
                      inputBinding:
                        position: 0
                        prefix: '-e2'
                        itemSeparator: ','
                        shellQuote: false
                    - id: binding_threshold
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '-b'
                        shellQuote: false
                    - id: percentile_threshold
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '--percentile-threshold'
                        shellQuote: false
                    - id: allele_specific_binding_thresholds
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '--allele-specific-binding-thresholds'
                        shellQuote: false
                    - id: iedb_retries
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '-r'
                        shellQuote: false
                    - id: keep_tmp_files
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '-k'
                        shellQuote: false
                    - id: normal_sample_name
                      type: string?
                      inputBinding:
                        position: 0
                        prefix: '--normal-sample-name'
                        shellQuote: false
                    - id: net_chop_method
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - cterm
                            - 20s
                          name: net_chop_method
                      inputBinding:
                        position: 0
                        prefix: '--net-chop-method'
                        shellQuote: false
                    - id: netmhc_stab
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '--netmhc-stab'
                        shellQuote: false
                    - id: run_reference_proteome_similarity
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '--run-reference-proteome-similarity'
                        shellQuote: false
                    - id: top_score_metric
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - lowest
                            - median
                          name: top_score_metric
                      inputBinding:
                        position: 0
                        prefix: '-m'
                        shellQuote: false
                    - id: net_chop_threshold
                      type: float?
                      inputBinding:
                        position: 0
                        prefix: '--net-chop-threshold'
                        shellQuote: false
                    - id: additional_report_columns
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - sample_name
                          name: additional_report_columns
                      inputBinding:
                        position: 0
                        prefix: '-a'
                        shellQuote: false
                    - id: fasta_size
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '-s'
                        shellQuote: false
                    - id: downstream_sequence_length
                      type: string?
                      inputBinding:
                        position: 0
                        prefix: '-d'
                        shellQuote: false
                    - id: exclude_nas
                      type: boolean?
                      inputBinding:
                        position: 0
                        prefix: '--exclude-NAs'
                        shellQuote: false
                    - id: phased_proximal_variants_vcf
                      type: File?
                      inputBinding:
                        position: 0
                        prefix: '-p'
                        shellQuote: false
                      secondaryFiles:
                        - .tbi
                    - id: minimum_fold_change
                      type: float?
                      inputBinding:
                        position: 0
                        prefix: '-c'
                        shellQuote: false
                    - id: normal_cov
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '--normal-cov'
                        shellQuote: false
                    - id: tdna_cov
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '--tdna-cov'
                        shellQuote: false
                    - id: trna_cov
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '--trna-cov'
                        shellQuote: false
                    - id: normal_vaf
                      type: float?
                      inputBinding:
                        position: 0
                        prefix: '--normal-vaf'
                        shellQuote: false
                    - id: tdna_vaf
                      type: float?
                      inputBinding:
                        position: 0
                        prefix: '--tdna-vaf'
                        shellQuote: false
                    - id: trna_vaf
                      type: float?
                      inputBinding:
                        position: 0
                        prefix: '--trna-vaf'
                        shellQuote: false
                    - id: expn_val
                      type: float?
                      inputBinding:
                        position: 0
                        prefix: '--expn-val'
                        shellQuote: false
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
                      inputBinding:
                        position: 0
                        prefix: '--maximum-transcript-support-level'
                        shellQuote: false
                    - default: 8
                      id: n_threads
                      type: int?
                      inputBinding:
                        position: 0
                        prefix: '--n-threads'
                        shellQuote: false
                  outputs:
                    - id: mhc_i_all_epitopes
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/MHC_Class_I/$(inputs.sample_name).all_epitopes.tsv
                    - id: mhc_i_aggregated_report
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/MHC_Class_I/$(inputs.sample_name).all_epitopes.aggregated.tsv
                    - id: mhc_i_filtered_epitopes
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/MHC_Class_I/$(inputs.sample_name).filtered.tsv
                    - id: mhc_ii_all_epitopes
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/MHC_Class_II/$(inputs.sample_name).all_epitopes.tsv
                    - id: mhc_ii_aggregated_report
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/MHC_Class_II/$(inputs.sample_name).all_epitopes.aggregated.tsv
                    - id: mhc_ii_filtered_epitopes
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/MHC_Class_II/$(inputs.sample_name).filtered.tsv
                    - id: combined_all_epitopes
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/combined/$(inputs.sample_name).all_epitopes.tsv
                    - id: combined_aggregated_report
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/combined/$(inputs.sample_name).all_epitopes.aggregated.tsv
                    - id: combined_filtered_epitopes
                      type: File?
                      outputBinding:
                        glob: >-
                          pvacseq_predictions/combined/$(inputs.sample_name).filtered.tsv
                    - id: pvacseq_predictions
                      type: Directory
                      outputBinding:
                        glob: pvacseq_predictions
                  label: pvacseq
                  arguments:
                    - position: 0
                      shellQuote: false
                      valueFrom: $TMPDIR
                    - /tmp/pvacseq
                    - position: 0
                      shellQuote: false
                      valueFrom: ' && '
                    - export
                    - TMPDIR=/tmp/pvacseq
                    - position: 0
                      shellQuote: false
                      valueFrom: ' && '
                    - /usr/local/bin/pvacseq
                    - run
                    - '--iedb-install-directory'
                    - /opt/iedb
                    - '--pass-only'
                    - position: 5
                      shellQuote: false
                      valueFrom: pvacseq_predictions
                    - position: 10
                      prefix: ''
                      shellQuote: false
                      valueFrom: '&& echo "--------" 1>&2 && ls 1>&2'
                    - position: 9
                      prefix: ''
                      shellQuote: false
                      valueFrom: '|| true'
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: 16000
                      coresMin: $(inputs.n_threads)
                    - class: DockerRequirement
                      dockerPull: 'griffithlab/pvactools:2.0.3'
                    - class: InlineJavascriptRequirement
                  'sbg:appVersion':
                    - v1.0
                  'sbg:content_hash': >-
                    a349a251986c6258ef34edff34b03ea30e2612957203caf08356acc94322e2944
                  'sbg:contributors':
                    - rbowen_james
                  'sbg:createdBy': rbowen_james
                  'sbg:createdOn': 1614652262
                  'sbg:id': mwonge/mwtest/pvacseq/21
                  'sbg:image_url': null
                  'sbg:latestRevision': 21
                  'sbg:modifiedBy': rbowen_james
                  'sbg:modifiedOn': 1628034717
                  'sbg:project': mwonge/mwtest
                  'sbg:projectName': zcc-cavatica
                  'sbg:publisher': sbg
                  'sbg:revision': 21
                  'sbg:revisionNotes': Updated pVACtools version to 2.0.3
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1614652262
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1614653222
                      'sbg:revision': 1
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: ./pvacseq.cwl
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1614653309
                      'sbg:revision': 2
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: ./pvacseq.cwl
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1625561275
                      'sbg:revision': 3
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626054029
                      'sbg:revision': 4
                      'sbg:revisionNotes': Update pVACseq version to 2.0.2
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626056407
                      'sbg:revision': 5
                      'sbg:revisionNotes': Removed '--pass-only'
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626056797
                      'sbg:revision': 6
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626056880
                      'sbg:revision': 7
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626063867
                      'sbg:revision': 8
                      'sbg:revisionNotes': changed allele string
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626064345
                      'sbg:revision': 9
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626064349
                      'sbg:revision': 10
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626064744
                      'sbg:revision': 11
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626064964
                      'sbg:revision': 12
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626065392
                      'sbg:revision': 13
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626067957
                      'sbg:revision': 14
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626068799
                      'sbg:revision': 15
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1626133601
                      'sbg:revision': 16
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1628034393
                      'sbg:revision': 17
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: ./pvacseq.cwl
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1628034406
                      'sbg:revision': 18
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: ./pvacfuse.cwl
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1628034410
                      'sbg:revision': 19
                      'sbg:revisionNotes': |-
                        Uploaded using sbpack v2020.09.13. 
                        Source: ./pvacseq.cwl
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1628034684
                      'sbg:revision': 20
                      'sbg:revisionNotes': Updated pVACtools version to 2.0.3
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1628034717
                      'sbg:revision': 21
                      'sbg:revisionNotes': Updated pVACtools version to 2.0.3
                  'sbg:sbgMaintained': false
                  'sbg:validationErrors': []
                label: pvacseq
                'sbg:x': 782.31494140625
                'sbg:y': 311.171875
            requirements: []
            'sbg:appVersion':
              - v1.0
            'sbg:content_hash': a430996f62b6d20e61b31cec5636968a23e5f363bf2f89046f3f8a746c5921839
            'sbg:contributors':
              - rbowen_james
            'sbg:createdBy': rbowen_james
            'sbg:createdOn': 1625804986
            'sbg:id': mwonge/mwtest/vep-pvactools/4
            'sbg:image_url': >-
              https://cavatica.sbgenomics.com/ns/brood/images/mwonge/mwtest/vep-pvactools/4.png
            'sbg:latestRevision': 4
            'sbg:modifiedBy': rbowen_james
            'sbg:modifiedOn': 1628034876
            'sbg:project': mwonge/mwtest
            'sbg:projectName': zcc-cavatica
            'sbg:publisher': sbg
            'sbg:revision': 4
            'sbg:revisionNotes': Updated pVACtools version to 2.0.3
            'sbg:revisionsInfo':
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1625804986
                'sbg:revision': 0
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626066538
                'sbg:revision': 1
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626066560
                'sbg:revision': 2
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626831039
                'sbg:revision': 3
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1628034876
                'sbg:revision': 4
                'sbg:revisionNotes': Updated pVACtools version to 2.0.3
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
          label: vep-pvactools
          'sbg:x': 313.03125
          'sbg:y': 1068.9609375
      requirements:
        - class: SubworkflowFeatureRequirement
      'sbg:image_url': >-
        https://cavatica.sbgenomics.com/ns/brood/images/mwonge/mwtest/pvacseq-pvacfuse/4.png
      'sbg:projectName': zcc-cavatica
      'sbg:revisionsInfo':
        - 'sbg:revision': 0
          'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627302527
          'sbg:revisionNotes': null
        - 'sbg:revision': 1
          'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627307019
          'sbg:revisionNotes': null
        - 'sbg:revision': 2
          'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627533435
          'sbg:revisionNotes': null
        - 'sbg:revision': 3
          'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1628034926
          'sbg:revisionNotes': Updated pVACtools version to 2.0.3
        - 'sbg:revision': 4
          'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1628035201
          'sbg:revisionNotes': Updated pVACtools version to 2.0.3
      'sbg:appVersion':
        - v1.0
      'sbg:id': mwonge/mwtest/pvacseq-pvacfuse/4
      'sbg:revision': 4
      'sbg:revisionNotes': Updated pVACtools version to 2.0.3
      'sbg:modifiedOn': 1628035201
      'sbg:modifiedBy': rbowen_james
      'sbg:createdOn': 1627302527
      'sbg:createdBy': rbowen_james
      'sbg:project': mwonge/mwtest
      'sbg:sbgMaintained': false
      'sbg:validationErrors': []
      'sbg:contributors':
        - rbowen_james
      'sbg:latestRevision': 4
      'sbg:publisher': sbg
      'sbg:content_hash': ae1130d6c90b45dbedb2434062fe5580dfffbf63eced1176e77d9f9a923227fa1
    label: pvacseq-pvacfuse
    'sbg:x': 1265.7939453125
    'sbg:y': 1004.953125
  - id: hlatyping_bowtie_samtools
    in:
      - id: normal_dna_reads
        source:
          - normal_dna_reads
      - id: bowtie_index_archive
        source: bowtie_index_archive
      - id: tumour_dna_reads
        source:
          - tumour_dna_reads
      - id: minimum_read_length
        default: 100
        source: minimum_read_length
      - id: sample_id
        source: sample_id
      - id: minimum_read_length_1
        default: 100
        source: minimum_read_length_1
      - id: minimum_read_length_2
        default: 100
        source: minimum_read_length_2
      - id: tumour_rna_reads_1
        source: tumour_rna_reads_1
      - id: tumour_rna_reads_2
        source: tumour_rna_reads_2
      - id: tumour_dna_outdir_name
        source: tumour_dna_outdir_name
      - id: normal_dna_outdir_name
        source: normal_dna_outdir_name
      - id: tumour_rna_outdir_name
        source: tumour_rna_outdir_name
    out:
      - id: tumour_rna_hla-hd_results
      - id: tumour_rna_hla-hd_final
      - id: hlad_reads1
      - id: hla_reads2
      - id: filtered_BAM
      - id: bowtie_sam
      - id: aligned_reads_only
      - id: unaligned_reads_only
      - id: tumour_dna_hla-hd_results
      - id: tumour_dna_hla-hd_final_result
      - id: hlad_reads2
      - id: hla_reads3
      - id: filtered_BAM_1
      - id: bowtie_sam_1
      - id: aligned_reads_only_1
      - id: unaligned_reads_only_1
      - id: normal_dna_hla-hd_results
      - id: normal_dna_hla-hd_final_result
    run:
      class: Workflow
      cwlVersion: v1.0
      id: mwonge/mwtest/hlatyping-bowtie-samtools/6
      label: hlatyping_bowtie_samtools
      $namespaces:
        sbg: 'https://sevenbridges.com'
      inputs:
        - id: normal_dna_reads
          'sbg:fileTypes': >-
            FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA, FQ, FASTQ.GZ,
            FQ.GZ, FASTQ.BZ2, FQ.BZ2
          type: 'File[]'
          label: Normal DNA Reads
          'sbg:x': 0
          'sbg:y': 1088.2265625
        - id: bowtie_index_archive
          'sbg:fileTypes': TAR
          type: File
          'sbg:x': 0
          'sbg:y': 1301.6796875
        - id: tumour_dna_reads
          'sbg:fileTypes': >-
            FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA, FQ, FASTQ.GZ,
            FQ.GZ, FASTQ.BZ2, FQ.BZ2
          type: 'File[]'
          label: Tumour DNA Reads
          'sbg:x': 0
          'sbg:y': 768.0546875
        - id: minimum_read_length
          type: int
          'sbg:exposed': true
        - id: sample_id
          type: string
          'sbg:x': 0
          'sbg:y': 981.5078125
        - id: minimum_read_length_1
          type: int
          'sbg:exposed': true
        - id: minimum_read_length_2
          type: int
          'sbg:exposed': true
        - id: tumour_rna_reads_1
          type: File
          label: Tumour RNA Reads 1
          'sbg:x': 0
          'sbg:y': 554.6015625
        - id: tumour_rna_reads_2
          type: File
          label: Tumour RNA Reads 2
          'sbg:x': 0
          'sbg:y': 447.8828125
        - id: tumour_dna_outdir_name
          type: string
          label: Tumour DNA Output Directory Name
          'sbg:x': 0
          'sbg:y': 874.78125
        - id: normal_dna_outdir_name
          type: string
          label: Normal DNA Output Directory Name
          'sbg:x': 0
          'sbg:y': 1194.953125
        - id: tumour_rna_outdir_name
          type: string
          label: Tumour RNA Output Directory Name
          'sbg:x': 0
          'sbg:y': 661.328125
      outputs:
        - id: tumour_rna_hla-hd_results
          outputSource:
            - hla_hd_tumour_rna/hlahd_results
          type: Directory
          label: Tumour RNA HLA-HD Results
          'sbg:x': 974.7388916015625
          'sbg:y': 821.421875
        - id: tumour_rna_hla-hd_final
          outputSource:
            - hla_hd_tumour_rna/hlahd_final_results
          type: File
          label: Tumour RNA HLA-HD Final Result
          'sbg:x': 974.7388916015625
          'sbg:y': 928.140625
        - id: hlad_reads1
          outputSource:
            - hla-hd_tumour_dna/hlad_reads1
          'sbg:fileTypes': FASTQ
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 747.0546875
        - id: hla_reads2
          outputSource:
            - hla-hd_tumour_dna/hla_reads2
          'sbg:fileTypes': FASTQ
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 960.4921875
        - id: filtered_BAM
          outputSource:
            - hla-hd_tumour_dna/filtered_BAM
          'sbg:fileTypes': 'BAM, SAM, CRAM'
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 1322.6484375
        - id: bowtie_sam
          outputSource:
            - hla-hd_tumour_dna/bowtie_sam
          'sbg:fileTypes': SAM
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 1536.0859375
        - id: aligned_reads_only
          outputSource:
            - hla-hd_tumour_dna/aligned_reads_only
          'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
          type: 'File[]?'
          'sbg:x': 659.8301391601562
          'sbg:y': 1749.546875
        - id: unaligned_reads_only
          outputSource:
            - hla-hd_tumour_dna/unaligned_reads_only
          'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
          type: 'File[]?'
          'sbg:x': 659.8301391601562
          'sbg:y': 106.734375
        - id: tumour_dna_hla-hd_results
          outputSource:
            - hla-hd_tumour_dna/hlahd_results
          type: Directory
          label: Tumour DNA HLA-HD Results
          'sbg:x': 659.8301391601562
          'sbg:y': 213.4609375
        - id: tumour_dna_hla-hd_final_result
          outputSource:
            - hla-hd_tumour_dna/hlahd_final_results
          type: File
          label: Tumour DNA HLA-HD Final Result
          'sbg:x': 659.8301391601562
          'sbg:y': 320.1796875
        - id: hlad_reads2
          outputSource:
            - hla-hd_normal_dna/hlad_reads1
          'sbg:fileTypes': FASTQ
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 640.3359375
        - id: hla_reads3
          outputSource:
            - hla-hd_normal_dna/hla_reads2
          'sbg:fileTypes': FASTQ
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 853.7734375
        - id: filtered_BAM_1
          outputSource:
            - hla-hd_normal_dna/filtered_BAM
          'sbg:fileTypes': 'BAM, SAM, CRAM'
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 1215.9296875
        - id: bowtie_sam_1
          outputSource:
            - hla-hd_normal_dna/bowtie_sam
          'sbg:fileTypes': SAM
          type: File?
          'sbg:x': 659.8301391601562
          'sbg:y': 1429.3671875
        - id: aligned_reads_only_1
          outputSource:
            - hla-hd_normal_dna/aligned_reads_only
          'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
          type: 'File[]?'
          'sbg:x': 659.8301391601562
          'sbg:y': 1642.8125
        - id: unaligned_reads_only_1
          outputSource:
            - hla-hd_normal_dna/unaligned_reads_only
          'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
          type: 'File[]?'
          'sbg:x': 659.8301391601562
          'sbg:y': 0
        - id: normal_dna_hla-hd_results
          outputSource:
            - hla-hd_normal_dna/hlahd_results
          type: Directory
          label: Normal DNA HLA-HD Results
          'sbg:x': 659.8301391601562
          'sbg:y': 426.8984375
        - id: normal_dna_hla-hd_final_result
          outputSource:
            - hla-hd_normal_dna/hlahd_final_results
          type: File
          label: Normal DNA HLA-HD Final Result
          'sbg:x': 659.8301391601562
          'sbg:y': 533.6171875
      steps:
        - id: hla-hd_tumour_dna
          in:
            - id: sample_id
              source: sample_id
            - id: minimum_read_length
              default: 0
              source: minimum_read_length_1
            - id: read_sequence
              source:
                - tumour_dna_reads
            - id: bowtie_index_archive
              source: bowtie_index_archive
            - id: output_dir_name
              source: tumour_dna_outdir_name
          out:
            - id: hlad_reads1
            - id: hla_reads2
            - id: filtered_BAM
            - id: bowtie_sam
            - id: aligned_reads_only
            - id: unaligned_reads_only
            - id: hlahd_results
            - id: hlahd_final_results
          run:
            class: Workflow
            cwlVersion: v1.0
            id: mwonge/mwtest/hlahd-bowtie-samtools/6
            label: hlahd_bowtie_samtools
            $namespaces:
              sbg: 'https://sevenbridges.com'
            inputs:
              - id: sample_id
                type: string
                label: Sample ID
                doc: ID of the sample being analysed.
                'sbg:x': 0
                'sbg:y': 14
              - id: minimum_read_length
                type: int
                'sbg:exposed': true
              - id: read_sequence
                'sbg:fileTypes': >-
                  FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA, FQ,
                  FASTQ.GZ, FQ.GZ, FASTQ.BZ2, FQ.BZ2
                type: 'File[]'
                'sbg:x': 0
                'sbg:y': 120.75
              - id: bowtie_index_archive
                'sbg:fileTypes': TAR
                type: File
                'sbg:x': 0
                'sbg:y': 334.25
              - id: output_dir_name
                type: string?
                'sbg:x': 0
                'sbg:y': 227.5
            outputs:
              - id: hlad_reads1
                outputSource:
                  - samtools_fastq/output_pe_1
                'sbg:fileTypes': FASTQ
                type: File?
                label: HLA Reads 1
                doc: >-
                  Reads 1 FASTQ containing HLA reads after Yara alignment,
                  Samtools View filtering and splitting into FASTQ files.
                'sbg:x': 1256.66357421875
                'sbg:y': 46.25
              - id: hla_reads2
                outputSource:
                  - samtools_fastq/output_pe_2
                'sbg:fileTypes': FASTQ
                type: File?
                label: HLA Reads 2
                doc: >-
                  Reads 2 FASTQ containing HLA reads after Yara alignment,
                  Samtools View filtering and splitting into FASTQ files.
                'sbg:x': 1256.66357421875
                'sbg:y': 153
              - id: filtered_BAM
                outputSource:
                  - samtools_view_1_9_cwl1_0/out_alignments
                'sbg:fileTypes': 'BAM, SAM, CRAM'
                type: File?
                label: Filtered BAM
                doc: Yara aligned BAM filtered by Samtools View.
                'sbg:x': 977.3401489257812
                'sbg:y': 227.5
              - id: bowtie_sam
                outputSource:
                  - bowtie2_aligner/result_sam_file
                'sbg:fileTypes': SAM
                type: File?
                'sbg:x': 547.8245239257812
                'sbg:y': 241.5
              - id: aligned_reads_only
                outputSource:
                  - bowtie2_aligner/aligned_reads_only
                'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                type: 'File[]?'
                'sbg:x': 547.8245239257812
                'sbg:y': 348.25
              - id: unaligned_reads_only
                outputSource:
                  - bowtie2_aligner/unaligned_reads_only
                'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                type: 'File[]?'
                'sbg:x': 547.8245239257812
                'sbg:y': 0
              - id: hlahd_results
                outputSource:
                  - hla_hd/hlahd_results
                type: Directory
                'sbg:x': 1571.572265625
                'sbg:y': 120.75
              - id: hlahd_final_results
                outputSource:
                  - hla_hd/hlahd_final_results
                type: File
                'sbg:x': 1571.572265625
                'sbg:y': 227.5
            steps:
              - id: samtools_view_1_9_cwl1_0
                in:
                  - id: output_format
                    default: BAM
                  - id: fast_bam_compression
                    default: true
                  - id: filter_exclude_any
                    default: 4
                  - id: in_alignments
                    source: bowtie2_aligner/result_sam_file
                out:
                  - id: out_alignments
                  - id: reads_not_selected_by_filters
                  - id: alignement_count
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: admin/sbg-public-data/samtools-view-1-9-cwl1-0/6
                  baseCommand:
                    - /opt/samtools-1.9/samtools
                    - view
                  inputs:
                    - 'sbg:category': File inputs
                      id: in_index
                      type: File?
                      label: Index file
                      doc: This tool requires index file for some use cases.
                      'sbg:fileTypes': 'BAI, CRAI, CSI'
                    - 'sbg:altPrefix': '-O'
                      'sbg:category': Config inputs
                      'sbg:toolDefaultValue': SAM
                      id: output_format
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - SAM
                            - BAM
                            - CRAM
                          name: output_format
                      inputBinding:
                        position: 1
                        prefix: '--output-fmt'
                        shellQuote: false
                      label: Output format
                      doc: Output file format
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: fast_bam_compression
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-1'
                        shellQuote: false
                      label: Fast BAM compression
                      doc: >-
                        Enable fast BAM compression (implies output in bam
                        format).
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: uncompressed_bam
                      type: boolean?
                      inputBinding:
                        position: 3
                        prefix: '-u'
                        shellQuote: false
                      label: Output uncompressed BAM
                      doc: >-
                        Output uncompressed BAM (implies output BAM format).
                        This option saves time spent on
                        compression/decompression and is thus preferred when the
                        output is piped to another SAMtools command.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: include_header
                      type: boolean?
                      inputBinding:
                        position: 4
                        prefix: '-h'
                        shellQuote: false
                      label: Include the header in the output
                      doc: Include the header in the output.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: output_header_only
                      type: boolean?
                      inputBinding:
                        position: 5
                        prefix: '-H'
                        shellQuote: false
                      label: Output the header only
                      doc: Output the header only.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: collapse_cigar
                      type: boolean?
                      inputBinding:
                        position: 6
                        prefix: '-B'
                        shellQuote: false
                      label: Collapse the backward CIGAR operation
                      doc: Collapse the backward CIGAR operation.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_include
                      type: int?
                      inputBinding:
                        position: 7
                        prefix: '-f'
                        shellQuote: false
                      label: Include reads with all of these flags
                      doc: >-
                        Only output alignments with all bits set in this integer
                        present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_any
                      type: int?
                      inputBinding:
                        position: 8
                        prefix: '-F'
                        shellQuote: false
                      label: Exclude reads with any of these flags
                      doc: >-
                        Do not output alignments with any bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_all
                      type: int?
                      inputBinding:
                        position: 9
                        prefix: '-G'
                        shellQuote: false
                      label: Exclude reads with all of these flags
                      doc: >-
                        Only exclude reads with all of the bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'null'
                      id: read_group
                      type: string?
                      inputBinding:
                        position: 10
                        prefix: '-r'
                        shellQuote: false
                      label: Read group
                      doc: Only output reads in the specified read group.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_mapq
                      type: int?
                      inputBinding:
                        position: 11
                        prefix: '-q'
                        shellQuote: false
                      label: Minimum mapping quality
                      doc: Skip alignments with MAPQ smaller than this value.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'null'
                      id: filter_library
                      type: string?
                      inputBinding:
                        position: 12
                        prefix: '-l'
                        shellQuote: false
                      label: Only include alignments in library
                      doc: Only output alignments in this library.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: min_cigar_operations
                      type: int?
                      inputBinding:
                        position: 13
                        prefix: '-m'
                        shellQuote: false
                      label: Minimum number of CIGAR bases consuming query sequence
                      doc: >-
                        Only output alignments with number of CIGAR bases
                        consuming query sequence  ≥ INT.
                    - 'sbg:category': Config Inputs
                      id: read_tag_to_strip
                      type: 'string[]?'
                      inputBinding:
                        position: 14
                        prefix: ''
                        itemSeparator: ' '
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (self)
                              {
                                  var cmd = [];
                                  for (var i = 0; i < self.length; i++) 
                                  {
                                      cmd.push('-x', self[i]);
                                      
                                  }
                                  return cmd.join(' ');
                              }
                          }
                      label: Read tags to strip
                      doc: Read tag to exclude from output (repeatable).
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: count_alignments
                      type: boolean?
                      inputBinding:
                        position: 15
                        prefix: '-c'
                        shellQuote: false
                      label: Output only count of matching records
                      doc: >-
                        Instead of outputing the alignments, only count them and
                        output the total number. All filter options, such as -f,
                        -F, and -q, are taken into account.
                    - 'sbg:category': Config Inputs
                      id: input_fmt_option
                      type: string?
                      inputBinding:
                        position: 16
                        prefix: '--input-fmt-option'
                        shellQuote: false
                      label: Input file format option
                      doc: >-
                        Specify a single input file format option in the form of
                        OPTION or OPTION=VALUE.
                    - 'sbg:category': Config Inputs
                      id: output_fmt_option
                      type: string?
                      inputBinding:
                        position: 17
                        prefix: '--output-fmt-option'
                        shellQuote: false
                      label: Output file format option
                      doc: >-
                        Specify a single output file format option in the form
                        of OPTION or OPTION=VALUE.
                    - 'sbg:category': Config Inputs
                      id: subsample_fraction
                      type: float?
                      inputBinding:
                        position: 18
                        prefix: '-s'
                        shellQuote: false
                      label: Subsample fraction
                      doc: >-
                        Output only a proportion of the input alignments. This
                        subsampling acts in the same way on all of the alignment
                        records in the same template or read pair, so it never
                        keeps a read but not its mate. The integer and
                        fractional parts of the INT.FRAC are used separately:
                        the part after the decimal point sets the fraction of
                        templates/pairs to be kept, while the integer part is
                        used as a seed that influences which subset of reads is
                        kept. When subsampling data that has previously been
                        subsampled, be sure to use a different seed value from
                        those used previously; otherwise more reads will be
                        retained than expected.
                    - 'sbg:altPrefix': '-@'
                      'sbg:category': Execution
                      'sbg:toolDefaultValue': '1'
                      id: threads
                      type: int?
                      inputBinding:
                        position: 19
                        prefix: '--threads'
                        shellQuote: false
                        valueFrom: |-
                          ${
                            if((inputs.threads)){
                              return (inputs.threads) - 1
                            }
                            else{
                              return
                            }
                          }
                      label: Number of threads
                      doc: >-
                        Number of threads. SAMtools uses argument --threads/-@
                        to specify number of additional threads. This parameter
                        sets total number of threads (and CPU cores). Command
                        line argument will be reduced by 1 to set number of
                        additional threads.
                    - 'sbg:category': Config Inputs
                      id: omitted_reads_filename
                      type: string?
                      inputBinding:
                        position: 20
                        prefix: '-U'
                        shellQuote: false
                      label: Filename for reads not selected by filters
                      doc: >-
                        Write alignments that are not selected by the various
                        filter options to this file. When this option is used,
                        all alignments (or all alignments intersecting the
                        regions specified) are written to either the output file
                        or this file, but never both.
                    - default: default_output_filename
                      'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': stdout
                      id: output_filename
                      type: string?
                      inputBinding:
                        position: 21
                        prefix: '-o'
                        shellQuote: false
                        valueFrom: |-
                          ${
                            if (inputs.output_filename!="default_output_filename"){
                              return (inputs.output_filename)
                            }
                            input_filename = [].concat(inputs.in_alignments)[0].path.split('/').pop()
                            input_name_base = input_filename.split('.').slice(0,-1).join('.')
                            ext = 'sam'
                            if (inputs.count_alignments){
                              return input_name_base + '.count.txt'
                            }
                            if ((inputs.uncompressed_bam) || (inputs.fast_bam_compression)){
                              ext = 'bam'
                            }
                            if (inputs.output_format){
                              ext = (inputs.output_format).toLowerCase()
                            }
                            if (inputs.output_header_only){
                              ext = 'header.' + ext
                            }
                            if (inputs.subsample_fraction){
                              ext = 'subsample.' + ext
                            }
                            if ((inputs.bed_file) || (inputs.read_group) || (inputs.read_group_list) ||
                                (inputs.filter_mapq) || (inputs.filter_library) || (inputs.min_cigar_operations) ||
                                (inputs.filter_include) || (inputs.filter_exclude_any) || 
                                (inputs.filter_exclude_all) || (inputs.regions_array)){
                              ext = 'filtered.' + ext
                            }
                              
                            return input_name_base + '.' + ext
                          }
                      label: Output filename
                      doc: Define a filename of the output.
                    - 'sbg:category': File Inputs
                      id: bed_file
                      type: File?
                      inputBinding:
                        position: 22
                        prefix: '-L'
                        shellQuote: false
                      label: BED region file
                      doc: Only output alignments overlapping the input BED file.
                      'sbg:fileTypes': BED
                    - 'sbg:category': File Inputs
                      id: read_group_list
                      type: File?
                      inputBinding:
                        position: 23
                        prefix: '-R'
                        shellQuote: false
                      label: Read group list
                      doc: Output alignments in read groups listed in this file.
                      'sbg:fileTypes': TXT
                    - 'sbg:altPrefix': '-T'
                      'sbg:category': File Inputs
                      id: in_reference
                      type: File?
                      inputBinding:
                        position: 24
                        prefix: '--reference'
                        shellQuote: false
                      label: Reference file
                      doc: >-
                        A FASTA format reference file, optionally compressed by
                        bgzip and ideally indexed by SAMtools Faidx. If an index
                        is not present, one will be generated for you. This file
                        is used for compression/decompression of CRAM files.
                        Please provide reference file when using CRAM
                        input/output file.
                      'sbg:fileTypes': 'FASTA, FA, FASTA.GZ, FA.GZ, GZ'
                    - 'sbg:category': File Inputs
                      id: reference_file_list
                      type: File?
                      inputBinding:
                        position: 25
                        prefix: '-t'
                        shellQuote: false
                      label: List of reference names and lengths
                      doc: >-
                        A tab-delimited file. Each line must contain the
                        reference name in the first column and the length of the
                        reference in the second column, with one line for each
                        distinct reference. Any additional fields beyond the
                        second column are ignored. This file also defines the
                        order of the reference sequences in sorting. If you run
                        SAMtools Faidx on reference FASTA file (<ref.fa>), the
                        resulting index file <ref.fa>.fai can be used as this
                        file.
                      'sbg:fileTypes': 'FAI, TSV, TXT'
                    - 'sbg:category': File Inputs
                      id: in_alignments
                      type: File
                      inputBinding:
                        position: 99
                        shellQuote: false
                      label: Input BAM/SAM/CRAM file
                      doc: Input BAM/SAM/CRAM file.
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - 'sbg:category': Config Inputs
                      id: regions_array
                      type: 'string[]?'
                      inputBinding:
                        position: 100
                        shellQuote: false
                      label: Regions array
                      doc: >-
                        With no options or regions specified, prints all
                        alignments in the specified input alignment file (in
                        SAM, BAM, or CRAM format) to output file in specified
                        format. Use of region specifications requires a
                        coordinate-sorted and indexed input file (in BAM or CRAM
                        format). Regions can be specified as:
                        RNAME[:STARTPOS[-ENDPOS]] and all position coordinates
                        are 1-based.  Important note: when multiple regions are
                        given, some alignments may be output multiple times if
                        they overlap more than one of the specified regions.
                        Examples of region specifications:  chr1 - Output all
                        alignments mapped to the reference sequence named `chr1'
                        (i.e. @SQ SN:chr1);  chr2:1000000 - The region on chr2
                        beginning at base position 1,000,000 and ending at the
                        end of the chromosome;  chr3:1000-2000 - The 1001bp
                        region on chr3 beginning at base position 1,000 and
                        ending at base position 2,000 (including both end
                        positions);  '*' - Output the unmapped reads at the end
                        of the file (this does not include any unmapped reads
                        placed on a reference sequence alongside their mapped
                        mates.);  . - Output all alignments (mostly unnecessary
                        as not specifying a region at all has the same effect).
                    - 'sbg:category': Config inputs
                      'sbg:toolDefaultValue': 'False'
                      id: multi_region_iterator
                      type: boolean?
                      inputBinding:
                        position: 22
                        prefix: '-M'
                        shellQuote: false
                      label: Use the multi-region iterator
                      doc: >-
                        Use the multi-region iterator on the union of the BED
                        file and command-line region arguments.
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1500'
                      id: mem_per_job
                      type: int?
                      label: Memory per job
                      doc: Memory per job in MB.
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1'
                      id: cpu_per_job
                      type: int?
                      label: CPU per job
                      doc: Number of CPUs per job.
                  outputs:
                    - id: out_alignments
                      doc: The output file.
                      label: 'Output BAM, SAM, or CRAM file'
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                            if ((inputs.output_filename!="default_output_filename")){
                              return (inputs.output_filename)
                            }
                            input_filename = [].concat((inputs.in_alignments))[0].path.split('/').pop()
                            input_name_base = input_filename.split('.').slice(0,-1). join('.')
                            ext = 'sam'
                            if ((inputs.count_alignments)){
                              return 
                            }
                            if ((inputs.uncompressed_bam) || (inputs.fast_bam_compression)){
                              ext = 'bam'
                            }
                            if ((inputs.output_format)){
                              ext = (inputs.output_format).toLowerCase()
                            }
                            if ((inputs.output_header_only)){
                              ext = 'header.' + ext
                            }
                            if ((inputs.subsample_fraction)){
                              ext = 'subsample.' + ext
                            }
                            if ((inputs.bed_file) || (inputs.read_group) || (inputs.read_group_list) ||
                                (inputs.filter_mapq) || (inputs.filter_library) || (inputs.min_cigar_operations) ||
                                (inputs.filter_include) || (inputs.filter_exclude_any) || 
                                (inputs.filter_exclude_all) || (inputs.regions_array)){
                              ext = 'filtered.' + ext
                            }
                              
                            return input_name_base + '.' + ext
                          }
                        outputEval: '$(inheritMetadata(self, inputs.in_alignments))'
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - id: reads_not_selected_by_filters
                      doc: File containing reads that are not selected by filters.
                      label: Reads not selected by filters
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                            if ((inputs.omitted_reads_filename)){
                              return (inputs.omitted_reads_filename)
                            }
                          }
                        outputEval: '$(inheritMetadata(self, inputs.in_alignments))'
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - id: alignement_count
                      doc: File containing number of alignments.
                      label: Alignment count
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                            input_filename = [].concat((inputs.in_alignments))[0].path.split('/').pop()
                            input_name_base = input_filename.split('.').slice(0,-1). join('.')
                            return input_name_base + '.count.txt'
                          }
                        outputEval: '$(inheritMetadata(self, inputs.in_alignments))'
                      'sbg:fileTypes': TXT
                  doc: >-
                    **SAMtools View** tool prints all alignments from a SAM,
                    BAM, or CRAM file to an output file in SAM format
                    (headerless). You may specify one or more space-separated
                    region specifications to restrict output to only those
                    alignments which overlap the specified region(s). Use of
                    region specifications requires a coordinate-sorted and
                    indexed input file (in BAM or CRAM format) [1].


                    *A list of **all inputs and parameters** with corresponding
                    descriptions can be found at the bottom of the page.*


                    ####Regions


                    Regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] and
                    all position coordinates are 1-based. 


                    **Important note:** when multiple regions are given, some
                    alignments may be output multiple times if they overlap more
                    than one of the specified regions.


                    Examples of region specifications:


                    - **chr1**  - Output all alignments mapped to the reference
                    sequence named `chr1' (i.e. @SQ SN:chr1).


                    - **chr2:1000000** - The region on chr2 beginning at base
                    position 1,000,000 and ending at the end of the chromosome.


                    - **chr3:1000-2000** - The 1001bp region on chr3 beginning
                    at base position 1,000 and ending at base position 2,000
                    (including both end positions).


                    - **'\*'** - Output the unmapped reads at the end of the
                    file. (This does not include any unmapped reads placed on a
                    reference sequence alongside their mapped mates.)


                    - **.** - Output all alignments. (Mostly unnecessary as not
                    specifying a region at all has the same effect.) [1]


                    ###Common Use Cases


                    This tool can be used for: 


                    - Filtering BAM/SAM/CRAM files - options set by the
                    following parameters and input files: **Include reads with
                    all of these flags** (`-f`), **Exclude reads with any of
                    these flags** (`-F`), **Exclude reads with all of these
                    flags** (`-G`), **Read group** (`-r`), **Minimum mapping
                    quality** (`-q`), **Only include alignments in library**
                    (`-l`), **Minimum number of CIGAR bases consuming query
                    sequence** (`-m`), **Subsample fraction** (`-s`), **Read
                    group list** (`-R`), **BED region file** (`-L`)

                    - Format conversion between SAM/BAM/CRAM formats - set by
                    the following parameters: **Output format**
                    (`--output-fmt/-O`), **Fast bam compression** (`-1`),
                    **Output uncompressed BAM** (`-u`)

                    - Modification of the data which is contained in each
                    alignment - set by the following parameters: **Collapse the
                    backward CIGAR operation** (`-B`), **Read tags to strip**
                    (`-x`)

                    - Counting number of alignments in SAM/BAM/CRAM file - set
                    by parameter **Output only count of matching records**
                    (`-c`)


                    ###Changes Introduced by Seven Bridges


                    - Parameters **Output BAM** (`-b`) and **Output CRAM**
                    (`-C`) were excluded from the wrapper since they are
                    redundant with parameter **Output format**
                    (`--output-fmt/-O`).

                    - Parameter **Input format** (`-S`) was excluded from
                    wrapper since it is ignored by the tool (input format is
                    auto-detected).

                    - Input file **Index file** was added to the wrapper to
                    enable operations that require an index file for BAM/CRAM
                    files.

                    - Parameter **Number of threads** (`--threads/-@`) specifies
                    the total number of threads instead of additional threads.
                    Command line argument (`--threads/-@`) will be reduced by 1
                    to set the number of additional threads.


                    ###Common Issues and Important Notes


                    - When multiple regions are given, some alignments may be
                    output multiple times if they overlap more than one of the
                    specified regions [1].

                    - Use of region specifications requires a coordinate-sorted
                    and indexed input file (in BAM or CRAM format) [1].

                    - Option **Output uncompressed BAM** (`-u`) saves time spent
                    on compression/decompression and is thus preferred when the
                    output is piped to another SAMtools command [1].


                    ###Performance Benchmarking


                    Multithreading can be enabled by setting parameter **Number
                    of threads** (`--threads/-@`). In the following table you
                    can find estimates of **SAMtools View** running time and
                    cost. 


                    *Cost can be significantly reduced by using **spot
                    instances**. Visit the [Knowledge
                    Center](https://docs.sevenbridges.com/docs/about-spot-instances)
                    for more details.*  


                    | Input type | Input size | # of reads | Read length |
                    Output format | # of threads | Duration | Cost | Instance
                    (AWS)|

                    |---------------|--------------|-----------------|---------------|------------------|-------------------|-----------------|-------------|--------|-------------|

                    | BAM | 5.26 GB | 71.5M | 76 | BAM | 1 | 13min. | \$0.12 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | 161.2M | 101 | BAM | 1 | 33min. | \$0.30
                    | c4.2xlarge |

                    | BAM | 18.36 GB | 179M | 76 | BAM | 1 | 60min. | \$0.54 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | 845.6M | 150 | BAM | 1 | 3h 25min. |
                    \$1.84 | c4.2xlarge |

                    | BAM | 5.26 GB | 71.5M | 76 | BAM | 8 | 5min. | \$0.04 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | 161.2M | 101 | BAM | 8 | 11min. | \$0.10
                    | c4.2xlarge |

                    | BAM | 18.36 GB | 179M | 76 | BAM | 8 | 19min. | \$0.17 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | 845.6M | 150 | BAM | 8 | 61min. | \$0.55
                    | c4.2xlarge |

                    | BAM | 5.26 GB | 71.5M | 76 | SAM | 8 | 14min. | \$0.13 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | 161.2M | 101 | SAM | 8 | 23min. | \$0.21
                    | c4.2xlarge |

                    | BAM | 18.36 GB | 179M | 76 | SAM | 8 | 35min. | \$0.31 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | 845.6M | 150 | SAM | 8 | 2h 29min. |
                    \$1.34 | c4.2xlarge |


                    ###References


                    [1] [SAMtools
                    documentation](http://www.htslib.org/doc/samtools-1.9.html)
                  label: Samtools View CWL1.0
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: |-
                        ${
                          if (inputs.mem_per_job) {
                              return inputs.mem_per_job
                          }    
                          else {
                          mem_offset = 1000
                          if((inputs.in_reference)){
                            mem_offset = mem_offset + 3000
                          }
                          if((inputs.threads)){
                            threads = (inputs.threads)
                          }
                          else{
                            threads = 1
                          }
                          return mem_offset + threads * 500
                          }
                        }
                      coresMin: |-
                        ${
                          if (inputs.cpu_per_job) {
                              return inputs.cpu_per_job
                          }
                          else {
                          if((inputs.threads)){
                            return (inputs.threads)
                          }
                          else{
                            return 1
                          }
                          }
                        }
                    - class: DockerRequirement
                      dockerPull: 'images.sbgenomics.com/jrandjelovic/samtools-1-9:1'
                    - class: InitialWorkDirRequirement
                      listing:
                        - $(inputs.in_reference)
                        - $(inputs.reference_file_list)
                        - $(inputs.in_index)
                        - $(inputs.in_alignments)
                    - class: InlineJavascriptRequirement
                      expressionLib:
                        - |-

                          var setMetadata = function(file, metadata) {
                              if (!('metadata' in file))
                                  file['metadata'] = metadata;
                              else {
                                  for (var key in metadata) {
                                      file['metadata'][key] = metadata[key];
                                  }
                              }
                              return file
                          };

                          var inheritMetadata = function(o1, o2) {
                              var commonMetadata = {};
                              if (!Array.isArray(o2)) {
                                  o2 = [o2]
                              }
                              for (var i = 0; i < o2.length; i++) {
                                  var example = o2[i]['metadata'];
                                  for (var key in example) {
                                      if (i == 0)
                                          commonMetadata[key] = example[key];
                                      else {
                                          if (!(commonMetadata[key] == example[key])) {
                                              delete commonMetadata[key]
                                          }
                                      }
                                  }
                              }
                              if (!Array.isArray(o1)) {
                                  o1 = setMetadata(o1, commonMetadata)
                              } else {
                                  for (var i = 0; i < o1.length; i++) {
                                      o1[i] = setMetadata(o1[i], commonMetadata)
                                  }
                              }
                              return o1;
                          };
                  'sbg:appVersion':
                    - v1.0
                  'sbg:categories':
                    - Utilities
                    - BAM Processing
                    - CWL1.0
                  'sbg:content_hash': >-
                    aa82916613444b2d378befd3fc8666677b6b22c3fb84f9dd8985aa73494c63afa
                  'sbg:contributors':
                    - admin
                  'sbg:createdBy': admin
                  'sbg:createdOn': 1576244128
                  'sbg:id': admin/sbg-public-data/samtools-view-1-9-cwl1-0/6
                  'sbg:image_url': null
                  'sbg:latestRevision': 6
                  'sbg:license': MIT License
                  'sbg:links':
                    - id: 'http://www.htslib.org/'
                      label: Homepage
                    - id: 'https://github.com/samtools/samtools'
                      label: Source Code
                    - id: 'https://github.com/samtools/samtools/wiki'
                      label: Wiki
                    - id: >-
                        https://sourceforge.net/projects/samtools/files/samtools/
                      label: Download
                    - id: 'http://www.ncbi.nlm.nih.gov/pubmed/19505943'
                      label: Publication
                    - id: 'http://www.htslib.org/doc/samtools-1.9.html'
                      label: Documentation
                  'sbg:modifiedBy': admin
                  'sbg:modifiedOn': 1578576084
                  'sbg:project': admin/sbg-public-data
                  'sbg:projectName': SBG Public Data
                  'sbg:publisher': sbg
                  'sbg:revision': 6
                  'sbg:revisionNotes': Added file requirements for in_index and in_alignments
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244128
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 1
                      'sbg:revisionNotes': Final version
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 2
                      'sbg:revisionNotes': 'Edited description, tag, default values.'
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 3
                      'sbg:revisionNotes': mem_per_job default value set
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 4
                      'sbg:revisionNotes': Description edited - references put before full stop
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 5
                      'sbg:revisionNotes': Categories edited
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1578576084
                      'sbg:revision': 6
                      'sbg:revisionNotes': Added file requirements for in_index and in_alignments
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': >-
                    Heng Li (Sanger Institute), Bob Handsaker (Broad Institute),
                    Jue Ruan (Beijing Genome Institute), Colin Hercus, Petr
                    Danecek
                  'sbg:toolkit': samtools
                  'sbg:toolkitVersion': '1.9'
                  'sbg:validationErrors': []
                label: Samtools View CWL1.0
                'sbg:x': 547.8245239257812
                'sbg:y': 120.75
              - id: samtools_fastq
                in:
                  - id: in_alignments
                    source: samtools_view_1_9_cwl1_0/out_alignments
                out:
                  - id: output_pe_1
                  - id: output_pe_2
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: admin_sbg_public_data_samtools_fastq_1_9_cwl1_0_4
                  baseCommand: []
                  inputs:
                    - 'sbg:category': File Inputs
                      id: in_alignments
                      type: File
                      inputBinding:
                        position: 101
                        shellQuote: false
                      label: BAM/SAM/CRAM file
                      doc: Input BAM/SAM/CRAM file.
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - default: 0
                      'sbg:category': Config Inputs
                      id: read_0_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-0'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.read_0_filename = null
                              };


                              if (inputs.single_fastq_file == true) {
                                  return
                              } else {
                                  if (inputs.read_0_filename) {
                                      return inputs.read_0_filename
                                  } else {
                                      var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                      var ext = bamname.split('.').pop().toLowerCase()
                                      if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                          return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.pe_0.fastq'
                                      } else {
                                          return bamname + '.pe_0.fastq'
                                      }
                                  }
                              }
                          }
                      label: Unflagged reads filename
                      doc: >-
                        Write reads with both or neither of the BAM_READ1 and
                        BAM_READ2 flags set to this file.
                    - default: 0
                      'sbg:category': Config Inputs
                      id: read_1_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-1'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.read_1_filename = null
                              };


                              if (inputs.single_fastq_file == true) {
                                  return
                              } else {
                                  if (inputs.read_1_filename) {
                                      return inputs.read_1_filename
                                  } else {
                                      var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                      var ext = bamname.split('.').pop().toLowerCase()
                                      if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                          return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.pe_1.fastq'
                                      } else {
                                          return bamname + '.pe_1.fastq'
                                      }
                                  }
                              }
                          }
                      label: Filename for BAM_READ1 flaged reads
                      doc: Write reads with the BAM_READ1 flag set to this file.
                    - default: 0
                      'sbg:category': Config Inputs
                      id: read_2_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-2'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.read_2_filename = null
                              };


                              if (inputs.single_fastq_file == true) {
                                  return
                              } else {
                                  if (inputs.read_2_filename) {
                                      return inputs.read_2_filename
                                  } else {
                                      var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                      var ext = bamname.split('.').pop().toLowerCase()
                                      if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                          return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.pe_2.fastq'
                                      } else {
                                          return bamname + '.pe_2.fastq'
                                      }
                                  }
                              }
                          }
                      label: Filename for BAM_READ2 flaged reads
                      doc: Write reads with the BAM_READ2 flag set to this file.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_include
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-f'
                        shellQuote: false
                      label: Include reads with all of these flags
                      doc: >-
                        Only output alignments with all bits set in this integer
                        present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_any
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-F'
                        shellQuote: false
                      label: Exclude reads with any of these flags
                      doc: >-
                        Do not output alignments with any bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_all
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-G'
                        shellQuote: false
                      label: Exclude reads with all of these flags
                      doc: >-
                        Only exclude reads with all of the bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: do_not_append_read_number_to_name
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-n'
                        shellQuote: false
                      label: Don't append /1 and /2 to the read name
                      doc: >-
                        By default, either '/1' or '/2' is added to the end of
                        read names where the corresponding BAM_READ1 or
                        BAM_READ2 flag is set. Setting this parameter to True
                        causes read names to be left as they are.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: always_append_read_number_to_name
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-N'
                        shellQuote: false
                      label: Always append /1 and /2 to the read name
                      doc: >-
                        Always add either '/1' or '/2' to the end of read names
                        even when put into different files.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: output_quality_in_OQ_tag
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-O'
                        shellQuote: false
                      label: Output quality in the OQ tag if present
                      doc: >-
                        Use quality values from OQ tags in preference to
                        standard quality string if available.
                    - 'sbg:category': Config Inputs
                      id: singleton_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-s'
                        shellQuote: false
                      label: Singleton reads filename
                      doc: Write singleton reads in FASTQ format to this file.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: copy_tags
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-t'
                        shellQuote: false
                      label: 'Copy RG, BC and QT tags to the FASTQ header line'
                      doc: >-
                        Copy RG, BC and QT tags to the FASTQ header line, if
                        they exist.
                    - default: 0
                      'sbg:category': Config Inputs
                      id: copy_taglist
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-T'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.copy_taglist = null
                              };


                              if (inputs.copy_taglist) {
                                  return inputs.copy_taglist.replace(/ /g, '')
                              } else {
                                  return
                              }
                          }
                      label: Taglist to copy to the FASTQ header line
                      doc: >-
                        Specify a comma-separated list of tags to copy to the
                        FASTQ header line, if they exist.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '1'
                      id: default_quality_score
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-v'
                        shellQuote: false
                      label: Default quality score if not given in file
                      doc: Default quality score if not given in file.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: add_Illumina_casava_format_entry
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-i'
                        shellQuote: false
                      label: Add Illumina Casava 1.8 format entry to header
                      doc: >-
                        Add Illumina Casava 1.8 format entry to header (eg
                        1:N:0:ATCACG).
                    - 'sbg:category': Config Inputs
                      id: compression_level
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-c'
                        shellQuote: false
                      label: 'Compression level [0..9]'
                      doc: >-
                        Compression level [0..9] to use when creating GZ or BGZF
                        FASTQ files.
                    - 'sbg:category': Config Inputs
                      id: first_index_read_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--i1'
                        shellQuote: false
                      label: First index reads filename
                      doc: >-
                        Specify filename to which the first index reads will be
                        written.
                    - 'sbg:category': Config Inputs
                      id: second_index_read_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--i2'
                        shellQuote: false
                      label: Second index reads filename
                      doc: >-
                        Specify filename to which the second index reads will be
                        written.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': BC
                      id: barcode_tag
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--barcode-tag'
                        shellQuote: false
                      label: Barcode tag
                      doc: 'Aux tag to find index reads in [default: BC].'
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': QT
                      id: quality_tag
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--quality-tag'
                        shellQuote: false
                      label: Quality tag
                      doc: 'Aux tag to find index quality in [default: QT].'
                    - 'sbg:category': Config Inputs
                      id: index_format
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--index-format'
                        shellQuote: false
                      label: Index format for parsing barcode and quality tags
                      doc: >-
                        String to describe how to parse the barcode and quality
                        tags. For example:  i14i8 - the first 14 characters are
                        index 1, the next 8 characters are index 2; n8i14 ignore
                        the first 8 characters, and use the next 14 characters
                        for index 1; If the tag contains a separator, then the
                        numeric part can be replaced with '*' to mean 'read
                        until the separator or end of tag', for example:  n*i* -
                        ignore the left part of the tag until the separator,
                        then use the second part.
                    - 'sbg:category': Config Inputs
                      id: single_fastq_filename
                      type: string?
                      label: Single FASTQ filename
                      doc: >-
                        Filename of an output FASTQ file if only one file is
                        required.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: single_fastq_file
                      type: boolean?
                      label: Single FASTQ file
                      doc: True if only one FASTQ file is required as output.
                    - 'sbg:category': Config Inputs
                      id: input_fmt_option
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--input-fmt-option'
                        shellQuote: false
                      label: Input file format option
                      doc: >-
                        Specify a single input file format option in the form of
                        OPTION or OPTION=VALUE.
                    - default: 0
                      'sbg:altPrefix': '-@'
                      'sbg:category': Execution
                      'sbg:toolDefaultValue': '1'
                      id: threads
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '--threads'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.threads = null
                              };


                              if (inputs.threads) {
                                  return inputs.threads - 1
                              } else
                                  return
                          }
                      label: Number of threads
                      doc: >-
                        Number of threads. SAMtools uses argument --threads/-@
                        to specify number of additional threads. This parameter
                        sets total number of threads (and CPU cores). Command
                        line argument will be reduced by 1 to set number of
                        additional threads.
                    - 'sbg:category': File Inputs
                      id: in_reference
                      type: File?
                      inputBinding:
                        position: 2
                        prefix: '--reference'
                        shellQuote: false
                      label: Reference file
                      doc: >-
                        Reference file. This file is used for
                        compression/decompression of CRAM files. Please provide
                        reference file when using CRAM input/output file.
                      'sbg:fileTypes': 'FASTA, FASTA.GZ, FASTA.BGZF, GZ, BGZF'
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1500'
                      id: mem_per_job
                      type: int?
                      label: Memory per job
                      doc: Memory per job in MB.
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1'
                      id: cpu_per_job
                      type: int?
                      label: CPU per job
                      doc: Number of CPUs per job.
                  outputs:
                    - id: output_pe_1
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                              if (inputs.read_1_filename) {
                                  var read1 = inputs.read_1_filename
                              } else {
                                  read1 = '*pe_1.fastq*'
                              }
                              return read1
                          }
                    - id: output_pe_2
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                              if (inputs.read_2_filename) {
                                  var read2 = inputs.read_2_filename
                              } else {
                                  read2 = '*pe_2.fastq*'
                              }
                              return read2
                          }
                  doc: >-
                    **SAMtools FASTQ** tool converts a BAM or CRAM into FASTQ
                    format. The FASTQ files will be automatically compressed if
                    the filenames have a .gz or .bgzf extension [1].


                    **SAMtools FASTQ does not perform any preprocessing of the
                    input.** If you want to use output FASTQ files with
                    alignment tools, please make sure that the **BAM/SAM/CRAM
                    file** is sorted by read name (or collated). Otherwise,
                    alignment tools will fail. Supplementary and secondary
                    alignments are filtered out regardless of filtering
                    parameters. 


                    Parameter **Index format for parsing barcode and quality
                    tags** (`--index-format`) is a string used to describe how
                    to parse the barcode and quality tags. For example:  


                    - i14i8 - the first 14 characters are index 1, the next 8
                    characters are index 2

                    - n8i14 - ignore the first 8 characters, and use the next 14
                    characters for index 1  

                    - If the tag contains a separator, then the numeric part can
                    be replaced with '\*' to mean 'read until the separator or
                    end of tag', for example:  

                    n\*i\* - ignore the left part of the tag until the
                    separator, then use the second part [1].


                    *A list of **all inputs and parameters** with corresponding
                    descriptions can be found at the bottom of the page.*


                    ###Common Use Cases


                    - Default use case provides two FASTQ files as outputs
                    (**Paired-end FASTQ files**). If no parameter is set, this
                    will be the case. The tool will output file with reads that
                    are not properly flagged (**Unflagged reads FASTQ file**)
                    only in case this file is not empty. **SAMtools FASTQ does
                    not perform any preprocessing of the input.** If you want to
                    use output FASTQ files with alignment tools, please make
                    sure that the **BAM/SAM/CRAM file** is sorted by read name
                    (or collated). Otherwise, alignment tools will fail.
                    Supplementary and secondary alignments are filtered out
                    regardless of filtering parameters.   

                    - If a single FASTQ file (with both paired end reads) is
                    required, it should be specified by setting the boolean
                    parameter **Single FASTQ file** to True. 


                    ###Changes Introduced by Seven Bridges


                    - Parameter **Single FASTQ file** was added to parameter
                    list to provide the option for outputting a single FASTQ
                    file with all the reads.

                    - Parameter **Single FASTQ filename** was added to parameter
                    list to specify the filename when **Single FASTQ file** is
                    set to True. This parameter is not mandatory. If **Single
                    FASTQ file** is set to True and **Single FASTQ filename** is
                    not specified, default value will be used (*input.fastq* for
                    input **BAM/SAM/CRAM file** named *input.bam*).

                    - Parameter **Number of threads** (`--threads/-@`) specifies
                    the total number of threads instead of additional threads.
                    Command line argument (`--threads/-@`) will be reduced by 1
                    to set number of additional threads.


                    ###Common Issues and Important Notes


                    - If parameters **First index reads filename** (`--i1`)
                    and/or **Second index reads filename** (`--i2`) are
                    specified, **Index format for parsing barcode and quality
                    tags** (`--index-format`) should be specified too and this
                    format should match the number of index files required. 

                    - When specifying output filenames, complete names should be
                    used (including extensions). If the extension is .fastq.gz
                    or .fastq.bgzf, the output will be compressed. The tool does
                    not validate extensions. If the extension is not valid, the
                    task will not fail.

                    - Parameter **Number of threads** (`--threads/-@`) does not
                    decrease running time significantly (not more than 10% with
                    8 threads on a c4.2xlarge instance (AWS)).

                    - **SAMtools FASTQ** does not perform any preprocessing of
                    the input. If you want to use output FASTQ files with
                    alignment tools, please make sure that **BAM/SAM/CRAM file**
                    is sorted by read name (or collated). Otherwise, alignment
                    tools will fail. Supplementary and secondary alignments are
                    filtered out regardless of filtering parameters. 


                    ###Performance Benchmarking


                    In the following table you can find estimates of **SAMtools
                    FASTQ** running time and cost. Adding additional threads
                    does not decrease running time significantly (not more than
                    10% with 8 threads on a c4.2xlarge instance (AWS)).


                    *Cost can be significantly reduced by using **spot
                    instances**. Visit the [Knowledge
                    Center](https://docs.sevenbridges.com/docs/about-spot-instances)
                    for more details.*  


                    | Input type | Input size | Paired-end | # of reads | Read
                    length |  # of threads | Duration | Cost | Instance (AWS) |

                    |---------------|--------------|-----------------|---------------|------------------|-----------------|-------------|--------|-------------|

                    | BAM | 5.26 GB | Yes | 71.5M | 76 | 1 | 7min. | \$0.06 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | Yes | 161.2M | 101 | 1 | 18min. | \$0.16
                    | c4.2xlarge |

                    | BAM | 18.36 GB | Yes | 179M | 76 | 1 | 21min. | \$0.19 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | Yes | 845.6M | 150 | 1 | 2h 7min. |
                    \$1.14 | c4.2xlarge |


                    ###References


                    [1] [SAMtools
                    documentation](http://www.htslib.org/doc/samtools-1.9.html)
                  label: Samtools FASTQ CWL1.0
                  arguments:
                    - position: 0
                      prefix: ''
                      shellQuote: false
                      valueFrom: /opt/samtools-1.9/samtools
                    - position: 1
                      shellQuote: false
                      valueFrom: fastq
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: |-
                        ${
                            if (inputs.mem_per_job) {
                                return inputs.mem_per_job
                            }
                            else {
                            if (inputs.threads) {
                                var threads = inputs.threads
                            } else {
                                threads = 1
                            }
                            return 1000 + 500 * threads
                            }
                        }
                      coresMin: |-
                        ${
                            if (inputs.cpu_per_job) {
                                return inputs.cpu_per_job
                            }
                            else {
                            if (inputs.threads) {
                                return inputs.threads
                            } else {
                                return 1
                            }
                            }
                        }
                    - class: DockerRequirement
                      dockerPull: 'images.sbgenomics.com/jrandjelovic/samtools-1-9:1'
                    - class: InitialWorkDirRequirement
                      listing:
                        - $(inputs.in_reference)
                    - class: InlineJavascriptRequirement
                      expressionLib:
                        - |
                          var setMetadata = function(file, metadata) {
                              if (!('metadata' in file)) {
                                  file['metadata'] = {}
                              }
                              for (var key in metadata) {
                                  file['metadata'][key] = metadata[key];
                              }
                              return file
                          };

                          var inheritMetadata = function(o1, o2) {
                              var commonMetadata = {};
                              if (!o2) {
                                  return o1;
                              };
                              if (!Array.isArray(o2)) {
                                  o2 = [o2]
                              }
                              for (var i = 0; i < o2.length; i++) {
                                  var example = o2[i]['metadata'];
                                  for (var key in example) {
                                      if (i == 0)
                                          commonMetadata[key] = example[key];
                                      else {
                                          if (!(commonMetadata[key] == example[key])) {
                                              delete commonMetadata[key]
                                          }
                                      }
                                  }
                                  for (var key in commonMetadata) {
                                      if (!(key in example)) {
                                          delete commonMetadata[key]
                                      }
                                  }
                              }
                              if (!Array.isArray(o1)) {
                                  o1 = setMetadata(o1, commonMetadata)
                              } else {
                                  for (var i = 0; i < o1.length; i++) {
                                      o1[i] = setMetadata(o1[i], commonMetadata)
                                  }
                              }
                              return o1;
                          };
                  stdout: |-
                    ${
                        if (inputs.single_fastq_file == true) {
                            if (inputs.single_fastq_filename) {
                                return inputs.single_fastq_filename
                            } else {
                                var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                var ext = bamname.split('.').pop().toLowerCase()
                                if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                    return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.fastq'
                                } else {
                                    return bamname + '.fastq'
                                }
                            }
                        } else {
                            return
                        }
                    }
                  'sbg:appVersion':
                    - v1.0
                  'sbg:categories':
                    - Utilities
                    - BAM Processing
                    - CWL1.0
                  'sbg:cmdPreview': /opt/samtools-1.9/samtools fastq  /path/to/input.bam
                  'sbg:content_hash': >-
                    a7f356d044f84a927ee18a8c03a269f0d6cad2ff36f6c6f10c57d4ed5df3e5e04
                  'sbg:contributors':
                    - admin
                    - Rachel Bowen-James <rbowen-james@ccia.org.au>
                  'sbg:createdBy': admin
                  'sbg:image_url': null
                  'sbg:latestRevision': 4
                  'sbg:license': MIT License
                  'sbg:links':
                    - id: 'http://www.htslib.org/'
                      label: Homepage
                    - id: 'https://github.com/samtools/samtools'
                      label: Source Code
                    - id: 'https://github.com/samtools/samtools/wiki'
                      label: Wiki
                    - id: >-
                        https://sourceforge.net/projects/samtools/files/samtools/
                      label: Download
                    - id: 'http://www.ncbi.nlm.nih.gov/pubmed/19505943'
                      label: Publication
                    - id: 'http://www.htslib.org/doc/samtools-1.9.html'
                      label: Documentation
                  'sbg:modifiedBy': admin
                  'sbg:modifiedOn': 1576244287
                  'sbg:project': admin/sbg-public-data
                  'sbg:projectName': SBG Public Data
                  'sbg:publisher': sbg
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': >-
                    Heng Li (Sanger Institute), Bob Handsaker (Broad Institute),
                    Jue Ruan (Beijing Genome Institute), Colin Hercus, Petr
                    Danecek
                  'sbg:toolkit': SAMtools
                  'sbg:toolkitVersion': '1.9'
                  'sbg:validationErrors': []
                  'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
                label: Samtools FASTQ CWL1.0
                'sbg:x': 977.3401489257812
                'sbg:y': 113.75
              - id: hla_hd
                in:
                  - id: threads
                    default: 2
                  - id: minimum_read_length
                    default: 0
                    source: minimum_read_length
                  - id: fastq_reads1
                    source: samtools_fastq/output_pe_1
                  - id: fastq_reads2
                    source: samtools_fastq/output_pe_2
                  - id: sample_id
                    source: sample_id
                  - id: output_dir_name
                    source: output_dir_name
                out:
                  - id: hlahd_results
                  - id: hlahd_final_results
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/hla-hd/30
                  baseCommand: []
                  inputs:
                    - id: threads
                      type: int
                      inputBinding:
                        position: 1
                        prefix: '-t'
                        shellQuote: false
                      label: Threads
                      doc: Number of cores used to execute the program.
                    - 'sbg:toolDefaultValue': '100'
                      id: minimum_read_length
                      type: int
                      inputBinding:
                        position: 1
                        prefix: '-m'
                        shellQuote: false
                      label: Minimum read length
                      doc: >-
                        A read whose length is shorter than this parameter is
                        ignored.
                    - id: fastq_reads1
                      type: File
                      inputBinding:
                        position: 2
                        separate: false
                        shellQuote: false
                      label: FASTQ Reads 1
                      doc: Paired-end reads 1 in FASTQ format.
                      'sbg:fileTypes': FASTQ
                    - id: fastq_reads2
                      type: File
                      inputBinding:
                        position: 2
                        separate: false
                        shellQuote: false
                      label: FASTQ Reads 2
                      doc: Paired-end reads 2 in FASTQ format.
                      'sbg:fileTypes': FASTQ
                    - id: sample_id
                      type: string
                      inputBinding:
                        position: 4
                        separate: false
                        shellQuote: false
                      label: Sample ID
                      doc: >-
                        Sample ID for the input FASTQs. This will be used as the
                        name of the output directory.
                    - id: output_dir_name
                      type: string?
                  outputs:
                    - id: hlahd_results
                      doc: Directory containing results of the HLA-HD run.
                      label: Output directory
                      type: Directory
                      outputBinding:
                        glob: |-
                          ${
                              if (!inputs.output_dir_name) {
                                  return inputs.sample_id + "_hlahd"
                              } else {
                                  return inputs.output_dir_name
                              }
                          }
                    - id: hlahd_final_results
                      type: File
                      outputBinding:
                        glob: |-
                          ${
                              if (!inputs.output_dir_name) {
                                  return inputs.sample_id + "_hlahd/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                              } else {
                                  return inputs.output_dir_name + "/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                              }
                          }
                  doc: >-
                    ## About HLA-HD

                    HLA-HD documentation and release notes can be found
                    [here](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/).

                    HLA-HD (HLA typing from High-quality Dictionary) can
                    accurately determine HLA alleles with 6-digit precision from
                    NGS data (FASTQ format). RNA-Seq data can also be applied.


                    ## About this CWL tool

                    This CWL tool runs HLA-HD version 1.4.0 on paired-end FASTQ
                    files to determine HLA type.


                    ### Inputs and parameters

                    - The input paired-end read files can be from **WGS/WES or
                    RNA-seq**.

                    - The input paired-end read files must be in FASTQ format
                    (**not zipped**).

                    - The default minimum read length is 100, however this is
                    often too strict. Choose a lower threshold to include more
                    reads.


                    ### Output

                    - HLA-HD results are output to a directory named using the
                    input sample id.

                    - The final summary of HLA typing results can be found at
                    the following path:
                    `<output_dir_name>/result/<sample_id>_final.result.txt`.


                    ### Other notes

                    - This tool uses the HLA dictionary created from release
                    3.15.0 of the
                    [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA) database.

                    - This tool by default uses HLA allele frequency data
                    included with the HLA-HD release 1.4.0.
                  label: HLA-HD
                  arguments:
                    - position: 2
                      prefix: '-f'
                      shellQuote: false
                      valueFrom: /app/hlahd.1.4.0/freq_data
                    - position: 3
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: /app/hlahd.1.4.0/HLA_gene.split.txt
                    - position: 3
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: /app/hlahd.1.4.0/dictionary
                    - position: 1
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: hlahd.sh
                    - position: 5
                      prefix: ''
                      valueFrom: ./
                    - position: 6
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: |-
                        ${
                            if (!inputs.output_dir_name) {
                                return "&& mv ./" + inputs.sample_id + " ./" + inputs.sample_id + "_hlahd"
                            } else {
                                return "&& mv ./" + inputs.sample_id + " ./" + inputs.output_dir_name
                            }
                        }
                    - position: 0
                      prefix: ''
                      shellQuote: false
                      valueFrom: |-
                        ${
                            if (!inputs.output_dir_name) {
                                return "mkdir ./" + inputs.sample_id + "_hlahd &&"
                            } else {
                                return "mkdir ./" + inputs.output_dir_name + " &&"
                            }
                        }
                  requirements:
                    - class: ShellCommandRequirement
                    - class: DockerRequirement
                      dockerPull: pgc-images.sbgenomics.com/rbowen_james/hla-hd
                    - class: InlineJavascriptRequirement
                  'sbg:appVersion':
                    - v1.0
                  'sbg:categories':
                    - WGS
                    - RNA
                    - HLA Typing
                    - HLA
                    - MHC
                    - WES (WXS)
                  'sbg:content_hash': >-
                    a4adfcaaafaaa6b9ac47bbe7e09595dae02a8f11b9f16fdaa0cfcb21948bf91bd
                  'sbg:contributors':
                    - rbowen_james
                  'sbg:createdBy': rbowen_james
                  'sbg:createdOn': 1622961333
                  'sbg:id': mwonge/mwtest/hla-hd/30
                  'sbg:image_url': null
                  'sbg:latestRevision': 30
                  'sbg:modifiedBy': rbowen_james
                  'sbg:modifiedOn': 1627370832
                  'sbg:project': mwonge/mwtest
                  'sbg:projectName': zcc-cavatica
                  'sbg:publisher': sbg
                  'sbg:revision': 30
                  'sbg:revisionNotes': output final result file
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622961333
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622961864
                      'sbg:revision': 1
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622962039
                      'sbg:revision': 2
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622962892
                      'sbg:revision': 3
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623026790
                      'sbg:revision': 4
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623027275
                      'sbg:revision': 5
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623043041
                      'sbg:revision': 6
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623043454
                      'sbg:revision': 7
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623044879
                      'sbg:revision': 8
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623046596
                      'sbg:revision': 9
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623048682
                      'sbg:revision': 10
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624408084
                      'sbg:revision': 11
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624411148
                      'sbg:revision': 12
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624580200
                      'sbg:revision': 13
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624581085
                      'sbg:revision': 14
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624581597
                      'sbg:revision': 15
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624583201
                      'sbg:revision': 16
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624583825
                      'sbg:revision': 17
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624586750
                      'sbg:revision': 18
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624588045
                      'sbg:revision': 19
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624593624
                      'sbg:revision': 20
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624606248
                      'sbg:revision': 21
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624606452
                      'sbg:revision': 22
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624608065
                      'sbg:revision': 23
                      'sbg:revisionNotes': 'Fixed output dir issue, added docs.'
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624779790
                      'sbg:revision': 24
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627028760
                      'sbg:revision': 25
                      'sbg:revisionNotes': Added output dir option
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627028962
                      'sbg:revision': 26
                      'sbg:revisionNotes': make outdir
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627029487
                      'sbg:revision': 27
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627029828
                      'sbg:revision': 28
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627104094
                      'sbg:revision': 29
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627370832
                      'sbg:revision': 30
                      'sbg:revisionNotes': output final result file
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': Shuji Kawaguchi <shuji@genome.med.kyoto-u.ac.jp>
                  'sbg:toolkit': HLA-HD
                  'sbg:toolkitVersion': 1.4.0
                  'sbg:validationErrors': []
                  'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
                label: HLA-HD
                'sbg:x': 1256.66357421875
                'sbg:y': 280.875
              - id: bowtie2_aligner
                in:
                  - id: bowtie_index_archive
                    source: bowtie_index_archive
                  - id: paired_aligned_reads
                    default: raw
                  - id: read_sequence
                    source:
                      - read_sequence
                  - id: suppress_sam_records
                    default: true
                out:
                  - id: result_sam_file
                  - id: aligned_reads_only
                  - id: unaligned_reads_only
                run:
                  cwlVersion: 'sbg:draft-2'
                  class: CommandLineTool
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: admin/sbg-public-data/bowtie2-aligner/18
                  label: Bowtie2 Aligner
                  description: >-
                    Bowtie 2 is an ultrafast and memory-efficient tool for
                    aligning sequencing reads to long reference sequences. It is
                    particularly good at aligning reads of about 50 up to 100s
                    or 1,000s of characters to relatively long (e.g. mammalian)
                    genomes. Bowtie 2 indexes the genome with an [FM
                    Index](http://portal.acm.org/citation.cfm?id=796543) (based
                    on the [Burrows-Wheeler
                    Transform](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform)
                    or
                    [BWT](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform))
                    to keep its memory footprint small: for the human genome,
                    its memory footprint is typically around 3.2 gigabytes of
                    RAM. In order to create needed index files, you should run
                    [Bowtie2
                    Indexer](https://igor.sbgenomics.com/public/apps#tool/admin/sbg-public-data/bowtie2-indexer),
                    which produces archived index files (containing 6 files with
                    suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and
                    .rev.2.bt2).


                    Bowtie 2 supports gapped, local, and paired-end alignment
                    modes. Bowtie 2 outputs alignments in SAM format, enabling
                    interoperation with a large number of other tools (e.g.
                    [SAMtools](http://samtools.sourceforge.net/),
                    [GATK](http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit))
                    that use SAM.


                    ###Common issues###

                    No issues have been reported.


                    **Q&A:**


                    ***Q: What should I do if I already have Bowtie2 index
                    files, not archived as tar bundle?***


                    ***A***: You can provide your *.bt2 files to [SBG
                    Compressor](https://igor.sbgenomics.com/public/apps#admin/sbg-public-data/sbg-compressor-1-0/)
                    app from our public apps and set "TAR" as your output
                    format. After the task is finished, **you should assign
                    common prefix of the index files to the `Reference genome`
                    metadata field** and your TAR is ready for use.


                    ***Example:***

                    Indexed files: chr20.1.bt2, chr20.2.bt2, chr20.3.bt2,
                    chr20.4.bt2, chr20.rev.1.bt2, chr20.rev.2.bt2


                    Metadata - `Reference genome`: **chr20**


                    __Important note: In case of paired-end alignment it is
                    crucial to set metadata 'paired-end' field to 1/2. Sequences
                    specified as mate 1s must correspond file-for-file and
                    read-for-read with those specified for mate 2s. Reads may be
                    a mix of different lengths. In case of unpaired reads, the
                    same metadata field should be set to '-'. Only one type of
                    alignment can be performed at once, so all specified reads
                    should be either paired or unpaired.__
                  baseCommand:
                    - class: Expression
                      engine: '#cwl-js-engine'
                      script: |-
                        {
                          var archive_name = $job.inputs.bowtie_index_archive.path.split("/").pop()
                          return "tar -xvf ".concat(archive_name, " && rm -rf ", archive_name, " && ")
                        }
                    - /opt/bowtie2-2.2.6/bowtie2
                  inputs:
                    - 'sbg:category': Input files
                      'sbg:stageInput': link
                      type:
                        - File
                      label: Bowtie index archive
                      description: Archive file produced by Bowtie2 Indexer.
                      'sbg:fileTypes': TAR
                      id: '#bowtie_index_archive'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': End-to-end
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - End-to-end
                            - Local
                          name: alignment_mode
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              if ($job.inputs.alignment_mode == "Local") {
                                return "--local"
                              }
                            }
                        'sbg:cmdInclude': true
                      label: Alignment mode
                      description: >-
                        Alignment mode. End-to-end: entire read must align; no
                        clipping. Local: local alignment; ends might be soft
                        clipped.
                      id: '#alignment_mode'
                    - 'sbg:altPrefix': '-s'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': '-'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--skip'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Skip reads
                      description: >-
                        Skip (i.e. do not align) the first given number of reads
                        or pairs in the input.
                      id: '#skip_reads'
                    - 'sbg:altPrefix': '-u'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': No limit
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--upto'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Align next n reads
                      description: >-
                        Align the first given number of reads or read pairs from
                        the input (after the <int> reads or pairs have been
                        skipped with "Skip reads"), then stop.
                      id: '#align_next_n_reads'
                    - 'sbg:altPrefix': '-5'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--trim5'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Trim from 5'
                      description: >-
                        Trim given number of bases from 5' (left) end of each
                        read before alignment.
                      id: '#trim_from_5'
                    - 'sbg:altPrefix': '-3'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--trim3'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Trim from 3'
                      description: >-
                        Trim given number of bases from 3' (right) end of each
                        read before alignment.
                      id: '#trim_from_3'
                    - 'sbg:category': Input
                      'sbg:toolDefaultValue': Phred+33
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Auto-detect
                            - Phred+33
                            - Phred+64
                            - Solexa
                          name: quality_scale
                      label: Quality scale
                      description: Set quality scale.
                      id: '#quality_scale'
                    - 'sbg:category': Input
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--int-quals'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Integer qualities
                      description: >-
                        Quality values are represented in the read input file as
                        space-separated ASCII integers, e.g., 40 40 30 40...,
                        rather than ASCII characters, e.g., II?I....
                      id: '#integer_qualities'
                    - 'sbg:category': Alignment
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - '0'
                            - '1'
                          name: allowed_mismatch_number
                      inputBinding:
                        position: 0
                        prefix: '-N'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Allowed mismatch number
                      description: >-
                        Sets the number of mismatches to allowed in a seed
                        alignment during multiseed alignment. Can be set to 0 or
                        1. Setting this higher makes alignment slower (often
                        much slower) but increases sensitivity.
                      id: '#allowed_mismatch_number'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 22 or 20 (depending on preset type and alignment mode)
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-L'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Seed substring length
                      description: >-
                        Sets the length of the seed substrings to align during
                        multiseed alignment. Smaller values make alignment
                        slower but more senstive. Must be > 3 and < 32. The
                        "Sensitive" preset is used by default, which sets this
                        option to 22 in "End-to-end" mode and to 20 in "Local"
                        mode.
                      id: '#seed_substring_length'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': '15'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--dpad'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Dynamic padding
                      description: >-
                        "Pads" dynamic programming problems by the given number
                        of columns on either side to allow gaps.
                      id: '#dynamic_padding'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': '4'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--gbar'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disallow gaps
                      description: >-
                        Disallow gaps within the given number of positions of
                        the beginning or end of the read.
                      id: '#disallow_gaps'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--ignore-quals'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Ignore qualities
                      description: >-
                        When calculating a mismatch penalty, always consider the
                        quality value at the mismatched position to be the
                        highest possible, regardless of the actual value. I.e.
                        treat all quality values as 30 on Phred scale. This is
                        also the default behavior when the input doesn't specify
                        quality values (e.g. when processing .fasta reads).
                      id: '#ignore_qualities'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--nofw'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Don't align forward
                      description: >-
                        If this option is specified, Bowtie2 will not attempt to
                        align unpaired reads to the forward (Watson) reference
                        strand. In paired-end mode, "Don't align forward" and
                        "Don't align reverse complement" pertain to the
                        fragments; i.e. specifying "Don't align forward" causes
                        Bowtie2 to explore only those paired-end configurations
                        corresponding to fragments from the reverse-complement
                        (Crick) strand. Default: both strands enabled.
                      id: '#dont_align_forward'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--norc'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Don't align reverse complement
                      description: >-
                        If this option is specified, Bowtie2 will not attempt to
                        align unpaired reads against the reverse-complement
                        (Crick) reference strand. In paired-end mode, "Don't
                        align forward" and "Don't align reverse complement"
                        pertain to the fragments; i.e. specifying "Don't align
                        forward" causes Bowtie2 to explore only those paired-end
                        configurations corresponding to fragments from the
                        reverse-complement (Crick) strand. Default: both strands
                        enabled.
                      id: '#dont_align_reverse_complement'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-1mm-upfront'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable 1 mismatch alignments
                      description: >-
                        By default, Bowtie2 will attempt to find either an exact
                        or a 1-mismatch end-to-end alignment for the read before
                        trying the multiseed heuristic. Such alignments can be
                        found very quickly, and many short read alignments have
                        exact or near-exact end-to-end alignments. However, this
                        can lead to unexpected alignments when the user also
                        sets options governing the multiseed heuristic, like
                        "Seed substring length" (-L) and "Allowed mismatch
                        number" (-N). For instance, if the user specifies 0 for
                        "Allowed mismatch number" and "Seed substring length"
                        equal to the length of the read, the user will be
                        surprised to find 1-mismatch alignments reported. This
                        option prevents Bowtie2 from searching for 1-mismatch
                        end-to-end alignments before using the multiseed
                        heuristic, which leads to the expected behavior when
                        combined with options such as "Seed substring length"
                        and "Allowed mismatch number". This comes at the expense
                        of speed.
                      id: '#disable_1_mismatch_alignments'
                    - 'sbg:category': Presets
                      'sbg:toolDefaultValue': Sensitive
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Very fast
                            - Fast
                            - Sensitive
                            - Very sensitive
                          name: preset_option
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              var preset_option = $job.inputs.preset_option
                              var alignment_mode = $job.inputs.alignment_mode
                              
                              var presets = {
                                "Very fast": "--very-fast",
                                "Fast": "--fast",
                                "Sensitive": "--sensitive",
                                "Very sensitive": "--very-sensitive"
                              }
                              if (alignment_mode == "Local" && preset_option) {
                                return presets[preset_option].concat("-local")
                              }
                              else if (preset_option){
                                return presets[preset_option]
                              }
                            }
                        'sbg:cmdInclude': true
                      label: Preset
                      description: >-
                        Preset options for "Seed extension attempts" (-D), "Max
                        number of re-seed" (-R), "Allowed mismatch number" (-N),
                        "Seed substring length" (-L) and "Interval function"
                        (-i) parameters. Values for these options vary depending
                        on whether the "Local" or "End-to-end" mode is selected
                        under "Alignment mode".
                      id: '#preset_option'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '0 for "End-to-end" mode, 2 for "Local" mode'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--ma'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Set match bonus
                      description: >-
                        Sets the match bonus. The given number is added to the
                        alignment score for each position where a read character
                        aligns to a reference character and the characters
                        match.
                      id: '#set_match_bonus'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '6'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--mp'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Maximum mismatch penalty
                      description: >-
                        Sets the maximum penalty for mismatch. Lower quality =
                        lower penalty.
                      id: '#maximum_mismatch_penalty'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '1'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--np'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Ambiguous character penalty
                      description: >-
                        Sets penalty for positions where the read, reference, or
                        both, contain an ambiguous character such as N.
                      id: '#ambiguous_character_penalty'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '5,3'
                      type:
                        - 'null'
                        - type: array
                          items: int
                      inputBinding:
                        position: 0
                        prefix: '--rdg'
                        separate: true
                        itemSeparator: ','
                        'sbg:cmdInclude': true
                      label: Read gap penalties
                      description: >-
                        Sets the read gap open (first value) and extend (second
                        value) penalty, respectively. A read gap of length N
                        gets a penalty of <gap-open-penalty> + N *
                        <gap-extend-penalty>.
                      id: '#read_gap_penalties'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '5,3'
                      type:
                        - 'null'
                        - type: array
                          items: int
                      inputBinding:
                        position: 0
                        prefix: '--rfg'
                        separate: true
                        itemSeparator: ','
                        'sbg:cmdInclude': true
                      label: Reference gap penalties
                      description: >-
                        Sets the reference gap open (first value) and extend
                        (second value) penalty, respectively. A reference gap of
                        length N gets a penalty of <gap-open-penalty> + N *
                        <gap-extend-penalty>.
                      id: '#reference_gap_penalties'
                    - 'sbg:category': Reporting
                      'sbg:toolDefaultValue': '-'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-k'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Report k alignments
                      description: >-
                        By default, Bowtie2 searches for distinct, valid
                        alignments for each read. When it finds a valid
                        alignment, it continues looking for alignments that are
                        nearly as good or better. The best alignment found is
                        reported (randomly selected from among best if tied).
                        Information about the best alignments is used to
                        estimate mapping quality and to set SAM optional fields,
                        such as AS:i and XS:i. When "Report k alignments" is
                        specified, however, Bowtie2 behaves differently.
                        Instead, it searches for at most <given-number>
                        distinct, valid alignments for each read. The search
                        terminates when it can't find more distinct valid
                        alignments, or when it finds <given-number>, whichever
                        happens first. All alignments found are reported in
                        descending order by alignment score. The alignment score
                        for a paired-end alignment equals the sum of the
                        alignment scores of the individual mates. Each reported
                        read or pair alignment beyond the first has the SAM
                        'secondary' bit (which equals 256) set in its FLAGS
                        field. For reads that have more than <given-number>
                        distinct, valid alignments, Bowtie2 does not gaurantee
                        that the <given-number> alignments reported are the best
                        possible in terms of alignment score. "Report k
                        alignments" is mutually exclusive with "Report all
                        alignments". Note: Bowtie 2 is not designed with large
                        values for "Report k alignments" in mind, and when
                        aligning reads to long, repetitive genomes alignment can
                        be very, very slow.
                      id: '#report_k_alignments'
                    - 'sbg:altPrefix': '-a'
                      'sbg:category': Reporting
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--all'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Report all alignments
                      description: >-
                        Like "Report k alignments" but with no upper limit on
                        number of alignments to search for. "Report all
                        alignments" is mutually exclusive with "Report k
                        alignments".
                      id: '#report_all_alignments'
                    - 'sbg:category': Effort
                      'sbg:toolDefaultValue': '15'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-D'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Seed extension attempts
                      description: >-
                        Maximum number of to consecutive seed extension attempts
                        that can "fail" before Bowtie2 moves on, using the
                        alignments found so far. A seed extension "fails" if it
                        does not yield a new best or a new second-best
                        alignment.
                      id: '#seed_extension_attempts'
                    - 'sbg:category': Effort
                      'sbg:toolDefaultValue': '2'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-R'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Max number of re-seed
                      description: >-
                        Given number is the maximum number of times Bowtie2 will
                        're-seed' reads with repetitive seeds. When
                        're-seeding', Bowtie2 simply chooses a new set of reads
                        (same length, same number of mismatches allowed) at
                        different offsets and searches for more alignments. A
                        read is considered to have repetitive seeds if the total
                        number of seed hits divided by the number of seeds that
                        aligned at least once is greater than 300.
                      id: '#max_number_of_re_seed'
                    - 'sbg:altPrefix': '-I'
                      'sbg:category': Paired-end
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--minins'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Minimum fragment length
                      description: >-
                        The minimum fragment length for valid paired-end
                        alignments. E.g. if 60 is specified for "Minimum
                        fragment length" (-I) and a paired-end alignment
                        consists of two 20-bp alignments in the appropriate
                        orientation with a 20-bp gap between them, that
                        alignment is considered valid (as long as "Maximum
                        fragment length" (-X) is also satisfied). A 19-bp gap
                        would not be valid in that case. If trimming options -3
                        or -5 are also used, the "Minimum fragment length"
                        constraint is applied with respect to the untrimmed
                        mates. The larger the difference between "Minimum
                        fragment length" and "Maximum fragment length", the
                        slower Bowtie2 will run. This is because larger
                        differences bewteen those two require that Bowtie2 scan
                        a larger window to determine if a concordant alignment
                        exists. For typical fragment length ranges (200 to 400
                        nucleotides), Bowtie2 is very efficient.
                      id: '#minimum_fragment_length'
                    - 'sbg:altPrefix': '-X'
                      'sbg:category': Paired-end
                      'sbg:toolDefaultValue': '500'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--maxins'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Maximum fragment length
                      description: >-
                        The maximum fragment length for valid paired-end
                        alignments. E.g. if "Maximum fragment length" (-X) 100
                        is specified and a paired-end alignment consists of two
                        20-bp alignments in the proper orientation with a 60-bp
                        gap between them, that alignment is considered valid (as
                        long as "Minimum fragment length" (-I) is also
                        satisfied). A 61-bp gap would not be valid in that case.
                        If trimming options -3 or -5 are also used, the "Maximum
                        fragment length" constraint is applied with respect to
                        the untrimmed mates, not the trimmed mates. The larger
                        the difference between "Minimum fragment length" and
                        "Maximum fragment length", the slower Bowtie2 will run.
                        This is because larger differences bewteen those two
                        require that Bowtie2 scan a larger window to determine
                        if a concordant alignment exists. For typical fragment
                        length ranges (200 to 400 nucleotides), Bowtie2 is very
                        efficient.
                      id: '#maximum_fragment_length'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': '--fr'
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - '--fr'
                            - '--rf'
                            - '--ff'
                          name: mates_alignment_orientation
                      inputBinding:
                        position: 0
                        separate: true
                        'sbg:cmdInclude': true
                      label: Mates alignment orientation
                      description: >-
                        The upstream/downstream mate orientations for a valid
                        paired-end alignment against the forward reference
                        strand. E.g., if --fr is specified and there is a
                        candidate paired-end alignment where mate 1 appears
                        upstream of the reverse complement of mate 2 and the
                        fragment length constraints ("Minimum fragment length"
                        (-I) and "Maximum fragment length" (-X)) are met, that
                        alignment is valid. Also, if mate 2 appears upstream of
                        the reverse complement of mate 1 and all other
                        constraints are met, that too is valid. --rf likewise
                        requires that an upstream mate1 be reverse-complemented
                        and a downstream mate2 be forward-oriented. --ff
                        requires both an upstream mate 1 and a downstream mate 2
                        to be forward-oriented. Default orientation --fr is
                        appropriate for Illumina's Paired-end Sequencing Assay.
                      id: '#mates_alignment_orientation'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-mixed'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable unpaired alignments
                      description: >-
                        By default, when Bowtie2 cannot find a concordant or
                        discordant alignment for a pair, it then tries to find
                        alignments for the individual mates. This option
                        disables that behavior.
                      id: '#disable_unpaired_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-discordant'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable discordant alignments
                      description: >-
                        By default, Bowtie2 looks for discordant alignments if
                        it cannot find any concordant alignments. A discordant
                        alignment is an alignment where both mates align
                        uniquely, but that does not satisfy the paired-end
                        constraints (--fr/--rf/--ff, -I, -X). This option
                        disables that behavior.
                      id: '#disable_discordant_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-dovetail'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable dovetail alignments
                      description: >-
                        If the mates "dovetail", that is if one mate alignment
                        extends past the beginning of the other such that the
                        wrong mate begins upstream, consider that to be
                        non-concordant.
                      id: '#disable_dovetail_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-contain'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable containing alignments
                      description: >-
                        If one mate alignment contains the other, consider that
                        to be non-concordant.
                      id: '#disable_containing_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-overlap'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable overlapping alignments
                      description: >-
                        If one mate alignment overlaps the other at all,
                        consider that to be non-concordant.
                      id: '#disable_overlapping_alignments'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-head'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Suppress header lines
                      description: Suppress SAM header lines (starting with @).
                      id: '#suppress_header_lines'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-sq'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Suppress SQ header lines
                      description: Suppress @SQ SAM header lines.
                      id: '#suppress_sq_header_lines'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': id
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        prefix: '--rg-id'
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              return ($job.inputs.read_group_id || "id") 
                            }
                        'sbg:cmdInclude': true
                      label: Set the read group ID
                      description: >-
                        Set the read group ID text. This causes the SAM @RG
                        header line to be printed, with the given text as the
                        value associated with the ID: tag. It also causes the
                        RG:Z: extra field to be attached to each SAM output
                        record, with value set to this text.
                      id: '#read_group_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - LS 454
                            - Helicos
                            - Illumina
                            - ABI SOLiD
                            - Ion Torrent PGM
                            - PacBio
                          name: platform
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              if($job.inputs.platform)
                                return "--rg PL:" +$job.inputs.platform.replace(/ /g,"_")
                              else if([].concat($job.inputs.read_sequence)[0].metadata){
                                if ([].concat($job.inputs.read_sequence)[0].metadata.platform) {
                                  return "--rg PL:" +[].concat($job.inputs.read_sequence)[0].metadata.platform.replace(/ /g,"_")
                                }
                              }
                              else {
                                return ""
                              }
                            }
                        'sbg:cmdInclude': true
                      label: Platform
                      description: >-
                        Specify the version of the technology that was used for
                        sequencing or assaying. Default: inferred from metadata.
                      id: '#platform'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.sample_id)\n    return \"--rg SM:\" +$job.inputs.sample_id\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.sample_id)\n  \t\treturn \"--rg SM:\" +[].concat($job.inputs.read_sequence)[0].metadata.sample_id\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Sample
                      description: Specify the sample ID for RG line.
                      id: '#sample_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.library_id)\n    return \"--rg LB:\" +$job.inputs.library_id\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.library_id)\n  \t\treturn \"--rg LB:\" +[].concat($job.inputs.read_sequence)[0].metadata.library_id\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Library
                      description: Specify the library ID for RG line.
                      id: '#library_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.platform_unit_id)\n    return \"--rg PU:\" +$job.inputs.platform_unit_id\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.platform_unit_id)\n  \t\treturn \"--rg PU:\" +[].concat($job.inputs.read_sequence)[0].metadata.platform_unit_id\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Platform unit
                      description: Specify the platform unit ID for RG line.
                      id: '#platform_unit_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.sequencing_center)\n    return \"--rg CN:\" +$job.inputs.sequencing_center\n    else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.seq_center)\n  \t\treturn \"--rg CN:\" +[].concat($job.inputs.read_sequence)[0].metadata.seq_center\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Sequencing center
                      description: Specify the sequencing center for RG line.
                      id: '#sequencing_center'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.median_fragment_length)\n    return \"--rg PI:\" +$job.inputs.median_fragment_length\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.median_fragment_length)\n  \t\treturn \"--rg PI:\" +[].concat($job.inputs.read_sequence)[0].metadata.median_fragment_length\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Median fragment length
                      description: Specify the median fragment length for RG line.
                      id: '#median_fragment_length'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--omit-sec-seq'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Omit SEQ and QUAL
                      description: >-
                        When printing secondary alignments, Bowtie 2 by default
                        will write out the SEQ and QUAL strings. Specifying this
                        option causes Bowtie 2 to print an asterisk ('*') in
                        those fields instead.
                      id: '#omit_seq_and_qual'
                    - 'sbg:category': Performance
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--reorder'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Reorder output
                      description: >-
                        Guarantees that output SAM records are printed in an
                        order corresponding to the order of the reads in the
                        original input file. Specifying "Reorder output" causes
                        Bowtie2 to run somewhat slower and use somewhat more
                        memory.
                      id: '#reorder_output'
                    - 'sbg:category': Other
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--seed'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Set seed
                      description: Set the seed for pseudo-random number generator.
                      id: '#set_seed'
                    - 'sbg:category': Other
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--non-deterministic'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Non deterministic
                      description: >-
                        Normally, Bowtie2 re-initializes its pseudo-random
                        generator for each read. It seeds the generator with a
                        number derived from (a) the read name, (b) the
                        nucleotide sequence, (c) the quality sequence, (d) the
                        value of the "Set seed" option. This means that if two
                        reads are identical (same name, same nucleotides, same
                        qualities) Bowtie2 will find and report the same
                        alignment(s) for both, even if there was ambiguity. When
                        "Non deterministic" is specified, Bowtie2 re-initializes
                        its pseudo-random generator for each read using the
                        current time. This means that Bowtie2 will not
                        necessarily report the same alignment for two identical
                        reads. This is counter-intuitive for some users, but
                        might be more appropriate in situations where the input
                        consists of many identical reads.
                      id: '#non_deterministic'
                    - 'sbg:category': Interval function
                      'sbg:toolDefaultValue': Square-root
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Constant
                            - Linear
                            - Square-root
                            - Natural log
                          name: function_i
                      label: Interval function
                      description: >-
                        Sets a function type F in function f governing the
                        interval between seed substrings, to use during
                        multiseed alignment. The interval function f is f(x) = A
                        + B * F(x), where x is the read length. By default,
                        function F is set to 'Square-root', Constant A to 1 and
                        Coefficient B to 1.15 or 0.75 for "End-to-end" and
                        "Local" mode respectively.
                      id: '#function_i'
                    - 'sbg:category': Interval function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '1'
                      type:
                        - 'null'
                        - string
                      label: 'Constant A [ interval function ]'
                      description: >-
                        Sets a constant A in function governing the interval
                        between seed substrings to use during multiseed
                        alignment. The interval function f is f(x) = A + B *
                        F(x), where x is the read length.
                      id: '#constant_i_a'
                    - 'sbg:category': Interval function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': 1.15 or 0.75 (depending on "Alignment mode")
                      type:
                        - 'null'
                        - string
                      label: 'Coefficient B [ interval function ]'
                      description: >-
                        Sets a coefficient B in function governing the interval
                        between seed substrings to use during multiseed
                        alignment. The interval function f is f(x) = A + B *
                        F(x), where x is the read length. Default: 1.15 in
                        "End-to-end" mode and 0.75 in "Local" mode.
                      id: '#coefficient_i_b'
                    - 'sbg:category': Ambiguous chars function
                      'sbg:toolDefaultValue': Linear
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Constant
                            - Linear
                            - Square-root
                            - Natural log
                          name: function_n_ceil
                      label: Ambiguous chars function
                      description: >-
                        Sets a function type F in function governing the maximum
                        number of ambiguous characters (usually Ns and/or .s)
                        allowed in a read as a function of read length. The
                        N-ceiling function f is f(x) = A + B * F(x), where x is
                        the read length. Reads exceeding this ceiling are
                        filtered out.
                      id: '#function_n_ceil'
                    - 'sbg:category': Ambiguous chars function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - string
                      label: 'Constant A [ ambiguous chars function ]'
                      description: >-
                        Sets a constant A in function governing the maximum
                        number of ambiguous characters (usually Ns and/or .s)
                        allowed in a read as a function of read length. The
                        N-ceiling function f is f(x) = A + B * F(x), where x is
                        the read length. Reads exceeding this ceiling are
                        filtered out.
                      id: '#constant_nceil_a'
                    - 'sbg:category': Ambiguous chars function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '0.15'
                      type:
                        - 'null'
                        - string
                      label: 'Coefficient B [ ambiguous chars function ]'
                      description: >-
                        Sets a coefficient B in function governing the maximum
                        number of ambiguous characters (usually Ns and/or .s)
                        allowed in a read as a function of read length. The
                        N-ceiling function f is f(x) = A + B * F(x), where x is
                        the read length. Reads exceeding this ceiling are
                        filtered out.
                      id: '#coefficient_nceil_b'
                    - 'sbg:category': Alignment score function
                      'sbg:toolDefaultValue': Natural log or Linear (depending on "Alignment mode")
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Constant
                            - Linear
                            - Square-root
                            - Natural log
                          name: function_score_min
                      label: Alignment score function
                      description: >-
                        Sets a function type F in function governing the minimum
                        alignment score needed for an alignment to be considered
                        "valid" (i.e. good enough to report). This is a function
                        of read length. The minimum-score function f is f(x) = A
                        + B * F(x), where x is the read length. By default,
                        function F is set to "Natural log" or "Linear", Constant
                        A to 20 or -0.6 and Coefficient B to 8 or -0.6 depending
                        on the "Alignment mode": "End-to-end" or "Local"
                        respectively.
                      id: '#function_score_min'
                    - 'sbg:category': Alignment score function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': 20 or -0.6 (depending on "Alignment mode")
                      type:
                        - 'null'
                        - string
                      label: 'Constant A [ alignment score function ]'
                      description: >-
                        Sets a constant A in function governing the minimum
                        alignment score needed for an alignment to be considered
                        'valid' (i.e. good enough to report). This is a function
                        of read length. The minimum-score function f is f(x) = A
                        + B * F(x), where x is the read length. Default: 20 in
                        "End-to-end" mode and -0.6 in "Local" mode.
                      id: '#constant_scoremin_a'
                    - 'sbg:category': Alignment score function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': 8 or -0.6 (depending on "Alignment mode")
                      type:
                        - 'null'
                        - string
                      label: 'Coefficient B [ alignment score function ]'
                      description: >-
                        Sets a coefficient B in function governing the minimum
                        alignment score needed for an alignment to be considered
                        'valid' (i.e. good enough to report). This is a function
                        of read length. The minimum-score function f is f(x) = A
                        + B * F(x), where x is the read length. Default: 8 in
                        "End-to-end" mode and -0.6 in "Local" mode.
                      id: '#coefficient_scoremin_b'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: unpaired_unaligned_reads
                      label: Unpaired unaligned reads
                      description: >-
                        Output unpaired reads that fail to align. These reads
                        correspond to the SAM records with the FLAGS 0x4 bit set
                        and neither the 0x40 nor 0x80 bits set. If "gzip
                        compressed" is specified, output will be gzip
                        compressed. If "bzip2 compressed" is specified, output
                        will be bzip2 compressed. Reads written in this way will
                        appear exactly as they did in the input file, without
                        any modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#unpaired_unaligned_reads'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: unpaired_aligned_reads
                      label: Unpaired aligned reads
                      description: >-
                        Output unpaired reads that align at least once. These
                        reads correspond to the SAM records with the FLAGS 0x4,
                        0x40, and 0x80 bits unset. If "gzip compressed" is
                        specified, output will be gzip compressed. If "bzip2
                        compressed" is specified, output will be bzip2
                        compressed. Reads written in this way will appear
                        exactly as they did in the input file, without any
                        modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#unpaired_aligned_reads'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: paired_unaligned_reads
                      label: Paired unaligned reads
                      description: >-
                        Output paired-end reads that fail to align concordantly.
                        These reads correspond to the SAM records with the FLAGS
                        0x4 bit set and either the 0x40 or 0x80 bit set
                        (depending on whether it's mate #1 or #2). .1 and .2
                        strings are added to the filename to distinguish which
                        file contains mate #1 and mate #2. If "gzip compressed"
                        is specified, output will be gzip compressed. If "bzip2
                        compressed" is specified, output will be bzip2
                        compressed. Reads written in this way will appear
                        exactly as they did in the input file, without any
                        modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#paired_unaligned_reads'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: paired_aligned_reads
                      label: Paired aligned reads
                      description: >-
                        Output paired-end reads that align concordantly at least
                        once. These reads correspond to the SAM records with the
                        FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set
                        (depending on whether it's mate #1 or #2). .1 and .2
                        strings are added to the filename to distinguish which
                        file contains mate #1 and mate #2. If "gzip compressed"
                        is specified, output will be gzip compressed. If "bzip2
                        compressed" is specified, output will be bzip2
                        compressed. Reads written in this way will appear
                        exactly as they did in the input file, without any
                        modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#paired_aligned_reads'
                    - 'sbg:category': Input files
                      type:
                        - type: array
                          items: File
                      label: Read sequence
                      description: >-
                        Read sequence in FASTQ or FASTA format. COuld be also
                        gzip'ed (extension .gz) or bzip2'ed (extension .bz2). In
                        case of paired-end alignment it is crucial to set
                        metadata 'paired-end' field to 1/2.
                      'sbg:fileTypes': >-
                        FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA,
                        FQ, FASTQ.GZ, FQ.GZ, FASTQ.BZ2, FQ.BZ2
                      id: '#read_sequence'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-unal'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Suppress SAM records for unaligned reads
                      description: Suppress SAM records for reads that failed to align.
                      id: '#suppress_sam_records'
                    - 'sbg:category': Input
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '-f'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Input FASTA files
                      description: >-
                        Reads (specified with <m1>, <m2>, <s>) are FASTA files.
                        FASTA files usually have extension .fa, .fasta, .mfa,
                        .fna or similar. FASTA files do not have a way of
                        specifying quality values, so when -f is set, the result
                        is as if --ignore-quals is also set.
                      id: '#input_fasta_files'
                    - 'sbg:altPrefix': '-threads'
                      'sbg:category': Performance
                      'sbg:toolDefaultValue': '8'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n\tif($job.inputs.threads)\n    {\n    \treturn \" -p \" + $job.inputs.threads\n    }\n  \telse\n    {\n    \treturn \" -p 8 \"\n    }\n}"
                        'sbg:cmdInclude': true
                      label: Number of threads
                      description: >-
                        Launch NTHREADS parallel search threads (default: 1).
                        Threads will run on separate processors/cores and
                        synchronize when parsing reads and outputting
                        alignments. Searching for alignments is highly parallel,
                        and speedup is close to linear. Increasing -p increases
                        Bowtie 2's memory footprint. E.g. when aligning to a
                        human genome index, increasing -p from 1 to 8 increases
                        the memory footprint by a few hundred megabytes. This
                        option is only available if bowtie is linked with the
                        pthreads library (i.e. if BOWTIE_PTHREADS=0 is not
                        specified at build time).
                      id: '#threads'
                  outputs:
                    - type:
                        - 'null'
                        - File
                      label: Result SAM file
                      description: >-
                        SAM file containing the results of the alignment. It
                        contains both aligned and unaligned reads.
                      'sbg:fileTypes': SAM
                      outputBinding:
                        streamable: false
                        glob: '*.sam'
                        'sbg:metadata':
                          __inherit__: fasta_reference
                        'sbg:inheritMetadataFrom': '#read_sequence'
                      id: '#result_sam_file'
                    - type:
                        - 'null'
                        - type: array
                          items: File
                      label: Aligned reads only
                      description: FASTQ file with reads that align at least once.
                      'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                      outputBinding:
                        streamable: false
                        glob: '*_aligned*'
                        'sbg:metadata':
                          __inherit__: fasta_reference
                        'sbg:inheritMetadataFrom': '#read_sequence'
                      id: '#aligned_reads_only'
                    - type:
                        - 'null'
                        - type: array
                          items: File
                      label: Unaligned reads only
                      description: FASTQ file with reads that failed to align.
                      'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                      outputBinding:
                        streamable: false
                        glob: '*_unaligned*'
                        'sbg:metadata':
                          __inherit__: fasta_reference
                        'sbg:inheritMetadataFrom': '#read_sequence'
                      id: '#unaligned_reads_only'
                  requirements:
                    - class: ExpressionEngineRequirement
                      engineCommand: cwl-engine.js
                      id: '#cwl-js-engine'
                      requirements:
                        - class: DockerRequirement
                          dockerPull: rabix/js-engine
                  hints:
                    - class: 'sbg:CPURequirement'
                      value: 8
                    - class: 'sbg:MemRequirement'
                      value: 6000
                    - class: DockerRequirement
                      dockerImageId: 029d3a264215
                      dockerPull: 'images.sbgenomics.com/ana_d/bowtie2:2.2.6'
                  arguments:
                    - position: 100
                      prefix: ''
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: "cmd = \"\"\nreads = [].concat($job.inputs.read_sequence)\nreads1 = [];\nreads2 = [];\nu_reads = [];\nfor (var i = 0; i < reads.length; i++){\n    if (reads[i].metadata.paired_end == 1){\n      reads1.push(reads[i].path);\n    }\n    else if (reads[i].metadata.paired_end == 2){\n      reads2.push(reads[i].path);\n    }\n  else {\n  \tu_reads.push(reads[i].path);\n   }\n  }\nif (reads1.length > 0 & reads1.length == reads2.length){\n\tcmd = \"-1 \" + reads1.join(\",\") + \" -2 \" + reads2.join(\",\");\n}\nif (u_reads.length > 0){\n\tcmd = \" -U \" + u_reads.join(\",\");\n}\ncmd\n"
                    - position: 101
                      prefix: '-S'
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                            var A= array.concat().sort(), 
                                a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                            while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                            return a1.substring(0, i);
                            }
                            path_list = []
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            return "./".concat(common_prefix.replace( /\-$|\_$|\.$/, '' ), ".", "sam")
                          }
                    - position: 0
                      prefix: '-x'
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |
                          {
                            var index_prefix = $job.inputs.bowtie_index_archive.metadata.reference_genome
                            return index_prefix
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.unpaired_unaligned_reads && $job.inputs.unpaired_unaligned_reads != "None") {
                              if ($job.inputs.unpaired_unaligned_reads == "raw") {
                                return "--un ".concat(common_prefix, ".unpaired_unaligned.fastq")
                              }
                              else if ($job.inputs.unpaired_unaligned_reads == "gzip compressed") {
                                return "--un ".concat(common_prefix, ".unpaired_unaligned.fastq.gz")
                              }
                              else if ($job.inputs.unpaired_unaligned_reads == "bzip2 compressed") {
                                return "--un ".concat(common_prefix, ".unpaired_unaligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.unpaired_aligned_reads && $job.inputs.unpaired_aligned_reads != "None") {
                              if ($job.inputs.unpaired_aligned_reads == "raw") {
                                return "--al ".concat(common_prefix, ".unpaired_aligned.fastq")
                              }
                              else if ($job.inputs.unpaired_aligned_reads == "gzip compressed") {
                                return "--al ".concat(common_prefix, ".unpaired_aligned.fastq.gz")
                              }
                              else if ($job.inputs.unpaired_aligned_reads == "bzip2 compressed") {
                                return "--al ".concat(common_prefix, ".unpaired_aligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.paired_unaligned_reads && $job.inputs.paired_unaligned_reads != "None") {
                              if ($job.inputs.paired_unaligned_reads == "raw") {
                                return "--un-conc ".concat(common_prefix, ".paired_unaligned.fastq")
                              }
                              else if ($job.inputs.paired_unaligned_reads == "gzip compressed") {
                                return "--un-conc ".concat(common_prefix, ".paired_unaligned.fastq.gz")
                              }
                              else if ($job.inputs.paired_unaligned_reads == "bzip2 compressed") {
                                return "--un-conc ".concat(common_prefix, ".paired_unaligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.paired_aligned_reads && $job.inputs.paired_aligned_reads != "None") {
                              if ($job.inputs.paired_aligned_reads == "raw") {
                                return "--al-conc ".concat(common_prefix, ".paired_aligned.fastq")
                              }
                              else if ($job.inputs.paired_aligned_reads == "gzip compressed") {
                                return "--al-conc ".concat(common_prefix, ".paired_aligned.fastq.gz")
                              }
                              else if ($job.inputs.paired_aligned_reads == "bzip2 compressed") {
                                return "--al-conc ".concat(common_prefix, ".paired_aligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      prefix: ''
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            var functions = {
                              "Constant": "C",
                              "Linear": "L",
                              "Square-root": "S",
                              "Natural log": "G"
                            }
                            function_type = $job.inputs.function_i
                            value_list = [functions[function_type], $job.inputs.constant_i_a, $job.inputs.coefficient_i_b]
                            if (functions[function_type] && $job.inputs.constant_i_a && $job.inputs.coefficient_i_b) {
                              return "-i ".concat(value_list.join(","))
                            }
                          }
                    - position: 0
                      prefix: ''
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            var functions = {
                              "Constant": "C",
                              "Linear": "L",
                              "Square-root": "S",
                              "Natural log": "G"
                            }
                            function_type = $job.inputs.function_n_ceil
                            value_list = [functions[function_type], $job.inputs.constant_nceil_a, $job.inputs.coefficient_nceil_b]
                            if (functions[function_type] && $job.inputs.constant_nceil_a && $job.inputs.coefficient_nceil_b) {
                              return "--n-ceil ".concat(value_list.join(","))
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            var functions = {
                              "Constant": "C",
                              "Linear": "L",
                              "Square-root": "S",
                              "Natural log": "G"
                            }
                            function_type = $job.inputs.function_score_min
                            
                            value_list = [functions[function_type], $job.inputs.constant_scoremin_a, $job.inputs.coefficient_scoremin_b]
                            if (functions[function_type] && $job.inputs.constant_scoremin_a && $job.inputs.coefficient_scoremin_b) {
                              return "--score-min ".concat(value_list.join(","))
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            meta_qual="";
                            if ([].concat($job.inputs.read_sequence)[0].metadata){
                              if ([].concat($job.inputs.read_sequence)[0].metadata.quality_scale){
                                meta_qual = [].concat($job.inputs.read_sequence)[0].metadata.quality_scale
                              }
                            }
                            
                            if ($job.inputs.quality_scale == "Phred+33") {
                              return "--phred33"
                            }
                            else if ($job.inputs.quality_scale == "Phred+64") {
                              return "--phred64"
                            }
                            else if ($job.inputs.quality_scale == "Solexa") {
                              return "--solexa-quals"
                            }
                            else if ($job.inputs.quality_scale == "Auto-detect") {
                              if (meta_qual == "solexa") {
                                return "--solexa-quals"
                              }
                              else if (meta_qual == "illumina13" || meta_qual == "illumina15") {
                                return "--phred64"
                              }
                            }
                          }
                    - position: 102
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            common_prefix = sharedStart(path_list)
                            
                            gzip = "gzip compressed"
                            bzip = "bzip2 compressed"
                            paired_aligned = $job.inputs.paired_aligned_reads
                            paired_unaligned = $job.inputs.paired_unaligned_reads
                            aligned_first = common_prefix.concat(".paired_aligned.fastq.1")
                            aligned_second = common_prefix.concat(".paired_aligned.fastq.2")
                            aligned_first_mv = common_prefix.concat(".paired_aligned.1.fastq")
                            aligned_second_mv = common_prefix.concat(".paired_aligned.2.fastq")
                            
                            unaligned_first = common_prefix.concat(".paired_unaligned.fastq.1")
                            unaligned_second = common_prefix.concat(".paired_unaligned.fastq.2")
                            unaligned_first_mv = common_prefix.concat(".paired_unaligned.1.fastq")
                            unaligned_second_mv = common_prefix.concat(".paired_unaligned.2.fastq")
                            
                            aligned = ""
                            unaligned = ""
                            
                            if (paired_aligned && paired_aligned == gzip) {
                              aligned = "&& mv ".concat(aligned_first, ".gz ", aligned_first_mv, ".gz && mv ", aligned_second, ".gz ", aligned_second_mv, ".gz ") 
                            }
                            else if (paired_aligned && paired_aligned == bzip) {
                              aligned = "&& mv ".concat(aligned_first, ".bz2 ", aligned_first_mv, ".bz2 && mv ", aligned_second, ".bz2 ", aligned_second_mv, ".bz2 ")
                            }
                            if (paired_unaligned && paired_unaligned == gzip) {
                              unaligned = "&& mv ".concat(unaligned_first, ".gz ", unaligned_first_mv, ".gz && mv ", unaligned_second, ".gz ", unaligned_second_mv, ".gz")
                            }
                            else if (paired_unaligned && paired_unaligned == bzip) {
                              unaligned = "&& mv ".concat(unaligned_first, ".bz2 ", unaligned_first_mv, ".bz2 && mv ", unaligned_second, ".bz2 ", unaligned_second_mv, ".bz2")
                            }
                            
                            return aligned.concat(unaligned)
                          }
                  'sbg:appVersion':
                    - 'sbg:draft-2'
                  'sbg:categories':
                    - Alignment
                  'sbg:cmdPreview': >-
                    tar -xvf chr20_bowtie2-2.2.6.tar && rm -rf
                    chr20_bowtie2-2.2.6.tar &&  /opt/bowtie2-2.2.6/bowtie2 -x
                    chr20  --un mate.unpaired_unaligned.fastq.gz  --al
                    mate.unpaired_aligned.fastq.gz  --un-conc
                    mate.paired_unaligned.fastq  --al-conc
                    mate.paired_aligned.fastq.gz  -i
                    L,constant_i_a-string-value,coefficient_i_b-string-value 
                    --n-ceil
                    S,constant_nceil_a-string-value,coefficient_nceil_b-string-value 
                    --score-min
                    L,constant_scoremin_a-string-value,coefficient_scoremin_b-string-value 
                    --phred33  -1 /demo/test-data/mate1.fastq -2
                    /demo/test-data/mate2.fastq -S ./mate.sam  && mv
                    mate.paired_aligned.fastq.1.gz
                    mate.paired_aligned.1.fastq.gz && mv
                    mate.paired_aligned.fastq.2.gz
                    mate.paired_aligned.2.fastq.gz
                  'sbg:content_hash': null
                  'sbg:contributors':
                    - sevenbridges-pgc
                    - admin
                  'sbg:createdBy': sevenbridges-pgc
                  'sbg:createdOn': 1454426457
                  'sbg:id': admin/sbg-public-data/bowtie2-aligner/18
                  'sbg:image_url': null
                  'sbg:job':
                    allocatedResources:
                      cpu: 8
                      mem: 6000
                    inputs:
                      alignment_mode: Local
                      allowed_mismatch_number: '0'
                      bowtie_index_archive:
                        class: File
                        metadata:
                          reference_genome: chr20
                        path: /demo/test-data/chr20_bowtie2-2.2.6.tar
                        secondaryFiles: []
                        size: 0
                      coefficient_i_b: coefficient_i_b-string-value
                      coefficient_nceil_b: coefficient_nceil_b-string-value
                      coefficient_scoremin_b: coefficient_scoremin_b-string-value
                      constant_i_a: constant_i_a-string-value
                      constant_nceil_a: constant_nceil_a-string-value
                      constant_scoremin_a: constant_scoremin_a-string-value
                      disable_overlapping_alignments: false
                      disable_unpaired_alignments: false
                      function_i: Linear
                      function_n_ceil: Square-root
                      function_score_min: Linear
                      input_fasta_files: true
                      mates_alignment_orientation: '--rf'
                      paired_aligned_reads: gzip compressed
                      paired_unaligned_reads: raw
                      platform: ABI SOLiD
                      preset_option: Very fast
                      quality_scale: Phred+33
                      read_gap_penalties:
                        - 0
                      read_sequence:
                        - class: File
                          metadata:
                            file_format: fastq
                            paired_end: '1'
                            quality_scale: illumina15
                          path: /demo/test-data/mate1.fastq
                          secondaryFiles: []
                          size: 0
                        - metadata:
                            file_format: fastq
                            paired_end: '2'
                            qual_scale: illumina15
                          path: /demo/test-data/mate2.fastq
                          secondaryFiles: []
                      reference_gap_penalties:
                        - 0
                      sample_id: nn
                      suppress_sam_records: true
                      threads: null
                      unpaired_aligned_reads: gzip compressed
                      unpaired_unaligned_reads: gzip compressed
                  'sbg:latestRevision': 18
                  'sbg:license': GNU General Public License v3.0 only
                  'sbg:links':
                    - id: 'http://bowtie-bio.sourceforge.net/bowtie2/index.shtml'
                      label: Homepage
                    - id: >-
                        http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/
                      label: Download
                    - id: >-
                        http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html
                      label: Publication
                    - id: 'http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml'
                      label: Manual
                    - id: 'https://github.com/BenLangmead/bowtie2'
                      label: Source Code
                  'sbg:modifiedBy': admin
                  'sbg:modifiedOn': 1533302766
                  'sbg:project': admin/sbg-public-data
                  'sbg:projectName': SBG Public Data
                  'sbg:publisher': sbg
                  'sbg:revision': 18
                  'sbg:revisionNotes': |-
                    - changed labels for Function type
                    - changed type to string for coefficients and constants
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1454426457
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1454426458
                      'sbg:revision': 1
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1462903844
                      'sbg:revision': 2
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1465231336
                      'sbg:revision': 3
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1471881576
                      'sbg:revision': 4
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1471881576
                      'sbg:revision': 5
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1472054036
                      'sbg:revision': 6
                      'sbg:revisionNotes': >-
                        Redesigned to accept archive with Bowtie2 index files on
                        the input.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1472208214
                      'sbg:revision': 7
                      'sbg:revisionNotes': Additional information updated.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1472208214
                      'sbg:revision': 8
                      'sbg:revisionNotes': Additional information updated.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1477482292
                      'sbg:revision': 9
                      'sbg:revisionNotes': >-
                        Javascript for reads written to check if metadata is
                        null.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1497287393
                      'sbg:revision': 10
                      'sbg:revisionNotes': >-
                        expressions for single end reads and scatter mode
                        updated
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1497287393
                      'sbg:revision': 11
                      'sbg:revisionNotes': |-
                        changed Docker file
                        added select output_format option
                        changed expression for handling SE reads
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1497287393
                      'sbg:revision': 12
                      'sbg:revisionNotes': added  piping and pipe status check
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 13
                      'sbg:revisionNotes': Works with fasta files and with only one input file
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 14
                      'sbg:revisionNotes': Works with fasta files and with only one input file
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 15
                      'sbg:revisionNotes': Number of threads added
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 16
                      'sbg:revisionNotes': Number of threads added
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501775358
                      'sbg:revision': 17
                      'sbg:revisionNotes': 'Threads: default 8'
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1533302766
                      'sbg:revision': 18
                      'sbg:revisionNotes': |-
                        - changed labels for Function type
                        - changed type to string for coefficients and constants
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': Ben Langmead/John Hopkins University
                  'sbg:toolkit': Bowtie2
                  'sbg:toolkitVersion': 2.2.6
                  'sbg:validationErrors': []
                label: Bowtie2 Aligner
                'sbg:x': 222.015625
                'sbg:y': 160.125
            requirements: []
            'sbg:appVersion':
              - v1.0
              - 'sbg:draft-2'
            'sbg:content_hash': aed84ae83d3a7c6b130e43ea33af7d1377c85e82c6ca358ff534111b7e6b030d7
            'sbg:contributors':
              - rbowen_james
            'sbg:createdBy': rbowen_james
            'sbg:createdOn': 1626440046
            'sbg:id': mwonge/mwtest/hlahd-bowtie-samtools/6
            'sbg:image_url': >-
              https://cavatica.sbgenomics.com/ns/brood/images/mwonge/mwtest/hlahd-bowtie-samtools/6.png
            'sbg:latestRevision': 6
            'sbg:modifiedBy': rbowen_james
            'sbg:modifiedOn': 1627371058
            'sbg:project': mwonge/mwtest
            'sbg:projectName': zcc-cavatica
            'sbg:publisher': sbg
            'sbg:revision': 6
            'sbg:revisionNotes': null
            'sbg:revisionsInfo':
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626440046
                'sbg:revision': 0
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029037
                'sbg:revision': 1
                'sbg:revisionNotes': added outdir
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029530
                'sbg:revision': 2
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029871
                'sbg:revision': 3
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627104164
                'sbg:revision': 4
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627370861
                'sbg:revision': 5
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627371058
                'sbg:revision': 6
                'sbg:revisionNotes': null
            'sbg:sbgMaintained': false
            'sbg:toolAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
            'sbg:validationErrors': []
            'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
          label: HLA-HD Tumour DNA
          'sbg:x': 286.875
          'sbg:y': 616.6875
        - id: hla-hd_normal_dna
          in:
            - id: sample_id
              source: sample_id
            - id: minimum_read_length
              default: 0
              source: minimum_read_length_2
            - id: read_sequence
              source:
                - normal_dna_reads
            - id: bowtie_index_archive
              source: bowtie_index_archive
            - id: output_dir_name
              source: normal_dna_outdir_name
          out:
            - id: hlad_reads1
            - id: hla_reads2
            - id: filtered_BAM
            - id: bowtie_sam
            - id: aligned_reads_only
            - id: unaligned_reads_only
            - id: hlahd_results
            - id: hlahd_final_results
          run:
            class: Workflow
            cwlVersion: v1.0
            id: mwonge/mwtest/hlahd-bowtie-samtools/6
            label: hlahd_bowtie_samtools
            $namespaces:
              sbg: 'https://sevenbridges.com'
            inputs:
              - id: sample_id
                type: string
                label: Sample ID
                doc: ID of the sample being analysed.
                'sbg:x': 0
                'sbg:y': 14
              - id: minimum_read_length
                type: int
                'sbg:exposed': true
              - id: read_sequence
                'sbg:fileTypes': >-
                  FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA, FQ,
                  FASTQ.GZ, FQ.GZ, FASTQ.BZ2, FQ.BZ2
                type: 'File[]'
                'sbg:x': 0
                'sbg:y': 120.75
              - id: bowtie_index_archive
                'sbg:fileTypes': TAR
                type: File
                'sbg:x': 0
                'sbg:y': 334.25
              - id: output_dir_name
                type: string?
                'sbg:x': 0
                'sbg:y': 227.5
            outputs:
              - id: hlad_reads1
                outputSource:
                  - samtools_fastq/output_pe_1
                'sbg:fileTypes': FASTQ
                type: File?
                label: HLA Reads 1
                doc: >-
                  Reads 1 FASTQ containing HLA reads after Yara alignment,
                  Samtools View filtering and splitting into FASTQ files.
                'sbg:x': 1256.66357421875
                'sbg:y': 46.25
              - id: hla_reads2
                outputSource:
                  - samtools_fastq/output_pe_2
                'sbg:fileTypes': FASTQ
                type: File?
                label: HLA Reads 2
                doc: >-
                  Reads 2 FASTQ containing HLA reads after Yara alignment,
                  Samtools View filtering and splitting into FASTQ files.
                'sbg:x': 1256.66357421875
                'sbg:y': 153
              - id: filtered_BAM
                outputSource:
                  - samtools_view_1_9_cwl1_0/out_alignments
                'sbg:fileTypes': 'BAM, SAM, CRAM'
                type: File?
                label: Filtered BAM
                doc: Yara aligned BAM filtered by Samtools View.
                'sbg:x': 977.3401489257812
                'sbg:y': 227.5
              - id: bowtie_sam
                outputSource:
                  - bowtie2_aligner/result_sam_file
                'sbg:fileTypes': SAM
                type: File?
                'sbg:x': 547.8245239257812
                'sbg:y': 241.5
              - id: aligned_reads_only
                outputSource:
                  - bowtie2_aligner/aligned_reads_only
                'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                type: 'File[]?'
                'sbg:x': 547.8245239257812
                'sbg:y': 348.25
              - id: unaligned_reads_only
                outputSource:
                  - bowtie2_aligner/unaligned_reads_only
                'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                type: 'File[]?'
                'sbg:x': 547.8245239257812
                'sbg:y': 0
              - id: hlahd_results
                outputSource:
                  - hla_hd/hlahd_results
                type: Directory
                'sbg:x': 1571.572265625
                'sbg:y': 120.75
              - id: hlahd_final_results
                outputSource:
                  - hla_hd/hlahd_final_results
                type: File
                'sbg:x': 1571.572265625
                'sbg:y': 227.5
            steps:
              - id: samtools_view_1_9_cwl1_0
                in:
                  - id: output_format
                    default: BAM
                  - id: fast_bam_compression
                    default: true
                  - id: filter_exclude_any
                    default: 4
                  - id: in_alignments
                    source: bowtie2_aligner/result_sam_file
                out:
                  - id: out_alignments
                  - id: reads_not_selected_by_filters
                  - id: alignement_count
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: admin/sbg-public-data/samtools-view-1-9-cwl1-0/6
                  baseCommand:
                    - /opt/samtools-1.9/samtools
                    - view
                  inputs:
                    - 'sbg:category': File inputs
                      id: in_index
                      type: File?
                      label: Index file
                      doc: This tool requires index file for some use cases.
                      'sbg:fileTypes': 'BAI, CRAI, CSI'
                    - 'sbg:altPrefix': '-O'
                      'sbg:category': Config inputs
                      'sbg:toolDefaultValue': SAM
                      id: output_format
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - SAM
                            - BAM
                            - CRAM
                          name: output_format
                      inputBinding:
                        position: 1
                        prefix: '--output-fmt'
                        shellQuote: false
                      label: Output format
                      doc: Output file format
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: fast_bam_compression
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-1'
                        shellQuote: false
                      label: Fast BAM compression
                      doc: >-
                        Enable fast BAM compression (implies output in bam
                        format).
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: uncompressed_bam
                      type: boolean?
                      inputBinding:
                        position: 3
                        prefix: '-u'
                        shellQuote: false
                      label: Output uncompressed BAM
                      doc: >-
                        Output uncompressed BAM (implies output BAM format).
                        This option saves time spent on
                        compression/decompression and is thus preferred when the
                        output is piped to another SAMtools command.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: include_header
                      type: boolean?
                      inputBinding:
                        position: 4
                        prefix: '-h'
                        shellQuote: false
                      label: Include the header in the output
                      doc: Include the header in the output.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: output_header_only
                      type: boolean?
                      inputBinding:
                        position: 5
                        prefix: '-H'
                        shellQuote: false
                      label: Output the header only
                      doc: Output the header only.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: collapse_cigar
                      type: boolean?
                      inputBinding:
                        position: 6
                        prefix: '-B'
                        shellQuote: false
                      label: Collapse the backward CIGAR operation
                      doc: Collapse the backward CIGAR operation.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_include
                      type: int?
                      inputBinding:
                        position: 7
                        prefix: '-f'
                        shellQuote: false
                      label: Include reads with all of these flags
                      doc: >-
                        Only output alignments with all bits set in this integer
                        present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_any
                      type: int?
                      inputBinding:
                        position: 8
                        prefix: '-F'
                        shellQuote: false
                      label: Exclude reads with any of these flags
                      doc: >-
                        Do not output alignments with any bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_all
                      type: int?
                      inputBinding:
                        position: 9
                        prefix: '-G'
                        shellQuote: false
                      label: Exclude reads with all of these flags
                      doc: >-
                        Only exclude reads with all of the bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'null'
                      id: read_group
                      type: string?
                      inputBinding:
                        position: 10
                        prefix: '-r'
                        shellQuote: false
                      label: Read group
                      doc: Only output reads in the specified read group.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_mapq
                      type: int?
                      inputBinding:
                        position: 11
                        prefix: '-q'
                        shellQuote: false
                      label: Minimum mapping quality
                      doc: Skip alignments with MAPQ smaller than this value.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'null'
                      id: filter_library
                      type: string?
                      inputBinding:
                        position: 12
                        prefix: '-l'
                        shellQuote: false
                      label: Only include alignments in library
                      doc: Only output alignments in this library.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: min_cigar_operations
                      type: int?
                      inputBinding:
                        position: 13
                        prefix: '-m'
                        shellQuote: false
                      label: Minimum number of CIGAR bases consuming query sequence
                      doc: >-
                        Only output alignments with number of CIGAR bases
                        consuming query sequence  ≥ INT.
                    - 'sbg:category': Config Inputs
                      id: read_tag_to_strip
                      type: 'string[]?'
                      inputBinding:
                        position: 14
                        prefix: ''
                        itemSeparator: ' '
                        shellQuote: false
                        valueFrom: |-
                          ${
                              if (self)
                              {
                                  var cmd = [];
                                  for (var i = 0; i < self.length; i++) 
                                  {
                                      cmd.push('-x', self[i]);
                                      
                                  }
                                  return cmd.join(' ');
                              }
                          }
                      label: Read tags to strip
                      doc: Read tag to exclude from output (repeatable).
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: count_alignments
                      type: boolean?
                      inputBinding:
                        position: 15
                        prefix: '-c'
                        shellQuote: false
                      label: Output only count of matching records
                      doc: >-
                        Instead of outputing the alignments, only count them and
                        output the total number. All filter options, such as -f,
                        -F, and -q, are taken into account.
                    - 'sbg:category': Config Inputs
                      id: input_fmt_option
                      type: string?
                      inputBinding:
                        position: 16
                        prefix: '--input-fmt-option'
                        shellQuote: false
                      label: Input file format option
                      doc: >-
                        Specify a single input file format option in the form of
                        OPTION or OPTION=VALUE.
                    - 'sbg:category': Config Inputs
                      id: output_fmt_option
                      type: string?
                      inputBinding:
                        position: 17
                        prefix: '--output-fmt-option'
                        shellQuote: false
                      label: Output file format option
                      doc: >-
                        Specify a single output file format option in the form
                        of OPTION or OPTION=VALUE.
                    - 'sbg:category': Config Inputs
                      id: subsample_fraction
                      type: float?
                      inputBinding:
                        position: 18
                        prefix: '-s'
                        shellQuote: false
                      label: Subsample fraction
                      doc: >-
                        Output only a proportion of the input alignments. This
                        subsampling acts in the same way on all of the alignment
                        records in the same template or read pair, so it never
                        keeps a read but not its mate. The integer and
                        fractional parts of the INT.FRAC are used separately:
                        the part after the decimal point sets the fraction of
                        templates/pairs to be kept, while the integer part is
                        used as a seed that influences which subset of reads is
                        kept. When subsampling data that has previously been
                        subsampled, be sure to use a different seed value from
                        those used previously; otherwise more reads will be
                        retained than expected.
                    - 'sbg:altPrefix': '-@'
                      'sbg:category': Execution
                      'sbg:toolDefaultValue': '1'
                      id: threads
                      type: int?
                      inputBinding:
                        position: 19
                        prefix: '--threads'
                        shellQuote: false
                        valueFrom: |-
                          ${
                            if((inputs.threads)){
                              return (inputs.threads) - 1
                            }
                            else{
                              return
                            }
                          }
                      label: Number of threads
                      doc: >-
                        Number of threads. SAMtools uses argument --threads/-@
                        to specify number of additional threads. This parameter
                        sets total number of threads (and CPU cores). Command
                        line argument will be reduced by 1 to set number of
                        additional threads.
                    - 'sbg:category': Config Inputs
                      id: omitted_reads_filename
                      type: string?
                      inputBinding:
                        position: 20
                        prefix: '-U'
                        shellQuote: false
                      label: Filename for reads not selected by filters
                      doc: >-
                        Write alignments that are not selected by the various
                        filter options to this file. When this option is used,
                        all alignments (or all alignments intersecting the
                        regions specified) are written to either the output file
                        or this file, but never both.
                    - default: default_output_filename
                      'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': stdout
                      id: output_filename
                      type: string?
                      inputBinding:
                        position: 21
                        prefix: '-o'
                        shellQuote: false
                        valueFrom: |-
                          ${
                            if (inputs.output_filename!="default_output_filename"){
                              return (inputs.output_filename)
                            }
                            input_filename = [].concat(inputs.in_alignments)[0].path.split('/').pop()
                            input_name_base = input_filename.split('.').slice(0,-1).join('.')
                            ext = 'sam'
                            if (inputs.count_alignments){
                              return input_name_base + '.count.txt'
                            }
                            if ((inputs.uncompressed_bam) || (inputs.fast_bam_compression)){
                              ext = 'bam'
                            }
                            if (inputs.output_format){
                              ext = (inputs.output_format).toLowerCase()
                            }
                            if (inputs.output_header_only){
                              ext = 'header.' + ext
                            }
                            if (inputs.subsample_fraction){
                              ext = 'subsample.' + ext
                            }
                            if ((inputs.bed_file) || (inputs.read_group) || (inputs.read_group_list) ||
                                (inputs.filter_mapq) || (inputs.filter_library) || (inputs.min_cigar_operations) ||
                                (inputs.filter_include) || (inputs.filter_exclude_any) || 
                                (inputs.filter_exclude_all) || (inputs.regions_array)){
                              ext = 'filtered.' + ext
                            }
                              
                            return input_name_base + '.' + ext
                          }
                      label: Output filename
                      doc: Define a filename of the output.
                    - 'sbg:category': File Inputs
                      id: bed_file
                      type: File?
                      inputBinding:
                        position: 22
                        prefix: '-L'
                        shellQuote: false
                      label: BED region file
                      doc: Only output alignments overlapping the input BED file.
                      'sbg:fileTypes': BED
                    - 'sbg:category': File Inputs
                      id: read_group_list
                      type: File?
                      inputBinding:
                        position: 23
                        prefix: '-R'
                        shellQuote: false
                      label: Read group list
                      doc: Output alignments in read groups listed in this file.
                      'sbg:fileTypes': TXT
                    - 'sbg:altPrefix': '-T'
                      'sbg:category': File Inputs
                      id: in_reference
                      type: File?
                      inputBinding:
                        position: 24
                        prefix: '--reference'
                        shellQuote: false
                      label: Reference file
                      doc: >-
                        A FASTA format reference file, optionally compressed by
                        bgzip and ideally indexed by SAMtools Faidx. If an index
                        is not present, one will be generated for you. This file
                        is used for compression/decompression of CRAM files.
                        Please provide reference file when using CRAM
                        input/output file.
                      'sbg:fileTypes': 'FASTA, FA, FASTA.GZ, FA.GZ, GZ'
                    - 'sbg:category': File Inputs
                      id: reference_file_list
                      type: File?
                      inputBinding:
                        position: 25
                        prefix: '-t'
                        shellQuote: false
                      label: List of reference names and lengths
                      doc: >-
                        A tab-delimited file. Each line must contain the
                        reference name in the first column and the length of the
                        reference in the second column, with one line for each
                        distinct reference. Any additional fields beyond the
                        second column are ignored. This file also defines the
                        order of the reference sequences in sorting. If you run
                        SAMtools Faidx on reference FASTA file (<ref.fa>), the
                        resulting index file <ref.fa>.fai can be used as this
                        file.
                      'sbg:fileTypes': 'FAI, TSV, TXT'
                    - 'sbg:category': File Inputs
                      id: in_alignments
                      type: File
                      inputBinding:
                        position: 99
                        shellQuote: false
                      label: Input BAM/SAM/CRAM file
                      doc: Input BAM/SAM/CRAM file.
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - 'sbg:category': Config Inputs
                      id: regions_array
                      type: 'string[]?'
                      inputBinding:
                        position: 100
                        shellQuote: false
                      label: Regions array
                      doc: >-
                        With no options or regions specified, prints all
                        alignments in the specified input alignment file (in
                        SAM, BAM, or CRAM format) to output file in specified
                        format. Use of region specifications requires a
                        coordinate-sorted and indexed input file (in BAM or CRAM
                        format). Regions can be specified as:
                        RNAME[:STARTPOS[-ENDPOS]] and all position coordinates
                        are 1-based.  Important note: when multiple regions are
                        given, some alignments may be output multiple times if
                        they overlap more than one of the specified regions.
                        Examples of region specifications:  chr1 - Output all
                        alignments mapped to the reference sequence named `chr1'
                        (i.e. @SQ SN:chr1);  chr2:1000000 - The region on chr2
                        beginning at base position 1,000,000 and ending at the
                        end of the chromosome;  chr3:1000-2000 - The 1001bp
                        region on chr3 beginning at base position 1,000 and
                        ending at base position 2,000 (including both end
                        positions);  '*' - Output the unmapped reads at the end
                        of the file (this does not include any unmapped reads
                        placed on a reference sequence alongside their mapped
                        mates.);  . - Output all alignments (mostly unnecessary
                        as not specifying a region at all has the same effect).
                    - 'sbg:category': Config inputs
                      'sbg:toolDefaultValue': 'False'
                      id: multi_region_iterator
                      type: boolean?
                      inputBinding:
                        position: 22
                        prefix: '-M'
                        shellQuote: false
                      label: Use the multi-region iterator
                      doc: >-
                        Use the multi-region iterator on the union of the BED
                        file and command-line region arguments.
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1500'
                      id: mem_per_job
                      type: int?
                      label: Memory per job
                      doc: Memory per job in MB.
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1'
                      id: cpu_per_job
                      type: int?
                      label: CPU per job
                      doc: Number of CPUs per job.
                  outputs:
                    - id: out_alignments
                      doc: The output file.
                      label: 'Output BAM, SAM, or CRAM file'
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                            if ((inputs.output_filename!="default_output_filename")){
                              return (inputs.output_filename)
                            }
                            input_filename = [].concat((inputs.in_alignments))[0].path.split('/').pop()
                            input_name_base = input_filename.split('.').slice(0,-1). join('.')
                            ext = 'sam'
                            if ((inputs.count_alignments)){
                              return 
                            }
                            if ((inputs.uncompressed_bam) || (inputs.fast_bam_compression)){
                              ext = 'bam'
                            }
                            if ((inputs.output_format)){
                              ext = (inputs.output_format).toLowerCase()
                            }
                            if ((inputs.output_header_only)){
                              ext = 'header.' + ext
                            }
                            if ((inputs.subsample_fraction)){
                              ext = 'subsample.' + ext
                            }
                            if ((inputs.bed_file) || (inputs.read_group) || (inputs.read_group_list) ||
                                (inputs.filter_mapq) || (inputs.filter_library) || (inputs.min_cigar_operations) ||
                                (inputs.filter_include) || (inputs.filter_exclude_any) || 
                                (inputs.filter_exclude_all) || (inputs.regions_array)){
                              ext = 'filtered.' + ext
                            }
                              
                            return input_name_base + '.' + ext
                          }
                        outputEval: '$(inheritMetadata(self, inputs.in_alignments))'
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - id: reads_not_selected_by_filters
                      doc: File containing reads that are not selected by filters.
                      label: Reads not selected by filters
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                            if ((inputs.omitted_reads_filename)){
                              return (inputs.omitted_reads_filename)
                            }
                          }
                        outputEval: '$(inheritMetadata(self, inputs.in_alignments))'
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - id: alignement_count
                      doc: File containing number of alignments.
                      label: Alignment count
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                            input_filename = [].concat((inputs.in_alignments))[0].path.split('/').pop()
                            input_name_base = input_filename.split('.').slice(0,-1). join('.')
                            return input_name_base + '.count.txt'
                          }
                        outputEval: '$(inheritMetadata(self, inputs.in_alignments))'
                      'sbg:fileTypes': TXT
                  doc: >-
                    **SAMtools View** tool prints all alignments from a SAM,
                    BAM, or CRAM file to an output file in SAM format
                    (headerless). You may specify one or more space-separated
                    region specifications to restrict output to only those
                    alignments which overlap the specified region(s). Use of
                    region specifications requires a coordinate-sorted and
                    indexed input file (in BAM or CRAM format) [1].


                    *A list of **all inputs and parameters** with corresponding
                    descriptions can be found at the bottom of the page.*


                    ####Regions


                    Regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] and
                    all position coordinates are 1-based. 


                    **Important note:** when multiple regions are given, some
                    alignments may be output multiple times if they overlap more
                    than one of the specified regions.


                    Examples of region specifications:


                    - **chr1**  - Output all alignments mapped to the reference
                    sequence named `chr1' (i.e. @SQ SN:chr1).


                    - **chr2:1000000** - The region on chr2 beginning at base
                    position 1,000,000 and ending at the end of the chromosome.


                    - **chr3:1000-2000** - The 1001bp region on chr3 beginning
                    at base position 1,000 and ending at base position 2,000
                    (including both end positions).


                    - **'\*'** - Output the unmapped reads at the end of the
                    file. (This does not include any unmapped reads placed on a
                    reference sequence alongside their mapped mates.)


                    - **.** - Output all alignments. (Mostly unnecessary as not
                    specifying a region at all has the same effect.) [1]


                    ###Common Use Cases


                    This tool can be used for: 


                    - Filtering BAM/SAM/CRAM files - options set by the
                    following parameters and input files: **Include reads with
                    all of these flags** (`-f`), **Exclude reads with any of
                    these flags** (`-F`), **Exclude reads with all of these
                    flags** (`-G`), **Read group** (`-r`), **Minimum mapping
                    quality** (`-q`), **Only include alignments in library**
                    (`-l`), **Minimum number of CIGAR bases consuming query
                    sequence** (`-m`), **Subsample fraction** (`-s`), **Read
                    group list** (`-R`), **BED region file** (`-L`)

                    - Format conversion between SAM/BAM/CRAM formats - set by
                    the following parameters: **Output format**
                    (`--output-fmt/-O`), **Fast bam compression** (`-1`),
                    **Output uncompressed BAM** (`-u`)

                    - Modification of the data which is contained in each
                    alignment - set by the following parameters: **Collapse the
                    backward CIGAR operation** (`-B`), **Read tags to strip**
                    (`-x`)

                    - Counting number of alignments in SAM/BAM/CRAM file - set
                    by parameter **Output only count of matching records**
                    (`-c`)


                    ###Changes Introduced by Seven Bridges


                    - Parameters **Output BAM** (`-b`) and **Output CRAM**
                    (`-C`) were excluded from the wrapper since they are
                    redundant with parameter **Output format**
                    (`--output-fmt/-O`).

                    - Parameter **Input format** (`-S`) was excluded from
                    wrapper since it is ignored by the tool (input format is
                    auto-detected).

                    - Input file **Index file** was added to the wrapper to
                    enable operations that require an index file for BAM/CRAM
                    files.

                    - Parameter **Number of threads** (`--threads/-@`) specifies
                    the total number of threads instead of additional threads.
                    Command line argument (`--threads/-@`) will be reduced by 1
                    to set the number of additional threads.


                    ###Common Issues and Important Notes


                    - When multiple regions are given, some alignments may be
                    output multiple times if they overlap more than one of the
                    specified regions [1].

                    - Use of region specifications requires a coordinate-sorted
                    and indexed input file (in BAM or CRAM format) [1].

                    - Option **Output uncompressed BAM** (`-u`) saves time spent
                    on compression/decompression and is thus preferred when the
                    output is piped to another SAMtools command [1].


                    ###Performance Benchmarking


                    Multithreading can be enabled by setting parameter **Number
                    of threads** (`--threads/-@`). In the following table you
                    can find estimates of **SAMtools View** running time and
                    cost. 


                    *Cost can be significantly reduced by using **spot
                    instances**. Visit the [Knowledge
                    Center](https://docs.sevenbridges.com/docs/about-spot-instances)
                    for more details.*  


                    | Input type | Input size | # of reads | Read length |
                    Output format | # of threads | Duration | Cost | Instance
                    (AWS)|

                    |---------------|--------------|-----------------|---------------|------------------|-------------------|-----------------|-------------|--------|-------------|

                    | BAM | 5.26 GB | 71.5M | 76 | BAM | 1 | 13min. | \$0.12 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | 161.2M | 101 | BAM | 1 | 33min. | \$0.30
                    | c4.2xlarge |

                    | BAM | 18.36 GB | 179M | 76 | BAM | 1 | 60min. | \$0.54 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | 845.6M | 150 | BAM | 1 | 3h 25min. |
                    \$1.84 | c4.2xlarge |

                    | BAM | 5.26 GB | 71.5M | 76 | BAM | 8 | 5min. | \$0.04 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | 161.2M | 101 | BAM | 8 | 11min. | \$0.10
                    | c4.2xlarge |

                    | BAM | 18.36 GB | 179M | 76 | BAM | 8 | 19min. | \$0.17 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | 845.6M | 150 | BAM | 8 | 61min. | \$0.55
                    | c4.2xlarge |

                    | BAM | 5.26 GB | 71.5M | 76 | SAM | 8 | 14min. | \$0.13 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | 161.2M | 101 | SAM | 8 | 23min. | \$0.21
                    | c4.2xlarge |

                    | BAM | 18.36 GB | 179M | 76 | SAM | 8 | 35min. | \$0.31 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | 845.6M | 150 | SAM | 8 | 2h 29min. |
                    \$1.34 | c4.2xlarge |


                    ###References


                    [1] [SAMtools
                    documentation](http://www.htslib.org/doc/samtools-1.9.html)
                  label: Samtools View CWL1.0
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: |-
                        ${
                          if (inputs.mem_per_job) {
                              return inputs.mem_per_job
                          }    
                          else {
                          mem_offset = 1000
                          if((inputs.in_reference)){
                            mem_offset = mem_offset + 3000
                          }
                          if((inputs.threads)){
                            threads = (inputs.threads)
                          }
                          else{
                            threads = 1
                          }
                          return mem_offset + threads * 500
                          }
                        }
                      coresMin: |-
                        ${
                          if (inputs.cpu_per_job) {
                              return inputs.cpu_per_job
                          }
                          else {
                          if((inputs.threads)){
                            return (inputs.threads)
                          }
                          else{
                            return 1
                          }
                          }
                        }
                    - class: DockerRequirement
                      dockerPull: 'images.sbgenomics.com/jrandjelovic/samtools-1-9:1'
                    - class: InitialWorkDirRequirement
                      listing:
                        - $(inputs.in_reference)
                        - $(inputs.reference_file_list)
                        - $(inputs.in_index)
                        - $(inputs.in_alignments)
                    - class: InlineJavascriptRequirement
                      expressionLib:
                        - |-

                          var setMetadata = function(file, metadata) {
                              if (!('metadata' in file))
                                  file['metadata'] = metadata;
                              else {
                                  for (var key in metadata) {
                                      file['metadata'][key] = metadata[key];
                                  }
                              }
                              return file
                          };

                          var inheritMetadata = function(o1, o2) {
                              var commonMetadata = {};
                              if (!Array.isArray(o2)) {
                                  o2 = [o2]
                              }
                              for (var i = 0; i < o2.length; i++) {
                                  var example = o2[i]['metadata'];
                                  for (var key in example) {
                                      if (i == 0)
                                          commonMetadata[key] = example[key];
                                      else {
                                          if (!(commonMetadata[key] == example[key])) {
                                              delete commonMetadata[key]
                                          }
                                      }
                                  }
                              }
                              if (!Array.isArray(o1)) {
                                  o1 = setMetadata(o1, commonMetadata)
                              } else {
                                  for (var i = 0; i < o1.length; i++) {
                                      o1[i] = setMetadata(o1[i], commonMetadata)
                                  }
                              }
                              return o1;
                          };
                  'sbg:appVersion':
                    - v1.0
                  'sbg:categories':
                    - Utilities
                    - BAM Processing
                    - CWL1.0
                  'sbg:content_hash': >-
                    aa82916613444b2d378befd3fc8666677b6b22c3fb84f9dd8985aa73494c63afa
                  'sbg:contributors':
                    - admin
                  'sbg:createdBy': admin
                  'sbg:createdOn': 1576244128
                  'sbg:id': admin/sbg-public-data/samtools-view-1-9-cwl1-0/6
                  'sbg:image_url': null
                  'sbg:latestRevision': 6
                  'sbg:license': MIT License
                  'sbg:links':
                    - id: 'http://www.htslib.org/'
                      label: Homepage
                    - id: 'https://github.com/samtools/samtools'
                      label: Source Code
                    - id: 'https://github.com/samtools/samtools/wiki'
                      label: Wiki
                    - id: >-
                        https://sourceforge.net/projects/samtools/files/samtools/
                      label: Download
                    - id: 'http://www.ncbi.nlm.nih.gov/pubmed/19505943'
                      label: Publication
                    - id: 'http://www.htslib.org/doc/samtools-1.9.html'
                      label: Documentation
                  'sbg:modifiedBy': admin
                  'sbg:modifiedOn': 1578576084
                  'sbg:project': admin/sbg-public-data
                  'sbg:projectName': SBG Public Data
                  'sbg:publisher': sbg
                  'sbg:revision': 6
                  'sbg:revisionNotes': Added file requirements for in_index and in_alignments
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244128
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 1
                      'sbg:revisionNotes': Final version
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 2
                      'sbg:revisionNotes': 'Edited description, tag, default values.'
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 3
                      'sbg:revisionNotes': mem_per_job default value set
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 4
                      'sbg:revisionNotes': Description edited - references put before full stop
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1576244287
                      'sbg:revision': 5
                      'sbg:revisionNotes': Categories edited
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1578576084
                      'sbg:revision': 6
                      'sbg:revisionNotes': Added file requirements for in_index and in_alignments
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': >-
                    Heng Li (Sanger Institute), Bob Handsaker (Broad Institute),
                    Jue Ruan (Beijing Genome Institute), Colin Hercus, Petr
                    Danecek
                  'sbg:toolkit': samtools
                  'sbg:toolkitVersion': '1.9'
                  'sbg:validationErrors': []
                label: Samtools View CWL1.0
                'sbg:x': 547.8245239257812
                'sbg:y': 120.75
              - id: samtools_fastq
                in:
                  - id: in_alignments
                    source: samtools_view_1_9_cwl1_0/out_alignments
                out:
                  - id: output_pe_1
                  - id: output_pe_2
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: admin_sbg_public_data_samtools_fastq_1_9_cwl1_0_4
                  baseCommand: []
                  inputs:
                    - 'sbg:category': File Inputs
                      id: in_alignments
                      type: File
                      inputBinding:
                        position: 101
                        shellQuote: false
                      label: BAM/SAM/CRAM file
                      doc: Input BAM/SAM/CRAM file.
                      'sbg:fileTypes': 'BAM, SAM, CRAM'
                    - default: 0
                      'sbg:category': Config Inputs
                      id: read_0_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-0'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.read_0_filename = null
                              };


                              if (inputs.single_fastq_file == true) {
                                  return
                              } else {
                                  if (inputs.read_0_filename) {
                                      return inputs.read_0_filename
                                  } else {
                                      var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                      var ext = bamname.split('.').pop().toLowerCase()
                                      if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                          return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.pe_0.fastq'
                                      } else {
                                          return bamname + '.pe_0.fastq'
                                      }
                                  }
                              }
                          }
                      label: Unflagged reads filename
                      doc: >-
                        Write reads with both or neither of the BAM_READ1 and
                        BAM_READ2 flags set to this file.
                    - default: 0
                      'sbg:category': Config Inputs
                      id: read_1_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-1'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.read_1_filename = null
                              };


                              if (inputs.single_fastq_file == true) {
                                  return
                              } else {
                                  if (inputs.read_1_filename) {
                                      return inputs.read_1_filename
                                  } else {
                                      var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                      var ext = bamname.split('.').pop().toLowerCase()
                                      if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                          return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.pe_1.fastq'
                                      } else {
                                          return bamname + '.pe_1.fastq'
                                      }
                                  }
                              }
                          }
                      label: Filename for BAM_READ1 flaged reads
                      doc: Write reads with the BAM_READ1 flag set to this file.
                    - default: 0
                      'sbg:category': Config Inputs
                      id: read_2_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-2'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.read_2_filename = null
                              };


                              if (inputs.single_fastq_file == true) {
                                  return
                              } else {
                                  if (inputs.read_2_filename) {
                                      return inputs.read_2_filename
                                  } else {
                                      var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                      var ext = bamname.split('.').pop().toLowerCase()
                                      if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                          return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.pe_2.fastq'
                                      } else {
                                          return bamname + '.pe_2.fastq'
                                      }
                                  }
                              }
                          }
                      label: Filename for BAM_READ2 flaged reads
                      doc: Write reads with the BAM_READ2 flag set to this file.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_include
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-f'
                        shellQuote: false
                      label: Include reads with all of these flags
                      doc: >-
                        Only output alignments with all bits set in this integer
                        present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_any
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-F'
                        shellQuote: false
                      label: Exclude reads with any of these flags
                      doc: >-
                        Do not output alignments with any bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '0'
                      id: filter_exclude_all
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-G'
                        shellQuote: false
                      label: Exclude reads with all of these flags
                      doc: >-
                        Only exclude reads with all of the bits set in this
                        integer present in the FLAG field.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: do_not_append_read_number_to_name
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-n'
                        shellQuote: false
                      label: Don't append /1 and /2 to the read name
                      doc: >-
                        By default, either '/1' or '/2' is added to the end of
                        read names where the corresponding BAM_READ1 or
                        BAM_READ2 flag is set. Setting this parameter to True
                        causes read names to be left as they are.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: always_append_read_number_to_name
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-N'
                        shellQuote: false
                      label: Always append /1 and /2 to the read name
                      doc: >-
                        Always add either '/1' or '/2' to the end of read names
                        even when put into different files.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: output_quality_in_OQ_tag
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-O'
                        shellQuote: false
                      label: Output quality in the OQ tag if present
                      doc: >-
                        Use quality values from OQ tags in preference to
                        standard quality string if available.
                    - 'sbg:category': Config Inputs
                      id: singleton_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-s'
                        shellQuote: false
                      label: Singleton reads filename
                      doc: Write singleton reads in FASTQ format to this file.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: copy_tags
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-t'
                        shellQuote: false
                      label: 'Copy RG, BC and QT tags to the FASTQ header line'
                      doc: >-
                        Copy RG, BC and QT tags to the FASTQ header line, if
                        they exist.
                    - default: 0
                      'sbg:category': Config Inputs
                      id: copy_taglist
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '-T'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.copy_taglist = null
                              };


                              if (inputs.copy_taglist) {
                                  return inputs.copy_taglist.replace(/ /g, '')
                              } else {
                                  return
                              }
                          }
                      label: Taglist to copy to the FASTQ header line
                      doc: >-
                        Specify a comma-separated list of tags to copy to the
                        FASTQ header line, if they exist.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': '1'
                      id: default_quality_score
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-v'
                        shellQuote: false
                      label: Default quality score if not given in file
                      doc: Default quality score if not given in file.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: add_Illumina_casava_format_entry
                      type: boolean?
                      inputBinding:
                        position: 2
                        prefix: '-i'
                        shellQuote: false
                      label: Add Illumina Casava 1.8 format entry to header
                      doc: >-
                        Add Illumina Casava 1.8 format entry to header (eg
                        1:N:0:ATCACG).
                    - 'sbg:category': Config Inputs
                      id: compression_level
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '-c'
                        shellQuote: false
                      label: 'Compression level [0..9]'
                      doc: >-
                        Compression level [0..9] to use when creating GZ or BGZF
                        FASTQ files.
                    - 'sbg:category': Config Inputs
                      id: first_index_read_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--i1'
                        shellQuote: false
                      label: First index reads filename
                      doc: >-
                        Specify filename to which the first index reads will be
                        written.
                    - 'sbg:category': Config Inputs
                      id: second_index_read_filename
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--i2'
                        shellQuote: false
                      label: Second index reads filename
                      doc: >-
                        Specify filename to which the second index reads will be
                        written.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': BC
                      id: barcode_tag
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--barcode-tag'
                        shellQuote: false
                      label: Barcode tag
                      doc: 'Aux tag to find index reads in [default: BC].'
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': QT
                      id: quality_tag
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--quality-tag'
                        shellQuote: false
                      label: Quality tag
                      doc: 'Aux tag to find index quality in [default: QT].'
                    - 'sbg:category': Config Inputs
                      id: index_format
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--index-format'
                        shellQuote: false
                      label: Index format for parsing barcode and quality tags
                      doc: >-
                        String to describe how to parse the barcode and quality
                        tags. For example:  i14i8 - the first 14 characters are
                        index 1, the next 8 characters are index 2; n8i14 ignore
                        the first 8 characters, and use the next 14 characters
                        for index 1; If the tag contains a separator, then the
                        numeric part can be replaced with '*' to mean 'read
                        until the separator or end of tag', for example:  n*i* -
                        ignore the left part of the tag until the separator,
                        then use the second part.
                    - 'sbg:category': Config Inputs
                      id: single_fastq_filename
                      type: string?
                      label: Single FASTQ filename
                      doc: >-
                        Filename of an output FASTQ file if only one file is
                        required.
                    - 'sbg:category': Config Inputs
                      'sbg:toolDefaultValue': 'False'
                      id: single_fastq_file
                      type: boolean?
                      label: Single FASTQ file
                      doc: True if only one FASTQ file is required as output.
                    - 'sbg:category': Config Inputs
                      id: input_fmt_option
                      type: string?
                      inputBinding:
                        position: 2
                        prefix: '--input-fmt-option'
                        shellQuote: false
                      label: Input file format option
                      doc: >-
                        Specify a single input file format option in the form of
                        OPTION or OPTION=VALUE.
                    - default: 0
                      'sbg:altPrefix': '-@'
                      'sbg:category': Execution
                      'sbg:toolDefaultValue': '1'
                      id: threads
                      type: int?
                      inputBinding:
                        position: 2
                        prefix: '--threads'
                        shellQuote: false
                        valueFrom: |-
                          ${
                              var self
                              if (self == 0) {
                                  self = null;
                                  inputs.threads = null
                              };


                              if (inputs.threads) {
                                  return inputs.threads - 1
                              } else
                                  return
                          }
                      label: Number of threads
                      doc: >-
                        Number of threads. SAMtools uses argument --threads/-@
                        to specify number of additional threads. This parameter
                        sets total number of threads (and CPU cores). Command
                        line argument will be reduced by 1 to set number of
                        additional threads.
                    - 'sbg:category': File Inputs
                      id: in_reference
                      type: File?
                      inputBinding:
                        position: 2
                        prefix: '--reference'
                        shellQuote: false
                      label: Reference file
                      doc: >-
                        Reference file. This file is used for
                        compression/decompression of CRAM files. Please provide
                        reference file when using CRAM input/output file.
                      'sbg:fileTypes': 'FASTA, FASTA.GZ, FASTA.BGZF, GZ, BGZF'
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1500'
                      id: mem_per_job
                      type: int?
                      label: Memory per job
                      doc: Memory per job in MB.
                    - 'sbg:category': Platform Options
                      'sbg:toolDefaultValue': '1'
                      id: cpu_per_job
                      type: int?
                      label: CPU per job
                      doc: Number of CPUs per job.
                  outputs:
                    - id: output_pe_1
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                              if (inputs.read_1_filename) {
                                  var read1 = inputs.read_1_filename
                              } else {
                                  read1 = '*pe_1.fastq*'
                              }
                              return read1
                          }
                    - id: output_pe_2
                      type: File?
                      outputBinding:
                        glob: |-
                          ${
                              if (inputs.read_2_filename) {
                                  var read2 = inputs.read_2_filename
                              } else {
                                  read2 = '*pe_2.fastq*'
                              }
                              return read2
                          }
                  doc: >-
                    **SAMtools FASTQ** tool converts a BAM or CRAM into FASTQ
                    format. The FASTQ files will be automatically compressed if
                    the filenames have a .gz or .bgzf extension [1].


                    **SAMtools FASTQ does not perform any preprocessing of the
                    input.** If you want to use output FASTQ files with
                    alignment tools, please make sure that the **BAM/SAM/CRAM
                    file** is sorted by read name (or collated). Otherwise,
                    alignment tools will fail. Supplementary and secondary
                    alignments are filtered out regardless of filtering
                    parameters. 


                    Parameter **Index format for parsing barcode and quality
                    tags** (`--index-format`) is a string used to describe how
                    to parse the barcode and quality tags. For example:  


                    - i14i8 - the first 14 characters are index 1, the next 8
                    characters are index 2

                    - n8i14 - ignore the first 8 characters, and use the next 14
                    characters for index 1  

                    - If the tag contains a separator, then the numeric part can
                    be replaced with '\*' to mean 'read until the separator or
                    end of tag', for example:  

                    n\*i\* - ignore the left part of the tag until the
                    separator, then use the second part [1].


                    *A list of **all inputs and parameters** with corresponding
                    descriptions can be found at the bottom of the page.*


                    ###Common Use Cases


                    - Default use case provides two FASTQ files as outputs
                    (**Paired-end FASTQ files**). If no parameter is set, this
                    will be the case. The tool will output file with reads that
                    are not properly flagged (**Unflagged reads FASTQ file**)
                    only in case this file is not empty. **SAMtools FASTQ does
                    not perform any preprocessing of the input.** If you want to
                    use output FASTQ files with alignment tools, please make
                    sure that the **BAM/SAM/CRAM file** is sorted by read name
                    (or collated). Otherwise, alignment tools will fail.
                    Supplementary and secondary alignments are filtered out
                    regardless of filtering parameters.   

                    - If a single FASTQ file (with both paired end reads) is
                    required, it should be specified by setting the boolean
                    parameter **Single FASTQ file** to True. 


                    ###Changes Introduced by Seven Bridges


                    - Parameter **Single FASTQ file** was added to parameter
                    list to provide the option for outputting a single FASTQ
                    file with all the reads.

                    - Parameter **Single FASTQ filename** was added to parameter
                    list to specify the filename when **Single FASTQ file** is
                    set to True. This parameter is not mandatory. If **Single
                    FASTQ file** is set to True and **Single FASTQ filename** is
                    not specified, default value will be used (*input.fastq* for
                    input **BAM/SAM/CRAM file** named *input.bam*).

                    - Parameter **Number of threads** (`--threads/-@`) specifies
                    the total number of threads instead of additional threads.
                    Command line argument (`--threads/-@`) will be reduced by 1
                    to set number of additional threads.


                    ###Common Issues and Important Notes


                    - If parameters **First index reads filename** (`--i1`)
                    and/or **Second index reads filename** (`--i2`) are
                    specified, **Index format for parsing barcode and quality
                    tags** (`--index-format`) should be specified too and this
                    format should match the number of index files required. 

                    - When specifying output filenames, complete names should be
                    used (including extensions). If the extension is .fastq.gz
                    or .fastq.bgzf, the output will be compressed. The tool does
                    not validate extensions. If the extension is not valid, the
                    task will not fail.

                    - Parameter **Number of threads** (`--threads/-@`) does not
                    decrease running time significantly (not more than 10% with
                    8 threads on a c4.2xlarge instance (AWS)).

                    - **SAMtools FASTQ** does not perform any preprocessing of
                    the input. If you want to use output FASTQ files with
                    alignment tools, please make sure that **BAM/SAM/CRAM file**
                    is sorted by read name (or collated). Otherwise, alignment
                    tools will fail. Supplementary and secondary alignments are
                    filtered out regardless of filtering parameters. 


                    ###Performance Benchmarking


                    In the following table you can find estimates of **SAMtools
                    FASTQ** running time and cost. Adding additional threads
                    does not decrease running time significantly (not more than
                    10% with 8 threads on a c4.2xlarge instance (AWS)).


                    *Cost can be significantly reduced by using **spot
                    instances**. Visit the [Knowledge
                    Center](https://docs.sevenbridges.com/docs/about-spot-instances)
                    for more details.*  


                    | Input type | Input size | Paired-end | # of reads | Read
                    length |  # of threads | Duration | Cost | Instance (AWS) |

                    |---------------|--------------|-----------------|---------------|------------------|-----------------|-------------|--------|-------------|

                    | BAM | 5.26 GB | Yes | 71.5M | 76 | 1 | 7min. | \$0.06 |
                    c4.2xlarge |

                    | BAM | 11.86 GB | Yes | 161.2M | 101 | 1 | 18min. | \$0.16
                    | c4.2xlarge |

                    | BAM | 18.36 GB | Yes | 179M | 76 | 1 | 21min. | \$0.19 |
                    c4.2xlarge |

                    | BAM | 58.61 GB | Yes | 845.6M | 150 | 1 | 2h 7min. |
                    \$1.14 | c4.2xlarge |


                    ###References


                    [1] [SAMtools
                    documentation](http://www.htslib.org/doc/samtools-1.9.html)
                  label: Samtools FASTQ CWL1.0
                  arguments:
                    - position: 0
                      prefix: ''
                      shellQuote: false
                      valueFrom: /opt/samtools-1.9/samtools
                    - position: 1
                      shellQuote: false
                      valueFrom: fastq
                  requirements:
                    - class: ShellCommandRequirement
                    - class: ResourceRequirement
                      ramMin: |-
                        ${
                            if (inputs.mem_per_job) {
                                return inputs.mem_per_job
                            }
                            else {
                            if (inputs.threads) {
                                var threads = inputs.threads
                            } else {
                                threads = 1
                            }
                            return 1000 + 500 * threads
                            }
                        }
                      coresMin: |-
                        ${
                            if (inputs.cpu_per_job) {
                                return inputs.cpu_per_job
                            }
                            else {
                            if (inputs.threads) {
                                return inputs.threads
                            } else {
                                return 1
                            }
                            }
                        }
                    - class: DockerRequirement
                      dockerPull: 'images.sbgenomics.com/jrandjelovic/samtools-1-9:1'
                    - class: InitialWorkDirRequirement
                      listing:
                        - $(inputs.in_reference)
                    - class: InlineJavascriptRequirement
                      expressionLib:
                        - |
                          var setMetadata = function(file, metadata) {
                              if (!('metadata' in file)) {
                                  file['metadata'] = {}
                              }
                              for (var key in metadata) {
                                  file['metadata'][key] = metadata[key];
                              }
                              return file
                          };

                          var inheritMetadata = function(o1, o2) {
                              var commonMetadata = {};
                              if (!o2) {
                                  return o1;
                              };
                              if (!Array.isArray(o2)) {
                                  o2 = [o2]
                              }
                              for (var i = 0; i < o2.length; i++) {
                                  var example = o2[i]['metadata'];
                                  for (var key in example) {
                                      if (i == 0)
                                          commonMetadata[key] = example[key];
                                      else {
                                          if (!(commonMetadata[key] == example[key])) {
                                              delete commonMetadata[key]
                                          }
                                      }
                                  }
                                  for (var key in commonMetadata) {
                                      if (!(key in example)) {
                                          delete commonMetadata[key]
                                      }
                                  }
                              }
                              if (!Array.isArray(o1)) {
                                  o1 = setMetadata(o1, commonMetadata)
                              } else {
                                  for (var i = 0; i < o1.length; i++) {
                                      o1[i] = setMetadata(o1[i], commonMetadata)
                                  }
                              }
                              return o1;
                          };
                  stdout: |-
                    ${
                        if (inputs.single_fastq_file == true) {
                            if (inputs.single_fastq_filename) {
                                return inputs.single_fastq_filename
                            } else {
                                var bamname = [].concat(inputs.in_alignments)[0].path.split('/')[[].concat(inputs.in_alignments)[0].path.split('/').length - 1]
                                var ext = bamname.split('.').pop().toLowerCase()
                                if ((ext == 'bam') || (ext == 'cram') || (ext == 'sam')) {
                                    return bamname.split('.').slice(0, bamname.split('.').length - 1).join('.') + '.fastq'
                                } else {
                                    return bamname + '.fastq'
                                }
                            }
                        } else {
                            return
                        }
                    }
                  'sbg:appVersion':
                    - v1.0
                  'sbg:categories':
                    - Utilities
                    - BAM Processing
                    - CWL1.0
                  'sbg:cmdPreview': /opt/samtools-1.9/samtools fastq  /path/to/input.bam
                  'sbg:content_hash': >-
                    a7f356d044f84a927ee18a8c03a269f0d6cad2ff36f6c6f10c57d4ed5df3e5e04
                  'sbg:contributors':
                    - admin
                    - Rachel Bowen-James <rbowen-james@ccia.org.au>
                  'sbg:createdBy': admin
                  'sbg:image_url': null
                  'sbg:latestRevision': 4
                  'sbg:license': MIT License
                  'sbg:links':
                    - id: 'http://www.htslib.org/'
                      label: Homepage
                    - id: 'https://github.com/samtools/samtools'
                      label: Source Code
                    - id: 'https://github.com/samtools/samtools/wiki'
                      label: Wiki
                    - id: >-
                        https://sourceforge.net/projects/samtools/files/samtools/
                      label: Download
                    - id: 'http://www.ncbi.nlm.nih.gov/pubmed/19505943'
                      label: Publication
                    - id: 'http://www.htslib.org/doc/samtools-1.9.html'
                      label: Documentation
                  'sbg:modifiedBy': admin
                  'sbg:modifiedOn': 1576244287
                  'sbg:project': admin/sbg-public-data
                  'sbg:projectName': SBG Public Data
                  'sbg:publisher': sbg
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': >-
                    Heng Li (Sanger Institute), Bob Handsaker (Broad Institute),
                    Jue Ruan (Beijing Genome Institute), Colin Hercus, Petr
                    Danecek
                  'sbg:toolkit': SAMtools
                  'sbg:toolkitVersion': '1.9'
                  'sbg:validationErrors': []
                  'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
                label: Samtools FASTQ CWL1.0
                'sbg:x': 977.3401489257812
                'sbg:y': 113.75
              - id: hla_hd
                in:
                  - id: threads
                    default: 2
                  - id: minimum_read_length
                    default: 0
                    source: minimum_read_length
                  - id: fastq_reads1
                    source: samtools_fastq/output_pe_1
                  - id: fastq_reads2
                    source: samtools_fastq/output_pe_2
                  - id: sample_id
                    source: sample_id
                  - id: output_dir_name
                    source: output_dir_name
                out:
                  - id: hlahd_results
                  - id: hlahd_final_results
                run:
                  class: CommandLineTool
                  cwlVersion: v1.0
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: mwonge/mwtest/hla-hd/30
                  baseCommand: []
                  inputs:
                    - id: threads
                      type: int
                      inputBinding:
                        position: 1
                        prefix: '-t'
                        shellQuote: false
                      label: Threads
                      doc: Number of cores used to execute the program.
                    - 'sbg:toolDefaultValue': '100'
                      id: minimum_read_length
                      type: int
                      inputBinding:
                        position: 1
                        prefix: '-m'
                        shellQuote: false
                      label: Minimum read length
                      doc: >-
                        A read whose length is shorter than this parameter is
                        ignored.
                    - id: fastq_reads1
                      type: File
                      inputBinding:
                        position: 2
                        separate: false
                        shellQuote: false
                      label: FASTQ Reads 1
                      doc: Paired-end reads 1 in FASTQ format.
                      'sbg:fileTypes': FASTQ
                    - id: fastq_reads2
                      type: File
                      inputBinding:
                        position: 2
                        separate: false
                        shellQuote: false
                      label: FASTQ Reads 2
                      doc: Paired-end reads 2 in FASTQ format.
                      'sbg:fileTypes': FASTQ
                    - id: sample_id
                      type: string
                      inputBinding:
                        position: 4
                        separate: false
                        shellQuote: false
                      label: Sample ID
                      doc: >-
                        Sample ID for the input FASTQs. This will be used as the
                        name of the output directory.
                    - id: output_dir_name
                      type: string?
                  outputs:
                    - id: hlahd_results
                      doc: Directory containing results of the HLA-HD run.
                      label: Output directory
                      type: Directory
                      outputBinding:
                        glob: |-
                          ${
                              if (!inputs.output_dir_name) {
                                  return inputs.sample_id + "_hlahd"
                              } else {
                                  return inputs.output_dir_name
                              }
                          }
                    - id: hlahd_final_results
                      type: File
                      outputBinding:
                        glob: |-
                          ${
                              if (!inputs.output_dir_name) {
                                  return inputs.sample_id + "_hlahd/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                              } else {
                                  return inputs.output_dir_name + "/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                              }
                          }
                  doc: >-
                    ## About HLA-HD

                    HLA-HD documentation and release notes can be found
                    [here](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/).

                    HLA-HD (HLA typing from High-quality Dictionary) can
                    accurately determine HLA alleles with 6-digit precision from
                    NGS data (FASTQ format). RNA-Seq data can also be applied.


                    ## About this CWL tool

                    This CWL tool runs HLA-HD version 1.4.0 on paired-end FASTQ
                    files to determine HLA type.


                    ### Inputs and parameters

                    - The input paired-end read files can be from **WGS/WES or
                    RNA-seq**.

                    - The input paired-end read files must be in FASTQ format
                    (**not zipped**).

                    - The default minimum read length is 100, however this is
                    often too strict. Choose a lower threshold to include more
                    reads.


                    ### Output

                    - HLA-HD results are output to a directory named using the
                    input sample id.

                    - The final summary of HLA typing results can be found at
                    the following path:
                    `<output_dir_name>/result/<sample_id>_final.result.txt`.


                    ### Other notes

                    - This tool uses the HLA dictionary created from release
                    3.15.0 of the
                    [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA) database.

                    - This tool by default uses HLA allele frequency data
                    included with the HLA-HD release 1.4.0.
                  label: HLA-HD
                  arguments:
                    - position: 2
                      prefix: '-f'
                      shellQuote: false
                      valueFrom: /app/hlahd.1.4.0/freq_data
                    - position: 3
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: /app/hlahd.1.4.0/HLA_gene.split.txt
                    - position: 3
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: /app/hlahd.1.4.0/dictionary
                    - position: 1
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: hlahd.sh
                    - position: 5
                      prefix: ''
                      valueFrom: ./
                    - position: 6
                      prefix: ''
                      separate: false
                      shellQuote: false
                      valueFrom: |-
                        ${
                            if (!inputs.output_dir_name) {
                                return "&& mv ./" + inputs.sample_id + " ./" + inputs.sample_id + "_hlahd"
                            } else {
                                return "&& mv ./" + inputs.sample_id + " ./" + inputs.output_dir_name
                            }
                        }
                    - position: 0
                      prefix: ''
                      shellQuote: false
                      valueFrom: |-
                        ${
                            if (!inputs.output_dir_name) {
                                return "mkdir ./" + inputs.sample_id + "_hlahd &&"
                            } else {
                                return "mkdir ./" + inputs.output_dir_name + " &&"
                            }
                        }
                  requirements:
                    - class: ShellCommandRequirement
                    - class: DockerRequirement
                      dockerPull: pgc-images.sbgenomics.com/rbowen_james/hla-hd
                    - class: InlineJavascriptRequirement
                  'sbg:appVersion':
                    - v1.0
                  'sbg:categories':
                    - WGS
                    - RNA
                    - HLA Typing
                    - HLA
                    - MHC
                    - WES (WXS)
                  'sbg:content_hash': >-
                    a4adfcaaafaaa6b9ac47bbe7e09595dae02a8f11b9f16fdaa0cfcb21948bf91bd
                  'sbg:contributors':
                    - rbowen_james
                  'sbg:createdBy': rbowen_james
                  'sbg:createdOn': 1622961333
                  'sbg:id': mwonge/mwtest/hla-hd/30
                  'sbg:image_url': null
                  'sbg:latestRevision': 30
                  'sbg:modifiedBy': rbowen_james
                  'sbg:modifiedOn': 1627370832
                  'sbg:project': mwonge/mwtest
                  'sbg:projectName': zcc-cavatica
                  'sbg:publisher': sbg
                  'sbg:revision': 30
                  'sbg:revisionNotes': output final result file
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622961333
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622961864
                      'sbg:revision': 1
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622962039
                      'sbg:revision': 2
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1622962892
                      'sbg:revision': 3
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623026790
                      'sbg:revision': 4
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623027275
                      'sbg:revision': 5
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623043041
                      'sbg:revision': 6
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623043454
                      'sbg:revision': 7
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623044879
                      'sbg:revision': 8
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623046596
                      'sbg:revision': 9
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1623048682
                      'sbg:revision': 10
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624408084
                      'sbg:revision': 11
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624411148
                      'sbg:revision': 12
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624580200
                      'sbg:revision': 13
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624581085
                      'sbg:revision': 14
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624581597
                      'sbg:revision': 15
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624583201
                      'sbg:revision': 16
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624583825
                      'sbg:revision': 17
                      'sbg:revisionNotes': ''
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624586750
                      'sbg:revision': 18
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624588045
                      'sbg:revision': 19
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624593624
                      'sbg:revision': 20
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624606248
                      'sbg:revision': 21
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624606452
                      'sbg:revision': 22
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624608065
                      'sbg:revision': 23
                      'sbg:revisionNotes': 'Fixed output dir issue, added docs.'
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1624779790
                      'sbg:revision': 24
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627028760
                      'sbg:revision': 25
                      'sbg:revisionNotes': Added output dir option
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627028962
                      'sbg:revision': 26
                      'sbg:revisionNotes': make outdir
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627029487
                      'sbg:revision': 27
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627029828
                      'sbg:revision': 28
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627104094
                      'sbg:revision': 29
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': rbowen_james
                      'sbg:modifiedOn': 1627370832
                      'sbg:revision': 30
                      'sbg:revisionNotes': output final result file
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': Shuji Kawaguchi <shuji@genome.med.kyoto-u.ac.jp>
                  'sbg:toolkit': HLA-HD
                  'sbg:toolkitVersion': 1.4.0
                  'sbg:validationErrors': []
                  'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
                label: HLA-HD
                'sbg:x': 1256.66357421875
                'sbg:y': 280.875
              - id: bowtie2_aligner
                in:
                  - id: bowtie_index_archive
                    source: bowtie_index_archive
                  - id: paired_aligned_reads
                    default: raw
                  - id: read_sequence
                    source:
                      - read_sequence
                  - id: suppress_sam_records
                    default: true
                out:
                  - id: result_sam_file
                  - id: aligned_reads_only
                  - id: unaligned_reads_only
                run:
                  cwlVersion: 'sbg:draft-2'
                  class: CommandLineTool
                  $namespaces:
                    sbg: 'https://sevenbridges.com'
                  id: admin/sbg-public-data/bowtie2-aligner/18
                  label: Bowtie2 Aligner
                  description: >-
                    Bowtie 2 is an ultrafast and memory-efficient tool for
                    aligning sequencing reads to long reference sequences. It is
                    particularly good at aligning reads of about 50 up to 100s
                    or 1,000s of characters to relatively long (e.g. mammalian)
                    genomes. Bowtie 2 indexes the genome with an [FM
                    Index](http://portal.acm.org/citation.cfm?id=796543) (based
                    on the [Burrows-Wheeler
                    Transform](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform)
                    or
                    [BWT](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform))
                    to keep its memory footprint small: for the human genome,
                    its memory footprint is typically around 3.2 gigabytes of
                    RAM. In order to create needed index files, you should run
                    [Bowtie2
                    Indexer](https://igor.sbgenomics.com/public/apps#tool/admin/sbg-public-data/bowtie2-indexer),
                    which produces archived index files (containing 6 files with
                    suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and
                    .rev.2.bt2).


                    Bowtie 2 supports gapped, local, and paired-end alignment
                    modes. Bowtie 2 outputs alignments in SAM format, enabling
                    interoperation with a large number of other tools (e.g.
                    [SAMtools](http://samtools.sourceforge.net/),
                    [GATK](http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit))
                    that use SAM.


                    ###Common issues###

                    No issues have been reported.


                    **Q&A:**


                    ***Q: What should I do if I already have Bowtie2 index
                    files, not archived as tar bundle?***


                    ***A***: You can provide your *.bt2 files to [SBG
                    Compressor](https://igor.sbgenomics.com/public/apps#admin/sbg-public-data/sbg-compressor-1-0/)
                    app from our public apps and set "TAR" as your output
                    format. After the task is finished, **you should assign
                    common prefix of the index files to the `Reference genome`
                    metadata field** and your TAR is ready for use.


                    ***Example:***

                    Indexed files: chr20.1.bt2, chr20.2.bt2, chr20.3.bt2,
                    chr20.4.bt2, chr20.rev.1.bt2, chr20.rev.2.bt2


                    Metadata - `Reference genome`: **chr20**


                    __Important note: In case of paired-end alignment it is
                    crucial to set metadata 'paired-end' field to 1/2. Sequences
                    specified as mate 1s must correspond file-for-file and
                    read-for-read with those specified for mate 2s. Reads may be
                    a mix of different lengths. In case of unpaired reads, the
                    same metadata field should be set to '-'. Only one type of
                    alignment can be performed at once, so all specified reads
                    should be either paired or unpaired.__
                  baseCommand:
                    - class: Expression
                      engine: '#cwl-js-engine'
                      script: |-
                        {
                          var archive_name = $job.inputs.bowtie_index_archive.path.split("/").pop()
                          return "tar -xvf ".concat(archive_name, " && rm -rf ", archive_name, " && ")
                        }
                    - /opt/bowtie2-2.2.6/bowtie2
                  inputs:
                    - 'sbg:category': Input files
                      'sbg:stageInput': link
                      type:
                        - File
                      label: Bowtie index archive
                      description: Archive file produced by Bowtie2 Indexer.
                      'sbg:fileTypes': TAR
                      id: '#bowtie_index_archive'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': End-to-end
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - End-to-end
                            - Local
                          name: alignment_mode
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              if ($job.inputs.alignment_mode == "Local") {
                                return "--local"
                              }
                            }
                        'sbg:cmdInclude': true
                      label: Alignment mode
                      description: >-
                        Alignment mode. End-to-end: entire read must align; no
                        clipping. Local: local alignment; ends might be soft
                        clipped.
                      id: '#alignment_mode'
                    - 'sbg:altPrefix': '-s'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': '-'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--skip'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Skip reads
                      description: >-
                        Skip (i.e. do not align) the first given number of reads
                        or pairs in the input.
                      id: '#skip_reads'
                    - 'sbg:altPrefix': '-u'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': No limit
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--upto'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Align next n reads
                      description: >-
                        Align the first given number of reads or read pairs from
                        the input (after the <int> reads or pairs have been
                        skipped with "Skip reads"), then stop.
                      id: '#align_next_n_reads'
                    - 'sbg:altPrefix': '-5'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--trim5'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Trim from 5'
                      description: >-
                        Trim given number of bases from 5' (left) end of each
                        read before alignment.
                      id: '#trim_from_5'
                    - 'sbg:altPrefix': '-3'
                      'sbg:category': Input
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--trim3'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Trim from 3'
                      description: >-
                        Trim given number of bases from 3' (right) end of each
                        read before alignment.
                      id: '#trim_from_3'
                    - 'sbg:category': Input
                      'sbg:toolDefaultValue': Phred+33
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Auto-detect
                            - Phred+33
                            - Phred+64
                            - Solexa
                          name: quality_scale
                      label: Quality scale
                      description: Set quality scale.
                      id: '#quality_scale'
                    - 'sbg:category': Input
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--int-quals'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Integer qualities
                      description: >-
                        Quality values are represented in the read input file as
                        space-separated ASCII integers, e.g., 40 40 30 40...,
                        rather than ASCII characters, e.g., II?I....
                      id: '#integer_qualities'
                    - 'sbg:category': Alignment
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - '0'
                            - '1'
                          name: allowed_mismatch_number
                      inputBinding:
                        position: 0
                        prefix: '-N'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Allowed mismatch number
                      description: >-
                        Sets the number of mismatches to allowed in a seed
                        alignment during multiseed alignment. Can be set to 0 or
                        1. Setting this higher makes alignment slower (often
                        much slower) but increases sensitivity.
                      id: '#allowed_mismatch_number'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 22 or 20 (depending on preset type and alignment mode)
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-L'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Seed substring length
                      description: >-
                        Sets the length of the seed substrings to align during
                        multiseed alignment. Smaller values make alignment
                        slower but more senstive. Must be > 3 and < 32. The
                        "Sensitive" preset is used by default, which sets this
                        option to 22 in "End-to-end" mode and to 20 in "Local"
                        mode.
                      id: '#seed_substring_length'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': '15'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--dpad'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Dynamic padding
                      description: >-
                        "Pads" dynamic programming problems by the given number
                        of columns on either side to allow gaps.
                      id: '#dynamic_padding'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': '4'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--gbar'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disallow gaps
                      description: >-
                        Disallow gaps within the given number of positions of
                        the beginning or end of the read.
                      id: '#disallow_gaps'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--ignore-quals'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Ignore qualities
                      description: >-
                        When calculating a mismatch penalty, always consider the
                        quality value at the mismatched position to be the
                        highest possible, regardless of the actual value. I.e.
                        treat all quality values as 30 on Phred scale. This is
                        also the default behavior when the input doesn't specify
                        quality values (e.g. when processing .fasta reads).
                      id: '#ignore_qualities'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--nofw'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Don't align forward
                      description: >-
                        If this option is specified, Bowtie2 will not attempt to
                        align unpaired reads to the forward (Watson) reference
                        strand. In paired-end mode, "Don't align forward" and
                        "Don't align reverse complement" pertain to the
                        fragments; i.e. specifying "Don't align forward" causes
                        Bowtie2 to explore only those paired-end configurations
                        corresponding to fragments from the reverse-complement
                        (Crick) strand. Default: both strands enabled.
                      id: '#dont_align_forward'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--norc'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Don't align reverse complement
                      description: >-
                        If this option is specified, Bowtie2 will not attempt to
                        align unpaired reads against the reverse-complement
                        (Crick) reference strand. In paired-end mode, "Don't
                        align forward" and "Don't align reverse complement"
                        pertain to the fragments; i.e. specifying "Don't align
                        forward" causes Bowtie2 to explore only those paired-end
                        configurations corresponding to fragments from the
                        reverse-complement (Crick) strand. Default: both strands
                        enabled.
                      id: '#dont_align_reverse_complement'
                    - 'sbg:category': Alignment
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-1mm-upfront'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable 1 mismatch alignments
                      description: >-
                        By default, Bowtie2 will attempt to find either an exact
                        or a 1-mismatch end-to-end alignment for the read before
                        trying the multiseed heuristic. Such alignments can be
                        found very quickly, and many short read alignments have
                        exact or near-exact end-to-end alignments. However, this
                        can lead to unexpected alignments when the user also
                        sets options governing the multiseed heuristic, like
                        "Seed substring length" (-L) and "Allowed mismatch
                        number" (-N). For instance, if the user specifies 0 for
                        "Allowed mismatch number" and "Seed substring length"
                        equal to the length of the read, the user will be
                        surprised to find 1-mismatch alignments reported. This
                        option prevents Bowtie2 from searching for 1-mismatch
                        end-to-end alignments before using the multiseed
                        heuristic, which leads to the expected behavior when
                        combined with options such as "Seed substring length"
                        and "Allowed mismatch number". This comes at the expense
                        of speed.
                      id: '#disable_1_mismatch_alignments'
                    - 'sbg:category': Presets
                      'sbg:toolDefaultValue': Sensitive
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Very fast
                            - Fast
                            - Sensitive
                            - Very sensitive
                          name: preset_option
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              var preset_option = $job.inputs.preset_option
                              var alignment_mode = $job.inputs.alignment_mode
                              
                              var presets = {
                                "Very fast": "--very-fast",
                                "Fast": "--fast",
                                "Sensitive": "--sensitive",
                                "Very sensitive": "--very-sensitive"
                              }
                              if (alignment_mode == "Local" && preset_option) {
                                return presets[preset_option].concat("-local")
                              }
                              else if (preset_option){
                                return presets[preset_option]
                              }
                            }
                        'sbg:cmdInclude': true
                      label: Preset
                      description: >-
                        Preset options for "Seed extension attempts" (-D), "Max
                        number of re-seed" (-R), "Allowed mismatch number" (-N),
                        "Seed substring length" (-L) and "Interval function"
                        (-i) parameters. Values for these options vary depending
                        on whether the "Local" or "End-to-end" mode is selected
                        under "Alignment mode".
                      id: '#preset_option'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '0 for "End-to-end" mode, 2 for "Local" mode'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--ma'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Set match bonus
                      description: >-
                        Sets the match bonus. The given number is added to the
                        alignment score for each position where a read character
                        aligns to a reference character and the characters
                        match.
                      id: '#set_match_bonus'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '6'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--mp'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Maximum mismatch penalty
                      description: >-
                        Sets the maximum penalty for mismatch. Lower quality =
                        lower penalty.
                      id: '#maximum_mismatch_penalty'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '1'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--np'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Ambiguous character penalty
                      description: >-
                        Sets penalty for positions where the read, reference, or
                        both, contain an ambiguous character such as N.
                      id: '#ambiguous_character_penalty'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '5,3'
                      type:
                        - 'null'
                        - type: array
                          items: int
                      inputBinding:
                        position: 0
                        prefix: '--rdg'
                        separate: true
                        itemSeparator: ','
                        'sbg:cmdInclude': true
                      label: Read gap penalties
                      description: >-
                        Sets the read gap open (first value) and extend (second
                        value) penalty, respectively. A read gap of length N
                        gets a penalty of <gap-open-penalty> + N *
                        <gap-extend-penalty>.
                      id: '#read_gap_penalties'
                    - 'sbg:category': Scoring
                      'sbg:toolDefaultValue': '5,3'
                      type:
                        - 'null'
                        - type: array
                          items: int
                      inputBinding:
                        position: 0
                        prefix: '--rfg'
                        separate: true
                        itemSeparator: ','
                        'sbg:cmdInclude': true
                      label: Reference gap penalties
                      description: >-
                        Sets the reference gap open (first value) and extend
                        (second value) penalty, respectively. A reference gap of
                        length N gets a penalty of <gap-open-penalty> + N *
                        <gap-extend-penalty>.
                      id: '#reference_gap_penalties'
                    - 'sbg:category': Reporting
                      'sbg:toolDefaultValue': '-'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-k'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Report k alignments
                      description: >-
                        By default, Bowtie2 searches for distinct, valid
                        alignments for each read. When it finds a valid
                        alignment, it continues looking for alignments that are
                        nearly as good or better. The best alignment found is
                        reported (randomly selected from among best if tied).
                        Information about the best alignments is used to
                        estimate mapping quality and to set SAM optional fields,
                        such as AS:i and XS:i. When "Report k alignments" is
                        specified, however, Bowtie2 behaves differently.
                        Instead, it searches for at most <given-number>
                        distinct, valid alignments for each read. The search
                        terminates when it can't find more distinct valid
                        alignments, or when it finds <given-number>, whichever
                        happens first. All alignments found are reported in
                        descending order by alignment score. The alignment score
                        for a paired-end alignment equals the sum of the
                        alignment scores of the individual mates. Each reported
                        read or pair alignment beyond the first has the SAM
                        'secondary' bit (which equals 256) set in its FLAGS
                        field. For reads that have more than <given-number>
                        distinct, valid alignments, Bowtie2 does not gaurantee
                        that the <given-number> alignments reported are the best
                        possible in terms of alignment score. "Report k
                        alignments" is mutually exclusive with "Report all
                        alignments". Note: Bowtie 2 is not designed with large
                        values for "Report k alignments" in mind, and when
                        aligning reads to long, repetitive genomes alignment can
                        be very, very slow.
                      id: '#report_k_alignments'
                    - 'sbg:altPrefix': '-a'
                      'sbg:category': Reporting
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--all'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Report all alignments
                      description: >-
                        Like "Report k alignments" but with no upper limit on
                        number of alignments to search for. "Report all
                        alignments" is mutually exclusive with "Report k
                        alignments".
                      id: '#report_all_alignments'
                    - 'sbg:category': Effort
                      'sbg:toolDefaultValue': '15'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-D'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Seed extension attempts
                      description: >-
                        Maximum number of to consecutive seed extension attempts
                        that can "fail" before Bowtie2 moves on, using the
                        alignments found so far. A seed extension "fails" if it
                        does not yield a new best or a new second-best
                        alignment.
                      id: '#seed_extension_attempts'
                    - 'sbg:category': Effort
                      'sbg:toolDefaultValue': '2'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '-R'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Max number of re-seed
                      description: >-
                        Given number is the maximum number of times Bowtie2 will
                        're-seed' reads with repetitive seeds. When
                        're-seeding', Bowtie2 simply chooses a new set of reads
                        (same length, same number of mismatches allowed) at
                        different offsets and searches for more alignments. A
                        read is considered to have repetitive seeds if the total
                        number of seed hits divided by the number of seeds that
                        aligned at least once is greater than 300.
                      id: '#max_number_of_re_seed'
                    - 'sbg:altPrefix': '-I'
                      'sbg:category': Paired-end
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--minins'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Minimum fragment length
                      description: >-
                        The minimum fragment length for valid paired-end
                        alignments. E.g. if 60 is specified for "Minimum
                        fragment length" (-I) and a paired-end alignment
                        consists of two 20-bp alignments in the appropriate
                        orientation with a 20-bp gap between them, that
                        alignment is considered valid (as long as "Maximum
                        fragment length" (-X) is also satisfied). A 19-bp gap
                        would not be valid in that case. If trimming options -3
                        or -5 are also used, the "Minimum fragment length"
                        constraint is applied with respect to the untrimmed
                        mates. The larger the difference between "Minimum
                        fragment length" and "Maximum fragment length", the
                        slower Bowtie2 will run. This is because larger
                        differences bewteen those two require that Bowtie2 scan
                        a larger window to determine if a concordant alignment
                        exists. For typical fragment length ranges (200 to 400
                        nucleotides), Bowtie2 is very efficient.
                      id: '#minimum_fragment_length'
                    - 'sbg:altPrefix': '-X'
                      'sbg:category': Paired-end
                      'sbg:toolDefaultValue': '500'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--maxins'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Maximum fragment length
                      description: >-
                        The maximum fragment length for valid paired-end
                        alignments. E.g. if "Maximum fragment length" (-X) 100
                        is specified and a paired-end alignment consists of two
                        20-bp alignments in the proper orientation with a 60-bp
                        gap between them, that alignment is considered valid (as
                        long as "Minimum fragment length" (-I) is also
                        satisfied). A 61-bp gap would not be valid in that case.
                        If trimming options -3 or -5 are also used, the "Maximum
                        fragment length" constraint is applied with respect to
                        the untrimmed mates, not the trimmed mates. The larger
                        the difference between "Minimum fragment length" and
                        "Maximum fragment length", the slower Bowtie2 will run.
                        This is because larger differences bewteen those two
                        require that Bowtie2 scan a larger window to determine
                        if a concordant alignment exists. For typical fragment
                        length ranges (200 to 400 nucleotides), Bowtie2 is very
                        efficient.
                      id: '#maximum_fragment_length'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': '--fr'
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - '--fr'
                            - '--rf'
                            - '--ff'
                          name: mates_alignment_orientation
                      inputBinding:
                        position: 0
                        separate: true
                        'sbg:cmdInclude': true
                      label: Mates alignment orientation
                      description: >-
                        The upstream/downstream mate orientations for a valid
                        paired-end alignment against the forward reference
                        strand. E.g., if --fr is specified and there is a
                        candidate paired-end alignment where mate 1 appears
                        upstream of the reverse complement of mate 2 and the
                        fragment length constraints ("Minimum fragment length"
                        (-I) and "Maximum fragment length" (-X)) are met, that
                        alignment is valid. Also, if mate 2 appears upstream of
                        the reverse complement of mate 1 and all other
                        constraints are met, that too is valid. --rf likewise
                        requires that an upstream mate1 be reverse-complemented
                        and a downstream mate2 be forward-oriented. --ff
                        requires both an upstream mate 1 and a downstream mate 2
                        to be forward-oriented. Default orientation --fr is
                        appropriate for Illumina's Paired-end Sequencing Assay.
                      id: '#mates_alignment_orientation'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-mixed'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable unpaired alignments
                      description: >-
                        By default, when Bowtie2 cannot find a concordant or
                        discordant alignment for a pair, it then tries to find
                        alignments for the individual mates. This option
                        disables that behavior.
                      id: '#disable_unpaired_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-discordant'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable discordant alignments
                      description: >-
                        By default, Bowtie2 looks for discordant alignments if
                        it cannot find any concordant alignments. A discordant
                        alignment is an alignment where both mates align
                        uniquely, but that does not satisfy the paired-end
                        constraints (--fr/--rf/--ff, -I, -X). This option
                        disables that behavior.
                      id: '#disable_discordant_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-dovetail'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable dovetail alignments
                      description: >-
                        If the mates "dovetail", that is if one mate alignment
                        extends past the beginning of the other such that the
                        wrong mate begins upstream, consider that to be
                        non-concordant.
                      id: '#disable_dovetail_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-contain'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable containing alignments
                      description: >-
                        If one mate alignment contains the other, consider that
                        to be non-concordant.
                      id: '#disable_containing_alignments'
                    - 'sbg:category': Paired-end
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-overlap'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Disable overlapping alignments
                      description: >-
                        If one mate alignment overlaps the other at all,
                        consider that to be non-concordant.
                      id: '#disable_overlapping_alignments'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-head'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Suppress header lines
                      description: Suppress SAM header lines (starting with @).
                      id: '#suppress_header_lines'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-sq'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Suppress SQ header lines
                      description: Suppress @SQ SAM header lines.
                      id: '#suppress_sq_header_lines'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': id
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        prefix: '--rg-id'
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              return ($job.inputs.read_group_id || "id") 
                            }
                        'sbg:cmdInclude': true
                      label: Set the read group ID
                      description: >-
                        Set the read group ID text. This causes the SAM @RG
                        header line to be printed, with the given text as the
                        value associated with the ID: tag. It also causes the
                        RG:Z: extra field to be attached to each SAM output
                        record, with value set to this text.
                      id: '#read_group_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - LS 454
                            - Helicos
                            - Illumina
                            - ABI SOLiD
                            - Ion Torrent PGM
                            - PacBio
                          name: platform
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: |-
                            {
                              if($job.inputs.platform)
                                return "--rg PL:" +$job.inputs.platform.replace(/ /g,"_")
                              else if([].concat($job.inputs.read_sequence)[0].metadata){
                                if ([].concat($job.inputs.read_sequence)[0].metadata.platform) {
                                  return "--rg PL:" +[].concat($job.inputs.read_sequence)[0].metadata.platform.replace(/ /g,"_")
                                }
                              }
                              else {
                                return ""
                              }
                            }
                        'sbg:cmdInclude': true
                      label: Platform
                      description: >-
                        Specify the version of the technology that was used for
                        sequencing or assaying. Default: inferred from metadata.
                      id: '#platform'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.sample_id)\n    return \"--rg SM:\" +$job.inputs.sample_id\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.sample_id)\n  \t\treturn \"--rg SM:\" +[].concat($job.inputs.read_sequence)[0].metadata.sample_id\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Sample
                      description: Specify the sample ID for RG line.
                      id: '#sample_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.library_id)\n    return \"--rg LB:\" +$job.inputs.library_id\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.library_id)\n  \t\treturn \"--rg LB:\" +[].concat($job.inputs.read_sequence)[0].metadata.library_id\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Library
                      description: Specify the library ID for RG line.
                      id: '#library_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.platform_unit_id)\n    return \"--rg PU:\" +$job.inputs.platform_unit_id\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.platform_unit_id)\n  \t\treturn \"--rg PU:\" +[].concat($job.inputs.read_sequence)[0].metadata.platform_unit_id\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Platform unit
                      description: Specify the platform unit ID for RG line.
                      id: '#platform_unit_id'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - string
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.sequencing_center)\n    return \"--rg CN:\" +$job.inputs.sequencing_center\n    else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.seq_center)\n  \t\treturn \"--rg CN:\" +[].concat($job.inputs.read_sequence)[0].metadata.seq_center\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Sequencing center
                      description: Specify the sequencing center for RG line.
                      id: '#sequencing_center'
                    - 'sbg:category': Read group
                      'sbg:toolDefaultValue': Infered from metadata
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n  if($job.inputs.median_fragment_length)\n    return \"--rg PI:\" +$job.inputs.median_fragment_length\n  else if([].concat($job.inputs.read_sequence)[0].metadata){\n      if([].concat($job.inputs.read_sequence)[0].metadata.median_fragment_length)\n  \t\treturn \"--rg PI:\" +[].concat($job.inputs.read_sequence)[0].metadata.median_fragment_length\n      else return \"\"\n    }\n  else\n    return \"\"\n}"
                        'sbg:cmdInclude': true
                      label: Median fragment length
                      description: Specify the median fragment length for RG line.
                      id: '#median_fragment_length'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--omit-sec-seq'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Omit SEQ and QUAL
                      description: >-
                        When printing secondary alignments, Bowtie 2 by default
                        will write out the SEQ and QUAL strings. Specifying this
                        option causes Bowtie 2 to print an asterisk ('*') in
                        those fields instead.
                      id: '#omit_seq_and_qual'
                    - 'sbg:category': Performance
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--reorder'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Reorder output
                      description: >-
                        Guarantees that output SAM records are printed in an
                        order corresponding to the order of the reads in the
                        original input file. Specifying "Reorder output" causes
                        Bowtie2 to run somewhat slower and use somewhat more
                        memory.
                      id: '#reorder_output'
                    - 'sbg:category': Other
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        prefix: '--seed'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Set seed
                      description: Set the seed for pseudo-random number generator.
                      id: '#set_seed'
                    - 'sbg:category': Other
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--non-deterministic'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Non deterministic
                      description: >-
                        Normally, Bowtie2 re-initializes its pseudo-random
                        generator for each read. It seeds the generator with a
                        number derived from (a) the read name, (b) the
                        nucleotide sequence, (c) the quality sequence, (d) the
                        value of the "Set seed" option. This means that if two
                        reads are identical (same name, same nucleotides, same
                        qualities) Bowtie2 will find and report the same
                        alignment(s) for both, even if there was ambiguity. When
                        "Non deterministic" is specified, Bowtie2 re-initializes
                        its pseudo-random generator for each read using the
                        current time. This means that Bowtie2 will not
                        necessarily report the same alignment for two identical
                        reads. This is counter-intuitive for some users, but
                        might be more appropriate in situations where the input
                        consists of many identical reads.
                      id: '#non_deterministic'
                    - 'sbg:category': Interval function
                      'sbg:toolDefaultValue': Square-root
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Constant
                            - Linear
                            - Square-root
                            - Natural log
                          name: function_i
                      label: Interval function
                      description: >-
                        Sets a function type F in function f governing the
                        interval between seed substrings, to use during
                        multiseed alignment. The interval function f is f(x) = A
                        + B * F(x), where x is the read length. By default,
                        function F is set to 'Square-root', Constant A to 1 and
                        Coefficient B to 1.15 or 0.75 for "End-to-end" and
                        "Local" mode respectively.
                      id: '#function_i'
                    - 'sbg:category': Interval function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '1'
                      type:
                        - 'null'
                        - string
                      label: 'Constant A [ interval function ]'
                      description: >-
                        Sets a constant A in function governing the interval
                        between seed substrings to use during multiseed
                        alignment. The interval function f is f(x) = A + B *
                        F(x), where x is the read length.
                      id: '#constant_i_a'
                    - 'sbg:category': Interval function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': 1.15 or 0.75 (depending on "Alignment mode")
                      type:
                        - 'null'
                        - string
                      label: 'Coefficient B [ interval function ]'
                      description: >-
                        Sets a coefficient B in function governing the interval
                        between seed substrings to use during multiseed
                        alignment. The interval function f is f(x) = A + B *
                        F(x), where x is the read length. Default: 1.15 in
                        "End-to-end" mode and 0.75 in "Local" mode.
                      id: '#coefficient_i_b'
                    - 'sbg:category': Ambiguous chars function
                      'sbg:toolDefaultValue': Linear
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Constant
                            - Linear
                            - Square-root
                            - Natural log
                          name: function_n_ceil
                      label: Ambiguous chars function
                      description: >-
                        Sets a function type F in function governing the maximum
                        number of ambiguous characters (usually Ns and/or .s)
                        allowed in a read as a function of read length. The
                        N-ceiling function f is f(x) = A + B * F(x), where x is
                        the read length. Reads exceeding this ceiling are
                        filtered out.
                      id: '#function_n_ceil'
                    - 'sbg:category': Ambiguous chars function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '0'
                      type:
                        - 'null'
                        - string
                      label: 'Constant A [ ambiguous chars function ]'
                      description: >-
                        Sets a constant A in function governing the maximum
                        number of ambiguous characters (usually Ns and/or .s)
                        allowed in a read as a function of read length. The
                        N-ceiling function f is f(x) = A + B * F(x), where x is
                        the read length. Reads exceeding this ceiling are
                        filtered out.
                      id: '#constant_nceil_a'
                    - 'sbg:category': Ambiguous chars function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': '0.15'
                      type:
                        - 'null'
                        - string
                      label: 'Coefficient B [ ambiguous chars function ]'
                      description: >-
                        Sets a coefficient B in function governing the maximum
                        number of ambiguous characters (usually Ns and/or .s)
                        allowed in a read as a function of read length. The
                        N-ceiling function f is f(x) = A + B * F(x), where x is
                        the read length. Reads exceeding this ceiling are
                        filtered out.
                      id: '#coefficient_nceil_b'
                    - 'sbg:category': Alignment score function
                      'sbg:toolDefaultValue': Natural log or Linear (depending on "Alignment mode")
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - Constant
                            - Linear
                            - Square-root
                            - Natural log
                          name: function_score_min
                      label: Alignment score function
                      description: >-
                        Sets a function type F in function governing the minimum
                        alignment score needed for an alignment to be considered
                        "valid" (i.e. good enough to report). This is a function
                        of read length. The minimum-score function f is f(x) = A
                        + B * F(x), where x is the read length. By default,
                        function F is set to "Natural log" or "Linear", Constant
                        A to 20 or -0.6 and Coefficient B to 8 or -0.6 depending
                        on the "Alignment mode": "End-to-end" or "Local"
                        respectively.
                      id: '#function_score_min'
                    - 'sbg:category': Alignment score function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': 20 or -0.6 (depending on "Alignment mode")
                      type:
                        - 'null'
                        - string
                      label: 'Constant A [ alignment score function ]'
                      description: >-
                        Sets a constant A in function governing the minimum
                        alignment score needed for an alignment to be considered
                        'valid' (i.e. good enough to report). This is a function
                        of read length. The minimum-score function f is f(x) = A
                        + B * F(x), where x is the read length. Default: 20 in
                        "End-to-end" mode and -0.6 in "Local" mode.
                      id: '#constant_scoremin_a'
                    - 'sbg:category': Alignment score function
                      'sbg:stageInput': null
                      'sbg:toolDefaultValue': 8 or -0.6 (depending on "Alignment mode")
                      type:
                        - 'null'
                        - string
                      label: 'Coefficient B [ alignment score function ]'
                      description: >-
                        Sets a coefficient B in function governing the minimum
                        alignment score needed for an alignment to be considered
                        'valid' (i.e. good enough to report). This is a function
                        of read length. The minimum-score function f is f(x) = A
                        + B * F(x), where x is the read length. Default: 8 in
                        "End-to-end" mode and -0.6 in "Local" mode.
                      id: '#coefficient_scoremin_b'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: unpaired_unaligned_reads
                      label: Unpaired unaligned reads
                      description: >-
                        Output unpaired reads that fail to align. These reads
                        correspond to the SAM records with the FLAGS 0x4 bit set
                        and neither the 0x40 nor 0x80 bits set. If "gzip
                        compressed" is specified, output will be gzip
                        compressed. If "bzip2 compressed" is specified, output
                        will be bzip2 compressed. Reads written in this way will
                        appear exactly as they did in the input file, without
                        any modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#unpaired_unaligned_reads'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: unpaired_aligned_reads
                      label: Unpaired aligned reads
                      description: >-
                        Output unpaired reads that align at least once. These
                        reads correspond to the SAM records with the FLAGS 0x4,
                        0x40, and 0x80 bits unset. If "gzip compressed" is
                        specified, output will be gzip compressed. If "bzip2
                        compressed" is specified, output will be bzip2
                        compressed. Reads written in this way will appear
                        exactly as they did in the input file, without any
                        modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#unpaired_aligned_reads'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: paired_unaligned_reads
                      label: Paired unaligned reads
                      description: >-
                        Output paired-end reads that fail to align concordantly.
                        These reads correspond to the SAM records with the FLAGS
                        0x4 bit set and either the 0x40 or 0x80 bit set
                        (depending on whether it's mate #1 or #2). .1 and .2
                        strings are added to the filename to distinguish which
                        file contains mate #1 and mate #2. If "gzip compressed"
                        is specified, output will be gzip compressed. If "bzip2
                        compressed" is specified, output will be bzip2
                        compressed. Reads written in this way will appear
                        exactly as they did in the input file, without any
                        modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#paired_unaligned_reads'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': None
                      type:
                        - 'null'
                        - type: enum
                          symbols:
                            - raw
                            - gzip compressed
                            - bzip2 compressed
                            - None
                          name: paired_aligned_reads
                      label: Paired aligned reads
                      description: >-
                        Output paired-end reads that align concordantly at least
                        once. These reads correspond to the SAM records with the
                        FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set
                        (depending on whether it's mate #1 or #2). .1 and .2
                        strings are added to the filename to distinguish which
                        file contains mate #1 and mate #2. If "gzip compressed"
                        is specified, output will be gzip compressed. If "bzip2
                        compressed" is specified, output will be bzip2
                        compressed. Reads written in this way will appear
                        exactly as they did in the input file, without any
                        modification (same sequence, same name, same quality
                        string, same quality encoding). Reads will not
                        necessarily appear in the same order as they did in the
                        input.
                      id: '#paired_aligned_reads'
                    - 'sbg:category': Input files
                      type:
                        - type: array
                          items: File
                      label: Read sequence
                      description: >-
                        Read sequence in FASTQ or FASTA format. COuld be also
                        gzip'ed (extension .gz) or bzip2'ed (extension .bz2). In
                        case of paired-end alignment it is crucial to set
                        metadata 'paired-end' field to 1/2.
                      'sbg:fileTypes': >-
                        FASTA, FASTA.GZ, FASTA.BZ2, FA.GZ, FA.BZ2, FASTQ, FA,
                        FQ, FASTQ.GZ, FQ.GZ, FASTQ.BZ2, FQ.BZ2
                      id: '#read_sequence'
                    - 'sbg:category': Output
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '--no-unal'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Suppress SAM records for unaligned reads
                      description: Suppress SAM records for reads that failed to align.
                      id: '#suppress_sam_records'
                    - 'sbg:category': Input
                      'sbg:toolDefaultValue': 'False'
                      type:
                        - 'null'
                        - boolean
                      inputBinding:
                        position: 0
                        prefix: '-f'
                        separate: true
                        'sbg:cmdInclude': true
                      label: Input FASTA files
                      description: >-
                        Reads (specified with <m1>, <m2>, <s>) are FASTA files.
                        FASTA files usually have extension .fa, .fasta, .mfa,
                        .fna or similar. FASTA files do not have a way of
                        specifying quality values, so when -f is set, the result
                        is as if --ignore-quals is also set.
                      id: '#input_fasta_files'
                    - 'sbg:altPrefix': '-threads'
                      'sbg:category': Performance
                      'sbg:toolDefaultValue': '8'
                      type:
                        - 'null'
                        - int
                      inputBinding:
                        position: 0
                        separate: true
                        valueFrom:
                          class: Expression
                          engine: '#cwl-js-engine'
                          script: "{\n\tif($job.inputs.threads)\n    {\n    \treturn \" -p \" + $job.inputs.threads\n    }\n  \telse\n    {\n    \treturn \" -p 8 \"\n    }\n}"
                        'sbg:cmdInclude': true
                      label: Number of threads
                      description: >-
                        Launch NTHREADS parallel search threads (default: 1).
                        Threads will run on separate processors/cores and
                        synchronize when parsing reads and outputting
                        alignments. Searching for alignments is highly parallel,
                        and speedup is close to linear. Increasing -p increases
                        Bowtie 2's memory footprint. E.g. when aligning to a
                        human genome index, increasing -p from 1 to 8 increases
                        the memory footprint by a few hundred megabytes. This
                        option is only available if bowtie is linked with the
                        pthreads library (i.e. if BOWTIE_PTHREADS=0 is not
                        specified at build time).
                      id: '#threads'
                  outputs:
                    - type:
                        - 'null'
                        - File
                      label: Result SAM file
                      description: >-
                        SAM file containing the results of the alignment. It
                        contains both aligned and unaligned reads.
                      'sbg:fileTypes': SAM
                      outputBinding:
                        streamable: false
                        glob: '*.sam'
                        'sbg:metadata':
                          __inherit__: fasta_reference
                        'sbg:inheritMetadataFrom': '#read_sequence'
                      id: '#result_sam_file'
                    - type:
                        - 'null'
                        - type: array
                          items: File
                      label: Aligned reads only
                      description: FASTQ file with reads that align at least once.
                      'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                      outputBinding:
                        streamable: false
                        glob: '*_aligned*'
                        'sbg:metadata':
                          __inherit__: fasta_reference
                        'sbg:inheritMetadataFrom': '#read_sequence'
                      id: '#aligned_reads_only'
                    - type:
                        - 'null'
                        - type: array
                          items: File
                      label: Unaligned reads only
                      description: FASTQ file with reads that failed to align.
                      'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTQ.BZ2'
                      outputBinding:
                        streamable: false
                        glob: '*_unaligned*'
                        'sbg:metadata':
                          __inherit__: fasta_reference
                        'sbg:inheritMetadataFrom': '#read_sequence'
                      id: '#unaligned_reads_only'
                  requirements:
                    - class: ExpressionEngineRequirement
                      engineCommand: cwl-engine.js
                      id: '#cwl-js-engine'
                      requirements:
                        - class: DockerRequirement
                          dockerPull: rabix/js-engine
                  hints:
                    - class: 'sbg:CPURequirement'
                      value: 8
                    - class: 'sbg:MemRequirement'
                      value: 6000
                    - class: DockerRequirement
                      dockerImageId: 029d3a264215
                      dockerPull: 'images.sbgenomics.com/ana_d/bowtie2:2.2.6'
                  arguments:
                    - position: 100
                      prefix: ''
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: "cmd = \"\"\nreads = [].concat($job.inputs.read_sequence)\nreads1 = [];\nreads2 = [];\nu_reads = [];\nfor (var i = 0; i < reads.length; i++){\n    if (reads[i].metadata.paired_end == 1){\n      reads1.push(reads[i].path);\n    }\n    else if (reads[i].metadata.paired_end == 2){\n      reads2.push(reads[i].path);\n    }\n  else {\n  \tu_reads.push(reads[i].path);\n   }\n  }\nif (reads1.length > 0 & reads1.length == reads2.length){\n\tcmd = \"-1 \" + reads1.join(\",\") + \" -2 \" + reads2.join(\",\");\n}\nif (u_reads.length > 0){\n\tcmd = \" -U \" + u_reads.join(\",\");\n}\ncmd\n"
                    - position: 101
                      prefix: '-S'
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                            var A= array.concat().sort(), 
                                a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                            while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                            return a1.substring(0, i);
                            }
                            path_list = []
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            return "./".concat(common_prefix.replace( /\-$|\_$|\.$/, '' ), ".", "sam")
                          }
                    - position: 0
                      prefix: '-x'
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |
                          {
                            var index_prefix = $job.inputs.bowtie_index_archive.metadata.reference_genome
                            return index_prefix
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.unpaired_unaligned_reads && $job.inputs.unpaired_unaligned_reads != "None") {
                              if ($job.inputs.unpaired_unaligned_reads == "raw") {
                                return "--un ".concat(common_prefix, ".unpaired_unaligned.fastq")
                              }
                              else if ($job.inputs.unpaired_unaligned_reads == "gzip compressed") {
                                return "--un ".concat(common_prefix, ".unpaired_unaligned.fastq.gz")
                              }
                              else if ($job.inputs.unpaired_unaligned_reads == "bzip2 compressed") {
                                return "--un ".concat(common_prefix, ".unpaired_unaligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.unpaired_aligned_reads && $job.inputs.unpaired_aligned_reads != "None") {
                              if ($job.inputs.unpaired_aligned_reads == "raw") {
                                return "--al ".concat(common_prefix, ".unpaired_aligned.fastq")
                              }
                              else if ($job.inputs.unpaired_aligned_reads == "gzip compressed") {
                                return "--al ".concat(common_prefix, ".unpaired_aligned.fastq.gz")
                              }
                              else if ($job.inputs.unpaired_aligned_reads == "bzip2 compressed") {
                                return "--al ".concat(common_prefix, ".unpaired_aligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.paired_unaligned_reads && $job.inputs.paired_unaligned_reads != "None") {
                              if ($job.inputs.paired_unaligned_reads == "raw") {
                                return "--un-conc ".concat(common_prefix, ".paired_unaligned.fastq")
                              }
                              else if ($job.inputs.paired_unaligned_reads == "gzip compressed") {
                                return "--un-conc ".concat(common_prefix, ".paired_unaligned.fastq.gz")
                              }
                              else if ($job.inputs.paired_unaligned_reads == "bzip2 compressed") {
                                return "--un-conc ".concat(common_prefix, ".paired_unaligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            
                            common_prefix = sharedStart(path_list)
                            
                            if ($job.inputs.paired_aligned_reads && $job.inputs.paired_aligned_reads != "None") {
                              if ($job.inputs.paired_aligned_reads == "raw") {
                                return "--al-conc ".concat(common_prefix, ".paired_aligned.fastq")
                              }
                              else if ($job.inputs.paired_aligned_reads == "gzip compressed") {
                                return "--al-conc ".concat(common_prefix, ".paired_aligned.fastq.gz")
                              }
                              else if ($job.inputs.paired_aligned_reads == "bzip2 compressed") {
                                return "--al-conc ".concat(common_prefix, ".paired_aligned.fastq.bz2")
                              }
                            }
                          }
                    - position: 0
                      prefix: ''
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            var functions = {
                              "Constant": "C",
                              "Linear": "L",
                              "Square-root": "S",
                              "Natural log": "G"
                            }
                            function_type = $job.inputs.function_i
                            value_list = [functions[function_type], $job.inputs.constant_i_a, $job.inputs.coefficient_i_b]
                            if (functions[function_type] && $job.inputs.constant_i_a && $job.inputs.coefficient_i_b) {
                              return "-i ".concat(value_list.join(","))
                            }
                          }
                    - position: 0
                      prefix: ''
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            var functions = {
                              "Constant": "C",
                              "Linear": "L",
                              "Square-root": "S",
                              "Natural log": "G"
                            }
                            function_type = $job.inputs.function_n_ceil
                            value_list = [functions[function_type], $job.inputs.constant_nceil_a, $job.inputs.coefficient_nceil_b]
                            if (functions[function_type] && $job.inputs.constant_nceil_a && $job.inputs.coefficient_nceil_b) {
                              return "--n-ceil ".concat(value_list.join(","))
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            var functions = {
                              "Constant": "C",
                              "Linear": "L",
                              "Square-root": "S",
                              "Natural log": "G"
                            }
                            function_type = $job.inputs.function_score_min
                            
                            value_list = [functions[function_type], $job.inputs.constant_scoremin_a, $job.inputs.coefficient_scoremin_b]
                            if (functions[function_type] && $job.inputs.constant_scoremin_a && $job.inputs.coefficient_scoremin_b) {
                              return "--score-min ".concat(value_list.join(","))
                            }
                          }
                    - position: 0
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            meta_qual="";
                            if ([].concat($job.inputs.read_sequence)[0].metadata){
                              if ([].concat($job.inputs.read_sequence)[0].metadata.quality_scale){
                                meta_qual = [].concat($job.inputs.read_sequence)[0].metadata.quality_scale
                              }
                            }
                            
                            if ($job.inputs.quality_scale == "Phred+33") {
                              return "--phred33"
                            }
                            else if ($job.inputs.quality_scale == "Phred+64") {
                              return "--phred64"
                            }
                            else if ($job.inputs.quality_scale == "Solexa") {
                              return "--solexa-quals"
                            }
                            else if ($job.inputs.quality_scale == "Auto-detect") {
                              if (meta_qual == "solexa") {
                                return "--solexa-quals"
                              }
                              else if (meta_qual == "illumina13" || meta_qual == "illumina15") {
                                return "--phred64"
                              }
                            }
                          }
                    - position: 102
                      separate: true
                      valueFrom:
                        class: Expression
                        engine: '#cwl-js-engine'
                        script: |-
                          {
                            function sharedStart(array){
                              var A= array.concat().sort(), 
                                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                              while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                              return a1.substring(0, i);
                            }
                            path_list = []
                            reads = [].concat($job.inputs.read_sequence)
                            reads.forEach(function(f){return path_list.push(f.path.replace(/\\/g,'/').replace( /.*\//, '' ))})
                            common_prefix = sharedStart(path_list)
                            
                            gzip = "gzip compressed"
                            bzip = "bzip2 compressed"
                            paired_aligned = $job.inputs.paired_aligned_reads
                            paired_unaligned = $job.inputs.paired_unaligned_reads
                            aligned_first = common_prefix.concat(".paired_aligned.fastq.1")
                            aligned_second = common_prefix.concat(".paired_aligned.fastq.2")
                            aligned_first_mv = common_prefix.concat(".paired_aligned.1.fastq")
                            aligned_second_mv = common_prefix.concat(".paired_aligned.2.fastq")
                            
                            unaligned_first = common_prefix.concat(".paired_unaligned.fastq.1")
                            unaligned_second = common_prefix.concat(".paired_unaligned.fastq.2")
                            unaligned_first_mv = common_prefix.concat(".paired_unaligned.1.fastq")
                            unaligned_second_mv = common_prefix.concat(".paired_unaligned.2.fastq")
                            
                            aligned = ""
                            unaligned = ""
                            
                            if (paired_aligned && paired_aligned == gzip) {
                              aligned = "&& mv ".concat(aligned_first, ".gz ", aligned_first_mv, ".gz && mv ", aligned_second, ".gz ", aligned_second_mv, ".gz ") 
                            }
                            else if (paired_aligned && paired_aligned == bzip) {
                              aligned = "&& mv ".concat(aligned_first, ".bz2 ", aligned_first_mv, ".bz2 && mv ", aligned_second, ".bz2 ", aligned_second_mv, ".bz2 ")
                            }
                            if (paired_unaligned && paired_unaligned == gzip) {
                              unaligned = "&& mv ".concat(unaligned_first, ".gz ", unaligned_first_mv, ".gz && mv ", unaligned_second, ".gz ", unaligned_second_mv, ".gz")
                            }
                            else if (paired_unaligned && paired_unaligned == bzip) {
                              unaligned = "&& mv ".concat(unaligned_first, ".bz2 ", unaligned_first_mv, ".bz2 && mv ", unaligned_second, ".bz2 ", unaligned_second_mv, ".bz2")
                            }
                            
                            return aligned.concat(unaligned)
                          }
                  'sbg:appVersion':
                    - 'sbg:draft-2'
                  'sbg:categories':
                    - Alignment
                  'sbg:cmdPreview': >-
                    tar -xvf chr20_bowtie2-2.2.6.tar && rm -rf
                    chr20_bowtie2-2.2.6.tar &&  /opt/bowtie2-2.2.6/bowtie2 -x
                    chr20  --un mate.unpaired_unaligned.fastq.gz  --al
                    mate.unpaired_aligned.fastq.gz  --un-conc
                    mate.paired_unaligned.fastq  --al-conc
                    mate.paired_aligned.fastq.gz  -i
                    L,constant_i_a-string-value,coefficient_i_b-string-value 
                    --n-ceil
                    S,constant_nceil_a-string-value,coefficient_nceil_b-string-value 
                    --score-min
                    L,constant_scoremin_a-string-value,coefficient_scoremin_b-string-value 
                    --phred33  -1 /demo/test-data/mate1.fastq -2
                    /demo/test-data/mate2.fastq -S ./mate.sam  && mv
                    mate.paired_aligned.fastq.1.gz
                    mate.paired_aligned.1.fastq.gz && mv
                    mate.paired_aligned.fastq.2.gz
                    mate.paired_aligned.2.fastq.gz
                  'sbg:content_hash': null
                  'sbg:contributors':
                    - sevenbridges-pgc
                    - admin
                  'sbg:createdBy': sevenbridges-pgc
                  'sbg:createdOn': 1454426457
                  'sbg:id': admin/sbg-public-data/bowtie2-aligner/18
                  'sbg:image_url': null
                  'sbg:job':
                    allocatedResources:
                      cpu: 8
                      mem: 6000
                    inputs:
                      alignment_mode: Local
                      allowed_mismatch_number: '0'
                      bowtie_index_archive:
                        class: File
                        metadata:
                          reference_genome: chr20
                        path: /demo/test-data/chr20_bowtie2-2.2.6.tar
                        secondaryFiles: []
                        size: 0
                      coefficient_i_b: coefficient_i_b-string-value
                      coefficient_nceil_b: coefficient_nceil_b-string-value
                      coefficient_scoremin_b: coefficient_scoremin_b-string-value
                      constant_i_a: constant_i_a-string-value
                      constant_nceil_a: constant_nceil_a-string-value
                      constant_scoremin_a: constant_scoremin_a-string-value
                      disable_overlapping_alignments: false
                      disable_unpaired_alignments: false
                      function_i: Linear
                      function_n_ceil: Square-root
                      function_score_min: Linear
                      input_fasta_files: true
                      mates_alignment_orientation: '--rf'
                      paired_aligned_reads: gzip compressed
                      paired_unaligned_reads: raw
                      platform: ABI SOLiD
                      preset_option: Very fast
                      quality_scale: Phred+33
                      read_gap_penalties:
                        - 0
                      read_sequence:
                        - class: File
                          metadata:
                            file_format: fastq
                            paired_end: '1'
                            quality_scale: illumina15
                          path: /demo/test-data/mate1.fastq
                          secondaryFiles: []
                          size: 0
                        - metadata:
                            file_format: fastq
                            paired_end: '2'
                            qual_scale: illumina15
                          path: /demo/test-data/mate2.fastq
                          secondaryFiles: []
                      reference_gap_penalties:
                        - 0
                      sample_id: nn
                      suppress_sam_records: true
                      threads: null
                      unpaired_aligned_reads: gzip compressed
                      unpaired_unaligned_reads: gzip compressed
                  'sbg:latestRevision': 18
                  'sbg:license': GNU General Public License v3.0 only
                  'sbg:links':
                    - id: 'http://bowtie-bio.sourceforge.net/bowtie2/index.shtml'
                      label: Homepage
                    - id: >-
                        http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/
                      label: Download
                    - id: >-
                        http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html
                      label: Publication
                    - id: 'http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml'
                      label: Manual
                    - id: 'https://github.com/BenLangmead/bowtie2'
                      label: Source Code
                  'sbg:modifiedBy': admin
                  'sbg:modifiedOn': 1533302766
                  'sbg:project': admin/sbg-public-data
                  'sbg:projectName': SBG Public Data
                  'sbg:publisher': sbg
                  'sbg:revision': 18
                  'sbg:revisionNotes': |-
                    - changed labels for Function type
                    - changed type to string for coefficients and constants
                  'sbg:revisionsInfo':
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1454426457
                      'sbg:revision': 0
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1454426458
                      'sbg:revision': 1
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1462903844
                      'sbg:revision': 2
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': sevenbridges-pgc
                      'sbg:modifiedOn': 1465231336
                      'sbg:revision': 3
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1471881576
                      'sbg:revision': 4
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1471881576
                      'sbg:revision': 5
                      'sbg:revisionNotes': null
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1472054036
                      'sbg:revision': 6
                      'sbg:revisionNotes': >-
                        Redesigned to accept archive with Bowtie2 index files on
                        the input.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1472208214
                      'sbg:revision': 7
                      'sbg:revisionNotes': Additional information updated.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1472208214
                      'sbg:revision': 8
                      'sbg:revisionNotes': Additional information updated.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1477482292
                      'sbg:revision': 9
                      'sbg:revisionNotes': >-
                        Javascript for reads written to check if metadata is
                        null.
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1497287393
                      'sbg:revision': 10
                      'sbg:revisionNotes': >-
                        expressions for single end reads and scatter mode
                        updated
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1497287393
                      'sbg:revision': 11
                      'sbg:revisionNotes': |-
                        changed Docker file
                        added select output_format option
                        changed expression for handling SE reads
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1497287393
                      'sbg:revision': 12
                      'sbg:revisionNotes': added  piping and pipe status check
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 13
                      'sbg:revisionNotes': Works with fasta files and with only one input file
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 14
                      'sbg:revisionNotes': Works with fasta files and with only one input file
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 15
                      'sbg:revisionNotes': Number of threads added
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501773135
                      'sbg:revision': 16
                      'sbg:revisionNotes': Number of threads added
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1501775358
                      'sbg:revision': 17
                      'sbg:revisionNotes': 'Threads: default 8'
                    - 'sbg:modifiedBy': admin
                      'sbg:modifiedOn': 1533302766
                      'sbg:revision': 18
                      'sbg:revisionNotes': |-
                        - changed labels for Function type
                        - changed type to string for coefficients and constants
                  'sbg:sbgMaintained': false
                  'sbg:toolAuthor': Ben Langmead/John Hopkins University
                  'sbg:toolkit': Bowtie2
                  'sbg:toolkitVersion': 2.2.6
                  'sbg:validationErrors': []
                label: Bowtie2 Aligner
                'sbg:x': 222.015625
                'sbg:y': 160.125
            requirements: []
            'sbg:appVersion':
              - v1.0
              - 'sbg:draft-2'
            'sbg:content_hash': aed84ae83d3a7c6b130e43ea33af7d1377c85e82c6ca358ff534111b7e6b030d7
            'sbg:contributors':
              - rbowen_james
            'sbg:createdBy': rbowen_james
            'sbg:createdOn': 1626440046
            'sbg:id': mwonge/mwtest/hlahd-bowtie-samtools/6
            'sbg:image_url': >-
              https://cavatica.sbgenomics.com/ns/brood/images/mwonge/mwtest/hlahd-bowtie-samtools/6.png
            'sbg:latestRevision': 6
            'sbg:modifiedBy': rbowen_james
            'sbg:modifiedOn': 1627371058
            'sbg:project': mwonge/mwtest
            'sbg:projectName': zcc-cavatica
            'sbg:publisher': sbg
            'sbg:revision': 6
            'sbg:revisionNotes': null
            'sbg:revisionsInfo':
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626440046
                'sbg:revision': 0
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029037
                'sbg:revision': 1
                'sbg:revisionNotes': added outdir
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029530
                'sbg:revision': 2
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029871
                'sbg:revision': 3
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627104164
                'sbg:revision': 4
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627370861
                'sbg:revision': 5
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627371058
                'sbg:revision': 6
                'sbg:revisionNotes': null
            'sbg:sbgMaintained': false
            'sbg:toolAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
            'sbg:validationErrors': []
            'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
          label: HLA-HD Normal DNA
          'sbg:x': 286.875
          'sbg:y': 821.40625
        - id: hla_hd_tumour_rna
          in:
            - id: threads
              default: 2
            - id: minimum_read_length
              default: 0
              source: minimum_read_length
            - id: fastq_reads1
              source: gunzip/unzipped_file
            - id: fastq_reads2
              source: gunzip_1/unzipped_file
            - id: sample_id
              source: sample_id
            - id: output_dir_name
              source: tumour_rna_outdir_name
          out:
            - id: hlahd_results
            - id: hlahd_final_results
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: mwonge/mwtest/hla-hd/30
            baseCommand: []
            inputs:
              - id: threads
                type: int
                inputBinding:
                  position: 1
                  prefix: '-t'
                  shellQuote: false
                label: Threads
                doc: Number of cores used to execute the program.
              - 'sbg:toolDefaultValue': '100'
                id: minimum_read_length
                type: int
                inputBinding:
                  position: 1
                  prefix: '-m'
                  shellQuote: false
                label: Minimum read length
                doc: A read whose length is shorter than this parameter is ignored.
              - id: fastq_reads1
                type: File
                inputBinding:
                  position: 2
                  separate: false
                  shellQuote: false
                label: FASTQ Reads 1
                doc: Paired-end reads 1 in FASTQ format.
                'sbg:fileTypes': FASTQ
              - id: fastq_reads2
                type: File
                inputBinding:
                  position: 2
                  separate: false
                  shellQuote: false
                label: FASTQ Reads 2
                doc: Paired-end reads 2 in FASTQ format.
                'sbg:fileTypes': FASTQ
              - id: sample_id
                type: string
                inputBinding:
                  position: 4
                  separate: false
                  shellQuote: false
                label: Sample ID
                doc: >-
                  Sample ID for the input FASTQs. This will be used as the name
                  of the output directory.
              - id: output_dir_name
                type: string?
            outputs:
              - id: hlahd_results
                doc: Directory containing results of the HLA-HD run.
                label: Output directory
                type: Directory
                outputBinding:
                  glob: |-
                    ${
                        if (!inputs.output_dir_name) {
                            return inputs.sample_id + "_hlahd"
                        } else {
                            return inputs.output_dir_name
                        }
                    }
              - id: hlahd_final_results
                type: File
                outputBinding:
                  glob: |-
                    ${
                        if (!inputs.output_dir_name) {
                            return inputs.sample_id + "_hlahd/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                        } else {
                            return inputs.output_dir_name + "/" + inputs.sample_id + "/result/" + inputs.sample_id + "_final.result.txt"
                        }
                    }
            doc: >-
              ## About HLA-HD

              HLA-HD documentation and release notes can be found
              [here](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/).

              HLA-HD (HLA typing from High-quality Dictionary) can accurately
              determine HLA alleles with 6-digit precision from NGS data (FASTQ
              format). RNA-Seq data can also be applied.


              ## About this CWL tool

              This CWL tool runs HLA-HD version 1.4.0 on paired-end FASTQ files
              to determine HLA type.


              ### Inputs and parameters

              - The input paired-end read files can be from **WGS/WES or
              RNA-seq**.

              - The input paired-end read files must be in FASTQ format (**not
              zipped**).

              - The default minimum read length is 100, however this is often
              too strict. Choose a lower threshold to include more reads.


              ### Output

              - HLA-HD results are output to a directory named using the input
              sample id.

              - The final summary of HLA typing results can be found at the
              following path:
              `<output_dir_name>/result/<sample_id>_final.result.txt`.


              ### Other notes

              - This tool uses the HLA dictionary created from release 3.15.0 of
              the [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA) database.

              - This tool by default uses HLA allele frequency data included
              with the HLA-HD release 1.4.0.
            label: HLA-HD
            arguments:
              - position: 2
                prefix: '-f'
                shellQuote: false
                valueFrom: /app/hlahd.1.4.0/freq_data
              - position: 3
                prefix: ''
                separate: false
                shellQuote: false
                valueFrom: /app/hlahd.1.4.0/HLA_gene.split.txt
              - position: 3
                prefix: ''
                separate: false
                shellQuote: false
                valueFrom: /app/hlahd.1.4.0/dictionary
              - position: 1
                prefix: ''
                separate: false
                shellQuote: false
                valueFrom: hlahd.sh
              - position: 5
                prefix: ''
                valueFrom: ./
              - position: 6
                prefix: ''
                separate: false
                shellQuote: false
                valueFrom: |-
                  ${
                      if (!inputs.output_dir_name) {
                          return "&& mv ./" + inputs.sample_id + " ./" + inputs.sample_id + "_hlahd"
                      } else {
                          return "&& mv ./" + inputs.sample_id + " ./" + inputs.output_dir_name
                      }
                  }
              - position: 0
                prefix: ''
                shellQuote: false
                valueFrom: |-
                  ${
                      if (!inputs.output_dir_name) {
                          return "mkdir ./" + inputs.sample_id + "_hlahd &&"
                      } else {
                          return "mkdir ./" + inputs.output_dir_name + " &&"
                      }
                  }
            requirements:
              - class: ShellCommandRequirement
              - class: DockerRequirement
                dockerPull: pgc-images.sbgenomics.com/rbowen_james/hla-hd
              - class: InlineJavascriptRequirement
            'sbg:appVersion':
              - v1.0
            'sbg:categories':
              - WGS
              - RNA
              - HLA Typing
              - HLA
              - MHC
              - WES (WXS)
            'sbg:content_hash': a4adfcaaafaaa6b9ac47bbe7e09595dae02a8f11b9f16fdaa0cfcb21948bf91bd
            'sbg:contributors':
              - rbowen_james
            'sbg:createdBy': rbowen_james
            'sbg:createdOn': 1622961333
            'sbg:id': mwonge/mwtest/hla-hd/30
            'sbg:image_url': null
            'sbg:latestRevision': 30
            'sbg:modifiedBy': rbowen_james
            'sbg:modifiedOn': 1627370832
            'sbg:project': mwonge/mwtest
            'sbg:projectName': zcc-cavatica
            'sbg:publisher': sbg
            'sbg:revision': 30
            'sbg:revisionNotes': output final result file
            'sbg:revisionsInfo':
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1622961333
                'sbg:revision': 0
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1622961864
                'sbg:revision': 1
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1622962039
                'sbg:revision': 2
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1622962892
                'sbg:revision': 3
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1623026790
                'sbg:revision': 4
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1623027275
                'sbg:revision': 5
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1623043041
                'sbg:revision': 6
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1623043454
                'sbg:revision': 7
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1623044879
                'sbg:revision': 8
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1623046596
                'sbg:revision': 9
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1623048682
                'sbg:revision': 10
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624408084
                'sbg:revision': 11
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624411148
                'sbg:revision': 12
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624580200
                'sbg:revision': 13
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624581085
                'sbg:revision': 14
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624581597
                'sbg:revision': 15
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624583201
                'sbg:revision': 16
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624583825
                'sbg:revision': 17
                'sbg:revisionNotes': ''
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624586750
                'sbg:revision': 18
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624588045
                'sbg:revision': 19
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624593624
                'sbg:revision': 20
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624606248
                'sbg:revision': 21
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624606452
                'sbg:revision': 22
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624608065
                'sbg:revision': 23
                'sbg:revisionNotes': 'Fixed output dir issue, added docs.'
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1624779790
                'sbg:revision': 24
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627028760
                'sbg:revision': 25
                'sbg:revisionNotes': Added output dir option
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627028962
                'sbg:revision': 26
                'sbg:revisionNotes': make outdir
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029487
                'sbg:revision': 27
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627029828
                'sbg:revision': 28
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627104094
                'sbg:revision': 29
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1627370832
                'sbg:revision': 30
                'sbg:revisionNotes': output final result file
            'sbg:sbgMaintained': false
            'sbg:toolAuthor': Shuji Kawaguchi <shuji@genome.med.kyoto-u.ac.jp>
            'sbg:toolkit': HLA-HD
            'sbg:toolkitVersion': 1.4.0
            'sbg:validationErrors': []
            'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
          label: HLA-HD Tumour RNA
          'sbg:x': 659.8301391601562
          'sbg:y': 1088.2109375
        - id: gunzip
          in:
            - id: vcf_gz
              source: tumour_rna_reads_1
          out:
            - id: unzipped_file
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: mwonge/mwtest/gunzip/175
            baseCommand: []
            inputs:
              - id: vcf_gz
                type: File
                label: bgzipped vcf file
            outputs:
              - id: unzipped_file
                label: unzipped file
                type: File
                outputBinding:
                  glob: $(inputs.vcf_gz.nameroot)
            doc: >

              # Samtools, Bcftools, HTSlib


              These apps all utilise a samtools docker (v1.10) in order to
              perform a variety of

              samtools transforms and functions.


              The docker contains (v1.10) of samtools, htslib and bcftools


              ## Documentation


              (Release Page)[http://www.htslib.org/download/]


              - (BCFtools)[http://www.htslib.org/doc/bcftools.html]

              - (bgzip)[http://www.htslib.org/doc/bgzip.html]

              - (HTSfile)[http://www.htslib.org/doc/htsfile.html]

              - (Samtools)[http://www.htslib.org/doc/samtools.html]

              - (Tabix)[http://www.htslib.org/doc/tabix.html]


              ## CWL Tools


              ### Bamname


              Uses samtools view to extract the sample name from the BAM header


              ### Gunzip


              Unzips a \*.gz file that was zipped by bgzip


              ### Tabix


              bgzips and tabix indexes a vcf using htslib tools bgzip and tabix


              ## Dependencies


              Access to the docker:

              `pgc-images.sbgenomics.com/syan/samtools:1.10`


              The docker contains the following:


              - wget

              - htslib (1.10.2)
                - bgzip
                - tabix
              - bcftools (1.10.2)

              - samtools (1.10)


              ## Instance Recommendations


              Notes:


              - Requirements will depend entirely on your inputs, sizes, and
              what tools you want to use

              - The following is just a suggestion


              - **ramMin** 1000

              - **coresMin** 1

              - **tmpdirMin** 10000

              - **outdirMin** 1000
            label: gunzip
            arguments:
              - position: 0
                valueFrom: gunzip $(inputs.vcf_gz.path)
            requirements:
              - class: ResourceRequirement
                ramMin: 100
                coresMin: 1
              - class: DockerRequirement
                dockerPull: 'pgc-images.sbgenomics.com/syan/samtools:production'
              - class: InitialWorkDirRequirement
                listing:
                  - $(inputs.vcf_gz)
              - class: InlineJavascriptRequirement
            'sbg:appVersion':
              - v1.0
            'sbg:content_hash': a83e312d670c3189e9725b5f990a11cbc3af131fd672a6587aed574bdd4849db1
            'sbg:contributors':
              - syan
              - rbowen_james
              - lcui
            'sbg:createdBy': syan
            'sbg:createdOn': 1582514245
            'sbg:id': mwonge/mwtest/gunzip/175
            'sbg:image_url': null
            'sbg:latestRevision': 175
            'sbg:license': MIT License (Expat)
            'sbg:links':
              - id: 'http://www.htslib.org/doc/'
                label: documentation
              - id: 'http://www.htslib.org'
                label: release
              - id: 'https://tldrlegal.com/license/mit-license#summary'
                label: license
            'sbg:modifiedBy': rbowen_james
            'sbg:modifiedOn': 1626442154
            'sbg:project': mwonge/mwtest
            'sbg:projectName': zcc-cavatica
            'sbg:publisher': sbg
            'sbg:revision': 175
            'sbg:revisionNotes': Input any file type
            'sbg:revisionsInfo':
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1582514245
                'sbg:revision': 0
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1582517956
                'sbg:revision': 1
                'sbg:revisionNotes': 'WORKING: small bug fix'
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599194096
                'sbg:revision': 2
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: vcfanno/gunzip.cwl
                  commit: 265529d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599543651
                'sbg:revision': 3
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599544004
                'sbg:revision': 4
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599544453
                'sbg:revision': 5
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599546459
                'sbg:revision': 6
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599546679
                'sbg:revision': 7
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599547207
                'sbg:revision': 8
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599547948
                'sbg:revision': 9
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599548541
                'sbg:revision': 10
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599549135
                'sbg:revision': 11
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625191
                'sbg:revision': 12
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 72c8c27
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625459
                'sbg:revision': 13
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625609
                'sbg:revision': 14
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625695
                'sbg:revision': 15
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625768
                'sbg:revision': 16
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600386605
                'sbg:revision': 17
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9f6ef41
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600386840
                'sbg:revision': 18
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9f6ef41
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600386989
                'sbg:revision': 19
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9a3adea
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600387750
                'sbg:revision': 20
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9a3adea
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600905338
                'sbg:revision': 21
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 43ce2ce
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600909753
                'sbg:revision': 22
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600909849
                'sbg:revision': 23
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600910147
                'sbg:revision': 24
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913171
                'sbg:revision': 25
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913382
                'sbg:revision': 26
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913550
                'sbg:revision': 27
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913828
                'sbg:revision': 28
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600914819
                'sbg:revision': 29
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600915127
                'sbg:revision': 30
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600915689
                'sbg:revision': 31
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600916424
                'sbg:revision': 32
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600916879
                'sbg:revision': 33
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600917077
                'sbg:revision': 34
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600917172
                'sbg:revision': 35
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600918777
                'sbg:revision': 36
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600922924
                'sbg:revision': 37
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600923324
                'sbg:revision': 38
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600924209
                'sbg:revision': 39
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601275270
                'sbg:revision': 40
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601357602
                'sbg:revision': 41
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601426727
                'sbg:revision': 42
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601427710
                'sbg:revision': 43
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601429396
                'sbg:revision': 44
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: caf3daa
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601429820
                'sbg:revision': 45
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b74fa94
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601430189
                'sbg:revision': 46
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8e036b4
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601430397
                'sbg:revision': 47
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9e6732f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601430532
                'sbg:revision': 48
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 09ebe6b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601431178
                'sbg:revision': 49
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c9cad2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601431458
                'sbg:revision': 50
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 47f3f1d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601512587
                'sbg:revision': 51
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e443cf6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601523880
                'sbg:revision': 52
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: c53e03d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601533362
                'sbg:revision': 53
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 1059db7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601597551
                'sbg:revision': 54
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a164da8
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601598695
                'sbg:revision': 55
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601598973
                'sbg:revision': 56
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601937742
                'sbg:revision': 57
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 39cf8a0
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601952251
                'sbg:revision': 58
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b4fc80e
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601956650
                'sbg:revision': 59
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ee9d914
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601958278
                'sbg:revision': 60
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cca9253
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601958378
                'sbg:revision': 61
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4bf226f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601958766
                'sbg:revision': 62
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 762b48c
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601959478
                'sbg:revision': 63
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4eedb5d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601959656
                'sbg:revision': 64
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b814705
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601960404
                'sbg:revision': 65
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 2fdfb35
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601961439
                'sbg:revision': 66
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8ea1f28
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601962570
                'sbg:revision': 67
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f9f7d03
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601964984
                'sbg:revision': 68
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 5602fe0
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602026318
                'sbg:revision': 69
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3fdd7fa
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602026419
                'sbg:revision': 70
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a5b345f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602030805
                'sbg:revision': 71
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f24c4f6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602031143
                'sbg:revision': 72
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 95f2bc2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602031769
                'sbg:revision': 73
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 99fed98
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602035405
                'sbg:revision': 74
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9da71e7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602039870
                'sbg:revision': 75
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f6ffa7c
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602044768
                'sbg:revision': 76
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b023f61
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602113318
                'sbg:revision': 77
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 33b12f2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602128487
                'sbg:revision': 78
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 316d9e9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602136128
                'sbg:revision': 79
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4a39cbc
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602211323
                'sbg:revision': 80
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 81b0723
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602556905
                'sbg:revision': 81
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7fd45ff
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602561458
                'sbg:revision': 82
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: af351ea
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603234936
                'sbg:revision': 83
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cfa7cd7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603235081
                'sbg:revision': 84
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ff066af
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603245673
                'sbg:revision': 85
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8fac7db
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603318097
                'sbg:revision': 86
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a167b08
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603344742
                'sbg:revision': 87
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 895c8ce
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603421167
                'sbg:revision': 88
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cf3371f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603421465
                'sbg:revision': 89
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8c1512c
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603673485
                'sbg:revision': 90
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e8f60ab
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603678826
                'sbg:revision': 91
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: d6e2f0b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603679109
                'sbg:revision': 92
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e952c85
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603681128
                'sbg:revision': 93
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8d8ea57
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603685680
                'sbg:revision': 94
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c3335a
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603857129
                'sbg:revision': 95
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e8f85c0
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603865020
                'sbg:revision': 96
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7003aab
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603942516
                'sbg:revision': 97
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: fd0b5ae
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603947521
                'sbg:revision': 98
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3e69899
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603949308
                'sbg:revision': 99
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: fd7fbb5
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603952569
                'sbg:revision': 100
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 775f506
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604277332
                'sbg:revision': 101
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c06064
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604277488
                'sbg:revision': 102
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c06064
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604280961
                'sbg:revision': 103
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6d53c85
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604365472
                'sbg:revision': 104
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 61b1a85
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604552841
                'sbg:revision': 105
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e63e5c4
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604964304
                'sbg:revision': 106
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 824c601
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605049096
                'sbg:revision': 107
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7862097
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605082174
                'sbg:revision': 108
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: bc4fdf8
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605085488
                'sbg:revision': 109
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7d7ee70
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605133898
                'sbg:revision': 110
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7d7ee70
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605484203
                'sbg:revision': 111
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 541cebe
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605572338
                'sbg:revision': 112
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 47dbb50
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605576089
                'sbg:revision': 113
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 38786f7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605665284
                'sbg:revision': 114
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7718959
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605665517
                'sbg:revision': 115
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7718959
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1605678242
                'sbg:revision': 116
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b47363e
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1605678518
                'sbg:revision': 117
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6da6c0f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605681519
                'sbg:revision': 118
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8b652ac
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605741716
                'sbg:revision': 119
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8b652ac
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605742265
                'sbg:revision': 120
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 0b4cf01
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1606864102
                'sbg:revision': 121
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 85328e2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1606864393
                'sbg:revision': 122
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a0d4aff
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1606971686
                'sbg:revision': 123
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 684d1ed
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607035060
                'sbg:revision': 124
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6d018e4
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607040615
                'sbg:revision': 125
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a98fa70
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607049957
                'sbg:revision': 126
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 578b0de
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607316131
                'sbg:revision': 127
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ef0c2a5
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607319792
                'sbg:revision': 128
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7a8240d
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1607902135
                'sbg:revision': 129
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: d59f4db
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607906092
                'sbg:revision': 130
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 5836465
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607906828
                'sbg:revision': 131
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ef6ed2b
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1607907895
                'sbg:revision': 132
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7f0bdb4
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1607911308
                'sbg:revision': 133
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: d04ae7d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607914411
                'sbg:revision': 134
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 18dfb49
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1608075467
                'sbg:revision': 135
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 18dfb49
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1608077674
                'sbg:revision': 136
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 56beeb1
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1610424226
                'sbg:revision': 137
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ecdc528
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1611035956
                'sbg:revision': 138
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 373907c
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1611098612
                'sbg:revision': 139
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9e52f26
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613085601
                'sbg:revision': 140
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6e0b99b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613525101
                'sbg:revision': 141
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3e6f0d4
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613624658
                'sbg:revision': 142
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 773b440
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613627544
                'sbg:revision': 143
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e5bdc2b
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613687495
                'sbg:revision': 144
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 32bb0c3
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613704590
                'sbg:revision': 145
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e670ba9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613704909
                'sbg:revision': 146
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e670ba9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613705150
                'sbg:revision': 147
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e6a7d65
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613705290
                'sbg:revision': 148
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e6a7d65
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613705408
                'sbg:revision': 149
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: befd8b1
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613706024
                'sbg:revision': 150
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 48965cf
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613718345
                'sbg:revision': 151
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 87d785c
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613968202
                'sbg:revision': 152
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 57998bf
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613969225
                'sbg:revision': 153
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: faf8958
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1614034790
                'sbg:revision': 154
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: abe98a3
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614048486
                'sbg:revision': 155
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e706c86
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614049815
                'sbg:revision': 156
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ddf9af6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614298550
                'sbg:revision': 157
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ae3a543
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614304902
                'sbg:revision': 158
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b48195d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614305593
                'sbg:revision': 159
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b1bc92f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614557665
                'sbg:revision': 160
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: c6d3eec
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614655443
                'sbg:revision': 161
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4edb2c7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614657425
                'sbg:revision': 162
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e8821d7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614728741
                'sbg:revision': 163
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8aaf9f9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615166236
                'sbg:revision': 164
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b6911aa
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615177336
                'sbg:revision': 165
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cf2acf6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615256213
                'sbg:revision': 166
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9637fb6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615335999
                'sbg:revision': 167
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 19905f3
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616024279
                'sbg:revision': 168
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 489075b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616369899
                'sbg:revision': 169
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: fcbfa76
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616470580
                'sbg:revision': 170
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a46e093
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616538078
                'sbg:revision': 171
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4ee52e6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1617254163
                'sbg:revision': 172
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4726466
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1618536828
                'sbg:revision': 173
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f3f5fd6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1618553641
                'sbg:revision': 174
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 5615bb2
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626442154
                'sbg:revision': 175
                'sbg:revisionNotes': Input any file type
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
          label: gunzip
          'sbg:x': 286.875
          'sbg:y': 1083.8671875
        - id: gunzip_1
          in:
            - id: vcf_gz
              source: tumour_rna_reads_2
          out:
            - id: unzipped_file
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              sbg: 'https://sevenbridges.com'
            id: mwonge/mwtest/gunzip/175
            baseCommand: []
            inputs:
              - id: vcf_gz
                type: File
                label: bgzipped vcf file
            outputs:
              - id: unzipped_file
                label: unzipped file
                type: File
                outputBinding:
                  glob: $(inputs.vcf_gz.nameroot)
            doc: >

              # Samtools, Bcftools, HTSlib


              These apps all utilise a samtools docker (v1.10) in order to
              perform a variety of

              samtools transforms and functions.


              The docker contains (v1.10) of samtools, htslib and bcftools


              ## Documentation


              (Release Page)[http://www.htslib.org/download/]


              - (BCFtools)[http://www.htslib.org/doc/bcftools.html]

              - (bgzip)[http://www.htslib.org/doc/bgzip.html]

              - (HTSfile)[http://www.htslib.org/doc/htsfile.html]

              - (Samtools)[http://www.htslib.org/doc/samtools.html]

              - (Tabix)[http://www.htslib.org/doc/tabix.html]


              ## CWL Tools


              ### Bamname


              Uses samtools view to extract the sample name from the BAM header


              ### Gunzip


              Unzips a \*.gz file that was zipped by bgzip


              ### Tabix


              bgzips and tabix indexes a vcf using htslib tools bgzip and tabix


              ## Dependencies


              Access to the docker:

              `pgc-images.sbgenomics.com/syan/samtools:1.10`


              The docker contains the following:


              - wget

              - htslib (1.10.2)
                - bgzip
                - tabix
              - bcftools (1.10.2)

              - samtools (1.10)


              ## Instance Recommendations


              Notes:


              - Requirements will depend entirely on your inputs, sizes, and
              what tools you want to use

              - The following is just a suggestion


              - **ramMin** 1000

              - **coresMin** 1

              - **tmpdirMin** 10000

              - **outdirMin** 1000
            label: gunzip
            arguments:
              - position: 0
                valueFrom: gunzip $(inputs.vcf_gz.path)
            requirements:
              - class: ResourceRequirement
                ramMin: 100
                coresMin: 1
              - class: DockerRequirement
                dockerPull: 'pgc-images.sbgenomics.com/syan/samtools:production'
              - class: InitialWorkDirRequirement
                listing:
                  - $(inputs.vcf_gz)
              - class: InlineJavascriptRequirement
            'sbg:appVersion':
              - v1.0
            'sbg:content_hash': a83e312d670c3189e9725b5f990a11cbc3af131fd672a6587aed574bdd4849db1
            'sbg:contributors':
              - syan
              - rbowen_james
              - lcui
            'sbg:createdBy': syan
            'sbg:createdOn': 1582514245
            'sbg:id': mwonge/mwtest/gunzip/175
            'sbg:image_url': null
            'sbg:latestRevision': 175
            'sbg:license': MIT License (Expat)
            'sbg:links':
              - id: 'http://www.htslib.org/doc/'
                label: documentation
              - id: 'http://www.htslib.org'
                label: release
              - id: 'https://tldrlegal.com/license/mit-license#summary'
                label: license
            'sbg:modifiedBy': rbowen_james
            'sbg:modifiedOn': 1626442154
            'sbg:project': mwonge/mwtest
            'sbg:projectName': zcc-cavatica
            'sbg:publisher': sbg
            'sbg:revision': 175
            'sbg:revisionNotes': Input any file type
            'sbg:revisionsInfo':
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1582514245
                'sbg:revision': 0
                'sbg:revisionNotes': null
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1582517956
                'sbg:revision': 1
                'sbg:revisionNotes': 'WORKING: small bug fix'
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599194096
                'sbg:revision': 2
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: vcfanno/gunzip.cwl
                  commit: 265529d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599543651
                'sbg:revision': 3
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599544004
                'sbg:revision': 4
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599544453
                'sbg:revision': 5
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599546459
                'sbg:revision': 6
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599546679
                'sbg:revision': 7
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599547207
                'sbg:revision': 8
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599547948
                'sbg:revision': 9
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599548541
                'sbg:revision': 10
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599549135
                'sbg:revision': 11
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625191
                'sbg:revision': 12
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 72c8c27
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625459
                'sbg:revision': 13
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625609
                'sbg:revision': 14
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625695
                'sbg:revision': 15
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1599625768
                'sbg:revision': 16
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: c1d2c23
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600386605
                'sbg:revision': 17
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9f6ef41
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600386840
                'sbg:revision': 18
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9f6ef41
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600386989
                'sbg:revision': 19
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9a3adea
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600387750
                'sbg:revision': 20
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 9a3adea
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600905338
                'sbg:revision': 21
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: samtools/gunzip.cwl
                  commit: 43ce2ce
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600909753
                'sbg:revision': 22
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600909849
                'sbg:revision': 23
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600910147
                'sbg:revision': 24
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913171
                'sbg:revision': 25
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913382
                'sbg:revision': 26
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913550
                'sbg:revision': 27
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600913828
                'sbg:revision': 28
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600914819
                'sbg:revision': 29
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600915127
                'sbg:revision': 30
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600915689
                'sbg:revision': 31
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600916424
                'sbg:revision': 32
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600916879
                'sbg:revision': 33
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600917077
                'sbg:revision': 34
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600917172
                'sbg:revision': 35
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600918777
                'sbg:revision': 36
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600922924
                'sbg:revision': 37
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600923324
                'sbg:revision': 38
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1600924209
                'sbg:revision': 39
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601275270
                'sbg:revision': 40
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601357602
                'sbg:revision': 41
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601426727
                'sbg:revision': 42
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601427710
                'sbg:revision': 43
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: 
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601429396
                'sbg:revision': 44
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: caf3daa
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601429820
                'sbg:revision': 45
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b74fa94
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601430189
                'sbg:revision': 46
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8e036b4
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601430397
                'sbg:revision': 47
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9e6732f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601430532
                'sbg:revision': 48
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 09ebe6b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601431178
                'sbg:revision': 49
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c9cad2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601431458
                'sbg:revision': 50
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 47f3f1d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601512587
                'sbg:revision': 51
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e443cf6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601523880
                'sbg:revision': 52
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: c53e03d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601533362
                'sbg:revision': 53
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 1059db7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601597551
                'sbg:revision': 54
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a164da8
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601598695
                'sbg:revision': 55
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601598973
                'sbg:revision': 56
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: (uncommitted file)
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601937742
                'sbg:revision': 57
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 39cf8a0
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601952251
                'sbg:revision': 58
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b4fc80e
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601956650
                'sbg:revision': 59
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ee9d914
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601958278
                'sbg:revision': 60
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cca9253
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601958378
                'sbg:revision': 61
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4bf226f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601958766
                'sbg:revision': 62
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 762b48c
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601959478
                'sbg:revision': 63
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4eedb5d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1601959656
                'sbg:revision': 64
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b814705
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601960404
                'sbg:revision': 65
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 2fdfb35
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601961439
                'sbg:revision': 66
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8ea1f28
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601962570
                'sbg:revision': 67
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f9f7d03
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1601964984
                'sbg:revision': 68
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 5602fe0
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602026318
                'sbg:revision': 69
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3fdd7fa
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602026419
                'sbg:revision': 70
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a5b345f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602030805
                'sbg:revision': 71
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f24c4f6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602031143
                'sbg:revision': 72
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 95f2bc2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602031769
                'sbg:revision': 73
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 99fed98
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602035405
                'sbg:revision': 74
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9da71e7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602039870
                'sbg:revision': 75
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f6ffa7c
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602044768
                'sbg:revision': 76
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b023f61
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602113318
                'sbg:revision': 77
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 33b12f2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602128487
                'sbg:revision': 78
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 316d9e9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602136128
                'sbg:revision': 79
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4a39cbc
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602211323
                'sbg:revision': 80
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 81b0723
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602556905
                'sbg:revision': 81
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7fd45ff
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1602561458
                'sbg:revision': 82
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: af351ea
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603234936
                'sbg:revision': 83
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cfa7cd7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603235081
                'sbg:revision': 84
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ff066af
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603245673
                'sbg:revision': 85
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8fac7db
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603318097
                'sbg:revision': 86
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a167b08
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603344742
                'sbg:revision': 87
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 895c8ce
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603421167
                'sbg:revision': 88
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cf3371f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603421465
                'sbg:revision': 89
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8c1512c
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603673485
                'sbg:revision': 90
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e8f60ab
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603678826
                'sbg:revision': 91
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: d6e2f0b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603679109
                'sbg:revision': 92
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e952c85
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603681128
                'sbg:revision': 93
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8d8ea57
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603685680
                'sbg:revision': 94
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c3335a
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603857129
                'sbg:revision': 95
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e8f85c0
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603865020
                'sbg:revision': 96
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7003aab
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603942516
                'sbg:revision': 97
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: fd0b5ae
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603947521
                'sbg:revision': 98
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3e69899
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603949308
                'sbg:revision': 99
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: fd7fbb5
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1603952569
                'sbg:revision': 100
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 775f506
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604277332
                'sbg:revision': 101
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c06064
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604277488
                'sbg:revision': 102
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3c06064
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604280961
                'sbg:revision': 103
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6d53c85
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604365472
                'sbg:revision': 104
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 61b1a85
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604552841
                'sbg:revision': 105
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e63e5c4
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1604964304
                'sbg:revision': 106
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 824c601
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605049096
                'sbg:revision': 107
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7862097
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605082174
                'sbg:revision': 108
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: bc4fdf8
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605085488
                'sbg:revision': 109
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7d7ee70
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605133898
                'sbg:revision': 110
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7d7ee70
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605484203
                'sbg:revision': 111
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 541cebe
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605572338
                'sbg:revision': 112
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 47dbb50
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605576089
                'sbg:revision': 113
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 38786f7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605665284
                'sbg:revision': 114
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7718959
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605665517
                'sbg:revision': 115
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7718959
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1605678242
                'sbg:revision': 116
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b47363e
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1605678518
                'sbg:revision': 117
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6da6c0f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605681519
                'sbg:revision': 118
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8b652ac
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605741716
                'sbg:revision': 119
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8b652ac
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1605742265
                'sbg:revision': 120
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 0b4cf01
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1606864102
                'sbg:revision': 121
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 85328e2
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1606864393
                'sbg:revision': 122
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a0d4aff
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1606971686
                'sbg:revision': 123
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 684d1ed
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607035060
                'sbg:revision': 124
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6d018e4
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607040615
                'sbg:revision': 125
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a98fa70
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607049957
                'sbg:revision': 126
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 578b0de
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607316131
                'sbg:revision': 127
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ef0c2a5
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607319792
                'sbg:revision': 128
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7a8240d
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1607902135
                'sbg:revision': 129
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: d59f4db
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607906092
                'sbg:revision': 130
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 5836465
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607906828
                'sbg:revision': 131
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ef6ed2b
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1607907895
                'sbg:revision': 132
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 7f0bdb4
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1607911308
                'sbg:revision': 133
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: d04ae7d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1607914411
                'sbg:revision': 134
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 18dfb49
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1608075467
                'sbg:revision': 135
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 18dfb49
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1608077674
                'sbg:revision': 136
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 56beeb1
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1610424226
                'sbg:revision': 137
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ecdc528
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1611035956
                'sbg:revision': 138
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 373907c
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1611098612
                'sbg:revision': 139
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9e52f26
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613085601
                'sbg:revision': 140
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 6e0b99b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613525101
                'sbg:revision': 141
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 3e6f0d4
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613624658
                'sbg:revision': 142
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 773b440
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613627544
                'sbg:revision': 143
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e5bdc2b
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613687495
                'sbg:revision': 144
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 32bb0c3
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613704590
                'sbg:revision': 145
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e670ba9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613704909
                'sbg:revision': 146
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e670ba9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613705150
                'sbg:revision': 147
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e6a7d65
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613705290
                'sbg:revision': 148
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e6a7d65
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613705408
                'sbg:revision': 149
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: befd8b1
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1613706024
                'sbg:revision': 150
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 48965cf
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613718345
                'sbg:revision': 151
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 87d785c
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613968202
                'sbg:revision': 152
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 57998bf
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1613969225
                'sbg:revision': 153
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: faf8958
              - 'sbg:modifiedBy': lcui
                'sbg:modifiedOn': 1614034790
                'sbg:revision': 154
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.09.13. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: abe98a3
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614048486
                'sbg:revision': 155
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.10.05. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e706c86
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614049815
                'sbg:revision': 156
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ddf9af6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614298550
                'sbg:revision': 157
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: ae3a543
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614304902
                'sbg:revision': 158
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b48195d
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614305593
                'sbg:revision': 159
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b1bc92f
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614557665
                'sbg:revision': 160
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: c6d3eec
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614655443
                'sbg:revision': 161
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4edb2c7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614657425
                'sbg:revision': 162
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: e8821d7
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1614728741
                'sbg:revision': 163
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 8aaf9f9
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615166236
                'sbg:revision': 164
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: b6911aa
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615177336
                'sbg:revision': 165
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: cf2acf6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615256213
                'sbg:revision': 166
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 9637fb6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1615335999
                'sbg:revision': 167
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 19905f3
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616024279
                'sbg:revision': 168
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 489075b
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616369899
                'sbg:revision': 169
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: fcbfa76
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616470580
                'sbg:revision': 170
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: a46e093
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1616538078
                'sbg:revision': 171
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4ee52e6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1617254163
                'sbg:revision': 172
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 4726466
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1618536828
                'sbg:revision': 173
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: f3f5fd6
              - 'sbg:modifiedBy': syan
                'sbg:modifiedOn': 1618553641
                'sbg:revision': 174
                'sbg:revisionNotes': |-
                  Uploaded using sbpack v2020.06.18. 
                  Source: 
                  repo: git@bitbucket.org:cciacb/cwl.git
                  file: tools/samtools/gunzip.cwl
                  commit: 5615bb2
              - 'sbg:modifiedBy': rbowen_james
                'sbg:modifiedOn': 1626442154
                'sbg:revision': 175
                'sbg:revisionNotes': Input any file type
            'sbg:sbgMaintained': false
            'sbg:validationErrors': []
          label: gunzip
          'sbg:x': 286.875
          'sbg:y': 977.1328125
      requirements:
        - class: SubworkflowFeatureRequirement
      'sbg:appVersion':
        - v1.0
        - 'sbg:draft-2'
      'sbg:content_hash': ab7ff4b2b5797ba8d0dfb091f0fe138216e553f7cea6b738da0a2885fe71cb1bb
      'sbg:contributors':
        - rbowen_james
      'sbg:createdBy': rbowen_james
      'sbg:createdOn': 1626440450
      'sbg:id': mwonge/mwtest/hlatyping-bowtie-samtools/6
      'sbg:image_url': >-
        https://cavatica.sbgenomics.com/ns/brood/images/mwonge/mwtest/hlatyping-bowtie-samtools/6.png
      'sbg:latestRevision': 6
      'sbg:modifiedBy': rbowen_james
      'sbg:modifiedOn': 1627371234
      'sbg:project': mwonge/mwtest
      'sbg:projectName': zcc-cavatica
      'sbg:publisher': sbg
      'sbg:revision': 6
      'sbg:revisionNotes': null
      'sbg:revisionsInfo':
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1626440450
          'sbg:revision': 0
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1626442245
          'sbg:revision': 1
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627029420
          'sbg:revision': 2
          'sbg:revisionNotes': Added outdirs
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627029679
          'sbg:revision': 3
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627030104
          'sbg:revision': 4
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627104343
          'sbg:revision': 5
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627371234
          'sbg:revision': 6
          'sbg:revisionNotes': null
      'sbg:sbgMaintained': false
      'sbg:validationErrors': []
    label: hlatyping_bowtie_samtools
    'sbg:x': 255.421875
    'sbg:y': 1160.546875
  - id: hlahd_consensus
    in:
      - id: tumour_rna
        source: hlatyping_bowtie_samtools/tumour_rna_hla-hd_final
      - id: tumour_dna
        source: hlatyping_bowtie_samtools/tumour_dna_hla-hd_final_result
      - id: normal_dna
        source: hlatyping_bowtie_samtools/normal_dna_hla-hd_final_result
      - id: sample_name
        source: sample_id
    out:
      - id: consensus_json
      - id: consensus_txt
      - id: consensus_string
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        sbg: 'https://sevenbridges.com'
      id: mwonge/mwtest/hlahd-consensus/2
      baseCommand:
        - python3
        - /app/hlahd_consensus_parser.py
      inputs:
        - id: tumour_rna
          type: File
          inputBinding:
            position: 0
            prefix: '-trna'
            shellQuote: false
          'sbg:fileTypes': txt
        - id: tumour_dna
          type: File
          inputBinding:
            position: 0
            prefix: '-tdna'
            shellQuote: false
          'sbg:fileTypes': txt
        - id: normal_dna
          type: File
          inputBinding:
            position: 0
            prefix: '-ndna'
            shellQuote: false
          'sbg:fileTypes': txt
        - id: sample_name
          type: string
          inputBinding:
            position: 1
            separate: false
            shellQuote: false
      outputs:
        - id: consensus_json
          type: File
          outputBinding:
            glob: '*.consensus.json'
        - id: consensus_txt
          type: File
          outputBinding:
            glob: '*.consensus.txt'
        - id: consensus_string
          type: string?
          outputBinding:
            glob: hla
      label: hlahd-consensus
      arguments:
        - position: 2
          prefix: ''
          separate: false
          shellQuote: false
          valueFrom: '> hla'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'rachelbj/hlahd-consensus:1.0'
      'sbg:appVersion':
        - v1.0
      'sbg:content_hash': ae25fd94df539081b78cdfb0a88227603591753d746b324b1cae8047e917f8a02
      'sbg:contributors':
        - rbowen_james
      'sbg:createdBy': rbowen_james
      'sbg:createdOn': 1627368908
      'sbg:id': mwonge/mwtest/hlahd-consensus/2
      'sbg:image_url': null
      'sbg:latestRevision': 2
      'sbg:modifiedBy': rbowen_james
      'sbg:modifiedOn': 1627370220
      'sbg:project': mwonge/mwtest
      'sbg:projectName': zcc-cavatica
      'sbg:publisher': sbg
      'sbg:revision': 2
      'sbg:revisionNotes': null
      'sbg:revisionsInfo':
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627368908
          'sbg:revision': 0
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627369732
          'sbg:revision': 1
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1627370220
          'sbg:revision': 2
          'sbg:revisionNotes': null
      'sbg:sbgMaintained': false
      'sbg:validationErrors': []
    label: hlahd-consensus
    'sbg:x': 952.7314453125
    'sbg:y': 1226.3203125
requirements:
  - class: SubworkflowFeatureRequirement
