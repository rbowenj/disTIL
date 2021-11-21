class: Workflow
cwlVersion: v1.0
id: distil
doc: >-
  # disTIL

  This workflow runs the disTIL immunoprofiling toolkit with:

  - Three-sample consensus HLA typing

  - Somatic VCF annotation with: VEP, gene expression, transcript expression,
  DNA BAM readcounts, RNA BAM readcounts

  - Fusion annotation with AGfusion

  - pVACseq and pVACfuse neoepitope prediction

  - 


  ## Three-Sample Consensus HLA Typing

  This version of
label: distil
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: fusion_tsv
    'sbg:fileTypes': TSV
    type: File
    label: Fusions TSV
    'sbg:x': 281.765625
    'sbg:y': 1461.140625
  - id: gene_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Gene Expression TSV
    'sbg:x': 0
    'sbg:y': 1349.140625
  - id: input_bam_dna
    'sbg:fileTypes': BAM
    type: File
    label: DNA BAM
    'sbg:x': 281.765625
    'sbg:y': 1205.203125
  - id: input_bam_rna
    'sbg:fileTypes': BAM
    type: File
    label: RNA BAM
    'sbg:x': 281.765625
    'sbg:y': 1098.234375
  - id: input_vcf
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    label: VCF
    'sbg:x': 0
    'sbg:y': 1242.171875
  - id: ref_genome_dna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: DNA Reference Genome
    'sbg:x': 0
    'sbg:y': 1028.234375
  - id: ref_genome_rna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: RNA Reference Genome
    'sbg:x': 0
    'sbg:y': 921.265625
  - id: transcript_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Transcript Expression File
    'sbg:x': 0
    'sbg:y': 279.453125
  - id: vep_cache
    'sbg:fileTypes': TAR.GZ
    type: File
    label: VEP Cache
    'sbg:x': 0
    'sbg:y': 172.484375
  - id: vep_plugin_files
    'sbg:fileTypes': PM
    type: 'File[]'
    label: VEP Plugin Files
    'sbg:x': 0
    'sbg:y': 65.515625
  - id: bowtie2_index
    'sbg:fileTypes': TAR
    type: File
    label: HLA Bowtie2 Index Archive
    'sbg:x': 0
    'sbg:y': 1777.015625
  - id: patient_id
    type: string
    label: Patient ID
    'sbg:x': 0
    'sbg:y': 1135.203125
  - id: sample_3_read2_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 3 Read 2 Sequences
    'sbg:x': 0
    'sbg:y': 386.421875
  - id: sample_3_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 3 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 493.390625
  - id: sample_2_read2_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 2 Read 2 Sequences
    'sbg:x': 0
    'sbg:y': 600.359375
  - id: sample_2_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 2 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 707.328125
  - id: sample_1_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 1 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 814.296875
  - id: bowtie2_index_prefix
    type: string
    label: HLA Bowtie2 Index Prefix
    'sbg:x': 0
    'sbg:y': 1670.046875
  - id: gene_col
    type: string
    label: Gene ID Column Name
    'sbg:x': 0
    'sbg:y': 1456.109375
  - id: expr_col
    type: string
    label: Gene Expression Column Name
    'sbg:x': 0
    'sbg:y': 1563.078125
  - id: phased_proximal_variants_vcf
    type: File?
    label: Phased Proximal Variants VCF
    'sbg:x': 281.765625
    'sbg:y': 842.296875
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
  - id: n_junction_reads
    type: int?
    'sbg:exposed': true
  - id: n_spanning_frags
    type: int?
    'sbg:exposed': true
outputs:
  - id: pvacseq_predictions
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacseq_predictions
    type: Directory
    label: pVACseq Predictions
    'sbg:x': 1597.5888671875
    'sbg:y': 439.90625
  - id: pvacfuse_predictions
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacfuse_predictions
    type: Directory?
    label: pVACfuse Predictions
    'sbg:x': 1597.5888671875
    'sbg:y': 974.75
  - id: agfusion_output
    outputSource:
      - pvactools_neoepitope_prediction_1/agfusion_output
    type: Directory?
    label: AGFusion Output Directory
    'sbg:x': 1597.5888671875
    'sbg:y': 1509.59375
  - id: variants
    outputSource:
      - tmb/variants
    type: File
    label: Canonical Missense Coding Variants
    'sbg:x': 862.2869873046875
    'sbg:y': 0
  - id: tmb_score
    outputSource:
      - tmb/tmb
    type: File
    label: TMB Score
    'sbg:x': 862.2869873046875
    'sbg:y': 106.96875
  - id: ipass_score_file
    outputSource:
      - ipass_patient/ipass_score_file
    type: File
    label: IPASS Score File
    'sbg:x': 862.2869873046875
    'sbg:y': 1307.6875
  - id: quanitseq_deconv
    outputSource:
      - immunedeconv_1/quanitseq_deconv
    'sbg:fileTypes': TSV
    type: File
    label: quanTIseq Deconvolution Results
    'sbg:x': 862.2869873046875
    'sbg:y': 855.75
  - id: epic_deconv
    outputSource:
      - immunedeconv_1/epic_deconv
    'sbg:fileTypes': TSV
    type: File
    label: EPIC Deconvolution Results
    'sbg:x': 862.2869873046875
    'sbg:y': 1414.65625
  - id: sample_1_results
    outputSource:
      - three_sample_hlatyping/tumour_rna_results
    type: Directory
    label: HLA-HD Results for Sample 1
    'sbg:x': 862.2869873046875
    'sbg:y': 748.78125
  - id: sample_3_results
    outputSource:
      - three_sample_hlatyping/tumour_dna_results
    type: Directory
    label: HLA-HD Results for Sample 3
    'sbg:x': 862.2869873046875
    'sbg:y': 534.84375
  - id: sample_2_results
    outputSource:
      - three_sample_hlatyping/normal_dna_results
    type: Directory
    label: HLA-HD Results for Sample 2
    'sbg:x': 862.2869873046875
    'sbg:y': 641.8125
  - id: sample3_json
    outputSource:
      - three_sample_hlatyping/sample3_json
    'sbg:fileTypes': JSON
    type: File?
    label: Sample 3 HLA-HD Results JSON
    'sbg:x': 862.2869873046875
    'sbg:y': 213.9375
  - id: sample2_json_1
    outputSource:
      - three_sample_hlatyping/sample2_json
    'sbg:fileTypes': JSON
    type: File
    label: Sample 2 HLA-HD Results JSON
    'sbg:x': 862.2869873046875
    'sbg:y': 320.90625
  - id: sample1_json_1
    outputSource:
      - three_sample_hlatyping/sample1_json
    'sbg:fileTypes': JSON
    type: File
    label: Sample 1 HLA-HD Results JSON
    'sbg:x': 862.2869873046875
    'sbg:y': 427.875
  - id: consensus_txt
    outputSource:
      - three_sample_hlatyping/consensus_txt
    'sbg:fileTypes': TXT
    type: File
    label: HLA Consensus Text File
    'sbg:x': 862.2869873046875
    'sbg:y': 1521.625
  - id: clin_sig_txt
    outputSource:
      - three_sample_hlatyping/clin_sig_txt
    'sbg:fileTypes': TXT
    type: File
    label: Cinically Significant HLA Consensus Text File
    'sbg:x': 862.2869873046875
    'sbg:y': 1735.5625
  - id: clin_sig_json
    outputSource:
      - three_sample_hlatyping/clin_sig_json
    'sbg:fileTypes': JSON
    type: File
    label: Cinically Significant HLA Consensus JSON
    'sbg:x': 862.2869873046875
    'sbg:y': 1842.53125
  - id: consensus_json
    outputSource:
      - three_sample_hlatyping/consensus_json
    'sbg:fileTypes': JSON
    type: File
    label: HLA Consensus JSON
    'sbg:x': 862.2869873046875
    'sbg:y': 1628.59375
  - id: pvacseq_mhc_ii_aggregated_report_1
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacseq_mhc_ii_aggregated_report
    type: File?
    label: pVACseq MHC-II Aggregated Report
    'sbg:x': 1597.5888671875
    'sbg:y': 546.875
  - id: pvacseq_mhc_i_aggregated_report_1
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacseq_mhc_i_aggregated_report
    type: File?
    label: pVACseq MHC-I Aggregated Report
    'sbg:x': 1597.5888671875
    'sbg:y': 653.84375
  - id: pvacfuse_star_filtered_mhc_ii_1
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacfuse_star_filtered_mhc_ii
    'sbg:fileTypes': TSV
    type: File
    label: pVACfuse MHC-II Star-Annotated Filtered Report
    'sbg:x': 1597.5888671875
    'sbg:y': 760.8125
  - id: pvacfuse_star_filtered_mhc_i_1
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacfuse_star_filtered_mhc_i
    'sbg:fileTypes': TSV
    type: File
    label: pVACfuse MHC-I Star-Annotated Filtered Report
    'sbg:x': 1597.5888671875
    'sbg:y': 867.78125
  - id: pvacfuse_mhc_ii_aggregated_report_1
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacfuse_mhc_ii_aggregated_report
    type: File?
    label: pVACfuse MHC-II Aggregated Report
    'sbg:x': 1597.5888671875
    'sbg:y': 1081.71875
  - id: pvacfuse_mhc_i_aggregated_report_1
    outputSource:
      - pvactools_neoepitope_prediction_1/pvacfuse_mhc_i_aggregated_report
    type: File?
    label: pVACfuse MHC-I Aggregated Report
    'sbg:x': 1597.5888671875
    'sbg:y': 1188.6875
  - id: mhc_ii_filtered_epitopes_1
    outputSource:
      - pvactools_neoepitope_prediction_1/mhc_ii_filtered_epitopes
    type: File?
    label: pVACseq MHC-II Filtered Neoepitope Report
    'sbg:x': 1597.5888671875
    'sbg:y': 1295.65625
  - id: mhc_i_filtered_epitopes_1
    outputSource:
      - pvactools_neoepitope_prediction_1/mhc_i_filtered_epitopes
    type: File?
    label: pVACseq MHC-I Filtered Neoepitope Report
    'sbg:x': 1597.5888671875
    'sbg:y': 1402.625
  - id: rc_exp_annotated_vcf_1
    outputSource:
      - pvactools_neoepitope_prediction_1/rc_exp_annotated_vcf
    'sbg:fileTypes': VCF.GZ
    type: File
    label: Annotated VCF
    'sbg:x': 1597.5888671875
    'sbg:y': 332.9375
steps:
  - id: tmb
    in:
      - id: vcf
        source: input_vcf
      - id: reference_build
        default: hg19
    out:
      - id: tmb
      - id: variants
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        sbg: 'https://sevenbridges.com'
      id: mwonge/mwtest/tmb/1
      baseCommand: []
      inputs:
        - id: vcf
          type: File
          inputBinding:
            position: 1
            separate: false
            shellQuote: false
            loadContents: true
          label: Input VCF
          doc: VCF from which to calculate TMB.
          'sbg:fileTypes': 'VCF, VCF.GZ'
        - id: reference_build
          type:
            type: enum
            symbols:
              - hg19
              - hg38
            name: reference_build
          label: Reference Build
          doc: >-
            Reference genome build. This determines the denominator of the TMB
            calculation.
      outputs:
        - id: tmb
          type: File
          outputBinding:
            loadContents: true
            glob: tmb
            outputEval: '$(String(self[0].contents))'
        - id: variants
          doc: Number of missense canonical variants in with filter 'PASS'.
          label: Canonical Missense Coding Variants
          type: File
          outputBinding:
            loadContents: true
            glob: variants
            outputEval: '$(String(self[0].contents))'
      doc: >-
        # About this tool

        This tool calculates Tumour Mutational Burden (TMB) from a somatic VCF. 


        ## Calculation

        TMB is calculated as the number of canonical somatic missense variants
        divided by the coding exon size.

        The coding exon size is calculated from a reference genome as follows:

        - Download the coding exon BED file for GENCODE/Ensembl genes from the
        UCSC Table Browser

        - Remove alternative contigs

        - Collpase overlapping regions

        - Remove problematic regions by subtracting the Boyle Lab blacklist -
        this BED is used for filtering VCFs

        - Count bases in the coding regions remaining in the BED file


        The input VCF is filtered for 'PASS' variants and intersected with the
        coding exon BED (indicated above) prior to counting the missense
        canonical variants.


        ## Outputs

        - Number of variants: the number of missense canonical variants in the
        input VCF (with 'PASS' filter and intersected with coding exons BED)

        - TMB: calculated by dividing the number of variants by the coding exon
        size for the genome build selected


        ## Docker

        This tool uses the Docker image: `rachelbj/tmb:1.0`


        ## Documentation

        - [Boyle Lab Blacklists](https://github.com/Boyle-Lab/Blacklist)

        - [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables)
      label: tmb
      arguments:
        - position: 0
          prefix: ''
          separate: false
          shellQuote: false
          valueFrom: 'e=`bedtools intersect -header -a '
        - position: 2
          prefix: ''
          separate: false
          shellQuote: false
          valueFrom: |-
            ${
                var denom = 0
                if (inputs.reference_build == "hg38") {
                    denom = 34.5
                    bed = "/app/gencode_basic_v38.coding.exons.collapsed.filtered.grch38.bed"
                } else {
                    denom = 34.9
                    bed = "/app/ensembl_genes.coding.exons.collapsed.filtered.grch37.bed"
                }
                var cmd = '-b ' + bed + '| bcftools view -H -f PASS | grep -E "missense[^,]+protein_coding[^,]+\\|YES\\|[^,]+Ensembl" | wc -l`'
                cmd = cmd + "; perl -E \"say sprintf('%.2f',$e/" + denom + ")\" > tmb;"
                return cmd
            }
        - position: 3
          prefix: ''
          separate: false
          shellQuote: false
          valueFrom: echo $e > variants;
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'rachelbj/tmb:1.0'
        - class: InlineJavascriptRequirement
      'sbg:appVersion':
        - v1.0
      'sbg:content_hash': a9e092b996bf30e387ebfe81c250dea670db591941cc5bb1eccd2ab4088b06b9b
      'sbg:contributors':
        - rbowen_james
      'sbg:createdBy': rbowen_james
      'sbg:createdOn': 1634349796
      'sbg:id': mwonge/mwtest/tmb/1
      'sbg:image_url': null
      'sbg:latestRevision': 1
      'sbg:modifiedBy': rbowen_james
      'sbg:modifiedOn': 1634592566
      'sbg:project': mwonge/mwtest
      'sbg:projectName': zcc-cavatica
      'sbg:publisher': sbg
      'sbg:revision': 1
      'sbg:revisionNotes': null
      'sbg:revisionsInfo':
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1634349796
          'sbg:revision': 0
          'sbg:revisionNotes': null
        - 'sbg:modifiedBy': rbowen_james
          'sbg:modifiedOn': 1634592566
          'sbg:revision': 1
          'sbg:revisionNotes': null
      'sbg:sbgMaintained': false
      'sbg:validationErrors': []
    label: tmb
    'sbg:x': 281.765625
    'sbg:y': 495.359375
  - id: three_sample_hlatyping
    in:
      - id: sample_1_read1_sequences
        source: sample_1_read1_sequences
      - id: bowtie2_index
        source: bowtie2_index
      - id: sample_3_read2_sequences
        source: sample_3_read2_sequences
      - id: sample_3_read1_sequences
        source: sample_3_read1_sequences
      - id: sample_2_read2_sequences
        source: sample_2_read2_sequences
      - id: sample_2_read1_sequences
        source: sample_2_read1_sequences
      - id: patient_id
        source: patient_id
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
    out:
      - id: tumour_rna_results
      - id: tumour_dna_results
      - id: normal_dna_results
      - id: clin_sig_json
      - id: consensus_txt
      - id: sample1_json
      - id: consensus_json
      - id: sample2_json
      - id: sample3_json
      - id: clin_sig_txt
    run: ./three-sample-hlatyping.cwl
    label: three-sample-hlatyping
    'sbg:x': 281.765625
    'sbg:y': 672.328125
  - id: ipass_patient
    in:
      - id: patient_id
        source: patient_id
      - id: gene_expr_file
        source: gene_expression_file
      - id: gene_col
        source: gene_col
      - id: expr_col
        source: expr_col
    out:
      - id: ipass_score_file
    run: ../tools/ipass-patient.cwl
    label: ipass-patient
    'sbg:x': 281.765625
    'sbg:y': 970.265625
  - id: immunedeconv_1
    in:
      - id: gene_expr
        source: gene_expression_file
      - id: patient_id
        source: patient_id
      - id: gene_col
        source: gene_col
      - id: expr_col
        source: expr_col
    out:
      - id: epic_deconv
      - id: quanitseq_deconv
    run: ../tools/immunedeconv.cwl
    label: immunedeconv
    'sbg:x': 281.765625
    'sbg:y': 1333.171875
  - id: pvactools_neoepitope_prediction_1
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
        source: three_sample_hlatyping/clin_sig_txt
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
        source: vcf_sample_names/tumour_sample_name
      - id: normal_sample_name
        source: vcf_sample_names/normal_sample_name
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
      - id: n_junction_reads
        source: n_junction_reads
      - id: n_spanning_frags
        source: n_spanning_frags
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
      - id: pvacfuse_star_filtered_mhc_ii
      - id: pvacfuse_star_filtered_mhc_i
    run: ./pvactools-neoepitope-prediction.cwl
    label: pvactools-neoepitope-prediction
    'sbg:x': 862.2869873046875
    'sbg:y': 1081.71875
  - id: vcf_sample_names
    in:
      - id: input_vcf
        source: input_vcf
    out:
      - id: tumour_sample_name
      - id: normal_sample_name
    run: ../tools/vcf-sample-names.cwl
    label: vcf-sample-names
    'sbg:x': 281.765625
    'sbg:y': 374.390625
requirements:
  - class: SubworkflowFeatureRequirement
