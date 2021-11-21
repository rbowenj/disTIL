class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - ln
  - '-s'
inputs:
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
  - id: allele_specific_binding_thresholds
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--allele-specific-binding-thresholds'
      shellQuote: false
  - id: alleles
    type: File
    inputBinding:
      position: 3
      prefix: ''
      separate: false
      shellQuote: false
      loadContents: true
      valueFrom: $(String(self.contents))
    label: HLA Alleles
    doc: >-
      A text file containing the HLA alleles for this sample as a comma
      separated string.
    'sbg:fileTypes': TXT
  - id: binding_threshold
    type: int?
    inputBinding:
      position: 0
      prefix: '-b'
      shellQuote: false
  - id: downstream_sequence_length
    type: string?
    inputBinding:
      position: 0
      prefix: '-d'
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
  - id: exclude_nas
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--exclude-NAs'
      shellQuote: false
  - id: expn_val
    type: float?
    inputBinding:
      position: 0
      prefix: '--expn-val'
      shellQuote: false
  - id: fasta_size
    type: int?
    inputBinding:
      position: 0
      prefix: '-s'
      shellQuote: false
  - id: iedb_retries
    type: int?
    inputBinding:
      position: 0
      prefix: '-r'
      shellQuote: false
  - id: input_vcf
    type: File
    inputBinding:
      position: 1
      shellQuote: false
    secondaryFiles:
      - .tbi
  - id: keep_tmp_files
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-k'
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
  - id: minimum_fold_change
    type: float?
    inputBinding:
      position: 0
      prefix: '-c'
      shellQuote: false
  - default: 8
    id: n_threads
    type: int?
    inputBinding:
      position: 0
      prefix: '--n-threads'
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
  - id: net_chop_threshold
    type: float?
    inputBinding:
      position: 0
      prefix: '--net-chop-threshold'
      shellQuote: false
  - id: netmhc_stab
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--netmhc-stab'
      shellQuote: false
  - id: normal_cov
    type: int?
    inputBinding:
      position: 0
      prefix: '--normal-cov'
      shellQuote: false
  - id: normal_sample_name
    type: string?
    inputBinding:
      position: 0
      prefix: '--normal-sample-name'
      shellQuote: false
  - id: normal_vaf
    type: float?
    inputBinding:
      position: 0
      prefix: '--normal-vaf'
      shellQuote: false
  - id: percentile_threshold
    type: int?
    inputBinding:
      position: 0
      prefix: '--percentile-threshold'
      shellQuote: false
  - id: phased_proximal_variants_vcf
    type: File?
    inputBinding:
      position: 0
      prefix: '-p'
      shellQuote: false
    secondaryFiles:
      - .tbi
  - id: prediction_algorithms
    type: 'string[]'
    inputBinding:
      position: 4
      shellQuote: false
  - id: run_reference_proteome_similarity
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--run-reference-proteome-similarity'
      shellQuote: false
  - id: sample_name
    type: string
    inputBinding:
      position: 2
      shellQuote: false
  - id: tdna_cov
    type: int?
    inputBinding:
      position: 0
      prefix: '--tdna-cov'
      shellQuote: false
  - id: tdna_vaf
    type: float?
    inputBinding:
      position: 0
      prefix: '--tdna-vaf'
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
  - id: trna_cov
    type: int?
    inputBinding:
      position: 0
      prefix: '--trna-cov'
      shellQuote: false
  - id: trna_vaf
    type: float?
    inputBinding:
      position: 0
      prefix: '--trna-vaf'
      shellQuote: false
outputs:
  - id: combined_aggregated_report
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/combined/$(inputs.sample_name).all_epitopes.aggregated.tsv
  - id: combined_all_epitopes
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/combined/$(inputs.sample_name).all_epitopes.tsv
  - id: combined_filtered_epitopes
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/combined/$(inputs.sample_name).filtered.tsv
  - id: mhc_i_aggregated_report
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/MHC_Class_I/$(inputs.sample_name).all_epitopes.aggregated.tsv
  - id: mhc_i_all_epitopes
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/MHC_Class_I/$(inputs.sample_name).all_epitopes.tsv
  - id: mhc_i_filtered_epitopes
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/MHC_Class_I/$(inputs.sample_name).filtered.tsv
  - id: mhc_ii_aggregated_report
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/MHC_Class_II/$(inputs.sample_name).all_epitopes.aggregated.tsv
  - id: mhc_ii_all_epitopes
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/MHC_Class_II/$(inputs.sample_name).all_epitopes.tsv
  - id: mhc_ii_filtered_epitopes
    type: File?
    outputBinding:
      glob: >-
        pvacseq_$(inputs.sample_name)/MHC_Class_II/$(inputs.sample_name).filtered.tsv
  - id: pvacseq_predictions
    type: Directory
    outputBinding:
      glob: pvacseq_$(inputs.sample_name)
doc: >-
  # About this tool

  This tool contains CWL from [the pVACtools
  team](https://github.com/genome/analysis-workflows/blob/master/definitions/tools/pvacseq.cwl),
  using pVACtools version 2.0.3.


  ## Edits

  Changes made to the original CWL provided by the pVACtools team include:

  - The allele input is now a text file (rather than an array). The text file
  should contain a string of comma-separated, single-quoted alleles, e.g.
  `'HLA-A*02:01','HLA-A*11:01','HLA-B*01:02','HLA-B*01:02''`

  - The output directory is now named `pvacseq_<sample_name>`.

  - The Docker image has been updated to use pVACtools 2.0.3 (rather than 2.0.1)


  ## Docker

  This tool uses the Docker image: `griffithlab/pvactools:2.0.3`.
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
    prefix: ''
    valueFrom: pvacseq_$(inputs.sample_name)
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: $(inputs.n_threads)
  - class: DockerRequirement
    dockerPull: 'griffithlab/pvactools:2.0.4'
  - class: InlineJavascriptRequirement
'sbg:toolAuthor': GriffithLab
'sbg:toolkit': pvactools
'sbg:toolkitVersion': 2.0.3
'sbg:wrapperAuthor': GriffithLab
