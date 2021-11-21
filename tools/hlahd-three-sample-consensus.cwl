class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: hlahd_three_sample_consensus_parser
baseCommand:
  - python3
  - /app/hlahd_consensus_parser_v2.py
inputs:
  - id: tumour_rna
    type: File
    inputBinding:
      position: 0
      prefix: '-s1'
      shellQuote: false
    label: Tumour RNA HLA-HD Results
    doc: >-
      The final results text file produced by running HLA-HD on the tumour RNA
      sample.
    'sbg:fileTypes': TXT
  - id: tumour_dna
    type: File?
    inputBinding:
      position: 0
      prefix: '-s3'
      shellQuote: false
    label: Tumour DNA HLA-HD Results
    doc: >-
      The final results text file produced by running HLA-HD on the tumour DNA
      sample.
    'sbg:fileTypes': TXT
  - id: normal_dna
    type: File
    inputBinding:
      position: 0
      prefix: '-s2'
      shellQuote: false
    label: Normal DNA HLA-HD Results
    doc: >-
      The final results text file produced by running HLA-HD on the normal DNA
      sample.
    'sbg:fileTypes': TXT
  - id: sample_name
    type: string
    inputBinding:
      position: 1
      shellQuote: false
    label: Sample Name
    doc: The sample name used to name the output files.
outputs:
  - id: research_grade_json
    doc: A JSON file containing the consensus HLA types.
    label: Consensus HLA Type JSON
    type: File
    outputBinding:
      glob: '*Sample_researchGrade_hla.json'
    'sbg:fileTypes': JSON
  - id: clinical_grade_txt
    doc: A text file containing the consensus HLA types.
    label: Consensus HLA Type Text File
    type: File
    outputBinding:
      glob: '*Sample_clinicalGrade_hla.txt'
    'sbg:fileTypes': TXT
  - id: consensus_string
    doc: >-
      A string containing the consensus HLA types in the format to be used by
      pVACseq.
    label: Consensus HLA String
    type: string
    outputBinding:
      glob: hla
  - id: sample1_json
    doc: >-
      A JSON file containing detailing which input samples supported each HLA
      allele consensus call.
    label: Consensus Support JSON
    type: File
    outputBinding:
      glob: '*sample1_hla.json'
    'sbg:fileTypes': JSON
  - id: clinical_grade_json
    type: File
    outputBinding:
      glob: '*Sample_clinicalGrade_hla.json'
  - id: sample2_json
    type: File
    outputBinding:
      glob: '*sample2_hla.json'
  - id: sample3_json
    type: File?
    outputBinding:
      glob: '*sample3_hla.json'
doc: >-
  # About this tool

  This tool produces a consensus HLA type based on three sets of results from
  HLA-HD. The three sets of results are derived from running HLA-HD on:

  - Tumour DNA

  - Tumour RNA

  - Normal DNA


  ## Before running this tool

  Before running this tool, you need to have obtained the HLA-HD results for
  each of the three aforementioned samples. You should use the 'final results'
  text files output by HLA-HD as inputs to this tool.

  **If you do not have tumour DNA, RNA and normal DNA data available, we
  recommend running HLA-HD (tool provided in this repo) on whichever samples you
  do have available and simply using those results.**


  ## Logic

  This tool compares the HLA types determined from three samples, and defines a
  'consensus' type for each allele according to the following rules:

  - An allele must have at least 2-field (4-digit) accuracy in order to be
  considered, otherwise it is considered 'non-typed'.

  - If a 3-field allele (e.g. HLA-A*01:01:01) is an exact match across at least
  two samples, it is taken as the consensus type.

  - If there is no 3-field match across at least two samples, then the allele is
  cut down to 2-fields (e.g. HLA-A*01:01) and again compared between samples. If
  this 2-field allele is an exact match across at least two samples, then it is
  taken as the consensus type.

  - If there is inadequate support for an allele, it is set as "Not typed".


  Using the above logic, this tool only provides a consensus type for an allele
  if there is adequate support for it (i.e. exact match across at least two
  samples). These consensus types may be alleles with 2- or 3-field accuracy.  

  Note that any allele set as "Not typed" may be a result of either HLA-HD not
  being able to type the allele, or there being inadequate concordance across
  the three samples.


  ## Outputs

  The consensus HLA types are output in three formats:

  - JSON format: this output contains **all predicted HLA alleles** and is
  useful for storing HLA data for patients, and useful for downstream analysis.

  - String format: this output **only contains the HLA alleles 'A', 'B', 'C',
  'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1'** (since these are the most comonly
  used HLA alleles for purposes such as donor matching). These alleles are
  formatted as a string ready to be passed to pVACseq.

  - Text format: this output is simply the string format, stored in a text file.


  Additionally a 'consensus support JSON' is produced, detailing which samples
  supported the HLA calls made for each allele. Note that these follow the order
  of alleles in the consensus HLA type JSON output.


  ## Docker

  This tool uses the Docker image: `rachelbj/hlahd-consensus:2.0.0`.


  ## Documentation

  - [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/)
label: hlahd-three-sample-consensus-parser
arguments:
  - position: 2
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: '> hla'
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/hlahd-consensus:2.0.0'
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
