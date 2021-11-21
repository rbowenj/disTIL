class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: hlahd_two_sample_consensus_parser
baseCommand:
  - python3
  - /app/hlahd_two_sample_consensus_parser.py
inputs:
  - id: hlahd_results_1
    type: File
    inputBinding:
      position: 0
      prefix: '-results_1'
      shellQuote: false
    label: HLA-HD Results File 1
    doc: The first HLA-HD final results file to use in consensus parsing.
    'sbg:fileTypes': TXT
  - id: hlahd_results_2
    type: File?
    inputBinding:
      position: 0
      prefix: '-results_2'
      shellQuote: false
    label: HLA-HD Results File 2
    doc: The second HLA-HD final results file to use in consensus parsing.
    'sbg:fileTypes': TXT
  - id: sample_name
    type: string
    inputBinding:
      position: 1
      shellQuote: false
    label: Sample Name
    doc: The sample name used to name the output files.
outputs:
  - id: consensus_json
    doc: A JSON file containing the consensus HLA types.
    label: Consensus HLA Type JSON
    type: File
    outputBinding:
      glob: '*.two_sample_consensus.json'
    'sbg:fileTypes': JSON
  - id: consensus_txt
    doc: A text file containing the consensus HLA types.
    label: Consensus HLA Type Text File
    type: File
    outputBinding:
      glob: '*.two_sample_consensus.txt'
    'sbg:fileTypes': TXT
  - id: support_json
    doc: >-
      A JSON file containing detailing which input samples supported each HLA
      allele consensus call.
    label: Consensus Support JSON
    type: File
    outputBinding:
      glob: '*.two_sample_support.json'
    'sbg:fileTypes': JSON
  - id: consensus_string
    doc: >-
      A string containing the consensus HLA types in the format to be used by
      pVACseq.
    label: Consensus HLA String
    type: string
    outputBinding:
      glob: hla
doc: >-
  # About this tool

  This tool calls consensus HLA alleles using two HLA-HD final results text
  files.  

  To make a consensus call, it is required that the allele be present in both
  results files with the same two- or three-field accuracy.


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
label: hlahd-two-sample-consensus-parser
arguments:
  - position: 2
    prefix: ''
    shellQuote: false
    valueFrom: '> hla'
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/hlahd-consensus:2.0.0'
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
