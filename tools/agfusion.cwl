class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: agfusion
baseCommand: []
inputs:
  - id: fusion_tsv
    type: File
    inputBinding:
      position: 1
      prefix: '-f'
      shellQuote: false
      valueFrom: infile
    label: Fusion TSV
    doc: A TSV file containing fusion variant calls.
    'sbg:fileTypes': TSV
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
    inputBinding:
      position: 2
      prefix: '-a'
      shellQuote: false
    label: Fusion Caller
    doc: The name of the algorithm used to generate the input Fusion TSV.
  - id: ref_genome
    type:
      type: enum
      symbols:
        - hg38
        - hg19
      name: ref_genome
    inputBinding:
      position: 0
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.ref_genome == "hg38") {
                return '&& agfusion download -s homo_sapiens -r 87 -d .'
            } else if (inputs.ref_genome == "hg19") {
                return '&& agfusion download -s homo_sapiens -r 75 -d .'
            }
        }
    label: Reference Genome Build
    doc: >-
      The referece genome build used to generate the fusion variant calls. This
      is used to determine which agfusion database to use (release 87 or 75).
outputs:
  - id: output_dir
    doc: AGFusion output directory.
    label: AGFusion Output Directory
    type: Directory
    outputBinding:
      glob: agfusion-output
doc: >-
  # About this tool

  This tool runs **agfusion v1.2** to annotate gene fusions. This analysis is
  necessary prior to running pVACfuse.


  ## Reference Genome Build

  agfusion requires a reference genome database as input. Two release versions
  are available:

  - Release 87 for hg38

  - Release 75 for hg19

  As part of the analysis, this tool downloads the required database release
  (using the command `agfusion download -g <genome_build>`) based on the input
  reference genome build selected.


  ## Docker

  This tool uses the Docker image: `zlskidmore/agfusion:1.2`


  ## Documentation

  - [agfusion](https://github.com/murphycj/AGFusion#basic-usage)

  - [agfusion publication](https://www.biorxiv.org/content/10.1101/080903v1)
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
      touch infile && if head -n 1 $(inputs.fusion_tsv.path) | grep 'est_J';
      then cut -f1-3,6- $(inputs.fusion_tsv.path) > infile ; else cat
      $(inputs.fusion_tsv.path) > infile ; fi && head infile 1>&2
  - position: 1
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: '&& agfusion batch'
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'zlskidmore/agfusion:1.2'
  - class: InlineJavascriptRequirement
'sbg:license': MIT license
'sbg:toolAuthor': 'Charlie Murphy, Olivier Elemento'
'sbg:toolkit': agfusion
'sbg:toolkitVersion': '1.2'
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
