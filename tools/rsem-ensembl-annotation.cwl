class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: rsem_ensembl_annotation
baseCommand: []
inputs:
  - id: rsem_results
    type: File
    inputBinding:
      position: 3
      separate: false
      shellQuote: false
outputs:
  - id: output
    type: File
    outputBinding:
      glob: '*_ensemblAnnot'
doc: >-
  # About this tool

  This tool uses R BioMart to convert HGNC symbols/RefSeq mRNA to Ensembl gene
  and transcript IDs. This is needed to prepare RSEM expression files prior to
  VCF expression annotation for pVACtools.


  ## Docker

  This tool uses the Docker image: `rachelbj/symbol_to_ensembl`
label: rsem-ensembl-annotation
arguments:
  - position: 2
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: Rscript /app/RSEM_symbolToEnsembl.R
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/symbol_to_ensembl:latest'
