class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: ipass_patient
baseCommand: []
inputs:
  - id: patient_id
    type: string
    inputBinding:
      position: 0
      separate: false
      shellQuote: false
    label: Patient ID
    doc: The patient ID for the sample being analysed.
  - id: gene_expr_file
    type: File
    inputBinding:
      position: 1
      separate: false
      shellQuote: false
    label: Gene Expression File
    doc: A TSV output by an expression quantification algorithm (e.g. RSEM).
    'sbg:fileTypes': TSV
  - id: gene_col
    type: string
    inputBinding:
      position: 3
      separate: false
      shellQuote: false
    label: Gene ID Column Name
    doc: >-
      Name of the column in the gene expression TSV containing gene IDs. E.g.
      gene_id
  - id: expr_col
    type: string
    inputBinding:
      position: 4
      separate: false
      shellQuote: false
    label: Gene Expression Column Name
    doc: >-
      Name of the column in the gene expression TSV containing expression
      values. E.g. TPM
outputs:
  - id: ipass_score_file
    type: File
    outputBinding:
      glob: $(inputs.patient_id)_IPASS.txt
doc: >-
  # About this tool

  This tool runs the IPASS gene expression classifier on a TSV of gene
  expression values output by an expression quantification algorithm (such as
  RSEM).


  ## Inputs

  - Patient ID: the ID of the patient being analysed. This is used to name
  output files.

  - Gene Expression TSV: a gene expression TSV produced by an expression
  quantification algorithm such as RSEM.

  - Gene ID Column Name: the name of the column containing gene IDs in the gene
  expression TSV.

  - Expression Value Column Name: the name of the column containing expression
  values (usually TPM) in the gene expression TSV.


  ## Docker

  This tool uses the Docker image `rachelbj/ipass-patient:1.0`
label: ipass-patient
arguments:
  - position: 0
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: Rscript /app/patientIPASS.R
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/ipass-patient:1.0'
  - class: InlineJavascriptRequirement
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
