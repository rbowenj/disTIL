class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: vatools_expression_annotation
baseCommand:
  - vcf-expression-annotator
inputs:
  - id: input_vcf
    type: File
    inputBinding:
      position: 0
      shellQuote: false
    label: Input VCF
    doc: The VCF to annotate with expression values.
    'sbg:fileTypes': VCF.GZ
    secondaryFiles:
      - .tbi
  - id: expression_file
    type: File
    inputBinding:
      position: 1
      shellQuote: false
    label: Expression File
    doc: >-
      The TSV file containing gene or transcript level expression values
      determined from RNA using the expression quantification algorithm.
    'sbg:fileTypes': TSV
  - id: quant_algo
    type:
      type: enum
      symbols:
        - kallisto
        - stringtie
        - cufflinks
        - custom
      name: quant_algo
    inputBinding:
      position: 2
      shellQuote: false
    label: Expression Quantification Algorithm
    doc: >-
      The expression quantification algorithm used to generate the expression
      file. Note that if 'custom' is selected, then the ID and expression
      columns must be provided as parameters.
  - id: expression_level
    type:
      type: enum
      symbols:
        - gene
        - transcript
      name: expression_level
    inputBinding:
      position: 3
      shellQuote: false
    label: Expression Level
    doc: >-
      Either 'gene' or 'transcript', depending on the level at which expression
      values are quantified in the expression file.
  - id: id_column
    type: string?
    inputBinding:
      position: 4
      prefix: '--id-column'
      shellQuote: false
    label: ID Column Name
    doc: >-
      The name of the column containing gene or transcript IDs in the expression
      file. Note that these IDs must be in the same format as IDs in the VCF
      (i.e. Ensembl, RefSeq)
  - id: expression_column
    type: string?
    inputBinding:
      position: 5
      prefix: '--expression-column'
      shellQuote: false
    label: Expression Column Name
    doc: >-
      The name of the column containing expression values (usually in TPM) in
      the expression file.
  - id: sample_name
    type: string
    inputBinding:
      position: 5
      prefix: '--sample-name'
      shellQuote: false
    label: Sample Name
    doc: Name of the sample to be annotated in the VCF file.
outputs:
  - id: exp_annotated_vcf
    type: File
    outputBinding:
      glob: |-
        ${
            if (inputs.expression_level == 'gene') {
                return "*.gx.vcf"
            } else {
                return "*.tx.vcf"
            }
        }
doc: >-
  # About this tool

  This tool runs VAtools expression annotation (v5.0.1) to add genome or
  transcript level expression values to a VCF file.  


  ## RNA expression values

  Genome and transcript level expression values in TSV format must be generated
  before running this analysis. Any expression quantification algorithm can be
  used to generate these values, but note that:

  - If you use a tool other than Kallisto, Stringtie or Cufflinks, you need to
  specify the names of the columns containing gene/transcript IDs and expression
  values.

  - The expression values must be in TSV file format.

  - The gene/transcript ID format used must match that in the VCF. **Given that
  this workflow is likely to be run after VEP annotation with Ensembl IDs (as an
  input for pVACseq), it is recommended that IDs are converted to Ensembl format
  prior to running this analysis (tool provided in this repo).**


  ## Output

  If the expression level is 'gene', then the output annotated VCF will have the
  added extension `gx`.  

  If the expression level is 'transcript', then the output annotated VCF will
  have the added extension `tx`.  



  ## Docker

  This tool uses the Docker image
  [griffithlab/vatools:5.0.1](https://hub.docker.com/r/griffithlab/vatools).


  ## Documentation

  -
  [VAtools](https://vatools.readthedocs.io/en/latest/vcf_expression_annotator.html)
label: vatools-expression-annotation
arguments:
  - position: 6
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: '--ignore-ensembl-id-version'
  - position: 4
    prefix: '--output-vcf'
    shellQuote: false
    valueFrom: |-
      ${
          var full = inputs.input_vcf.path
          var split_full = full.split('/')
          var base = split_full[split_full.length -1]
          var split_base = base.split('vcf')
          
          if (inputs.expression_level == 'gene') {
              return split_base[0] + 'gx.vcf'
          } else {
              return split_base[0] + 'tx.vcf'
          }
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'griffithlab/vatools:5.0.1'
  - class: InlineJavascriptRequirement
'sbg:toolAuthor': GriffithLab
'sbg:toolkit': vatools
'sbg:toolkitVersion': 5.0.1
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
