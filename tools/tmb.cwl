class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: tmb
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

  TMB is calculated as the number of canonical somatic missense variants divided
  by the coding exon size.

  The coding exon size is calculated from a reference genome as follows:

  - Download the coding exon BED file for GENCODE/Ensembl genes from the UCSC
  Table Browser

  - Remove alternative contigs

  - Collpase overlapping regions

  - Remove problematic regions by subtracting the Boyle Lab blacklist - this BED
  is used for filtering VCFs

  - Count bases in the coding regions remaining in the BED file


  The input VCF is filtered for 'PASS' variants and intersected with the coding
  exon BED (indicated above) prior to counting the missense canonical variants.


  ## Outputs

  - Number of variants: the number of missense canonical variants in the input
  VCF (with 'PASS' filter and intersected with coding exons BED)

  - TMB: calculated by dividing the number of variants by the coding exon size
  for the genome build selected


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
