class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: vatools_readcount_annotation
baseCommand:
  - vcf-readcount-annotator
inputs:
  - id: input_vcf
    type: File
    inputBinding:
      position: 2
      shellQuote: false
    label: Input VCF
    doc: VCF to annotate.
    'sbg:fileTypes': 'VCF.GZ, VCF'
    secondaryFiles:
      - .tbi
  - id: bam_readcount_file
    type: File
    inputBinding:
      position: 3
      shellQuote: false
    label: Bam-readcount File
    doc: TSV file containing readcounts from bam-readcount.
    'sbg:fileTypes': TSV
  - id: sample_name
    type: string
    inputBinding:
      position: 0
      prefix: '-s'
      shellQuote: false
    label: Sample Name
    doc: Name of the sample to annotate in the input VCF file.
  - id: variant_type
    type:
      type: enum
      symbols:
        - snv
        - indel
        - all
      name: variant_type
    inputBinding:
      position: 1
      prefix: '-t'
      shellQuote: false
    label: Variant Type
    doc: Which types of variants to annotate.
  - id: data_type
    type:
      - 'null'
      - type: enum
        symbols:
          - DNA
          - RNA
        name: data_type
    inputBinding:
      position: 4
      shellQuote: false
    label: Data Type
    doc: >-
      Either DNA or RNA, depending on whether the BAM used to generate
      readcounts  contained DNA or RNA data.
outputs:
  - id: annotated_vcf
    doc: VCF annotated with readcounts.
    label: Annotated VCF
    type: File
    outputBinding:
      glob: |-
        ${
            var split_vcf = inputs.input_vcf.path.split("/")
            var ext = ""
            if (inputs.variant_type == "snv") {
                if (inputs.data_type == "DNA") {
                    ext = "dsr"
                } else {
                    ext = "rsr"
                }
            } else if (inputs.variant_type == "indel") {
                if (inputs.data_type == "DNA") {
                    ext = "dir"
                } else {
                    ext = "rir"
                }
            } else {
                if (inputs.data_type == "DNA") {
                    ext = "dr"
                } else {
                    ext = "rr"
                }
            }
            return split_vcf[split_vcf.length - 1].replace(/\.vcf.*/g, "") + '.' + ext + '.vcf'
        }
      outputEval: '$(inheritMetadata(self, inputs.input_vcf))'
    'sbg:fileTypes': VCF
doc: >-
  # About this tool

  This tool runs VAtools readcount annotation (v5.0.1) to add readcounts to a
  VCF file.  


  ## Inputs and Parameters

  - The input VCF file should be decomposed with
  [vt-decompose](https://genome.sph.umich.edu/wiki/Vt) prior to running VAtools.

  - The bam-readcount input file can be generated using bam-readcount from a DNA
  or RNA BAM file.

  - The 'data type' parameter specifies whether the BAM used to generate
  readcounts contained DNA or RNA data.

  - The 'variant type' parameter specifies whether indels, SNVs or all should be
  annotated.

  - The sample name input must match the name of the sample being annotated in
  the VCF file. 


  ## Outputs

  Depending on the data type (DNA/RNA) and variant type (SNV/indel/all), the
  following extensions are added to the annotated output file name:

  - `dsr`: **D**NA **S**NV **r**eadcounts

  - `dir`: **D**NA **i**ndel **r**eadcounts

  - `dr`: **D**NA **r**eadcounts (SNVs and indels)

  - `rsr`: **R**NA **S**NV **r**eadcounts

  - `rir`: **R**NA **i**ndel **r**eadcounts

  - `rr`: **R**NA **r**eadcounts (SNVs and indels)


  ## Docker

  This tool uses the Docker image
  [griffithlab/vatools:5.0.1](https://hub.docker.com/r/griffithlab/vatools).


  ## Documentation

  -
  [VAtools](https://vatools.readthedocs.io/en/latest/vcf_readcount_annotator.html)

  - [bam-readcount](https://github.com/genome/bam-readcount)
label: vatools-readcount-annotation
arguments:
  - position: 1
    prefix: '-o'
    valueFrom: |-
      ${
          var split_vcf = inputs.input_vcf.path.split("/")
          var ext = ""
          if (inputs.variant_type == "snv") {
              if (inputs.data_type == "DNA") {
                  ext = "dsr"
              } else {
                  ext = "rsr"
              }
          } else if (inputs.variant_type == "indel") {
              if (inputs.data_type == "DNA") {
                  ext = "dir"
              } else {
                  ext = "rir"
              }
          } else {
              if (inputs.data_type == "DNA") {
                  ext = "dr"
              } else {
                  ext = "rr"
              }
          }
          return split_vcf[split_vcf.length - 1].replace(/\.vcf.*/g, "") + '.' + ext + '.vcf'
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'griffithlab/vatools:5.0.1'
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
'sbg:toolAuthor': GriffithLab
'sbg:toolkit': vatools
'sbg:toolkitVersion': 5.0.1
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
