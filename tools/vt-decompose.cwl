class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: vt_decompose
baseCommand:
  - /vt/vt
  - decompose
inputs:
  - id: seq_regions
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-s'
      separate: false
      shellQuote: false
    label: Sequential Region Selection
    doc: >-
      Sequential region selection as opposed to random access of regions
      specified by the i option.
  - id: input_vcf
    type: File
    inputBinding:
      position: 0
      shellQuote: false
    label: Input VCF
    doc: VCF of variants to decompose.
    'sbg:fileTypes': 'VCF, VCF.GZ'
    secondaryFiles:
      - .tbi
  - id: intervals_string
    type: string?
    inputBinding:
      position: 1
      prefix: ''
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (self) {
                return "-i " + self
            } else {
                return ""
            }
        }
    label: Intervals String
    doc: 'Multiple intervals in <seq>:<start>-<end> format delimited by commas'
  - id: intervals_file
    type: File?
    inputBinding:
      position: 1
      prefix: ''
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (self) {
                return "-I " + self.path
            } else {
                return ""
            }
        }
    label: Intervals File
    doc: >-
      Multiple intervals in <seq>:<start>-<end> format listed in a text file
      line by line
    'sbg:fileTypes': TXT
outputs:
  - id: decomposed_vcf
    doc: Input VCF after decomposition of variants.
    label: Decomposed VCF
    type: File
    outputBinding:
      glob: '*.decomp.vcf'
      outputEval: '$(inheritMetadata(self, inputs.input_vcf))'
    secondaryFiles:
      - .tbi
    'sbg:fileTypes': VCF
doc: >-
  # About this tool

  This tool split multi-allelic sites in a VCF using vt decompose (v0.57721). vt
  decompose should be run prior to VCF readcount annotation.


  ## Docker

  This tool uses the Docker image `rachelbj/vt-decompose:1.0.0`.


  ## Output

  The decomposed output VCF will have the extension `.decomp.vcf`.


  ## Documentation

  Documentation for vt decompose is available
  [here](https://genome.sph.umich.edu/wiki/Vt)
label: vt-decompose
arguments:
  - position: 2
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |-
      ${
          var split_vcf = inputs.input_vcf.path.split("/")
          return "-o " + split_vcf[split_vcf.length - 1].replace(/\.vcf.*/g, "") + '.decomp.vcf'
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/vt-decompose:1.0.0'
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
'sbg:toolAuthor': Adrian Tan <atks@umich.edu>
'sbg:toolkit': vt
'sbg:toolkitVersion': '0.57721'
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
