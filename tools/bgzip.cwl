class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: bgzip
baseCommand:
  - bgzip
inputs:
  - id: input_file
    type: File
    inputBinding:
      position: 0
      shellQuote: false
    label: Input File
    doc: File to run bgzip on.
outputs:
  - id: bgzipped_file
    doc: The bgzipped version of the input file.
    label: Bgzipped File
    type: File
    outputBinding:
      glob: $(inputs.input_file.basename).gz
      outputEval: '$(inheritMetadata(self, inputs.input_file))'
    'sbg:fileTypes': GZ
doc: |-
  # About this tool
  This app runs bgzip compression on a single input file.

  ## Docker
  This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
  - htslib v1.10.2
  - bcftools v1.10.2
  - samtools v1.10

  ## Documentation
  - [bgzip](http://www.htslib.org/doc/bgzip.html)
  - [bcftools](http://www.htslib.org/doc/bcftools.html)
  - [samtools](http://www.htslib.org/doc/samtools.html)
label: bgzip
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/samtools:1.10.0'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.input_file)
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
'sbg:toolkit': samtools
'sbg:toolkitVersion': '1.10'
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
