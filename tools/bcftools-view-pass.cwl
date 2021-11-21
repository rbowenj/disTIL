class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: bcftools_view_pass
baseCommand: []
inputs:
  - id: vcf
    type: File
    label: VCF
    doc: VCF to filter.
outputs:
  - id: pass_filtered_vcf
    doc: Input VCF filtered to contain only PASS variants.
    label: Pass Filtered VCF
    type: File
    outputBinding:
      glob: '*pass.vcf'
      outputEval: '$(inheritMetadata(self, inputs.vcf))'
label: bcftools-view-pass
arguments:
  - position: 1
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |-
      ${
          var full = inputs.vcf.path
          var split_full = full.split('/')
          var base = split_full[split_full.length -1]
          var split_base = base.split('vcf')
          var out = split_base[0] + 'pass.vcf'
          
          var cmd = "bcftools view -f PASS " + inputs.vcf.path + " > " + out
          return cmd
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/samtools:1.10.0'
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
