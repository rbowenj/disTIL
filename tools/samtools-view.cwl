class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: samtools_view
baseCommand:
  - samtools
  - view
inputs:
  - id: output_format
    type:
      - 'null'
      - type: enum
        symbols:
          - BAM
          - SAM
          - CRAM
        name: output_format
    inputBinding:
      position: 0
      prefix: '--output-fmt'
      shellQuote: false
    label: Output Format
    doc: 'Specifies BAM, SAM or CRAM output format.'
  - id: input_alignment
    type: File
    inputBinding:
      position: 10
      shellQuote: false
    label: Input Alignment
    doc: 'An input BAM, SAM or CRAM file.'
    'sbg:fileTypes': 'BAM, SAM, CRAM'
  - id: fast_bam_compression
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-1'
      shellQuote: false
    label: Fast BAM Compression
    doc: Whether the output file (which must be BAM format) should be compressed.
  - id: include_header
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-h'
      shellQuote: false
    label: Include Header
    doc: Whether the input alignment header should be included in the output file.
  - id: header_only
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-H'
      shellQuote: false
    label: Output Header Only
    doc: >-
      When this option is selected, the output file will only contain the header
      from the input alignment.
  - id: include_reads
    type: string?
    inputBinding:
      position: 0
      prefix: '-f'
      shellQuote: false
    label: Include reads with all of these flags
    doc: >-
      Only output alignments with all bits set in this integer present in the
      FLAG field.
  - id: exclude_reads_any
    type: string?
    inputBinding:
      position: 0
      prefix: '-F'
      shellQuote: false
    label: Exclude reads with any of these flags
    doc: >-
      Do not output alignments with any bits set in this integer present in the
      FLAG field.
  - id: exclude_reads_all
    type: string?
    inputBinding:
      position: 0
      prefix: '-G'
      shellQuote: false
    label: Exclude reads with all of these flags
    doc: >-
      Only exclude reads with all of the bits set in this integer present in the
      FLAG field.
  - id: bed_file
    type: File?
    inputBinding:
      position: 0
      prefix: '-L'
      shellQuote: false
    label: BED File
    doc: Only output alignments overlapping the regions specified in this BED file.
    'sbg:fileTypes': BED
outputs:
  - id: output_alignment
    doc: Output from samtools view.
    label: Output Alignment
    type: File
    outputBinding:
      glob: |-
        ${
            //Find the input basename (without the file extension)
            var input_split = inputs.input_alignment.path.split('/')
            var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
            
            var ext = ""
            //Determine the output file extension
            if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
                ext = ".bam"
            } else if (inputs.output_format == "SAM") {
                ext = ".sam"
            } else if (inputs.output_format == "CRAM") {
                ext = ".cram"
            } else {
                ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
            }
            
            //If only output header then add '.header'
            if (inputs.header_only) {
                return input_base + '.header' + ext
            }
            
            //If filtered on flags/bed file then add '.filtered'
            if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
                return input_base + '.filtered' + ext
            }
            
            return input_base + ext
        }
      outputEval: '$(inheritMetadata(self, inputs.input_alignment))'
    'sbg:fileTypes': 'BAM, SAM, CRAM'
doc: |-
  # About this tool
  This app runs samtools view on an input alignment.

  ## Docker
  This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:
  - htslib v1.10.2
  - bcftools v1.10.2
  - samtools v1.10

  ## Documentation
  - [samtools](http://www.htslib.org/doc/samtools.html)
label: samtools-view
arguments:
  - position: 0
    prefix: '-o'
    shellQuote: false
    valueFrom: |-
      ${
          //Find the input basename (without the file extension)
          var input_split = inputs.input_alignment.path.split('/')
          var input_base = input_split[input_split.length - 1].split('.').slice(0,-1). join('.')
          
          var ext = ""
          //Determine the output file extension
          if (inputs.fast_bam_compression || inputs.output_format == "BAM") {
              ext = ".bam"
          } else if (inputs.output_format == "SAM") {
              ext = ".sam"
          } else if (inputs.output_format == "CRAM") {
              ext = ".cram"
          } else {
              ext = "." + input_split[input_split.length - 1].split('.').slice(-1)
          }
          
          //If only output header then add '.header'
          if (inputs.header_only) {
              return input_base + '.header' + ext
          }
          
          //If filtered on flags/bed file then add '.filtered'
          if (inputs.include_reads || inputs.exclude_reads_any || inputs.exclude_reads_all || inputs.bed_file) {
              return input_base + '.filtered' + ext
          }
          
          return input_base + ext
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
'sbg:toolkit': samtools
'sbg:toolkitVersion': '1.10'
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
