class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: samtools_fastq
baseCommand:
  - samtools
  - fastq
inputs:
  - id: input_alignment
    type: File
    inputBinding:
      position: 10
      shellQuote: false
    label: Input Alignment
    doc: 'An input BAM, SAM or CRAM file.'
    'sbg:fileTypes': 'BAM, SAM, CRAM'
outputs:
  - id: output_fastq_1
    doc: Paired-end FASTQ 1 output by samtools fastq.
    label: Paired-End FASTQ 1
    type: File
    outputBinding:
      glob: |-
        ${
            var input_split = inputs.input_alignment.path.split('/')
            var input_base = input_split[input_split.length - 1]

            return input_base + '.pe_1.fastq'
        }
      outputEval: |-
        ${
          var out = self[0];
          out.metadata = {'paired_end' : ''}
          out.metadata['paired_end'] = 1;
          
          return out
        }
    'sbg:fileTypes': FASTQ
  - id: output_fastq_2
    doc: Paired-end FASTQ 2 output by samtools fastq.
    label: Paired-End FASTQ 2
    type: File
    outputBinding:
      glob: |-
        ${
            var input_split = inputs.input_alignment.path.split('/')
            var input_base = input_split[input_split.length - 1]

            return input_base + '.pe_2.fastq'
        }
      outputEval: |-
        ${
          var out = self[0];
          out.metadata = {'paired_end' : ''}
          out.metadata['paired_end'] = 1;
          
          return out
        }
    'sbg:fileTypes': FASTQ
doc: >-
  # About this tool

  This app runs samtools fastq on an input alignment.  

  Note that the output paired-end FASTQs have the 'paired_end' metadata field
  populated with '1' or '2'.


  ## Docker

  This CWL tool uses the Docker image `rachelbj/samtools:1.10.0` which contains:

  - htslib v1.10.2

  - bcftools v1.10.2

  - samtools v1.10


  ## Documentation

  - [samtools](http://www.htslib.org/doc/samtools.html)
label: samtools-fastq
arguments:
  - position: 0
    prefix: '-1'
    shellQuote: false
    valueFrom: |-
      ${
          var input_split = inputs.input_alignment.path.split('/')
          var input_base = input_split[input_split.length - 1]

          return input_base + '.pe_1.fastq'
      }
  - position: 0
    prefix: '-2'
    shellQuote: false
    valueFrom: |-
      ${
          var input_split = inputs.input_alignment.path.split('/')
          var input_base = input_split[input_split.length - 1]

          return input_base + '.pe_2.fastq'
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/samtools:1.10.0'
  - class: InlineJavascriptRequirement
'sbg:toolkit': samtools
'sbg:toolkitVersion': '1.10'
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
