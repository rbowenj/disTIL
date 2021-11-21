class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: sample_id_from_metadata
baseCommand:
  - bash
  - '-c'
  - '"echo -n'
inputs:
  - id: input_file
    type: File
    inputBinding:
      position: 0
      prefix: ''
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
          var input_files = [].concat(inputs.input_file)

          if (input_files[0].metadata && input_files[0].metadata.sample_id) {
            var sample = input_files[0].metadata.sample_id
          } else {
            var sample = "sample_unknown"
          }
         
          return sample
        }
    label: Input File
    doc: Input file from which to extract sample ID metadata.
outputs:
  - id: sample_id
    doc: The sample ID retrieved from the input file's metadata.
    label: Sample ID
    type: string
    outputBinding:
      loadContents: true
      glob: sample_id.txt
      outputEval: '$(String(self[0].contents))'
label: sample-id-from-metadata
arguments:
  - position: 1
    prefix: ''
    shellQuote: false
    valueFrom: '> sample_id.txt"'
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:18.04'
  - class: InlineJavascriptRequirement
