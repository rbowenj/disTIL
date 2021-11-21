class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: get_aggregate_report
baseCommand: []
inputs:
  - id: mhc_i
    type: File?
  - id: mhc_ii
    type: File?
  - id: combined
    type: File?
  - id: patient_id
    type: string
outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.patient_id).all_epitopes.aggregated.tsv
label: get-aggregate-report
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          if (inputs.combined) {
              return "cp " + inputs.combined.path + " ./" + inputs.patient_id + ".all_epitopes.aggregated.tsv"
          }
          if (inputs.mhc_i) {
              return "cp " + inputs.mhc_i.path + " ./" + inputs.patient_id + ".all_epitopes.aggregated.tsv"
          }
          if (inputs.mhc_ii) {
              return "cp " + inputs.mhc_ii.path + " ./" + inputs.patient_id + ".all_epitopes.aggregated.tsv"
          }
      }
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/vcf-sample-names:1.0'
  - class: InlineJavascriptRequirement
