class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: vcf_sample_names
baseCommand:
  - python3
  - /app/get_vcf_sample_names.py
inputs:
  - id: input_vcf
    type: File
    inputBinding:
      position: 0
      shellQuote: false
      loadContents: true
    label: Input VCF
    doc: VCF to get sample names from.
    'sbg:fileTypes': 'VCF, VCF.GZ'
outputs:
  - id: tumour_sample_name
    doc: >-
      The name of the tumour sample (sample not ending with "G") in the input
      VCF.
    label: Tumour Sample Name
    type: string
    outputBinding:
      loadContents: true
      glob: ./tumour.txt
      outputEval: '$(self[0].contents)'
  - id: normal_sample_name
    doc: The name of the normal sample (sample ending with "G") in the input VCF.
    label: Normal Sample Name
    type: string
    outputBinding:
      loadContents: true
      glob: ./normal.txt
      outputEval: '$(self[0].contents)'
label: vcf-sample-names
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/vcf-sample-names:1.0'
  - class: InlineJavascriptRequirement
