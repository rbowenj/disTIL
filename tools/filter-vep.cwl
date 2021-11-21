class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: filter_vep
baseCommand: []
inputs:
  - id: vep_cache
    type: File
    'sbg:fileTypes': TAR.GZ
  - id: vcf
    type: File
    inputBinding:
      position: 0
      prefix: '-i'
      shellQuote: false
outputs: []
label: filter-vep
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: mkdir ./cache && tar -xvf $(inputs.vep_cache.basename) -C ./cache && vep
  - position: 1
    prefix: ''
    shellQuote: false
    valueFrom: >-
      --cache --dir_cache ./cache --sift b --canonical --symbol --tab --fields
      Uploaded_variation,SYMBOL,CANONICAL,SIFT -o STDOUT | \

      ./filter_vep --filter "CANONICAL is YES" 1>&2
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ensemblorg/ensembl-vep:release_104.3'
