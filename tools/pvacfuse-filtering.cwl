class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: pvacfuse_filtering
baseCommand:
  - python
  - /app/pvacfuse_star_filter.py
inputs:
  - id: pvacfuse_report
    type: File
    inputBinding:
      position: 0
      shellQuote: false
      loadContents: true
    label: pVACfuse Filtered Report
    doc: A filtered report output by pVACfuse.
    'sbg:fileTypes': TSV
  - id: star_fusions
    type: File?
    inputBinding:
      position: 1
      shellQuote: false
    label: STAR-Fusion TSV
    doc: TSV of fusions called by STAR-Fusion.
    'sbg:fileTypes': TSV
outputs:
  - id: pvacfuse_star_filtered
    doc: >-
      pVACfuse filtered report further filtered using junction read counts
      annotated from STAR-Fusion TSV.
    label: Junction Read Count Filtered pVACfuse Report
    type: File
    outputBinding:
      glob: '*.junctionFiltered.tsv'
    'sbg:fileTypes': TSV
doc: >-
  # About this tool

  This tool filters a pVACfuse report using junction read counts annotated from
  the STAR-Fusion TSV used as input to pVACfuse.

  The left and right gene names and breakpoints are extracted from the pVAcfuse
  report and STAR-Fusion TSV. These fields are then used to join pVACfuse and
  STAR-Fusion records.

  The user can set a junction read count threshold - predicted neoepitopes with
  junction read counts lower than this threshold are discarded.


  ## Docker

  This tool uses the Docker image `rachelbj/pvactools-filters:1.0`.


  ## Documentation

  - [pVACfuse](https://pvactools.readthedocs.io/en/latest/pvacfuse.html)
label: pvacfuse-filtering
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/pvactools-filters:1.0'
  - class: InlineJavascriptRequirement
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
