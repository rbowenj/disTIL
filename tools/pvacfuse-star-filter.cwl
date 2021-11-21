class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: pvacfuse_star_filter
baseCommand: []
inputs:
  - id: pvacfuse_report
    type: File?
    inputBinding:
      position: 0
      shellQuote: false
      loadContents: true
      valueFrom: |-
        ${
            if (self) {
                var flags = ""
                if (inputs.n_junction_reads) {
                    flags = flags + "-n_junction_reads " + inputs.n_junction_reads
                }
                if (inputs.n_spanning_frags) {
                    flags = flags + " -n_spanning_frags " + inputs.n_spanning_frags
                }
                return "python /app/pvacfuse_star_filter.py " + self.path + ' ' + inputs.star_fusions.path + ' ' + flags
            } else {
                return "echo 'No inputs report' 1>&2"
            }
        }
    label: pVACfuse Filtered Report
    doc: A filtered report output by pVACfuse.
    'sbg:fileTypes': TSV
  - id: star_fusions
    type: File
    label: STAR-Fusion TSV
    doc: TSV of fusions called by STAR-Fusion.
    'sbg:fileTypes': TSV
  - 'sbg:toolDefaultValue': '0'
    id: n_junction_reads
    type: int?
    label: Junction Read Count Threshold
    doc: >-
      Predicted neoepitopes with junction read counts below this value will be
      discarded.
  - 'sbg:toolDefaultValue': '0'
    id: n_spanning_frags
    type: int?
    label: Spanning Fragment Count Threshold
    doc: >-
      Predicted neoepitopes with spanning fragment counts below this value will
      be discarded.
outputs:
  - id: pvacfuse_star_filtered
    doc: >-
      pVACfuse filtered report further filtered using junction read counts
      annotated from STAR-Fusion TSV.
    label: Junction Read Count Filtered pVACfuse Report
    type: File?
    outputBinding:
      glob: '*.junctionFiltered.tsv'
    'sbg:fileTypes': TSV
doc: >-
  # About this tool

  This tool annotates (and optionally filters) a pVACfuse report using junction
  read counts from the STAR-Fusion TSV used as input to pVACfuse.

  The left and right gene names and breakpoints are extracted from the pVAcfuse
  report and STAR-Fusion TSV. These fields are then used to join pVACfuse and
  STAR-Fusion records.

  Optionally, the user can set junction read count and spanning fragment count
  minimum thresholds - predicted neoepitopes with junction read counts/spanning
  fragment counts lower than these thresholds are discarded.


  ## Docker

  This tool uses the Docker image `rachelbj/pvactools-filters:1.0`.


  ## Documentation

  - [pVACfuse](https://pvactools.readthedocs.io/en/latest/pvacfuse.html)
label: pvacfuse-star-filter
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/pvacfuse-star-filter:1.0'
  - class: InlineJavascriptRequirement
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
