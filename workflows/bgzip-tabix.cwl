class: Workflow
cwlVersion: v1.0
id: bgzip_tabix
doc: |-
  # About this tool
  This tool runs:
  1. bgzip compression of an input file
  2. tabix indexing of the bgzipped file

  # Documentation
  - [bgzip](http://www.htslib.org/doc/bgzip.html)
  - [tabix](http://www.htslib.org/doc/tabix.html)
label: bgzip-tabix
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: input_file
    type: File
    label: Input File
    doc: An unzipped file to be bgzipped.
    'sbg:x': 0
    'sbg:y': 60.5
outputs:
  - id: index_file
    outputSource:
      - tabix/index_file
    'sbg:fileTypes': GZ.TBI
    type: File?
    label: Index file
    doc: Index for the bgzipped file.
    'sbg:x': 600.4171142578125
    'sbg:y': 114
  - id: bgzipped_file
    outputSource:
      - bgzip/bgzipped_file
    'sbg:fileTypes': GZ
    type: File
    label: Bgzipped File
    doc: Bgzipped version of the input file.
    'sbg:x': 343.84375
    'sbg:y': 114
  - id: zipped_with_index
    outputSource:
      - tabix/zipped_with_index
    'sbg:fileTypes': GZ
    type: File
    label: Zipped File with Index
    doc: The bgzipped input file with its newly created index as a secondary file.
    secondaryFiles:
      - .tbi
    'sbg:x': 600.5211181640625
    'sbg:y': -68.5
steps:
  - id: bgzip
    in:
      - id: input_file
        source: input_file
    out:
      - id: bgzipped_file
    run: ../tools/bgzip.cwl
    label: bgzip
    doc: Run bgzip compression on the input file.
    'sbg:x': 130.828125
    'sbg:y': 60.5
  - id: tabix
    in:
      - id: input_file
        source: bgzip/bgzipped_file
    out:
      - id: index_file
      - id: zipped_with_index
    run: ../tools/tabix.cwl
    label: tabix
    doc: Run tabix to create an index of the bgzipped file.
    'sbg:x': 343.84375
    'sbg:y': 0
requirements: []
