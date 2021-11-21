class: Workflow
cwlVersion: v1.0
id: bcftools_bgzip_tabix
doc: >-
  # About this workflow

  This workflow runs`bcftools view -f PASS` on an input VCF followed by bgzip
  compression and tabix to produce an index.
label: bcftools-bgzip-tabix
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: vcf
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    label: Input VCF
    doc: VCF to filter for PASS variants.
    'sbg:x': -536.6011352539062
    'sbg:y': -151
outputs:
  - id: zipped_with_index
    outputSource:
      - bgzip_tabix/zipped_with_index
    'sbg:fileTypes': GZ
    type: File
    label: Zipped Filtered VCF with Index
    'sbg:x': 93.39886474609375
    'sbg:y': -212
  - id: index_file
    outputSource:
      - bgzip_tabix/index_file
    'sbg:fileTypes': GZ.TBI
    type: File?
    label: Filtered VCF Index
    'sbg:x': 92.39886474609375
    'sbg:y': -82
  - id: bgzipped_file
    outputSource:
      - bgzip_tabix/bgzipped_file
    'sbg:fileTypes': GZ
    type: File
    label: ZIpped Filtered VCF
    'sbg:x': 47.39886474609375
    'sbg:y': 121
steps:
  - id: bcftools_view_pass
    in:
      - id: vcf
        source: vcf
    out:
      - id: pass_filtered_vcf
    run: ../tools/bcftools-view-pass.cwl
    label: bcftools-view-pass
    'sbg:x': -353
    'sbg:y': -150
  - id: bgzip_tabix
    in:
      - id: input_file
        source: bcftools_view_pass/pass_filtered_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': -109
    'sbg:y': -159
requirements:
  - class: SubworkflowFeatureRequirement
