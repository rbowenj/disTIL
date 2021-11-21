class: Workflow
cwlVersion: v1.0
id: pvactools_neoepitope_filtering
label: pvactools-neoepitope-filtering
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: pvacseq_report
    'sbg:fileTypes': TSV
    type: File
    label: pVACseq Aggregate Report
    doc: pVACseq aggregate report of predicted neoepitopes.
    'sbg:x': 0
    'sbg:y': 214
  - id: pvacfuse_report
    'sbg:fileTypes': TSV
    type: File
    label: pVACfuse Aggregate Report
    doc: pVACfuse aggregate report of predicted neoepitopes.
    'sbg:x': 0
    'sbg:y': 321
  - id: star_fusions
    'sbg:fileTypes': TSV
    type: File?
    label: STAR-Fusion TSV
    doc: Optional TSV containing fusions predicted by STAR-Fusion.
    'sbg:x': 0
    'sbg:y': 107
  - id: patient_id
    type: string
    label: Patient ID
    doc: Patient ID used for naming output files.
    'sbg:x': 0
    'sbg:y': 428
outputs:
  - id: pvacseq_snvs_filtered
    outputSource:
      - pvacseq_filtering/pvacseq_snvs_filtered
    'sbg:fileTypes': CSV
    type: File
    label: pVACseq Filtered SNV-Derived Neoepitopes
    doc: Filtered SNV-derived neoepitopes predicted by pVACseq.
    'sbg:x': 675.50244140625
    'sbg:y': 0
  - id: pvacseq_indels_filtered
    outputSource:
      - pvacseq_filtering/pvacseq_indels_filtered
    'sbg:fileTypes': CSV
    type: File
    label: pVACseq Filtered Indel-Derived Neoepitopes
    doc: Filtered indel-derived neoepitopes predicted by pVACseq.
    'sbg:x': 675.50244140625
    'sbg:y': 107
  - id: pvacseq_filtered_pass
    outputSource:
      - pvacseq_filtering/pvacseq_filtered_pass
    'sbg:fileTypes': CSV
    type: File
    label: pVACseq Filtered Pass Neoepitopes
    doc: >-
      Filtered SNV and indel-derived neoepitopes predicted by pVACseq with the
      'Pass' tier.
    'sbg:x': 675.50244140625
    'sbg:y': 214
  - id: pvacseq_filtered
    outputSource:
      - pvacseq_filtering/pvacseq_filtered
    'sbg:fileTypes': CSV
    type: File
    label: pVACseq Filtered Neoepitopes
    doc: Filtered SNV and indel-derived neoepitopes predicted by pVACseq.
    'sbg:x': 675.50244140625
    'sbg:y': 321
  - id: pvacfuse_filtered_strict
    outputSource:
      - pvacfuse_filtering/pvacfuse_filtered_strict
    'sbg:fileTypes': CSV
    type: File?
    label: pVACfuse Strict Filtered Neoepitopes
    doc: >-
      Filtered neoepitopes predicted by pVACfuse, restricted to those with
      junction read count of at least 5 (obtained from STAR-Fusion file).
    'sbg:x': 675.50244140625
    'sbg:y': 428
  - id: pvacfuse_filtered
    outputSource:
      - pvacfuse_filtering/pvacfuse_filtered
    'sbg:fileTypes': CSV
    type: File
    label: pVACfuse Filtered Neoepitopes
    doc: Filtered neoepitopes predicted by pVACfuse.
    'sbg:x': 675.50244140625
    'sbg:y': 535
steps:
  - id: pvacseq_filtering
    in:
      - id: pvacseq_report
        source: pvacseq_report
      - id: patient_id
        source: patient_id
    out:
      - id: pvacseq_filtered
      - id: pvacseq_filtered_pass
      - id: pvacseq_indels_filtered
      - id: pvacseq_snvs_filtered
    run: ../tools/pvacseq-filtering.cwl
    label: pvacseq-filtering
    'sbg:x': 212.78125
    'sbg:y': 179
  - id: pvacfuse_filtering
    in:
      - id: pvacfuse_report
        source: pvacfuse_report
      - id: patient_id
        source: patient_id
      - id: star_fusions
        source: star_fusions
    out:
      - id: pvacfuse_filtered
      - id: pvacfuse_filtered_strict
    run: ../tools/pvacfuse-filtering.cwl
    label: pvacfuse-filtering
    'sbg:x': 212.78125
    'sbg:y': 321
requirements: []
