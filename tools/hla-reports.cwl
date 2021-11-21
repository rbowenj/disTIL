class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: mwonge/mwtest/hla-reports/2
baseCommand: []
inputs:
  - id: full_hla
    type: File
    label: HLA Consensus JSON
    doc: >-
      A JSON file containing all HLA consensus alleles. All HLA genes typed by
      HLA-HD are included at two- or three-field accuracies.
    'sbg:fileTypes': JSON
  - id: clin_sig_hla
    type: File
    label: Clinically Significant HLA Consensus JSON
    doc: >-
      A JSON file containing consensus alleles for the clinically significant
      classical HLA genes. Only classical Class-I and Class-II genes are
      included and all alleles are truncated to two-field accuracy.
    'sbg:fileTypes': JSON
  - id: patient_id
    type: string
    label: Patient ID
    doc: ID for the patient.
outputs:
  - id: hla_report
    doc: A PDF report containing the HLA consensus results.
    label: HLA Report
    type: File
    outputBinding:
      glob: |-
        ${
            return inputs.patient_id + '_hlaReport.pdf'
        }
    'sbg:fileTypes': PDF
doc: >-
  # About this tool


  This tool generates an HLA report from the output of disTIL HLA consensus HLA
  typing.
label: hla-report
arguments:
  - position: 2
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          var cmd = 'Rscript -e "rmarkdown::render(\'hla_report_generator.Rmd\',params=list(full_hla_json=\'./' + inputs.full_hla.basename + '\', clin_hla_json=\'./' + inputs.clin_sig_hla.basename + '\', pid=\'' + inputs.patient_id + '\'), output_file=paste(\'' + inputs.patient_id + '\', \'_hlaReport.pdf\', sep=\'\'))\"'
          return cmd
      }
  - position: 10
    prefix: ''
    shellQuote: false
    valueFrom: 1>&2
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: ' cp /hla_report_generator.Rmd . && ls 1>&2 &&'
  - position: 1
    prefix: ''
    shellQuote: false
    valueFrom: >-
      echo "RES" 1>&2 && head $(inputs.full_hla.basename) 1>&2 && echo "CLIN"
      1>&2 && head $(inputs.clin_sig_hla.basename) 1>&2 &&
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/hla-report:1.0.0'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.clin_sig_hla)
      - $(inputs.full_hla)
  - class: InlineJavascriptRequirement
'sbg:appVersion':
  - v1.0
'sbg:content_hash': aa85380446ce78b154d4fc82224a7a8ca4d58e19fa5393ea4a8f4ea870b0e3a81
'sbg:contributors':
  - rbowen_james
'sbg:createdBy': rbowen_james
'sbg:createdOn': 1633932511
'sbg:id': mwonge/mwtest/hla-reports/2
'sbg:image_url': null
'sbg:latestRevision': 2
'sbg:modifiedBy': rbowen_james
'sbg:modifiedOn': 1633952116
'sbg:project': mwonge/mwtest
'sbg:projectName': zcc-cavatica
'sbg:publisher': sbg
'sbg:revision': 2
'sbg:revisionNotes': null
'sbg:revisionsInfo':
  - 'sbg:modifiedBy': rbowen_james
    'sbg:modifiedOn': 1633932511
    'sbg:revision': 0
    'sbg:revisionNotes': null
  - 'sbg:modifiedBy': rbowen_james
    'sbg:modifiedOn': 1633952040
    'sbg:revision': 1
    'sbg:revisionNotes': null
  - 'sbg:modifiedBy': rbowen_james
    'sbg:modifiedOn': 1633952116
    'sbg:revision': 2
    'sbg:revisionNotes': null
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
'sbg:toolkit': disTIL
'sbg:toolkitVersion': 1.0.0
