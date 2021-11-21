class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: distil_report
baseCommand: []
inputs:
  - id: hla_json
    type: File
    label: HLA Consensus JSON
    doc: >-
      A JSON file containing all HLA consensus alleles. All HLA genes typed by
      HLA-HD are included at two- or three-field accuracies.
    'sbg:fileTypes': JSON
  - id: patient_id
    type: string
    label: Patient ID
    doc: ID for the patient.
  - id: pvacseq_i
    type: File
  - id: pvacseq_ii
    type: File
  - id: pvacfuse_i
    type: File
  - id: pvacfuse_ii
    type: File
  - id: coding_missense_variants
    type: int
  - id: tmb
    type: float
  - id: ipass
    type: File
  - id: epic_deconv
    type: File
  - id: quant_deconv
    type: File
outputs:
  - id: hla_report
    doc: A PDF report containing the HLA consensus results.
    label: HLA Report
    type: File
    outputBinding:
      glob: |-
        ${
            return inputs.patient_id + '_distilReport.html'
        }
    'sbg:fileTypes': PDF
doc: >-
  # About this tool

  This tool generates a disTIL HTML report using the outputs of disTIL's
  analysis modules.


  ## Inputs

  - HLA JSON file

  - pVACseq MHC-I Filtered Report

  - pVACseq MHC-II Filtered Report

  - pVACfuse MHC-I STAR Annotated Filtered Report

  - pVACfuse MHC-II STAR Annotated Filtered Report

  - TMB (string)

  - Number of coding missense variants (string)

  - EPIC Deconvolution TSV

  - quanTIseq Deconvolution TSV

  - IPASS TSV

  - Patient ID
label: distil-report
arguments:
  - position: 2
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          var cmd = 'Rscript -e "rmarkdown::render(\'report_generator.Rmd\',params=list(hla_json=\'./' + inputs.hla_json.basename + '\', pvacseq_i=\'./' + inputs.pvacseq_i.basename + '\', pvacseq_ii=\'./' + inputs.pvacseq_ii.basename + '\', pvacfuse_i=\'./' + inputs.pvacfuse_i.basename + '\', pvacfuse_ii=\'./' + inputs.pvacfuse_ii.basename + '\', tmb=\'' + inputs.tmb + '\', coding_missense_variants=\'' + inputs.coding_missense_variants + '\', ipass=\'' + inputs.ipass.basename + '\',  epic_deconv=\'' + inputs.epic_deconv.basename + '\', quantiseq_deconv=\'' + inputs.quant_deconv.basename +  '\', pid=\'' + inputs.patient_id + '\'), output_file=paste(\'' + inputs.patient_id + '\', \'_distilReport.html\', sep=\'\'))\"'
          return cmd
      }
  - position: 10
    prefix: ''
    shellQuote: false
    valueFrom: 1>&2
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: ' cp /report_generator.Rmd . &&  cp /report_styles.css . && ls 1>&2 &&'
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/distil-report:1.0.0'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.hla_json)
      - $(inputs.pvacseq_i)
      - $(inputs.pvacseq_ii)
      - $(inputs.pvacfuse_i)
      - $(inputs.pvacfuse_ii)
      - $(inputs.ipass)
      - $(inputs.epic_deconv)
      - $(inputs.quant_deconv)
  - class: InlineJavascriptRequirement
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
'sbg:toolkit': disTIL
'sbg:toolkitVersion': 1.0.0
