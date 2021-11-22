class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: distil_report
baseCommand: []
inputs:
  - id: hla_json
    type: File?
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
    type: File?
    label: pVACseq MHC-I Filtered Report
  - id: pvacseq_ii
    type: File?
    label: pVACseq MHC-II Filtered Report
  - id: pvacfuse_i
    type: File?
    label: pVACfuse MHC-I Filtered Report
  - id: pvacfuse_ii
    type: File?
    label: pVACfuse MHC-II Filtered Report
  - id: coding_missense_variants
    type: File?
    label: Coding Missense Variants File
  - id: tmb
    type: File?
    label: TMB File
    'sbg:fileTypes': TXT
  - id: ipass
    type: File?
    label: IPASS Results
  - id: epic_deconv
    type: File?
    label: EPIC Deconvolution Results
  - id: quant_deconv
    type: File?
    label: quanTIseq Deconvolution Results
outputs:
  - id: disTIL_report
    doc: An HTML report containing the summarised results of immunoprofiling.
    label: disTIL Immunoprofiling Report
    type: File
    outputBinding:
      glob: |-
        ${
            return inputs.patient_id + '_distilReport.html'
        }
    'sbg:fileTypes': HTML
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
          var cmd = 'Rscript -e "rmarkdown::render(\'report_generator.Rmd\',params=list(pid=\'' + inputs.patient_id + '\''
          
          
          if (inputs.hla_json) {
              cmd = cmd + ', hla_json=\'./' + inputs.hla_json.basename + '\''
          }
          if (inputs.pvacseq_i) {
              cmd = cmd + ', pvacseq_i=\'./' + inputs.pvacseq_i.basename + '\''
          }
          if (inputs.pvacseq_ii) {
              cmd = cmd + ', pvacseq_ii=\'./' + inputs.pvacseq_ii.basename + '\''
          }
          if (inputs.pvacfuse_i) {
              cmd = cmd + ', pvacfuse_i=\'./' + inputs.pvacfuse_i.basename + '\''
          }
          if (inputs.pvacfuse_ii) {
              cmd = cmd + ', pvacfuse_ii=\'./' + inputs.pvacfuse_ii.basename + '\''
          }
          if (inputs.coding_missense_variants) {
              cmd = cmd + ', coding_missense_variants=\'./' + inputs.coding_missense_variants.basename + '\''
          }
          if (inputs.tmb) {
              cmd = cmd + ', tmb=\'./' + inputs.tmb.basename + '\''
          }
          if (inputs.ipass) {
              cmd = cmd + ', ipass=\'./' + inputs.ipass.basename + '\''
          }
          if (inputs.epic_deconv) {
              cmd = cmd + ', epic_deconv=\'./' + inputs.epic_deconv.basename + '\''
          }
          if (inputs.quantiseq_deconv) {
              cmd = cmd + ', quantiseq_deconv=\'./' + inputs.quantiseq_deconv.basename + '\''
          }
          
          cmd = cmd + '), output_file=paste(\'' + inputs.patient_id + '\', \'_distilReport.html\', sep=\'\'))\"'

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
      - $(inputs.tmb)
      - $(inputs.coding_missense_variants)
  - class: InlineJavascriptRequirement
'sbg:toolkit': disTIL
'sbg:toolkitVersion': 1.0.0
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
