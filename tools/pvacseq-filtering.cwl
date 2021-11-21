class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: pvacseq_filtering
baseCommand:
  - python
  - /app/pvacseq_filtering.py
inputs:
  - id: pvacseq_report
    type: File
    inputBinding:
      position: 0
      shellQuote: false
      loadContents: true
    label: pVACseq Aggregate Report
    doc: An aggregate report output by pVACseq.
    'sbg:fileTypes': TSV
  - id: patient_id
    type: string
    inputBinding:
      position: 1
      shellQuote: false
    label: Patient ID
    doc: The ID of the patient.
outputs:
  - id: pvacseq_filtered
    doc: A reformatted pVACseq report ranked by IC50.
    label: pVACseq Filtered Neoepitopes
    type: File
    outputBinding:
      glob: '*pvacseq.filtered.csv'
      outputEval: '$(inheritMetadata(self, inputs.pvacseq_report))'
    'sbg:fileTypes': CSV
  - id: pvacseq_filtered_pass
    doc: pVACseq neoepitopes with the 'Pass' tier.
    label: pVACseq Filtered Pass Neoepitopes
    type: File
    outputBinding:
      glob: '*pvacseq.filtered.pass.csv'
      outputEval: '$(inheritMetadata(self, inputs.pvacseq_report))'
    'sbg:fileTypes': CSV
  - id: pvacseq_indels_filtered
    doc: Filtered indel-derived neoepitopes.
    label: pVACseq Filtered Indel Neoepitopes
    type: File
    outputBinding:
      glob: '*pvacseq.indels.filtered.csv'
      outputEval: '$(inheritMetadata(self, inputs.pvacseq_report))'
    'sbg:fileTypes': CSV
  - id: pvacseq_snvs_filtered
    doc: Filtered SNV-derived neoepitopes.
    label: pVACseq Filtered SNV Neoepitopes
    type: File
    outputBinding:
      glob: '*pvacseq.snvs.filtered.csv'
      outputEval: '$(inheritMetadata(self, inputs.pvacseq_report))'
    'sbg:fileTypes': CSV
doc: >-
  # About this tool

  This tool runs filtering and ranking of neopepitopes predicted by pVACseq. The
  main input to this tool is a pVACseq aggregate report.


  ## Reformatting

  - Count the number of the patient's HLA Class I and II molecules each
  predicted neoepitope binds to. Summarise this as '# HLA-I alleles' and '#
  HLA-II alleles'.

  - Add a column called 'DAI' containing the Differential Agretopicity Index for
  each neoepitope, calculated as `ic50_MT/ic50_WT`.

  - Add a column called 'Percentile mutant/ wildtype' calculated as
  'percentile_MT/percentile_WT'.

  - Rename each column to be more readable.

  - Reorder the columns to prioritise gene fusion information.


  ## Filtering and Ranking Steps

  - Only keep neoepitopes which are predicted to bind to at least 1 of the
  patient's HLA Class I or II molecules.

  - Rank the neoepitopes by mutant IC50 ascending (smaller IC50 is better).


  ## Output

  There are 4 output files:

  - `<patient_id>_pvacseq.snvs.filtered.csv`: filtered SNV neoepitopes.

  - `<patient_id>_pvacseq.indels.filtered.csv`: filtered indel neoepitopes.

  - `<patient_id>_pvacseq.filtered.csv`: filtered SNV and indel neoepitopes.

  - `<patient_id>_pvacseq.filtered.pass.csv`: filtered SNV and indel neoepitopes
  **where the tier is 'Pass'**.


  ## Docker

  This tool uses the Docker image `rachelbj/pvactools-filters:1.0`.


  ## Documentation

  - [pVACseq](https://pvactools.readthedocs.io/en/latest/pvacseq.html)

  - [Explanation of pVACseq
  tiers](https://pvactools.readthedocs.io/en/latest/pvacseq/output_files.html)

  - [Differential Agretopicity
  Index](https://rupress.org/jem/article/211/11/2231/41489/Genomic-and-bioinformatic-profiling-of-mutational)
label: pvacseq-filtering
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'rachelbj/pvactools-filters:1.0'
  - class: InlineJavascriptRequirement
    expressionLib:
      - |-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };
