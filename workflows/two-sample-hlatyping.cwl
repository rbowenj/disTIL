class: Workflow
cwlVersion: v1.0
id: two_sample_hlatyping
doc: >-
  # About this workflow

  This workflow runs HLA-HD (with restricted reads for speed) on two sets of
  paired-end FASTQs (typically tumour RNA, and normal DNA) and generates a
  consensus type for each HLA allele. See the `hlahd-consensus-parser` tool
  documentation for more information on how a consensus HLA type is determined.


  ## Outputs

  There are 4 main outputs from this tool:

  - **HLA Consensus JSON:** a JSON file containing consensus alleles for **all
  typed HLA genes at the highest possible accuracy level** (two- or
  three-field).

  - **HLA Consensus Text File:** a text file containing all **unique consensus
  alleles** from all HLA genes which were successfully typed by HLA-HD for which
  a consensus could be determined (i.e. 'Not typed' and 'No consensus' alleles
  are removed). The alleles are truncated to two-field accuracy. This file can
  be used as input to pVACtools.

  - **Clinically Significant HLA Consensus JSON:** a JSON file containing
  consensus alleles for **a subset of clinically significant HLA genes (HLA-A,
  -B, -C, -DRA, -DRB1, -DRB3, -DRB4, -DRB5, -DQA1, -DQB1, -DPA1, -DPB1) at the
  highest possible accuracy level** (two- or three-field).

  - **HLA Consensus Text File:** a text file containing all **unique consensus
  alleles** from the clinically significant subset of HLA genes(defined above)
  which were successfully typed by HLA-HD for which a consensus could be
  determined (i.e. 'Not typed' and 'No consensus' alleles are removed). The
  alleles are truncated to two-field accuracy. This file can be used as input to
  pVACtools.


  In addition, a JSON file is generated containing the two sets of input HLA-HD
  results in a format consistent with that of the consensus allele JSONs.


  ## Documentation

  - [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/)
label: two-sample-hlatyping
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: output_prefix
    type: string?
    'sbg:exposed': true
  - id: output_prefix_1
    type: string?
    'sbg:exposed': true
  - id: patient_id
    type: string
    label: Patient ID
    'sbg:x': 0
    'sbg:y': 427.625
  - id: sample1_read2_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 1 Read 2 Sequences
    'sbg:x': 0
    'sbg:y': 213.8125
  - id: sample1_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 1 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 320.71875
  - id: sample2_read2_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 2 Read 2 Sequences
    'sbg:x': 0
    'sbg:y': 0
  - id: sample2_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 2 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 106.90625
  - id: bowtie2_index_prefix
    type: string
    label: HLA Bowtie2 Index Prefix
    doc: >-
      The prefix of every file in the HLA Bowtie2 Index Archive (e.g. `hla_gen`
      if using an index built from IMGT `hla_gen.fasta`).
    'sbg:x': 0
    'sbg:y': 534.53125
  - id: bowtie2_index
    'sbg:fileTypes': TAR
    type: File
    label: HLA Bowtie2 Index Archive
    doc: Bowtie2 Index Archive for an HLA reference.
    'sbg:x': 0
    'sbg:y': 641.4375
outputs:
  - id: sample1_hlahd_output
    outputSource:
      - sample1_hlahd/hlahd_output
    type: Directory
    label: HLA-HD Results from Sample 1
    'sbg:x': 661.2220458984375
    'sbg:y': 278.71875
  - id: sample2_hlahd_output
    outputSource:
      - sample2_hlahd/hlahd_output
    type: Directory
    label: HLA-HD Results from Sample 2
    'sbg:x': 661.2220458984375
    'sbg:y': 171.8125
  - id: sample2_json
    outputSource:
      - hlahd_consensus_parser/sample2_json
    'sbg:fileTypes': JSON
    type: File
    label: Sample 2 HLA-HD Results JSON
    doc: >-
      HLA-HD final results obtained from sample 2 and represented in JSON
      format.
    'sbg:x': 1224.8511962890625
    'sbg:y': 53.453125
  - id: sample1_json
    outputSource:
      - hlahd_consensus_parser/sample1_json
    'sbg:fileTypes': JSON
    type: File
    label: Sample 1 HLA-HD Results JSON
    doc: >-
      HLA-HD final results obtained from sample 1 and represented in JSON
      format.
    'sbg:x': 1224.8511962890625
    'sbg:y': 160.359375
  - id: filtered_txt
    outputSource:
      - hlahd_consensus_parser/filtered_txt
    'sbg:fileTypes': TXT
    type: File
    label: HLA Consensus Text File
    doc: >-
      A text file containing a comma-separated list of consensus alleles for all
      typed HLA genes. Only includes alleles for which a consensus could be
      determined. All alleles are truncated to two-field accuracy. This file
      should be used as input for pVACseq.
    'sbg:x': 1224.8511962890625
    'sbg:y': 267.265625
  - id: clin_sig_json
    outputSource:
      - hlahd_consensus_parser/clin_sig_json
    'sbg:fileTypes': JSON
    type: File
    label: Clinically Significant HLA Consensus JSON
    doc: >-
      A JSON file containing consensus alleles for the clinically significant
      classical HLA genes: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1.
    'sbg:x': 1224.8511962890625
    'sbg:y': 587.984375
  - id: consensus_json
    outputSource:
      - hlahd_consensus_parser/consensus_json
    'sbg:fileTypes': JSON
    type: File
    label: HLA Consensus JSON
    doc: A JSON file containing cosnensus alleles for all typed HLA genes.
    'sbg:x': 1224.8511962890625
    'sbg:y': 374.171875
  - id: clin_sig_txt
    outputSource:
      - hlahd_consensus_parser/clin_sig_txt
    'sbg:fileTypes': TXT
    type: File
    label: Cinically Significant HLA Consensus Text File
    doc: >-
      A text file containing a comma-separated list of consensus alleles for the
      clinically significant classical HLA genes: A, B, C, DRB1, DQA1, DQB1,
      DPA1, DPB1. Only includes alleles for which a consensus could be
      determined. All alleles are truncated to two-field accuracy. This file
      should be used as input for pVACseq.
    'sbg:x': 1224.8511962890625
    'sbg:y': 481.078125
steps:
  - id: sample2_hlahd
    in:
      - id: sample_name
        source: patient_id
      - id: read2_sequences
        source: sample1_read2_sequences
      - id: read1_sequences
        source: sample1_read1_sequences
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      - id: bowtie2_index
        source: bowtie2_index
      - id: output_prefix
        source: output_prefix
    out:
      - id: hla_fastq_reads1
      - id: hla_fastq_reads2
      - id: hlahd_output
      - id: hlahd_final
    run: ./restricted-reads-hlahd.cwl
    label: sample2-hlahd
    'sbg:x': 275.296875
    'sbg:y': 211.265625
  - id: sample1_hlahd
    in:
      - id: sample_name
        source: patient_id
      - id: read2_sequences
        source: sample2_read2_sequences
      - id: read1_sequences
        source: sample2_read1_sequences
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      - id: bowtie2_index
        source: bowtie2_index
      - id: output_prefix
        source: output_prefix_1
    out:
      - id: hla_fastq_reads1
      - id: hla_fastq_reads2
      - id: hlahd_output
      - id: hlahd_final
    run: ./restricted-reads-hlahd.cwl
    label: sample1-hlahd
    'sbg:x': 275.296875
    'sbg:y': 374.171875
  - id: hlahd_consensus_parser
    in:
      - id: tumour_rna
        source: sample1_hlahd/hlahd_final
      - id: normal_dna
        source: sample2_hlahd/hlahd_final
      - id: sample_name
        source: patient_id
    out:
      - id: clin_sig_json
      - id: filtered_txt
      - id: sample1_json
      - id: consensus_json
      - id: sample2_json
      - id: sample3_json
      - id: clin_sig_txt
    run: ../tools/hlahd-consensus-parser.cwl
    label: hlahd-consensus-parser
    'sbg:x': 661.2220458984375
    'sbg:y': 427.625
requirements:
  - class: SubworkflowFeatureRequirement
