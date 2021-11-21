class: Workflow
cwlVersion: v1.0
id: three_sample_hlatyping
doc: >-
  # About this workflow

  This workflow runs HLA-HD (with restricted reads for speed) on three sets of
  paired-end FASTQs (tumour DNA, tumour RNA, and normal DNA) and generates a
  consensus type for each HLA allele.  See the `hlahd-consensus-parser` tool
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


  In addition, a JSON file is generated containing the two/three sets of input
  HLA-HD results in a format consistent with that of the consensus allele JSONs.


  ## Documentation

  - [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/)
label: three-sample-hlatyping
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: sample_1_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 1 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 106.578125
  - id: bowtie2_index
    'sbg:fileTypes': TAR
    type: File
    label: HLA Bowtie2 Index Archive
    doc: Bowtie2 Index Archive for an HLA reference.
    'sbg:x': 0
    'sbg:y': 852.4296875
  - id: sample_3_read2_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 3 Read 2 Sequences
    'sbg:x': 0
    'sbg:y': 213.15625
  - id: sample_3_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 3 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 319.734375
  - id: sample_2_read2_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 2 Read 2 Sequences
    'sbg:x': 0
    'sbg:y': 532.8125
  - id: sample_2_read1_sequences
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
    type: File
    label: Sample 2 Read 1 Sequences
    'sbg:x': 0
    'sbg:y': 639.390625
  - id: patient_id
    type: string
    label: Patient ID
    'sbg:x': 0
    'sbg:y': 426.2734375
  - id: bowtie2_index_prefix
    type: string
    label: HLA Bowtie2 Index Prefix
    'sbg:x': 0
    'sbg:y': 745.9296875
outputs:
  - id: tumour_rna_results
    outputSource:
      - sample_1_hlahd/hlahd_output
    type: Directory
    label: HLA-HD Results from Sample 1
    doc: Output directory from running HLA-HD on sample 1.
    'sbg:x': 689.6439208984375
    'sbg:y': 224.234375
  - id: tumour_dna_results
    outputSource:
      - sample_3_hlahd/hlahd_output
    type: Directory
    label: HLA-HD Results from Sample 3
    doc: Output directory from running HLA-HD on sample 3.
    'sbg:x': 689.6439208984375
    'sbg:y': 330.734375
  - id: normal_dna_results
    outputSource:
      - sample_2_hlahd/hlahd_output
    type: Directory
    label: HLA-HD Results from Sample 2
    doc: Output directory from running HLA-HD on sample 2.
    'sbg:x': 689.6439208984375
    'sbg:y': 437.234375
  - id: clin_sig_json
    outputSource:
      - hlahd_consensus_parser/clin_sig_json
    'sbg:fileTypes': JSON
    type: File
    label: Clinically Significant HLA Consensus JSON
    'sbg:x': 1265.2135009765625
    'sbg:y': 746.2421875
  - id: consensus_txt
    outputSource:
      - hlahd_consensus_parser/filtered_txt
    'sbg:fileTypes': TXT
    type: File
    label: HLA Consensus Text File
    'sbg:x': 1265.2135009765625
    'sbg:y': 426.1171875
  - id: sample1_json
    outputSource:
      - hlahd_consensus_parser/sample1_json
    'sbg:fileTypes': JSON
    type: File
    label: Sample 1 HLA-HD Results JSON
    'sbg:x': 1265.2135009765625
    'sbg:y': 319.4609375
  - id: consensus_json
    outputSource:
      - hlahd_consensus_parser/consensus_json
    'sbg:fileTypes': JSON
    type: File
    label: HLA Consensus JSON
    'sbg:x': 1265.2135009765625
    'sbg:y': 532.7734375
  - id: sample2_json
    outputSource:
      - hlahd_consensus_parser/sample2_json
    'sbg:fileTypes': JSON
    type: File
    label: Sample 2 HLA-HD Results JSON
    'sbg:x': 1265.2135009765625
    'sbg:y': 212.7265625
  - id: sample3_json
    outputSource:
      - hlahd_consensus_parser/sample3_json
    'sbg:fileTypes': JSON
    type: File?
    label: Sample 3 HLA-HD Results JSON
    'sbg:x': 1265.2135009765625
    'sbg:y': 105.9921875
  - id: clin_sig_txt
    outputSource:
      - hlahd_consensus_parser/clin_sig_txt
    'sbg:fileTypes': TXT
    type: File
    label: Clinically Significant HLA Consensus Text File
    'sbg:x': 1265.2135009765625
    'sbg:y': 639.5078125
steps:
  - id: sample_1_hlahd
    in:
      - id: sample_name
        source: patient_id
      - id: read1_sequences
        source: sample_1_read1_sequences
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      - id: bowtie2_index
        source: bowtie2_index
      - id: output_prefix
        default: sample_1
    out:
      - id: hla_fastq_reads1
      - id: hla_fastq_reads2
      - id: hlahd_output
      - id: hlahd_final
    run: ./restricted-reads-hlahd.cwl
    label: sample-1-hlahd
    doc: Run restricted reads HLA-HD typing from tumour RNA.
    'sbg:x': 303.71875
    'sbg:y': 235.734375
  - id: sample_3_hlahd
    in:
      - id: sample_name
        source: patient_id
      - id: read2_sequences
        source: sample_3_read2_sequences
      - id: read1_sequences
        source: sample_3_read1_sequences
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      - id: bowtie2_index
        source: bowtie2_index
      - id: output_prefix
        default: sample_3
    out:
      - id: hla_fastq_reads1
      - id: hla_fastq_reads2
      - id: hlahd_output
      - id: hlahd_final
    run: ./restricted-reads-hlahd.cwl
    label: sample-3-hlahd
    doc: Run restricted reads HLA-HD typing from tumour DNA.
    'sbg:x': 303.71875
    'sbg:y': 398.234375
  - id: sample_2_hlahd
    in:
      - id: sample_name
        source: patient_id
      - id: read2_sequences
        source: sample_2_read2_sequences
      - id: read1_sequences
        source: sample_2_read1_sequences
      - id: bowtie2_index_prefix
        source: bowtie2_index_prefix
      - id: bowtie2_index
        source: bowtie2_index
      - id: output_prefix
        default: sample_2
    out:
      - id: hla_fastq_reads1
      - id: hla_fastq_reads2
      - id: hlahd_output
      - id: hlahd_final
    run: ./restricted-reads-hlahd.cwl
    label: sample-2-hlahd
    doc: Run restricted reads HLA-HD typing from normal DNA.
    'sbg:x': 303.71875
    'sbg:y': 560.734375
  - id: hlahd_consensus_parser
    in:
      - id: tumour_rna
        source: sample_1_hlahd/hlahd_final
      - id: tumour_dna
        source: sample_3_hlahd/hlahd_final
      - id: normal_dna
        source: sample_2_hlahd/hlahd_final
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
    'sbg:x': 689.6439208984375
    'sbg:y': 585.984375
requirements:
  - class: SubworkflowFeatureRequirement
'sbg:toolAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
