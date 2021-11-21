class: Workflow
cwlVersion: v1.0
id: vcf_readcount_annotation
doc: >-
  # About this workflow

  This workflow runs readcount annotation of a VCF. This is part of the file
  preparation process for pVACseq.  

  The following tool versions are used in this workflow:

  - bam-readcount-helper v1.1.1 (uses bam-readcount v1.0.0)

  - vatools vcf-readcount-annotator v5.0.1



  ## Inputs and parameters

  To run this workflow, you will need:

  - **Input VCF:** An input VCF to be annotated with readcounts. Note that this
  VCF should be decomposed prior to running this analysis (see vt-decompose).

  - **Sample Name:** The name of the sample to annotate in the input VCF. Note
  that this string must exactly match the sample name in the VCF.

  -  **Input BAM:** An input BAM containing DNA or RNA data for the sample to be
  annotated.

  - **Reference Genome:** The reference genome FASTA used to generate the input
  BAM.

  - **Genomic Data Type:** The data type contained in the intput BAM (either DNA
  or RNA). Note that this parameter determines whether the DNA or RNA readcount
  fields will be annotated in the VCF.


  Optionally, you can include minimum base quality, minimum mapping wuality, and
  an interval string for use in running bam-readcount. For more information, see
  the bam-readcount documentation (below).


  ## Steps

  This workflow runs the following steps:

  1. Generates SNV and indel readcounts from a BAM using **bam-readcount** via
  an adapted version of the bam-readcount-helper from mgibio.

  2. Annotates the VCF with SNV readcounts using **vatools
  vcf-readcount-annotator**.

  2. Annotates the VCF with indel readcounts using **vatools
  vcf-readcount-annotator**.


  ## Outputs

  Depending on the data type (DNA/RNA) and variant type (SNV/indel/all), the
  following extensions are added to the annotated output file name:

  - `dsr`: **D**NA **S**NV **r**eadcounts

  - `dir`: **D**NA **i**ndel **r**eadcounts

  - `dr`: **D**NA **r**eadcounts (SNVs and indels)

  - `rsr`: **R**NA **S**NV **r**eadcounts

  - `rir`: **R**NA **i**ndel **r**eadcounts

  - `rr`: **R**NA **r**eadcounts (SNVs and indels)


  ## Documentation

  - [bam-readcount](https://github.com/genome/bam-readcount)

  -
  [bam-readcount-helper](https://github.com/genome/docker-bam_readcount_helper-cwl)

  - [vatools
  vcf-readcount-annotator](https://vatools.readthedocs.io/en/latest/vcf_readcount_annotator.html)
label: vcf-readcount-annotation
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: ref_genome
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: Reference Genome
    doc: FASTA reference genome used for the generation of the input BAM.
    secondaryFiles:
      - .fai
    'sbg:x': 0
    'sbg:y': 102.4296875
  - id: sample_name
    type: string
    label: Sample Name
    doc: Name of the sample to annotate in the VCF.
    'sbg:x': 408.390625
    'sbg:y': 88.4296875
  - id: input_bam
    'sbg:fileTypes': BAM
    type: File
    label: Input BAM
    doc: BAM file to be used for annotation of readcounts.
    secondaryFiles:
      - .bai
    'sbg:x': 0
    'sbg:y': 316.1484375
  - id: data_type
    type:
      - 'null'
      - type: enum
        symbols:
          - DNA
          - RNA
        name: data_type
    label: Genomic Data Type
    doc: >-
      Either DNA or RNA, depending on whether the input BAM file contains DNA or
      RNA data.
    'sbg:x': 408.390625
    'sbg:y': 195.2890625
  - id: min_base_qual
    type: int?
    'sbg:exposed': true
  - id: min_mapping_qual
    type: int?
    'sbg:exposed': true
  - id: input_vcf
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    label: Input VCF
    doc: >-
      The VCF to annotate with readcounts. Note that this VCF should be
      decomposed prior to running this workflow (see vt-decompose).
    secondaryFiles:
      - .tbi
    'sbg:x': -88.7520751953125
    'sbg:y': 235.1882781982422
outputs:
  - id: snv_indel_annot_zipped
    outputSource:
      - bgzip_tabix_2/zipped_with_index
    'sbg:fileTypes': VCF.GZ
    type: File
    label: SNV and Indel Annotated VCF
    doc: VCF annotated with SNV and indel readcounts from the input BAM.
    secondaryFiles:
      - .tbi
    'sbg:x': 1859.51953125
    'sbg:y': 209.2890625
  - id: snv_readcount
    outputSource:
      - bam_readcount_helper/snv_readcount
    'sbg:fileTypes': TSV
    type: File
    label: SNV Readcount
    doc: SNV readcounts from the input BAM.
    'sbg:x': 1006.9417724609375
    'sbg:y': 88.4296875
  - id: indel_readcount
    outputSource:
      - bam_readcount_helper/indel_readcount
    'sbg:fileTypes': TSV
    type: File
    label: Indel Readcount
    doc: Indel readcounts from the input BAM.
    'sbg:x': 1006.9417724609375
    'sbg:y': 195.2890625
steps:
  - id: snv_readcount_annot
    in:
      - id: input_vcf
        source: input_vcf
      - id: bam_readcount_file
        source: bam_readcount_helper/snv_readcount
      - id: sample_name
        source: sample_name
      - id: variant_type
        default: snv
      - id: data_type
        source: data_type
    out:
      - id: annotated_vcf
    run: ../tools/vatools-readcount-annotation.cwl
    label: readcount-annot-snv
    doc: VAtools annotation of the input VCF with SNV readcounts.
    'sbg:x': 683.0892333984375
    'sbg:y': 0
  - id: bgzip_tabix_1
    in:
      - id: input_file
        source: snv_readcount_annot/annotated_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': 1006.9417724609375
    'sbg:y': 316.1484375
  - id: indel_readcount_annot
    in:
      - id: input_vcf
        source: bgzip_tabix_1/zipped_with_index
      - id: bam_readcount_file
        source: bam_readcount_helper/indel_readcount
      - id: sample_name
        source: sample_name
      - id: variant_type
        default: indel
      - id: data_type
        source: data_type
    out:
      - id: annotated_vcf
    run: ../tools/vatools-readcount-annotation.cwl
    label: readcount-annot-indel
    doc: VAtools annotation of the input VCF with indel readcounts.
    'sbg:x': 1281.640380859375
    'sbg:y': 262.71875
  - id: bgzip_tabix_2
    in:
      - id: input_file
        source: indel_readcount_annot/annotated_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': 1584.8209228515625
    'sbg:y': 195.2890625
  - id: bam_readcount_helper
    in:
      - id: input_vcf
        source: input_vcf
      - id: sample_name
        source: sample_name
      - id: reference_fasta
        source: ref_genome
      - id: input_bam
        source: input_bam
      - id: data_type
        source: data_type
      - id: min_base_qual
        source: min_base_qual
      - id: min_mapping_qual
        source: min_mapping_qual
    out:
      - id: indel_readcount
      - id: snv_readcount
    run: ../tools/bam-readcount-helper.cwl
    label: bam-readcount-helper
    doc: Run bam-readcount for SNVs and indels.
    'sbg:x': 683.0892333984375
    'sbg:y': 369.578125
requirements:
  - class: SubworkflowFeatureRequirement
'sbg:toolAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
