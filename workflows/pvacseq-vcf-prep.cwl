class: Workflow
cwlVersion: v1.0
id: pvacseq_vcf_prep
doc: >-
  # About this workflow

  This workflow runs the VCF preparation steps necessary for pVACseq. This
  includes decomposition of the VCF, followed by annotation with DNA and RNA BAM
  readcounts, and annotation with genome and transcript level expression values
  (from prior RNA analysis).  


  The following tool versions are used in this workflow:

  - vep v104.3

  - vt-decompose v0.57721

  - bam-readcount-helper v1.1.1 (uses bam-readcount v1.0.0)

  - vatools vcf-readcount-annotator v5.0.1

  - vatools vcf-expression-annotator v5.0.1


  ## Before you run this workflow

  Note that this workflow requires a range of DNA and RNA inputs, as well as RNA
  expression values.  

  To use this workflow, you will need the following:

  - Inputs for VEP annotation:

  1. **VEP Cache:** A vep cache TAR file.

  2. **VEP Cache Version:** The version of the vep cache TAR file.

  3.**VEP Plugin Files:** Plugin files to be used in VEP annotation (for
  pVACseq, you must provide Wildtype.pm and Frameshift.pm plugin files).

  - Inputs for readcount annotation

  1. **Input VCF:** An input VCF to be annotated with readcounts. Note that this
  VCF does **not** need to be decomposed prior to running this analysis.

  2. **Sample Name:** The name of the sample to annotate in the input VCF. Note
  that this string must exactly match the sample name in the VCF.

  3. **DNA BAM:** A BAM containing DNA data for the sample of interest.

  4. **DNA Reference Genome:** The reference genome FASTA used to align the DNA
  BAM.

  5. **RNA BAM:** A BAM containing RNA data for the sample of interest.

  6. **RNA Reference Genome:** The reference genome FASTA used to align the RNA
  BAM.

  - Inputs for expression annotation

  1. **Gene Level Expression File:** A TSV file containing gene level expression
  values.

  2. **Gene Expression Quntification Algorithm:** The algorithm used to quantify
  gene level expression (and generate the gene level expression file).

  3. **Gene ID Column:** The name of the column containing gene IDs in the gene
  expression file.

  4. **Gene Expression Column:** The name of the column containing expression
  values in the gene expression file.

  5. **Transcript Level Expression File:** A TSV file containing gene level
  expression values.

  6. **Transcript Expression Quntification Algorithm:** The algorithm used to
  quantify transcript level expression (and generate the transcript level
  expression file).

  7. **Transcript ID Column:** The name of the column containing transcript IDs
  in the transcript expression file.

  8. **Transcript Expression Column:** The name of the column containing
  expression values in the transcript expression file.



  **If you do not have all of these inputs available, you can run components of
  this analysis using the constituent CWL tools and subworkflows (all of which
  are provided in this repo).**


  ### RNA expression values

  Genome and transcript level expression values in TSV format must be generated
  before running this analysis. Any tool can be used to generate these
  expression values, but note that:

  - If you use a tool other than Kallisto, Stringtie or Cufflinks, you need to
  specify the names of the columns containing gene/transcript IDs and expression
  values.

  - The expression values must be in TSV file format.

  - The gene/transcript ID format used.


  ## Steps

  This workflow runs the following steps:

  1.  Filters the input VCF to retain only PASS variants.

  2. Annotates the input VCF using **vep**.

  3. Decomposes the input VCF using **vt-decompose**.

  4. Generates SNV and indel readcounts from a BAM using **bam-readcount** via
  an adapted version of the bam-readcount-helper from mgibio.

  5. Annotates the VCF with SNV readcounts using **vatools
  vcf-readcount-annotator**.

  6. Annotates the VCF with indel readcounts using **vatools
  vcf-readcount-annotator**.

  7. Annotates the VCF with gene level expression values using **vatools
  vcf-expression-annotator**.

  8. Annotates the VCF with transcript level expression values using **vatools
  vcf-expression-annotator**.


  ## Documentation

  - [vep](https://www.ensembl.org/info/docs/tools/vep/index.html)

  - [vt-decompose](https://genome.sph.umich.edu/wiki/Vt)

  - [bam-readcount](https://github.com/genome/bam-readcount)

  -
  [bam-readcount-helper](https://github.com/genome/docker-bam_readcount_helper-cwl)

  - [vatools
  vcf-readcount-annotator](https://vatools.readthedocs.io/en/latest/vcf_readcount_annotator.html)

  - [vatools
  vcf-expression-annotator](https://vatools.readthedocs.io/en/latest/vcf_expression_annotator.html)
label: pvacseq-vcf-prep
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: min_base_qual
    type: int?
    'sbg:exposed': true
  - id: min_mapping_qual
    type: int?
    'sbg:exposed': true
  - id: sample_name
    type: string
    label: Sample Name
    doc: Name of the sample to annotate in the input VCF file.
    'sbg:x': 1585.253173828125
    'sbg:y': 113.6875
  - id: ref_genome_dna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: DNA Reference Genome
    doc: Reference sequence used to align the DNA BAM.
    secondaryFiles:
      - .fai
    'sbg:x': 466.84375
    'sbg:y': 220.375
  - id: ref_genome_rna
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: RNA Reference Genome
    doc: Reference sequence used to align the RNA BAM.
    secondaryFiles:
      - .fai
    'sbg:x': 1859.95166015625
    'sbg:y': 234.3984375
  - id: input_vcf
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    label: Input VCF
    doc: >-
      The VCF to be decomposed then annotated with DNA and RNA BAM readcounts.
      Note that this VCF does not need to be decomposed prior to running this
      analysis.
    secondaryFiles:
      - .tbi
    'sbg:x': 0
    'sbg:y': 234.3984375
  - id: input_bam_rna
    'sbg:fileTypes': BAM
    type: File
    label: RNA BAM
    doc: The RNA BAM for the sample of interest.
    secondaryFiles:
      - .bai
    'sbg:x': 1859.95166015625
    'sbg:y': 341.0859375
  - id: input_bam_dna
    'sbg:fileTypes': BAM
    type: File
    label: DNA BAM
    doc: The DNA BAM for the sample of interest.
    secondaryFiles:
      - .bai
    'sbg:x': 1585.253173828125
    'sbg:y': 220.375
  - id: intervals_string
    type: string?
    'sbg:exposed': true
  - id: id_column
    type: string?
    'sbg:exposed': true
  - id: expression_column
    type: string?
    'sbg:exposed': true
  - id: gene_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Gene Expression File
    doc: Gene level expression file.
    'sbg:x': 0
    'sbg:y': 341.0859375
  - id: gene_quant_algo
    type:
      type: enum
      symbols:
        - kallisto
        - stringtie
        - cufflinks
        - custom
      name: gene_quant_algo
    label: Gene Expression Quantification Algorithm
    doc: >-
      The expression quantification algorithm used to generate gene level
      expression values.
    'sbg:x': 2255.010009765625
    'sbg:y': 341.0859375
  - id: transcript_expression_file
    'sbg:fileTypes': TSV
    type: File
    label: Transcript Expression File
    doc: Transcript level expression file.
    'sbg:x': 0
    'sbg:y': 127.7109375
  - id: transcript_id_column
    type: string?
    'sbg:exposed': true
  - id: transcript_expression_column
    type: string?
    'sbg:exposed': true
  - id: transcript_quant_algo
    type:
      type: enum
      symbols:
        - kallisto
        - stringtie
        - cufflinks
        - custom
      name: transcript_quant_algo
    label: Transcript Expression Quantification Algorithm
    doc: >-
      The expression quantification algorithm used to generate transcript level
      expression values.
    'sbg:x': 3068.954833984375
    'sbg:y': 167.0078125
  - id: vep_cache
    'sbg:fileTypes': TAR.GZ
    type: File
    label: VEP Cache
    doc: The VEP cache TAR file.
    'sbg:x': 466.84375
    'sbg:y': 113.6875
  - id: vep_plugin_files
    'sbg:fileTypes': PM
    type: 'File[]'
    label: VEP Plugin Files
    doc: >-
      Plugin files to use in VEP annotation (for pVACseq, must use Wildtype and
      Frameshift).
    'sbg:x': 466.84375
    'sbg:y': 6.9765625
  - id: cache_version
    type: int
    label: VEP Cache Version
    doc: VEP cache version.
    'sbg:x': 466.84375
    'sbg:y': 327.0625
outputs:
  - id: dna_annot_zipped
    outputSource:
      - vcf_readcount_annotation_dna/snv_indel_annot_zipped
    'sbg:fileTypes': VCF.GZ
    type: File
    label: DNA Readcount Annotated VCF
    doc: Decomposed VCF annotated with SNV and indel readcounts from the DNA BAM.
    secondaryFiles:
      - .tbi
    'sbg:x': 2255.010009765625
    'sbg:y': 447.796875
  - id: dna_rna_annot_zipped
    outputSource:
      - vcf_readcount_annotation_rna/snv_indel_annot_zipped
    'sbg:fileTypes': VCF.GZ
    type: File
    label: DNA and RNA Readcount Annotated VCF
    doc: >-
      Decomposed VCF annotated with SNV and indel readcounts from the DNA and
      RNA BAMs.
    secondaryFiles:
      - .tbi
    'sbg:x': 2650.068359375
    'sbg:y': 394.4296875
  - id: snv_readcount_rna
    outputSource:
      - vcf_readcount_annotation_rna/snv_readcount
    'sbg:fileTypes': TSV
    type: File
    label: RNA SNV Readcount
    doc: SNV readcounts generated from the RNA BAM.
    'sbg:x': 2650.068359375
    'sbg:y': 181.0546875
  - id: indel_readcount_rna
    outputSource:
      - vcf_readcount_annotation_rna/indel_readcount
    'sbg:fileTypes': TSV
    type: File
    label: RNA Indel Readcount
    doc: Indel readcounts generated from the RNA BAM.
    'sbg:x': 2650.068359375
    'sbg:y': 287.7421875
  - id: snv_readcount_dna
    outputSource:
      - vcf_readcount_annotation_dna/snv_readcount
    'sbg:fileTypes': TSV
    type: File
    label: DNA SNV Readcount
    doc: SNV readcounts generated from the DNA BAM.
    'sbg:x': 2255.010009765625
    'sbg:y': 127.6875
  - id: indel_readcount_dna
    outputSource:
      - vcf_readcount_annotation_dna/indel_readcount
    'sbg:fileTypes': TSV
    type: File
    label: DNA Indel Readcount
    doc: Indel readcounts generated from the DNA BAM.
    'sbg:x': 2255.010009765625
    'sbg:y': 234.375
  - id: rc_exp_annotated_vcf
    outputSource:
      - bgzip_tabix_2/zipped_with_index
    'sbg:fileTypes': VCF.GZ
    type: File
    label: Readcount and Expression Annotated VCF
    doc: >-
      Decomposed VCF annotated with DNA and RNA readcounts as well as gene and
      transcript level expression.
    secondaryFiles:
      - .tbi
    'sbg:x': 4055.4619140625
    'sbg:y': 234.3984375
  - id: vep_annot
    outputSource:
      - bgzip_tabix_3/zipped_with_index
    'sbg:fileTypes': VCF.GZ
    type: File
    label: VEP Annotated VCF
    doc: VEP annotated VCF.
    secondaryFiles:
      - .tbi
    'sbg:x': 1343.003173828125
    'sbg:y': 287.7421875
steps:
  - id: vcf_readcount_annotation_dna
    in:
      - id: ref_genome
        source: ref_genome_dna
      - id: sample_name
        source: sample_name
      - id: input_bam
        source: input_bam_dna
      - id: data_type
        default: DNA
      - id: min_base_qual
        source: min_base_qual
      - id: min_mapping_qual
        source: min_mapping_qual
      - id: input_vcf
        source: bgzip_tabix/zipped_with_index
    out:
      - id: snv_indel_annot_zipped
      - id: snv_readcount
      - id: indel_readcount
    run: ./vcf-readcount-annotation.cwl
    label: vcf-readcount-annotation-dna
    doc: Annotation of the VCF with readcounts from a DNA BAM.
    'sbg:x': 1859.95166015625
    'sbg:y': 106.7109375
  - id: vcf_readcount_annotation_rna
    in:
      - id: ref_genome
        source: ref_genome_rna
      - id: sample_name
        source: sample_name
      - id: input_bam
        source: input_bam_rna
      - id: data_type
        default: RNA
      - id: min_base_qual
        source: min_base_qual
      - id: min_mapping_qual
        source: min_mapping_qual
      - id: input_vcf
        source: vcf_readcount_annotation_dna/snv_indel_annot_zipped
    out:
      - id: snv_indel_annot_zipped
      - id: snv_readcount
      - id: indel_readcount
    run: ./vcf-readcount-annotation.cwl
    label: vcf-readcount-annotation-rna
    doc: Annotation of the VCF with readcounts from an RNA BAM.
    'sbg:x': 2255.010009765625
    'sbg:y': 0
  - id: vt_decompose
    in:
      - id: seq_regions
        default: true
      - id: input_vcf
        source: bgzip_tabix_3/zipped_with_index
      - id: intervals_string
        source: intervals_string
    out:
      - id: decomposed_vcf
    run: ../tools/vt-decompose.cwl
    label: vt-decompose
    doc: Decompose the VCF file prior to readcount annotation.
    'sbg:x': 1343.003173828125
    'sbg:y': 181.0546875
  - id: bgzip_tabix
    in:
      - id: input_file
        source: vt_decompose/decomposed_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': 1585.253173828125
    'sbg:y': 341.0859375
  - id: vatools_gene_expression_annotation
    in:
      - id: input_vcf
        source: vcf_readcount_annotation_rna/snv_indel_annot_zipped
      - id: expression_file
        source: rsem_ensembl_annotation_2/output
      - id: quant_algo
        source: gene_quant_algo
      - id: expression_level
        default: gene
      - id: id_column
        source: id_column
      - id: expression_column
        source: expression_column
      - id: sample_name
        source: sample_name
    out:
      - id: exp_annotated_vcf
    run: ../tools/vatools-expression-annotation.cwl
    label: vatools-gene-expression-annotation
    doc: VCF annotation with gene level expression values.
    'sbg:x': 2650.068359375
    'sbg:y': 53.34375
  - id: vatools_transcript_expression_annotation
    in:
      - id: input_vcf
        source: bgzip_tabix_1/zipped_with_index
      - id: expression_file
        source: rsem_ensembl_annotation/output
      - id: quant_algo
        source: transcript_quant_algo
      - id: expression_level
        default: transcript
      - id: id_column
        source: transcript_id_column
      - id: expression_column
        source: transcript_expression_column
      - id: sample_name
        source: sample_name
    out:
      - id: exp_annotated_vcf
    run: ../tools/vatools-expression-annotation.cwl
    label: vatools-transcript-expression-annotation
    doc: VCF annotation with transcript level expression values.
    'sbg:x': 3361.876708984375
    'sbg:y': 213.3984375
  - id: bgzip_tabix_1
    in:
      - id: input_file
        source: vatools_gene_expression_annotation/exp_annotated_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': 3068.954833984375
    'sbg:y': 287.7421875
  - id: bgzip_tabix_2
    in:
      - id: input_file
        source: vatools_transcript_expression_annotation/exp_annotated_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': 3780.76318359375
    'sbg:y': 220.375
  - id: vep_with_plugins
    in:
      - id: input_file
        source: bgzip_tabix_4/zipped_with_index
      - id: vep_cache
        source: vep_cache
      - id: ref_genome
        source: ref_genome_dna
      - id: vep_plugin_files
        source:
          - vep_plugin_files
      - id: cache_version
        source: cache_version
      - id: merged
        default: true
      - id: symbol
        default: true
      - id: biotype
        default: true
      - id: numbers
        default: true
      - id: canonical
        default: true
      - id: total_length
        default: true
      - id: sift
        default: b
      - id: polyphen
        default: b
      - id: terms
        default: SO
    out:
      - id: vep_vcf
      - id: vep_stats
    run: ../tools/vep-with-plugins.cwl
    label: vep-with-plugins
    doc: Run VEP annotation of the input VCF.
    'sbg:x': 741.5423583984375
    'sbg:y': 206.375
  - id: bgzip_tabix_3
    in:
      - id: input_file
        source: vep_with_plugins/vep_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': 1068.3045654296875
    'sbg:y': 220.375
  - id: rsem_ensembl_annotation_2
    in:
      - id: rsem_results
        source: gene_expression_file
    out:
      - id: output
    run: ../tools/rsem-ensembl-annotation.cwl
    label: rsem-ensembl-annotation
    'sbg:x': 254.53125
    'sbg:y': 127.7109375
  - id: rsem_ensembl_annotation
    in:
      - id: rsem_results
        source: transcript_expression_file
    out:
      - id: output
    run: ../tools/rsem-ensembl-annotation.cwl
    label: rsem-ensembl-annotation
    'sbg:x': 254.53125
    'sbg:y': 234.3984375
  - id: bcftools_view_pass
    in:
      - id: vcf
        source: input_vcf
    out:
      - id: pass_filtered_vcf
    run: ../tools/bcftools-view-pass.cwl
    label: bcftools-view-pass
    'sbg:x': 254.53125
    'sbg:y': 341.0859375
  - id: bgzip_tabix_4
    in:
      - id: input_file
        source: bcftools_view_pass/pass_filtered_vcf
    out:
      - id: index_file
      - id: bgzipped_file
      - id: zipped_with_index
    run: ./bgzip-tabix.cwl
    label: bgzip-tabix
    'sbg:x': 466.84375
    'sbg:y': 447.7734375
requirements:
  - class: SubworkflowFeatureRequirement
'sbg:wrapperAuthor': Rachel Bowen-James <rbowen-james@ccia.org.au>
