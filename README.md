# The *disTIL* immunoprofiling toolkit

*disTIL* is a bioinformatics toolkit that integrates multiple analytical tools to enable the efficient and comprehensive characterisation of the tumour immune landscape. By incorporating current best-in-class algorithms into a unified, flexible, and cloud-agnostic workflow, *disTIL* facilitates scalable, reproducible analyses of the tumour immune microenvironment (TIME) and neoantigen landscape of a tumour. *disTIL* supports reference genome builds GRCh37 and GRCh38.

## Implementation

*disTIL* is implemented as a collection of [Common Workflow Language](https://www.commonwl.org/) (CWL) v1.0 tools and workflows, each of which can be run as a standalone analysis. Included in this repository are all CWL tools (in the `tools` directory) and subworkflows used in the toolkit, as well as a selection of end-to-end *disTIL* analysis workflows. Since CWL v1.0 does not support conditional workflows, multiple versions are supplied in the `workflows` directory. The `workflows` directory also contains subworkflows used in the *disTIL* modules. Additionally, the `docker` directory contains the Dockerfiles for some of the CWL tools in the toolkit.

All *disTIL* tools and subworkflows can be run as standalone analyses or be integrated into other CWL workflows. If you do not have all necessary patient data to run the full *disITL* workflow, you can select a subset of the provided tools/modules and create a customised *disTIL* workflow (we recommend using [Rabix Composer](https://rabix.io/) to do this).

## License
The *disTIL* toolkit is licensed under the [Artistic 2.0 license](https://opensource.org/licenses/Artistic-2.0) (see `LICENSE` for further detail).

## Prerequisites 

To run the *disTIL* toolkit, you need a CWL runner. We recommend using [cwltool](https://github.com/common-workflow-language/cwltool).
Clone this repository in order to locally modify and run CWL tools and workflows.

## Modules

*disTIL* consists of five analysis modules and a reporting module, together enabling the end-to-end evaluation and presentation of a range of tumour immune parameters. The following five analysis modules are included:

- **Consensus HLA typing** using [HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) v1.4.0
- **Neoepitope prediction** using [pVACtools](https://github.com/griffithlab/pVACtools) v2.0.4
- **Tumour mutational burden (TMB) calculation**
- **Gene expression classification** using IPASS (unpublished) - note that this module is specific to the analysis of paediatric cancer
- **Immune cell type deconvolution** using [EPIC](https://github.com/GfellerLab/EPIC) and [quanTIseq](https://github.com/icbi-lab/quanTIseq) via the [immunedeconv](https://github.com/icbi-lab/immunedeconv) R package

To summarise the results of analysis in an interpretable format, *disTIL*'s reporting module generates an interactive HTML report using knitr and Rmarkdown. A sample HTML report is available in the `docker/distil_report` directory. The *disTIL* HLA typing module also outputs a PDF HLA report (sample available in the `docker/hla_report` directory).

**See *disTIL* module details below for more information about each module.**


## Provided configurations of the *disTIL* workflow
### `disTIL_full`

This workflow runs all *disTIL* modules with all current capabilities, requiring the following patient data:

- Two or three sets of paired-end FASTQs for any combination of tumour/normal WGS/WES/RNA-seq
- Tumour gene-level expression TSV
- Tumour transcript-level expression TSV
- Tumour WGS BAM
- Tumour RNA-seq BAM
- Somatic VCF for the tumour sample (and optionally matched normal sample)
- Gene fusion TSV (generated using [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion))

Additional required inputs include:

- Reference genome used to generate the tumour DNA BAM
- Reference genome used to generate the tumour RNA BAM
- An HLA region Bowtie2 index to be used in HLA typing - **a Bowtie2 index generated from the IMGT/HLA's [`hla_gen.fasta`](https://github.com/ANHIG/IMGTHLA/blob/Latest/hla_gen.fasta) reference is provided**
- VEP plugin files to use in somatic VCF annotation - **the plugins required for downstream neoepitope prediction (`Wildtype.pm` and `Frameshift.pm`) are provided**
- Indexed VEP cache (version is dependent on the reference build being used) - available for download [here](http://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/)

Optionally, a phased proximal variants VCF can also be provided for use in pVACseq neoepitope prediction

### `disTIL_no_geneExp`

This workflow runs all *disTIL* modules except for the gene expression classification module. This version of the toolkit can be used to analyse adult tumour data.
This version requires all the same inputs as the full *disTIL* workflow above.

## *disTIL* module details
### Consensus HLA typing

The disTIL consensus HLA typing module uses HLA-HD to conduct HLA typing from two or three NGS samples for a given patient, producing two or three sets of candidate alleles. **There is no restriction on which samples are input to the module provided they are compatible with HLA-HD, meaning any combination of tumour and normal WGS, WES, and RNA-seq samples can be used.** *disTIL* then uses a consensus algorithm to determine consensus alleles from the two or three candidate allele sets produced by running HLA-HD on the input samples. To be deemed a consensus allele, an allele must be present across at least two candidate sets. Consensus alleles are reported to the highest possible level of accuracy (two or three fields). Note that the alleles within an allele pair for a given gene are unordered. If no consensus allele can be determined, the allele is recorded as 'No consensus'.

We defined a subset of clinically significant HLA genes to be used in neoepitope prediction: HLA-A, -B, -C, -DRA, -DRB1, -DRB3, -DRB4, -DRB5, -DQA1, -DQB1, -DPA1 and -DPB1.

#### Inputs

Two or three sets of paired-end FASTQs for any combination of tumour/normal WGS/WES/RNA-seq

#### Outputs

- *Consensus HLA JSON:* a JSON containing consensus alleles for all typed HLA genes
- *Clinically Significant HLA JSON:* a JSON containing consensus alleles for the subset of clinically significant HLA genes
- *Consensus HLA Text File:* text file containing a comma-separated list of unique consensus alleles (to two-field accuracy) for all typed HLA genes
- *Clinically Significant HLA Text File:* text file containing a comma-separated list of unique consensus alleles (to two-field accuracy) for the subset of clinically significant HLA genes
- *Sample 1 HLA-HD Results JSON:* a JSON file representing the candidate alleles typed by HLA-HD from the first input sample
- *Sample 2 HLA-HD Results JSON:* a JSON file representing the candidate alleles typed by HLA-HD from the second input sample
- *Sample 3 HLA-HD Results JSON:* a JSON file representing the candidate alleles typed by HLA-HD from the third input sample
- *HLA report:* a PDF report summarising the results of consensus HLA typing

### Neoepitope prediction

The *disTIL* neoepitope prediction module uses pVACtools (v2.0.4) to conduct neoepitope enumeration, filtering, and ranking. pVACseq is used to predict SNV- and indel-derived neoepitopes, and pVACfuse is used to predict fusion-derived neoepitopes. AGFusion annotation of the input gene fusion TSV is included in the module. Likewise, the following somatic VCF pre-processing steps are included:

- VCF filtering for 'PASS' variants
- Annotation with VEP
- Variant decomposition
- Coverage annotation using tumour DNA BAM and tumour RNA BAM
- Gene and transcript level expression annotation - a tool is included to first convert gene symbols and RefSeq transcript IDs to Ensembl IDs (using [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html))

After pVACtools prediction, this module includes a CWL tool to annotate (and optionally filter) pVACfuse fusion-derived neoepitopes with junction read and spanning fragment counts from the input gene fusion file. **Currently, this feature requires that the fusion TSV be generated using STAR-Fusion.**

#### Inputs

- Somatic variant VCF
- Gene fusion TSV
- Text file containing a comma-separated list of HLA alleles to use in neoepitope prediction
- VEP cache
- VEP plugin files
- Reference genomes used in alignment of tumour DNA and RNA BAM files

#### Outputs

- *pVACseq MHC-I Filtered Report:* top tier HLA-I binding SNV- and indel-derived neoepitopes
- *pVACseq MHC-II Filtered Report:* top tier HLA-II binding SNV- and indel-derived neoepitopes
- *pVACseq MHC-I Aggregate Report:* second tier HLA-I binding SNV- and indel-derived neoepitopes
- *pVACseq MHC-II Aggregate Report:* second tier HLA-II binding SNV- and indel-derived neoepitopes
- *pVACfuse MHC-I STAR-Annotated Filtered Report:* top tier HLA-I binding fusion-derived neoepitopes (including junction read and spanning fragment counts)
- *pVACfuse MHC-II STAR-Annotated Filtered Report:* top tier HLA-II binding fusion-derived neoepitopes (including junction read and spanning fragment counts)
- *pVACfuse MHC-I Aggregate Report:* second tier HLA-I binding fusion-derived neoepitopes
- *pVACfuse MHC-II Aggregate Report:* second tier HLA-II binding fusion-derived neoepitopes

### TMB calculation

To calculate the number of coding missense variants, the *disTIL* TMB module filters a VEP-annotated VCF to retain only ‘PASS’ variants, intersects it with a coding exome BED for the selected genome assembly, and counts the number of variants which VEP predicted to have a missense effect on the overlapping Ensembl canonical transcript. The number of coding missense variants is divided by the coding exome size for the selected genome assembly to yield TMB as the number of missense canonical variants per megabase of the coding exome.

#### Inputs

- VEP annotated VCF

#### Outputs

- *Canonical Missense Coding Variants:* a text file containing the number of coding missense canonical variants
- *TMB Score:* a text file containing the calculated TMB

### Gene expression classification

The *disTIL* gene expression classification module uses the Immune Paediatric Signature Score (IPASS) (unpublished) to classify paediatric tumours as immune hot or cold based on the predicted level of CD8+ T cell infiltration. Tumours with an IPASS score between -0.15 and 0.45 are considered as having ‘intermediate’ phenotype (neither hot nor cold).

#### Inputs

- Gene expression TSV

#### Outputs

- *IPASS Score:* a CSV file containing the patient's IPASS score as well as the expression scores for the 20 genes making up the IPASS immune signature

### Immune cell type deconvolution

The *disTIL* cell type deconvolution module uses EPIC and quanTIseq to predict the cell fractions of multiple immune cell types. It takes as input a gene expression TSV file, and outputs two CSVs containing the results of EPIC and quanTIseq deconvolution. Note that the cell fractions produced by both tools can be interpreted as absolute cell fractions which are comparable between cell types and between patients.

#### Inputs

- Gene expression TSV

#### Outputs

- *EPIC Deconvolution Results:* a CSV file containing the cell fractions predicted by EPIC
- *quanTIseq Deconvolution Results:* a CSV file containing the cell fractions predicted by quanTIseq