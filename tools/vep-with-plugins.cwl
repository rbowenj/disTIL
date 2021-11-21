class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: vep_with_plugins
baseCommand: []
inputs:
  - id: input_file
    type: File
    inputBinding:
      position: 3
      prefix: '--input_file'
      shellQuote: false
    label: Input File
    doc: Input VCF to VEP annotate. Can be bgzipped.
    'sbg:fileTypes': 'VCF, VCF.GZ'
  - id: vep_cache
    type: File
    label: VEP Cache
    doc: >-
      VEP cache supplied as a tar.gz. Caches can be downloaded from [VEP release
      page](http://ftp.ensembl.org/pub/). It is recommended that the cache
      version used matches the VEP release number (note that this app uses the
      latest VEP release. See
      [here](https://hub.docker.com/r/ensemblorg/ensembl-vep/tags?page=1&ordering=last_updated)
      for more info).
    'sbg:fileTypes': TAR.GZ
  - id: ref_genome
    type: File
    inputBinding:
      position: 3
      prefix: '--fasta'
      shellQuote: false
    label: Reference Genome
    doc: The reference genome FASTA file to use to look up reference sequence.
    'sbg:fileTypes': 'FASTA, FA'
  - id: vep_plugin_files
    type: 'File[]?'
    label: VEP Plugin Files
    doc: Optional VEP plugin files to use when annotating the input VCF.
    'sbg:fileTypes': PM
  - id: cache_version
    type: int
    inputBinding:
      position: 3
      prefix: '--cache_version'
      shellQuote: false
    label: Cache Version
    doc: >-
      Version number of the cache used. It is recommended that the cache version
      used matches the VEP release number (note that this app uses the latest
      VEP release. See
      [here](https://hub.docker.com/r/ensemblorg/ensembl-vep/tags?page=1&ordering=last_updated)
      for more info).
  - id: merged
    type: boolean
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.merged == true) {
                return '--merged'
            } else {
                return ''
            }
            
        }
    label: Merged
    doc: >-
      Use the merged Ensembl and RefSeq cache. Consequences are flagged with the
      SOURCE of each transcript used.

      NOTE: This flag MUST be used if the cache is merged.
  - 'sbg:toolDefaultValue': 'No'
    id: symbol
    type: boolean?
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.symbol == true) {
                return '--symbol'
            } else {
                return ''
            }
        }
    label: Symbol
    doc: Adds the gene symbol (e.g. HGNC) (where available) to the output.
  - 'sbg:toolDefaultValue': 'No'
    id: biotype
    type: boolean?
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.biotype == true) {
                return '--biotype'
            } else {
                return ''
            }
            
        }
    label: Biotype
    doc: Adds the biotype of the transcript or regulatory feature.
  - 'sbg:toolDefaultValue': 'No'
    id: numbers
    type: boolean?
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.numbers == true) {
                return '--numbers'
            } else {
                return ''
            }
        }
    label: Numbers
    doc: >-
      Adds affected exon and intron numbering to to output. Format is
      Number/Total.
  - 'sbg:toolDefaultValue': 'No'
    id: canonical
    type: boolean?
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.canonical == true) {
                return '--canonical'
            } else {
                return ''
            }
        }
    label: Canonical
    doc: >-
      Adds a flag indicating if the transcript is the canonical transcript for
      the gene.
  - 'sbg:toolDefaultValue': 'No'
    id: total_length
    type: boolean?
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.total_length == true) {
                return '--total_length'
            } else {
                return ''
            }
        }
    label: Total Length
    doc: 'Give cDNA, CDS and protein positions as Position/Length.'
  - 'sbg:toolDefaultValue': 'No'
    id: sift
    type:
      - 'null'
      - type: enum
        symbols:
          - p
          - s
          - b
        name: sift
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.sift == null) {
                return ''
            } else {
                return '--sift ' + inputs.sift
            }
        }
    label: Sift
    doc: >-
      Species limited SIFT predicts whether an amino acid substitution affects
      protein function based on sequence homology and the physical properties of
      amino acids. VEP can output the prediction term, score or both.
  - 'sbg:toolDefaultValue': 'No'
    id: polyphen
    type:
      - 'null'
      - type: enum
        symbols:
          - p
          - s
          - b
        name: polyphen
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.polyphen == null) {
                return ''
            } else {
                return '--polyphen ' + inputs.polyphen
            }
        }
    label: Polyphen
    doc: >-
      Human only PolyPhen is a tool which predicts possible impact of an amino
      acid substitution on the structure and function of a human protein using
      straightforward physical and comparative considerations. VEP can output
      the prediction term, score or both.
  - 'sbg:toolDefaultValue': SO
    id: terms
    type:
      - 'null'
      - type: enum
        symbols:
          - SO
          - display
          - NCBI
        name: terms
    inputBinding:
      position: 6
      separate: false
      shellQuote: false
      valueFrom: |-
        ${
            if (inputs.terms == null) {
                return ''
            } else {
                return '--terms ' + inputs.terms
            }
        }
    label: Terms
    doc: >-
      The type of consequence terms to output. The Ensembl terms are described
      [here](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences).
      The Sequence Ontology is a joint effort by genome annotation centres to
      standardise descriptions of biological sequences.
outputs:
  - id: vep_vcf
    doc: VEP annotated VCF file.
    label: Annotated VCF
    type: File
    outputBinding:
      glob: '*.vep.vcf'
    'sbg:fileTypes': VCF
  - id: vep_stats
    doc: Stats file produced by VEP.
    label: VEP Stats
    type: File?
    outputBinding:
      glob: '*.vep.html'
    'sbg:fileTypes': HTML
doc: >-
  #About VEP

  VEP determines the effect of your variants (SNPs, insertions, deletions, CNVs
  or structural variants) on genes, transcripts, and protein sequence, as well
  as regulatory regions.


  Simply input the coordinates of your variants and the nucleotide changes to
  find out the:

  Genes and Transcripts affected by the variants

  Location of the variants (e.g. upstream of a transcript, in coding sequence,
  in non-coding RNA, in regulatory regions)

  Consequence of your variants on the protein sequence (e.g. stop gained,
  missense, stop lost, frameshift), see variant consequences

  Known variants that match yours, and associated minor allele frequencies from
  the 1000 Genomes Project

  SIFT and PolyPhen-2 scores for changes to protein sequence

  ... And more! See data types, versions.


  # About this CWL tool

  This tool is intended for VEP annotation of VCF files prior to neoantigen
  prediction using pVACseq
  ([info](https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html)).
  All VEP options required for pVACseq are exposed as app settings, plus some
  additional options. 


  ## Cache

  The VEP cache must be supplied as a `tar.gz`. Caches can be downloaded from
  [VEP release page](http://ftp.ensembl.org/pub/). It is recommended that the
  cache version used matches the VEP release number (note that this app uses the
  latest VEP release. See
  [here](https://hub.docker.com/r/ensemblorg/ensembl-vep/tags?page=1&ordering=last_updated)
  to find out the latest VEP release number). As of July 2021, the latest VEP
  release was version 104.


  ### Merged

  If the cache version used is merged, the `--merged` flag must be used. This is
  achieved by setting the 'merged' input option to 'Yes'. 


  ## Plugins

  pVACseq requires the use of the Frameshift and Wildtype plugins, available for
  download
  [here](https://github.com/griffithlab/pVACtools/tree/master/tools/pvacseq/VEP_plugins).
label: vep-with-plugins
arguments:
  - position: 4
    prefix: '--dir_cache'
    shellQuote: false
    valueFrom: |-
      ${
          return './cache'
      }
  - position: 4
    prefix: '--output_file'
    shellQuote: false
    valueFrom: |-
      ${
          var in_file = inputs.input_file.basename
          
          if (in_file.endsWith(".gz")) {
              var out_file = in_file.replace('.vcf.gz', '.vep.vcf')
          } else {
              var out_file = in_file.replace('.vcf', '.vep.vcf')
          }
          
          return out_file
      }
  - position: 5
    prefix: '--stats_file'
    shellQuote: false
    valueFrom: |-
      ${
          var in_file = inputs.input_file.basename
          
          if (in_file.endsWith(".gz")) {
              var out_file = in_file.replace('.vcf.gz', '.vep.html')
          } else {
              var out_file = in_file.replace('.vcf', '.vep.html')
          }
          
          return out_file
      }
  - position: 5
    prefix: '--species'
    shellQuote: false
    valueFrom: homo_sapiens
  - position: 3
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: >-
      --cache --format vcf --vcf --offline --fork 8 --no_progress --tsl --hgvs
      --shift_hgvs 1
  - position: 0
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |-
      mkdir ./plugins ${
        if (inputs.vep_plugin_files == null) {
          return ""
        }

        let mv_cmd = "";
        for (var i=0; i < inputs.vep_plugin_files.length; i++) {
          mv_cmd = mv_cmd + " && mv " + inputs.vep_plugin_files[i].path + " ./plugins/";
        }
        return mv_cmd;
      } && ls ./plugins 1>&2 &&
  - position: 2
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: vep
  - position: 7
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |2-
       ${
          let plugin_cmd = "";
          for (var i=0; i < inputs.vep_plugin_files.length; i++) {
              plugin_split = inputs.vep_plugin_files[i].path.split('/')
              plugin = plugin_split[plugin_split.length-1].split('.')[0]
              plugin_cmd = plugin_cmd + "--plugin " + plugin + " ";
          }
          return plugin_cmd;
      }
  - position: 1
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |-
      ${
        var cache_bundle = inputs.vep_cache.basename
        return 'mkdir ./cache' + ' && tar -xvf ' + cache_bundle + ' -C ./cache &&'
      }
  - position: 10
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |-
      ${
          if (inputs.vep_plugin_files.length != 0) {
              return '--dir_plugins ./plugins'
          } else {
              return ''
          }
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 52000
    coresMin: 2
  - class: DockerRequirement
    dockerPull: 'ensemblorg/ensembl-vep:release_104.3'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.vep_plugin_files)
      - $(inputs.vep_cache)
  - class: InlineJavascriptRequirement
'sbg:categories':
  - Annotation
  - VCF Processing
'sbg:license': Apache 2.0
'sbg:toolAuthor': Ensembl
'sbg:toolkit': vep
'sbg:toolkitVersion': '104.3'
'sbg:wrapperAuthor': ''
