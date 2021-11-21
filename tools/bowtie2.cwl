class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: bowtie2
baseCommand: []
inputs:
  - id: bowtie2_index
    type: File
    label: Bowtie2 Index Archive
    doc: A TAR archive containing Bowtie2 index files.
    'sbg:fileTypes': TAR
  - id: read1_sequences
    type: File
    inputBinding:
      position: 2
      prefix: '-1'
      shellQuote: false
    label: Read 1 Sequences
    doc: Read 1 sequences in FASTA or FASTQ format (may be bgzipped).
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
  - id: read2_sequences
    type: File
    inputBinding:
      position: 3
      prefix: '-2'
      shellQuote: false
    label: Read 2 Sequences
    doc: Read 2 sequences in FASTA or FASTQ format (may be bgzipped).
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ, FQ, FQ.GZ'
  - id: no_unaligned
    type: boolean?
    inputBinding:
      position: 1
      shellQuote: false
      valueFrom: |-
        ${
            if (self == true) {
                return "--no-unal"
            } else {
                return ""
            }
        }
    label: Suppress Unaligned Reads
    doc: Suppres SAM records for unaligned reads.
  - id: sample_name
    type: string
    label: Sample Name
    doc: Sample name used to name the output SAM file.
  - id: bowtie2_index_prefix
    type: string
    inputBinding:
      position: 1
      prefix: '-x'
      shellQuote: false
    label: Bowtie2 Index Prefix
    doc: >-
      The prefix of index files contained in the Bowtie2 index TAR. Note that
      all Bowtie2 nidex files in the TAR should have this prefix.
outputs:
  - id: aligned_sam
    type: File
    outputBinding:
      glob: '*.sam'
      outputEval: '$(inheritMetadata(self, inputs.read1_sequences))'
doc: >-
  # About this tool

  This tool runs Bowtie2 (**v2.4.1**) alignment of reads to a reference
  sequence.


  ## Inputs and parameters

  - **Bowtie2 Index Archive:** an archive (TAR) of index files generated from a
  reference genome (using Bowtie2 build command).

  - **Bowtie2 Index Prefix:** the basename of the index files in the Bowtie2
  Index Archive. The basename should be shared by all index files in the TAR,
  and is the name of the index files up to but not including the final `.1.bt2`
  etc.

  - **Read 1 Sequences:** file containing mate 1 reads (for paired-end
  sequencing).

  - **Read 2 Sequences:** file containing mate 2 reads (for paired-end
  sequencing).

  - **Sample Name:** the name of the sample being analysed (used as the name for
  the output SAM file).


  ## Docker

  This tool uses the Docker image: `biocontainers/bowtie2:v2.4.1_cv1`.


  ## Documentation

  - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
label: bowtie2
arguments:
  - position: 4
    prefix: '-S'
    shellQuote: false
    valueFrom: $(inputs.sample_name).sam
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          var ref = inputs.bowtie2_index.path.split('/').splice(-1)
          return "tar -xvf " + ref + " && rm " + ref + " &&"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: bowtie2
  - position: 2
    prefix: '-p'
    shellQuote: false
    valueFrom: '8'
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 6000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'biocontainers/bowtie2:v2.4.1_cv1'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.bowtie2_index)
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
'sbg:toolAuthor': 'Langmead B, Salzberg S'
'sbg:toolkit': bowtie2
'sbg:toolkitVersion': 2.4.1
