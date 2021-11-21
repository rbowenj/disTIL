FROM rachelbj/samtools:1.10.0

RUN apt-get update -y && apt-get install -y bedtools

COPY ./ensembl_genes.coding.exons.collapsed.filtered.grch37.bed /app/ensembl_genes.coding.exons.collapsed.filtered.grch37.bed
COPY ./gencode_basic_v38.coding.exons.collapsed.filtered.grch38.bed /app/gencode_basic_v38.coding.exons.collapsed.filtered.grch38.bed

CMD /bin/bash