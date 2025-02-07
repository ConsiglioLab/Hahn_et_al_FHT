#!/usr/bin/env bash
# from https://divingintogeneticsandgenomics.com/post/how-to-make-a-transcript-to-gene-mapping-file/
zless -S data/Macaca_mulatta_genome/Macaca_mulatta.Mmul_10.110.gtf.gz | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$6"\t"$2}' | sort | uniq |  sed 's/\"//g' > data/Macaca_mulatta_genome/tx2gene_ensembl_Macaca_mulatta_genome.txt
cat data/Macaca_mulatta_genome/tx2gene_ensembl_Macaca_mulatta_genome.txt | head
echo "Done!"