#### LASTAL+ Search

```bash
lastdb chicken rna.fa
lastal+ -K 1 -P 50 -o chicken.dove.blast ../../blast/chicken ../../../assembly/RockDove.HPG.v1.0.3.fasta
```

#### Mapping IDs

```bash
curl -LO ftp://ftp.ncbi.nih.gov/genomes/Gallus_gallus/GFF/ref_Gallus_gallus-5.0_top_level.gff3.gz
gzip -d *gz

for i in $(awk -F "|" '{print $4}' chicken.dove.blast); do echo $(grep $i ref_Gallus_gallus-5.0_top_level.gff3) | tee -a list13; done

awk -F ";" '{print $6}' list13 | awk -F "=" '{print $2}' > chicken.genes
awk -F ";" '{print $3}' list13 | awk -F ":" '{print $2}' | awk -F "," '{print $1}' > chicken.entrez
awk -F ";" '{print $4}' list13 | awk -F "=" '{print $2}' > chicken.genbank
cut -f1 chicken.dove.blast > chicken.name
paste -d, chicken.name chicken.entrez chicken.genbank chicken.genes > rockdove.chicken.mapping
```
