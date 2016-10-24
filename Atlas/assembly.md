#### RCorrector

```bash
perl run_rcorrector.pl -k 31 -t 30
-1 3te.1.fq.gz,3h.1.fq.gz,3p.1.fq.gz,12p.1.fq.gz,12h.1.fq.gz,12ov.1.fq.gz \
-2 3te.2.fq.gz,3h.2.fq.gz,3p.2.fq.gz,12p.2.fq.gz,12h.2.fq.gz,12ov.2.fq.gz
```

#### Trinity

```bash
Trinity --seqType fq --max_memory 300G --trimmomatic --SS_lib_type RF \
--left 3te.1.fq.gz,3h.1.fq.gz,3p.1.fq.gz,12p.1.fq.gz,12h.1.fq.gz,12ov.1.fq.gz \
--right 3te.2.fq.gz,3h.2.fq.gz,3p.2.fq.gz,12p.2.fq.gz,12h.2.fq.gz,12ov.2.fq.gz \
--CPU 10 --output trinity_HPG_atlas --inchworm_cpu 10 --full_cleanup \
--quality_trimming_params "ILLUMINACLIP:/share/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

#### Binpacker

```bash
BinPacker -q -d -s fq -p pair -m RF -k 31 -g 200 -o binpacker_atlas \
-l hpg.axis.1.fq \
-r hpg.axis.2.fq
```

#### Transfuse
```bash
transfuse -t 30 -i 0.98 -o transfuse_mouse \
-l hpg.axis.1.fq \
-r hpg.axis.2.fq \
-a BinPacker.fa,trinity_HPG_atlas.Trinity.fasta
```

```bash
kallisto index -i kallisto.idx transfuse_mouse.fasta
kallisto quant -t 32 -i kallisto.idx -o kallisto_orig hpg.axis.1.fq hpg.axis.2.fq
salmon index -t transfuse_mouse.fasta -i salmon.idx --type quasi -k 31
salmon quant --seqBias --gcBias -p 32 -i salmon.idx -l a -1 hpg.axis.1.fq hpg.axis.2.fq -o salmon_orig
awk '1>$5{next}1' kallisto_orig/abundance.tsv | awk '{print $1}' > kallist
awk '1>$4{next}1' salmon_orig/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list
python filter.py transfuse_mouse.fasta uniq_list > RockDove.HPG.v1.0.1.fa
awk '/^>/{$0=">'RockDoveHPG'_"++i}1' RockDove.HPG.v1.0.1.fa > RockDove.HPG.v1.0.2.fasta
```

```bash
dammit annotate ../assembly/RockDove.HPG.v1.0.2.fasta \
--user-databases dammit_databases/tcdb.fasta dammit_databases/Gallus_gallus.Galgal4.ncrna.fa \
--busco-group vertebrata \
--n_threads 15 \
--full \
--database-dir dammit_databases/

```
