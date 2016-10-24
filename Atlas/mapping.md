#### Salmon

```bash
salmon index --type quasi --threads 4 --index pigeon -t Rock Dove 1.0.3.fasta

for i in $(ls -d1 y*atlas/); do
    F=$(basename $i);
    mkdir salmon_denovo_$F
    salmon quant -i pigeon -l a  \
    -1 <(zcat $i/fastq/*R1*fastq.gz) \
    -2 <(zcat $i/fastq/*R2*fastq.gz) \
    -p 36 --gcBias --seqBias --output salmon_denovo_$F/$F
done
```
