#### TransRate

```bash
transrate -o rockdove103 -t 30 -a RockDove.HPG.v1.0.3.fasta --left hpg.axis.1.fq.gz --right hpg.axis.2.fq.gz
```

#### BUSCO

```bash
python BUSCO.py -i RockDove.HPG.v1.0.3.fasta -o busco_103 -m tran \
-l vertebrata_odb9 -c 40
```
