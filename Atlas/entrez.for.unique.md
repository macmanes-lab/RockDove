### Finding entrez IDs for genes expressed uniquely in 1 tissue.

```
for i in $(awk '{if ($2 == 1 && $3 == 0 && $4 == 0 && $5 == 0) print $0}' final.inset_med.txt | cut -f1); do  grep -w $i tximport.passed.genes4 | cut -d, -f2 | tee -a entrez.only.ovary.txt; done
for i in $(awk '{if ($2 == 0 && $3 == 1 && $4 == 0 && $5 == 0) print $0}' final.inset_med.txt | cut -f1); do  grep -w $i tximport.passed.genes4 | cut -d, -f2 | tee -a entrez.only.testes.txt; done
for i in $(awk '{if ($2 == 0 && $3 == 0 && $4 == 1 && $5 == 0) print $0}' final.inset_med.txt | cut -f1); do  grep -w $i tximport.passed.genes4 | cut -d, -f2 | tee -a entrez.only.pituitary.txt; done
for i in $(awk '{if ($2 == 0 && $3 == 0 && $4 == 0 && $5 == 1) print $0}' final.inset_med.txt | cut -f1); do  grep -w $i tximport.passed.genes4 | cut -d, -f2 | tee -a entrez.only.hypothalamus.txt; done
```
