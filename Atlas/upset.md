#### Extracting count data

```bash
hypomed <-  apply(cpm(atlasobject)[, c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,44,47,50,52,59,62)], 1, median)
pitmed <-  apply(cpm(atlasobject)[, c(2,5,9,11,14,17,20,23,27,30,33,36,38,42,46,49,51,53,55,57,60,63,65)], 1, median)
testesmed <-  apply(cpm(atlasobject)[, c(3,6,12,15,18,21,24,39,54,58,61,64,66)], 1, median)
ovarymed <-  apply(cpm(atlasobject)[, c(8,26,29,32,35,41,43,45,48,56)], 1, median)
write.table(ovarymed, file="ovarym.txt", sep = "," , row.names = T)
write.table(testesmed, file="testesm.txt", sep = "," , row.names = T)
write.table(pitmed, file="pitm.txt", sep = "," , row.names = T)
write.table(hypomed, file="hypom.txt", sep = "," , row.names = T)

hypomeans <- rowMeans(cpm(atlasobject)[, c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,44,47,50,52,59,62)])
pitmeans <- rowMeans(cpm(atlasobject)[, c(2,5,9,11,14,17,20,23,27,30,33,36,38,42,46,49,51,53,55,57,60,63,65)])
testesmeans <- rowMeans(cpm(atlasobject)[, c(3,6,12,15,18,21,24,39,54,58,61,64,66)])
ovarymeans <- rowMeans(cpm(atlasobject)[, c(8,26,29,32,35,41,43,45,48,56)])

write.table(ovarymeans, file="ovary.txt", sep = "," , row.names = T)
write.table(testesmeans, file="testes.txt", sep = "," , row.names = T)
write.table(pitmeans, file="pit.txt", sep = "," , row.names = T)
write.table(hyomeans, file="hypo.txt", sep = "," , row.names = T)

paste <(sort -t, -k1,1 ovarym.txt) <(sort -t, -k1,1 testesm.txt) <(sort -t, -k1,1 pitm.txt) <(sort -t, -k1,1 hypom.txt) | tr , " "  | awk '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8}' | sed 's_"_ _g' > for_inset_med.txt
awk '{if ($2 > 10) print "1"; else print "0"}' for_inset_med.txt > ovm
awk '{if ($3 > 10) print "1"; else print "0"}' for_inset_med.txt > tem
awk '{if ($4 > 10) print "1"; else print "0"}' for_inset_med.txt > pim
awk '{if ($5 > 10) print "1"; else print "0"}' for_inset_med.txt > hym
paste <(cut -f1 for_inset_med.txt) ovm tem pim hym > final.inset_med.txt
```

#### Making graph in R

```R
upset <- read.delim(“final.inset_med.txt”, header=T)
upset(upset, order.by = “freq”)
```
