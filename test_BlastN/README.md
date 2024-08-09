### Test samples used in:
Wagner et al., Microbiol Spectr. 2022 Apr 27;10(2):e0256421
https://pubmed.ncbi.nlm.nih.gov/35234489/

### Example of search method from Wagner et al. (2022):
`awk -F'\t' '{if($17~/virus/){print($1,$2,$3,$15,$18)}}' SRR1106548-Blast-NT.csv | grep "ACCESSIONX"`

<<<<<<< HEAD

### Alternate search by species name
`awk -F'\t' '{if($17~/virus/){print($1,$2,$3,$15,$18)}}' SRR1106548-Blast-NT.csv | grep "Zed Zombie Virus"`
=======
>>>>>>> e2dfb0828946f94c4af33c12eaf7f17ae067442c
