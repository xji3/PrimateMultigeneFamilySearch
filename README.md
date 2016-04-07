# PrimateMultigeneFamilySearch
Search paralog ortholog database to construct multigene families in primates

Two lists of 2-paralog multigene families:

1. A less careful list that combining [Duplicated Genes Database](http://dgd.genouest.org/) and [OrthoMAM Database](http://www.orthomam.univ-montp2.fr/orthomam/html/). This list has 113 pairs of genes.

2. List from above being selected with [metaPhors database](http://betaorthology.phylomedb.org/). This list has 8 pairs of genes left.

DNA CDS sequences were obtained from OrthoMAM v9 that contain 5 primate species: 

1.  Homo sapiens
2.  Pan troglodytes
3.  Gorilla gorilla
4.  Pongo abelii
5.  Macaca mulatta


Now, add Callithrix (NCBI taxonomy id 9583) as outgroup. The two lists after this step have:

1. Less accurate list has 104 pairs of genes left.
2. More careful list has 7 pairs of genes left.

Note:

In STEAP4_STEAP2.fasta, MacacaSTEAP4 was added G in the end of sequence to complete the stop codon manually.

CNTN6 CNTN4 pair was removed.

In IGSF3_CD101.fasta, PongoIGSF3 was added A in the end of sequence to complete the stop codon manually.