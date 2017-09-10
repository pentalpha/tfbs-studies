# tfbs-studies
Bioinformatics study on transcription factors binding sites (TFBS) and their association with specific tissues.

This repository needs another directory, at the same parent directory, called "input" (or "../input", relatively"). It is available at http://177.20.147.141/~pitagoras/TF-findings/input/.

The results are stored in the "../results/" directory. The scripts must be executed in order(1-..., 2-..., 3-... etc).

## Task 1. Finding TFBSs for each gene

First we need to know that **transcription** is the process where a gene's DNA sequence is copied (transcribed) into 
an RNA molecule. So when enzyme RNA polymerase, which makes a new RNA molecule from a DNA template, his must attach 
to the DNA of the gene. It attaches at a spot called the **promoter** and his can attach to the promoter region only
with the help of proteins called **transcription factors** inside their **binding sites**. So transcription factors
are proteins that regulate the transcription of genes (activating or repressing), that is their copying into RNA, on
the way to making a protein. So the binding sites are regions that are linked to transcriptional regulation of a gene. 

For this task, our input data is:
- Transcription factors file: TFBS-ENCODE.bed
- Genes reference file: coord_upstream5kb_refGene.bed

This command has shown to be efficient at intersecting the two inputs ([script](0-createIntersectionBed.sh)):
```sh
bedtools intersect -wa -wb -a TFBS-ENCODE.bed -b coord_upstream5kb_refGene.bed > bedIntersectWaWbTFBSinGenes.bed
```

Ouput size: +- 196MB

1.969.951 lines, almost 2 million transcription factor binding sites related to some gene on the reference input.

Available at: http://177.20.147.141/~pitagoras/TF-findings/results/bedIntersectWaWbTFBSinGenes.bed

## Task 2. Finding Transcription Factors related with tissues


## References
Khan Academy. Gene regulation. <https://www.khanacademy.org/science/biology/gene-regulation>.
