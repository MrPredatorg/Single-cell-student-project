##############Making Nanopore 190c matrix for analisys

library("R.utils")
library('Matrix')
library("dplyr")

gunzip("sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt.gz", remove=TRUE)

df <-read.table("sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt",header=TRUE,sep="")

df <- df %>%
  group_by(geneId) %>%
  summarize(across(CAACTAGAGCTGTTCA:TTCTTAGTCTGTTGAG, sum))

df2 <- df[, -1]

rownames(df2) <- df$geneId

# save sparse matrix

gene_matrix <- Matrix(as.matrix(df2) , sparse = T )

head(gene_matrix)

writeMM(obj = gene_matrix, file="sc_nanopore_190/matrix.mtx")

# save genes and cells names

write.table(x = rownames(gene_matrix), file = "sc_nanopore_190/genes.tsv", col.names = FALSE, sep = ",")

write.table(x = colnames(gene_matrix), file = "sc_nanopore_190/barcodes.tsv", col.names = FALSE, sep = ",")

file.remove("sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt")


###############Making Nanopore 190c 951c matrix for analisys

gunzip("sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt.gz", remove=TRUE)

df <-read.table("sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt",header=TRUE,sep="")

df <- df %>%
  group_by(geneId) %>%
  summarize(across(GGCGACTAGGCTCATT:GGCTCGAGTACCGAGA, sum))

df2 <- df[, -1]

rownames(df2) <- df$geneId

# save sparse matrix
gene_matrix <- Matrix(as.matrix(df2) , sparse = T )

head(gene_matrix)

writeMM(obj = gene_matrix, file="sc_nanopore_951/matrix.mtx")

# save genes and cells names

write.table(x = rownames(gene_matrix), file = "sc_nanopore_951/genes.tsv", col.names = FALSE, sep = ",")

write.table(x = colnames(gene_matrix), file = "sc_nanopore_951/barcodes.tsv", col.names = FALSE, sep = ",")

rm(list=ls())

file.remove("sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt")
