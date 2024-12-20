library("R.utils")

untar("GSE130708_RAW.tar")

dir.create("sc_illumina_190")
dir.create("sc_illumina_951")
dir.create("sc_nanopore_190")
dir.create("sc_nanopore_951")


file.copy(from = "GSM3748086_190c_barcodes.tsv.gz",
          to   = "sc_illumina_190/barcodes.tsv.gz")

file.copy(from = "GSM3748086_190c_genes.tsv.gz",
          to   = "sc_illumina_190/features.tsv.gz")

file.copy(from = "GSM3748086_190c_matrix.mtx.gz",
          to   = "sc_illumina_190/matrix.mtx.gz")

file.copy(from = "GSM3748088_951c_barcodes.tsv.gz",
          to   = "sc_illumina_951/barcodes.tsv.gz")

file.copy(from = "GSM3748088_951c_genes.tsv.gz",
          to   = "sc_illumina_951/features.tsv.gz")

file.copy(from = "GSM3748088_951c_matrix.mtx.gz",
          to   = "sc_illumina_951/matrix.mtx.gz")

file.copy(from = "GSM3748087_190c.isoforms.matrix.txt.gz",
          to   = "sc_nanopore_190/GSM3748087_190c.isoforms.matrix.txt.gz")

file.copy(from = "GSM3748089_951c.isoforms.matrix.txt.gz",
          to   = "sc_nanopore_951/GSM3748089_951c.isoforms.matrix.txt.gz")

files_list <- list.files(pattern = ".gz")
files_list
file.remove(files_list)

rm(list=ls())