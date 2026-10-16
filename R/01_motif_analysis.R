

rm(list = ls())
gc()


# libraries ---------

library(data.table)
library(stringr)
library(DESeq2)
library(genekitr)
library(biomaRt)
library(Biostrings)

# load data ---------

sample_metadata = "sample_metadata.xlsx" |> 
    readxl::read_xlsx(sheet = 1) |> 
    setDT()


df = "DESeq2/DESeq2_Condition_SCA1_vs_Control.xlsx" |> 
    readxl::read_xlsx() |>
    setDT()



# clean data --------

df = df[!is.na(padj)]

df = df[padj < 0.05 & log2FoldChange > 0]

df = df[order(padj, -log2FoldChange)]


df = head(df, 50)
