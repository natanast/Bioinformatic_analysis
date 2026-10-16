

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

# map to ensebl ids -----

gene_map = df$gene_name |> 
    transId(transTo = c("symbol", "entrez", "ensembl", "uniprot"), unique = TRUE) |> 
    setDT()


gene_map = gene_map[!is.na(ensembl)]


df = df |> merge(gene_map, by.x = "gene_name", by.y = "input_id")


df = df[order(gene_name, padj)]
df = df[, by = gene_name, head(.SD, 1)]
df = df[, by = ensembl,   head(.SD, 1)]



fwrite(
    df[, .(gene_name, ensembl, log2FoldChange, pvalue, padj)],
    "selected_genes.txt",
    sep = "\t", quote = FALSE
)


# retrieve sequences --------

mart = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# seqs = getSequence(
#     id      = df$ensembl,
#     type    = "ensembl_gene_id",
#     seqType = "gene_exon_intron",
#     mart    = mart
# ) |> setDT()

get_seq_safe = function(id, mart, retries = 3) {
    
    for(i in seq_len(retries)) {
        
        out = tryCatch(
            getSequence(
                id      = id,
                type    = "ensembl_gene_id",
                seqType = "gene_exon_intron",
                mart    = mart
            ) |> setDT(),
            error = function(e) NULL
        )
        
        if(!is.null(out) && nrow(out)) return(out)
        
        Sys.sleep(5)   # περίμενε και ξαναδοκίμασε
    }
    
    message("  failed: ", id)
    
    return(NULL)
}


seqs = list()

for(i in seq_len(nrow(df))) {
    
    id = df$ensembl[i]
    
    message(i, "/", nrow(df), "  ", df$gene_name[i], "  (", id, ")")
    
    seqs[[id]] = get_seq_safe(id, mart)
    
    Sys.sleep(1)
}


seqs = seqs |> rbindlist(use.names = TRUE, fill = TRUE)

seqs = seqs[gene_exon_intron != "Sequence unavailable"]

nrow(seqs)


seqs = seqs |> merge(
    df[, .(ensembl, gene_name)], 
    by.x = "ensembl_gene_id", by.y = "ensembl"
)


seqs[, id := paste0(gene_name, "|", ensembl_gene_id)]


# split into overlapping fragments --------

fragment = function(seq, id, width = 100, step = 50) {
    
    n = nchar(seq)
    
    if(n <= width) return(data.table(frag_id = paste0(id, "_1_", n), seq = seq))
    
    starts = seq(1, n - width + 1, by = step)
    
    if(tail(starts, 1) + width - 1 < n) starts = c(starts, n - width + 1)
    
    data.table(
        frag_id = paste0(id, "_", starts, "_", starts + width - 1),
        seq     = substring(seq, starts, starts + width - 1)
    )
}


frags = seqs[, by = id, fragment(gene_exon_intron, id)]

frags = frags[!seq %like% "N"]


# write fasta --------

out = frags$seq |> str_replace_all("T", "U") |> RNAStringSet()

names(out) = frags$frag_id

out |> writeXStringSet("fragments_100nt.fa")

 




