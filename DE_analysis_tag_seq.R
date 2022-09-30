# Preliminary analysis of 3' Tag Seq data
# 12 Sept 2022

# Step 2: Attempt differential expression analysis

library(ape) # for read.gff
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(readxl)
library(vegan)

samples = read_excel("sample_to_condition_table.xlsx")
tpm = read_csv("kallisto_tpm_output.csv")
est_counts = read_csv("kallisto_est_counts.csv")
# label_linker = read_tsv("../Pseudotsuga_annotations/Psme.1_0.gmap.gff3",header=F, comment.char="#")
label_linker = read.gff("../Pseudotsuga_annotations/Psme.1_0.gmap.gff3")
psme_gtf = read.gff("../Pseudotsuga_annotations/Psme.1_0.gtf", GFF3 = FALSE)

justHeb = subset(samples, samples$fungus != "Suillus")

tpm_matrix = as.matrix(tpm[,2:ncol(tpm)])
rownames(tpm_matrix) = tpm$target_id

justheb_tpm = tpm_matrix[,1:16]

count_matrix = as.matrix(est_counts[,2:ncol(est_counts)])
rownames(count_matrix) = est_counts$target_id

justheb_counts = count_matrix[,1:16]

justHeb$group = paste(justHeb$colonized, justHeb$harvest_water_potential, sep = "_")

# I am following along with this DESeq2 beginner's guide:
# https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf

# The Internet assures me DESeq2 can only use count data
# not using TPMs then, at least for now.
# ddsFullCountTable = DESeqDataSetFromMatrix(
#   countData = round(justheb_counts), # error if non-integer
#   colData = justHeb,
#   design = ~ colonized * harvest_water_potential
# )

ddsFullCountTable = DESeqDataSetFromMatrix(
  countData = round(justheb_counts), # error if non-integer
  colData = justHeb,
  design = ~ group
)

dds = DESeq(ddsFullCountTable)

resultsNames(dds)
# res = results(dds)
res_wet = results(dds, contrast = c("group", "uncolonized_wet", "colonized_wet"))
# This will show fold changes and p values for the 
# "last variable in the design formula," so harvest water potential here.

# How many genes were significantly differentially
# regulated in my treatments?
sum(res$padj[is.na(res$padj) == F] < 0.05) #74
sum(res_wet$padj[is.na(res$padj) == F] < 0.05) #51

# Looking at adjusted p values less that 0.05,
# the answer is 74.

resSig = res[which(res$padj < 0.05),]
head(resSig[order(resSig$log2FoldChange),])

reswetSig = res_wet[which(res_wet$padj < 0.05),]
head(resSig[order(resSig$log2FoldChange),])

# regularized-logarithm transformation

rld = rlog(dds) # 18 warnings, all about graphics state,
# first one said "display list redraw incomplete"
rld_wet = rlog(dds)
# maybe will still work?

sampleDists = dist(t(assay(rld))) # I believe this gives Euclidean distance

plotPCA(rld, intgroup = c("colonized", "harvest_water_potential"))
# hmmmm
# colonized and uncolonized do not look as different as I had hoped.
# they seem to overlap.
# colonized plants do seem to exhibit greater transcriptional
# variability, though, FWIW.

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 40 )
top40_headmap = pheatmap( assay(rld)[ topVarGenes, ], 
          scale = "row",
          cluster_cols = TRUE)

save_plot("plots/top40_heatmap.pdf", top40_headmap)

# Col/uncol contrast in wet only

topVarGenes <- head( order( rowVars( assay(rld_wet) ), decreasing=TRUE ), 40 )
top40_heatmap_wet = pheatmap( assay(rld_wet)[ topVarGenes, ], 
                          scale = "row",
                          cluster_cols = TRUE)

save_plot("plots/top40_heatmap.pdf", top40_headmap)

# What if I am more interested in colonization?

ddsFullCountTable_col = DESeqDataSetFromMatrix(
  countData = round(justheb_counts), # error if non-integer
  colData = justHeb,
  design = ~ harvest_water_potential * colonized
)

dds_col = DESeq(ddsFullCountTable_col)
res_col = results(dds_col)

count(res_col$padj[is.na(res_col$padj) == F] < 0.05)
# Looking at adjusted p values less that 0.05,
# the answer is 74 still.

resSig_col = res_col[which(res_col$padj < 0.05),]
head(resSig_col[order(resSig_col$log2FoldChange),])

# regularized-logarithm transformation

rld_col = rlog(dds_col) # No warnings this time

sampleDists = dist(t(assay(rld_col)))

plotPCA(rld_col, intgroup = c("colonized", "harvest_water_potential"))
# same exact PCA, but that figures.

topVarGenes <- head( order( rowVars( assay(rld_col) ), decreasing=TRUE ), 35 )
pheatmap( assay(rld_col)[ topVarGenes, ], cluster_cols = TRUE)
# aaand heatmap is ALSO unaffected by the order of the variables in the 
# model. Good to know!

# This is really not a glaringly obvious acse of differential expression.
# I don't see anything that is totally consistent in colonized plants
# compared to uncolonized ones.

# Regardless, at least it is something. Let's see what these genes might be.

### Applying more meaningful names to the genes ####

# GenBankIDs = sub('ID=(.+)', '\\1', label_linker$attributes)
# GenBankIDs_alone = unlist(strsplit(GenBankIDs, "(?<=\\.[0-9])", perl = TRUE))[1] # lookback so first element is GenBank ID
# label_linker$target_id = GenBankIDs_alone
# 
# targ_ids_tidy = sub('Name=(.+)', '\\1', label_linker$attributes)

# trying substring extraction
# testset = label_linker[1:10,]
# startandend = str_locate(testset$attributes, "Target=(.+)")
# str_sub(testset$attributes[1], startandend[1,1], end = startandend[1,2])

label_linker$target_id = as.vector(sub('Target=', '\\1', (str_extract_all(label_linker$attributes, "Target=(.{14})", "\\1"))))

newtest = label_linker[grep("jcf7190000020720", label_linker$seqid),]
# Looks like every seqid has many unique target_id entries.
# pretty sure only the seqids are linked to annotations,
# though. collapse this.

justlabels = select(label_linker, seqid, target_id)
justlabels = distinct(justlabels)

est_counts_linked = left_join(est_counts, justlabels)

# This one has no seqid: GAZW02000005.1
grep("GAZW02000005.1", label_linker$attributes) # nada
# How many transcripts do we have with no seqID?

count(is.na(est_counts_linked$seqid)) #131,384 with no seqid
nrow(est_counts_linked) # out of 331,745 total entries.
# This is nearly half my data! Dang it!
# I guess I will just have to proceed for now.
#### ^PLEASE REMEMBER TO ASK COLLABORATORS WHY THIS MIGHT BE^ ####

psme_gtf = as_tibble(psme_gtf)
psme_gtf = rename(psme_gtf, "attributes" = "gene_name", "seqname" = "seqid")
psme_gtf_merge = select(psme_gtf, seqid, gene_name) %>% distinct()

est_counts_linked = as_tibble(est_counts_linked)
linked_formerge = select(est_counts_linked, target_id, seqid)

final_lookup = left_join(linked_formerge, psme_gtf_merge)

write_csv(final_lookup, "intermediate_output/lookup_table_gene_target_seqid.csv")

# It appears that MANY seqid entries correspond to
# multiple PSME_#### genes. Why?
# This is feeling like a very poor use of time
# when I know the genome center could do a better job.

# Well, I will just do my best for now, with the understanding
# that I do need to get external help with this.

### Statistics in Vegan ####

mygroups = select(samples, sample_ID, colonized, harvest_water_potential)
mygroups = arrange(mygroups, sample_ID)
mygroups = mygroups[-grep("RS", mygroups$sample_ID),]

groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))


# Euclidean distances from earlier:
sampleDists


###
