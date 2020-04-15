#### Input ####

# windows of differentiation
gen_fn = "differentiation_output.csv"

# genome annotations
txgdict = "../metadata/Anogam_tx2ge.csv"
gomapfi = "../metadata/Anogam_genes_eggnog_diamond.emapper.annotations.GO"
panfile = "../metadata/Anogam_genes_Pfamscan.seqs"
gff_fn  = "../metadata/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3"

# where to store output?
outcode = "differentiation_output_post" # (folder + prefix)

# load libraries
library(topGO)
library(fdrtool)
library(GenomicRanges)
library(rtracklayer)
source("../scripts_other/geneSetAnalysis.R")

# chromosomes to include
chromlist = c("2R","2L","3R","3L","X")

# genewise differentiation
gen = read.table(gen_fn, sep="\t", header = T)

# transcript to gene dictionary
di           = read.table(txgdict)
colnames(di) = c("transcript_id","gene_id")
di           = di[,1:2]

# functional mappings: GO, pfam
gomap = readMappings(gomapfi)
panno = read.table(file = panfile)
colnames(panno) = c("gene","pstart","pend","pfamid","domain","domseq")
panno = merge(panno,di, by.x="gene", by.y="gene_id")
pannu = subset(panno, select=c("gene","pfamid","domain"))
pannu = pannu[!duplicated(pannu), ]
panno_ag = aggregate(domain ~ gene, data = panno, paste, collapse = ",")


# load GFF and build GenomicRanges with genes
gff_gr      = rtracklayer::import.gff(gff_fn)
gff_gr_gene = gff_gr[gff_gr$type=="gene"]
gff_df_gene = data.frame(gff_gr_gene)
gene_names = gff_df_gene[,c("ID","Name","description")]
names(gff_gr_gene) = gff_gr_gene$ID


# adjust pvalues (FDR)
gen$PBS_p_adj = fdrtool(gen$PBS_p, statistic = 'pvalue', plot = F)$qval
write.table(gen, file = sprintf("differentiation_output_wFDR.csv", outcode), 
            sep="\t", quote = F, row.names = F)


# find regions pval<0.001
gef = gen[gen$PBS_p_adj < 0.001 & gen$PBS > 0, ]
# GenomicRanges of peaks
gef_gr = makeGRangesFromDataFrame(gef,keep.extra.columns = T,start.field = "start",end.field = "end")

# find overlapping genes
ov_gef_gr = findOverlapPairs(query=gef_gr,subject=gff_gr_gene,type="any",select="all",ignore.strand=T)
ov_gef_df = data.frame(
  chrom   = seqnames(ov_gef_gr@first),
  window  = start(ov_gef_gr@first),
  gene    = names(ov_gef_gr@second),
  start_g = start(ov_gef_gr@second),
  end_g   = end(ov_gef_gr@second),
  pbs     = ov_gef_gr@first$PBS,
  FDR     = ov_gef_gr@first$PBS_p_adj,
  fst     = ov_gef_gr@first$Fst
)
ov_list = unique(as.vector(ov_gef_df$gene))

# output table
ov_gef_df = merge(ov_gef_df, gene_names, by.x="gene", by.y="ID")
ov_gef_df = merge(ov_gef_df, panno_ag, by.x="gene", by.y="gene")
write.table(ov_gef_df, file = sprintf("%s.PBStop_genes.csv", outcode), sep="\t", quote = F, row.names = F)


# highest pbs
# functional enrichments
hygeofun(list_interest=ov_list, 
         annotation=pannu, gene_col="gene", ano_col="domain",
         outputname=outcode,
         name_geneset="PBStop_genes",topnum = 20, padjbool = T)
suppressMessages(topgofun(list_interest=ov_list, 
                          gomap=gomap,
                          ontologyset=c("BP","MF","CC"),tg_test="fisher",tg_algorithm="elim",
                          outputname=outcode,
                          name_geneset="PBStop_genes",topnum=20))



