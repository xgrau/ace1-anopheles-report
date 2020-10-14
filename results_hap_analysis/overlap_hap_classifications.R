# libraries
source("../scripts_other/geneSetAnalysis.R")
library(stringr)

# input
c_umap_fn = "umap_classification.csv"
c_umap = read.table(c_umap_fn, sep = "\t", header = T, stringsAsFactors = F)
c_umap$ox_code_sample = stringr::str_replace(c_umap$ox_code, pattern = "[a|b]$", replacement = "")

# get lists of haps
c_umap_280S = c_umap[c_umap$p2_dtc_prediction == 1,"ox_code"]
c_umap_280S = c_umap[c_umap$p2_rfc_prediction == 1,"ox_code"]
c_umap_280S = c_umap[c_umap$p2_umap_0 >6 & c_umap$p2_umap_1 >20 & c_umap$p2_umap_1 < 40 ,"ox_code"]

for (ldi in c("out_Ace1.Hapnet_mst_var_3465693_ALT.haps.list",
              "out_Ace1.Hapnet_mst_var_3469441_ALT.haps.list",
              "out_Ace1.Hapnet_mst_var_3481632_ALT.haps.list",
              "out_Ace1.Hapnet_mst_var_3504796_ALT.haps.list")) {
    
    c_ldbased_280S = read.table(ldi, header = F, sep="\t")
    c_ldbased_280S = stringr::str_replace(c_ldbased_280S$V1, pattern = " ", replacement = "")
    
    
    venn.two(list1 = c_umap_280S, list2 = c_ldbased_280S,
             catname1 = "UMAP", catname2 = "LD", 
             col1 = "springgreen3", col2 = "orange3",
             main = sprintf("overlap haps UMAP v. %s", ldi),eulerbool = T)
    
}

# 
# stop("ara")
# for (ldi in c("out_Ace1.Hapnet_mst_var_3465693_ALT.samp.list",
#               "out_Ace1.Hapnet_mst_var_3469441_ALT.samp.list",
#               "out_Ace1.Hapnet_mst_var_3481632_ALT.samp.list",
#               "out_Ace1.Hapnet_mst_var_3504796_ALT.samp.list")) {
#     
#     c_ldbased_280S = read.table(ldi, header = F, sep="\t")
#     c_ldbased_280S = unique(stringr::str_replace(c_ldbased_280S$V1, pattern = " ", replacement = ""))
#     
#     # get lists of haps
#     c_umap_280S = unique(c_umap[c_umap$p2_dtc_prediction == 1,"ox_code_sample"])
#     
#     venn.two(list1 = c_umap_280S, list2 = c_ldbased_280S,
#              catname1 = "UMAP", catname2 = "LD", 
#              col1 = "springgreen3", col2 = "purple",
#              main = sprintf("overlap samples UMAP v. %s", ldi),eulerbool = T)
#     
#     
# }
