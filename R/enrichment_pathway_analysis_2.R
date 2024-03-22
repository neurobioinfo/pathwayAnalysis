
#############################################################################
# This script contains functions for an Enrichment pathway analysis workflow
# An important part of this workflow is inspired from:
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html
###########################################################################


library(tidyverse)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library("ggupset")

# 1. Load gene list FROM CSV
# use read.csv function

# 2. Go overpresentation analysis
# Extract genes IDs
# Input: DEG gene list
# Output: genes ID list
extract_genesID<- function(DEG_geneList){
  genes_ids_list <- as.character(DEG_geneList$gene)
	return(genes_ids_list)
}

#3. ## Run GO enrichment analysis 
# Input: Significant genes id list
#        All genes id list
#        Go cathegory to test
#        q value cutoff
#        output file name
# Output: enrichGo object
#        enriched go_terms csv file
go_enrichment_obj<-function(genes_ids_list, universe_ids_list, go_cathegory, q_cutoff, SS_cutoff, output_dir, annotationType){
  file_name_full=(paste(output_dir, "/", go_cathegory, "_go.csv", sep=""))
  file_name_simpl=(paste(output_dir, "/", go_cathegory, "_go_simpl.csv", sep=""))
  
  go <- enrichGO(gene = genes_ids_list, 
                 universe = universe_ids_list,
                 keyType = annotationType,
                 OrgDb = org.Hs.eg.db, 
                 ont = go_cathegory, 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = q_cutoff, 
                 readable = TRUE)
  go_simpl <- simplify(go, cutoff=SS_cutoff, by="p.adjust", select_fun=min)
  print(data.frame(go))
  print(data.frame(go_simpl))
  go_cluster_summary <- data.frame(go)
  go_cluster_summary_simpl <- data.frame(go_simpl)
  write.csv(go_cluster_summary, file_name_full)
  write.csv(go_cluster_summary_simpl, file_name_simpl)
  return(go_simpl)
}

# 4. GO Gene Set Enrichment Analysis
gseaGO_analysis<-function(foldchanges, ont_cat, SS_cutoff, output_dir, annotationType){
  
  gseaGO_obj <- gseGO(geneList= foldchanges,
                      OrgDb= org.Hs.eg.db,
                      ont=ont_cat,
                      keyType = annotationType)
  
  # Simplify
  go_simpl <- simplify(gseaGO_obj, cutoff=SS_cutoff, by="p.adjust", select_fun=min)
  
  ## Extract the GSEA results
  gseaGO_results <- gseaGO_obj@result
  gseaGO_simpl_results <- go_simpl@result
  
  ## Write GSEA results to file
  file_name_full=(paste(output_dir, "/", ont_cat, "_gsea_GO_.csv", sep=""))
  file_name_simpl=(paste(output_dir, "/", ont_cat, "_gsea_GO_simpl_005.csv", sep=""))
  # View(gseaGO_results)
  write.csv(gseaGO_results, file_name_full)
  write.csv(gseaGO_simpl_results,  file_name_simpl)
  
  return(go_simpl)
}


#4. generate dotplots
# Input: BP_go Biological process (BP) EnrichGo object
#        MF_go Molecular function (MF) EnrichGo object
#        CC_goCelular component (CC) EnrichGo object
# Output: Dotplots for each GO cathegory
go_dotPlot<-function(BP_go, MF_go, CC_go, ALL_go, output_dir, analysisType){
  
  go_objs=c(BP_go, MF_go, CC_go, ALL_go)
  go_cats=c("BP", "MF", "CC", "ALL")
  i=1
  for (go_obj in go_objs){
    dotPlot_file=(paste(output_dir, "/", go_cats[i], "_", analysisType, "_dotplot.png", sep=""))
    png(file=dotPlot_file, width=700, height=1000)
    dp=dotplot(go_obj, showCategory=50)
    #print(dp)
    try(print(dp), silent = TRUE)
    dev.off()
    i=i+1
  }
}


#5. Generate enrichment GO plot, which shows the relationship between the top 50 most
# significantly enriched GO terms (padj.), by grouping similar terms together.
# Input: BP_go BP EnrichGo object
#        MF_go MF EnrichGo object
#        CC_go CC EnrichGo object
#        output_dir output directory
#        BP_cats BP go terms list or Number of top go terms 
#        MF_cats MF go terms list or Number of top go terms
#        CC_cats CC go terms list or Number of top go terms
# Output: Enrichment GO plots for each cathegory
go_emaPlot<-function(go_BP, go_MF, go_CC, go_ALL, output_dir, BP_cats, MF_cats, CC_cats, ALL_cats, analysisType){
  
  go_objs=c(go_BP, go_MF, go_CC, go_ALL)
  go_terms=list(BP_cats, MF_cats, CC_cats, ALL_cats)
  go_cats=c("BP", "MF", "CC", "ALL")
  i=1
  for (go_obj in go_objs){
    print(paste("generate Emaplot for", go_cats[i]))
    
    tryCatch({
      x2 <- pairwise_termsim(go_obj)
      plot_file=(paste(output_dir, "/", go_cats[i], "_", analysisType, "_emaplot.png", sep=""))
      par(mfrow=c(1,1))
      png(file=plot_file, width=800, height=600)
      ep=emapplot(x2, showCategory=go_terms[[i]])
      print(ep)
      dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    i=i+1
  }
}

# 6. Generate category netplots:
# the category cnetplot shows the relationships between the genes associated with the
# significant GO terms and fold changes of the significant associated genes
# Input: DEG_list Significant DEG list including fold change
#        BP_go BP EnrichGo object
#        MF_go MF EnrichGo object
#        CC_go CC EnrichGo object
#        output_dir output directory
#        BP_cats BP go terms list or Number of top go terms 
#        MF_cats MF go terms list or Number of top go terms
#        CC_cats CC go terms list or Number of top go terms
# Output: Cathegory cnetplot for each GO cathegory
go_cnetPlot<-function(DEG_list, go_BP, go_MF, go_CC, go_ALL, output_dir, 
                      BP_cats, MF_cats, CC_cats, ALL_cats, analysisType){
  
  ## extract the log2 fold changes from our results table 
  foldchanges <- DEG_list$b
  names(foldchanges) <- DEG_list$gene
  
  ## Cnetplot details the genes associated with one or more terms - by default gives
  # the top 5 significant terms (by padj)
  go_objs=c(go_BP, go_MF, go_CC, go_ALL)
  go_terms=list(BP_cats, MF_cats, CC_cats, ALL_cats)
  go_cats=c("BP", "MF", "CC", "ALL")
  i=1
  for (go_obj in go_objs){
    
    tryCatch({
      
      plot_file=(paste(output_dir, "/", go_cats[i], "_", analysisType, "_cnetplot.png", sep=""))
      par(mfrow=c(1,1))
      png(file=plot_file, width=800, height=600)
      cp=cnetplot(go_obj, 
                  categorySize="padj", 
                  showCategory=go_terms[[i]], 
                  foldChange=foldchanges, 
                  vertex.label.font=6)
      print(cp)
      dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    i=i+1
  }
}

#7. Functional class scoring methods
#Functional class scoring (FCS) tools, such as GSEA, most often use the gene-level
# statistics or log2 fold changes for all genes from the differential expression results,
# then look to see whether gene sets for particular biological pathways are enriched
# among the large positive or negative fold changes.

# However, an exception exists for RNA-seq datasets where GSEA may benefit from the
# removal of extremely low count genes (i.e., genes with artifactual levels of expression
# such that they are likely not actually expressed in any of the samples in the dataset).
#https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html

# Gene set enrichment analysis uses clusterProfiler and Pathview
# To perform GSEA analysis of KEGG gene sets. clusterProfiler requires the genes to be
# identified using Entrez IDs for all genes in our results dataset. 

# 0. Rename gene_symbol, foldChange and pvalue column
rename_columns<-function(DEG_list, inputType){
  if (inputType=="Seurat") {
    colnames(DEG_list)[6] ="padj"
    colnames(DEG_list)[1] ="gene"
    colnames(DEG_list)[3] ="b"
  } else if (inputType=="Deseq2") {
    colnames(DEG_list)[7] ="padj"
    colnames(DEG_list)[1] ="gene"
    colnames(DEG_list)[3] ="b"
  } else {
    return(DEG_list)
  }
  return(DEG_list)
}


#Add ENTREZID notation to DEG table
# Input: DEG_list Significant DEG list including fold change
# Output: DEG list including Entrez id
add_entrezid<-function(DEG_list){
  
  ## Remove any NA values (reduces the data by quite a bit)
  nrow(DEG_list)
  res_entrez <- dplyr::filter(DEG_list, "gene" != "NA")
  print("Remove NAs")
  print(nrow(res_entrez))
  
  ## Remove any Entrez duplicates
  print("Remove any Entrez duplicates")
  res_entrez <- res_entrez[which(duplicated(res_entrez$gene) == F), ]
  print("Remove Entrez duplicates")
  print(nrow(res_entrez))
  
  ## Name each fold change with the corresponding gene name
  hs <- org.Hs.eg.db
  
  my.symbols <- c(res_entrez$gene)
  entrezid <- AnnotationDbi::select(hs, 
                                    keys = my.symbols,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
  
  entrezid <- entrezid[which(duplicated(entrezid$SYMBOL) == F), ]
  if (all(res_entrez$gene==entrezid$SYMBOL)){
    res_entrez$ENTREZID<-entrezid$ENTREZID
  }else{
    print("Error: The number of entrezids doesn't match the number of genes in list")
    return
  }
  res_entrez2<-res_entrez %>% filter(!is.na(ENTREZID))
  return(res_entrez2)
}


## Create a dataframe including fold changes and gene Entrez IDs (foldChanges-EntrezIDs)
# Input: significant DEG list including Entrez IDs
# Output: foldchanges-entrezid
name_foldchanges<-function(res_entrez, annotationType){
  foldchanges <- res_entrez$b
  if (annotationType=="SYMBOL"){
    names(foldchanges) <- res_entrez$gene
  }else if (annotationType=="ENTREZ"){
    names(foldchanges) <- res_entrez$ENTREZID
  } 
  ## Sort fold changes in decreasing order
  foldchanges <- sort(foldchanges, decreasing = TRUE)
  foldchanges <- foldchanges[!duplicated(names(foldchanges))]
  return(foldchanges)
}


## GSEA analysis
# Input: foldChanges-EntrezIDs
# Output: GSEA KEGG pathways csv file
#         gseaKEGG object
gseaKEGG_analysis<-function(foldchanges, output_dir, file_name){
  ## GSEA using gene sets 
  gseaKEGG_obj <- gseKEGG(foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                      organism = "hsa", # supported organisms listed below
                      minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                      eps=0,
                      pvalueCutoff = 0.05,
                      verbose = FALSE)
  
  ## Extract the GSEA results
  gseaKEGG_results <- gseaKEGG_obj@result
  
  ## Write GSEA results to file
  file_name=(paste(output_dir, "/", file_name, sep=""))
  View(gseaKEGG_results)
  write.csv(gseaKEGG_results, file_name)
  
  return(gseaKEGG_obj)
  
}

#8. Generate GSEA plot and KEGG pathway plot
#  Input: gseaKEGG  object
#         foldChanges-EntrezIDs
#         pathways to generate image list
#         output_dir output directory
#         prefix to add to KEGG image name file
go_gseaKEGGplot<-function(gseaKEGG_obj, foldchanges,pathways_list, output_dir, prefix_name){
  
  current_dir<-getwd()
  setwd(paste(current_dir, "/", output_dir, sep=""))
  
  for (path in pathways_list){
    plot_file=(paste(prefix_name, "_", path, "_gseaplot.png", sep=""))
    png(file=plot_file, width=800, height=600)
    gp<-gseaplot(gseaKEGG_obj, geneSetID = path)
    print(gp)
    dev.off()
    
    ## Use the Pathview R package to integrate the KEGG pathway data from clusterProfiler into pathway images:
    #detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts
    
    ## Output images for a single significant KEGG pathway
    pathview(gene.data = foldchanges,
             pathway.id = path,
             species = "hsa",
             limit = list(gene = 2, # value gives the max/min limit for foldchanges
                          cpd = 1))
  }
  setwd(current_dir)
}

# 8. KEGG pathway over-representation analysis
# This function performs an enrichment pathway analysis from a subset of genes 
# Input: Relevant DEG list including gene Entrez IDs
#        output_dir
# Output: enrich KEGG pathways csv file

enrichKEGG_analysis<-function(res_entrez, output_dir){
  genes=res_entrez$ENTREZID
  yy = enrichKEGG(genes, organism="hsa", pvalueCutoff=0.05)
  enrichKEGG_results<-yy@result
  head(enrichKEGG_results)
  #plot(yy)
  
  ## Write results to file
  file_name=(paste(output_dir, "/enrichKEGG_results.csv", sep=""))
  write.csv(enrichKEGG_results, file_name)
  return( enrichKEGG_results)
  
  gseaKEGG_results <- gseaKEGG_obj@result
  
  ## Write GSEA results to file
  file_name=(paste(output_dir, "/", file_name, sep=""))
  View(gseaKEGG_results)
  write.csv(gseaKEGG_results, file_name)
  
  return(gseaKEGG_obj)
}

## (. Extract KEGG pathways gene set
# Input: pathway enrichment csv file
# Output: Gene set txt file
# Warning: this function is in development

get_geneID<-function(gseaKEGG_obj, output_dir){
  library(annotate)
  library(stringr)
  library(purrr)
  
  ## Extract the GSEA results
  gsea_results <- transpose(gseaKEGG_obj@result, .names = NULL)
  
  for (pathway in gsea_results){
    path_id=pathway[1]
    print(path_id)
    entrezid_genes<-as.character(pathway[11])
    print(entrezid_genes)
    entrezid_genes<-strsplit(entrezid_genes, split = "/")
    entrezid_genes<-unlist(entrezid_genes)
    gene<-getSYMBOL(as.character(entrezid_genes), data='org.Hs.eg')
    file_name=paste(output_dir, "/", path_id, "_gene_set.txt", sep="")
    write.table(gene, file=file_name)
  }
}

# 10. Get geneset list from gsea_results
get_geneID_fromFile<-function(gsea_results_file){
  library(annotate)
  library(stringr)
  gsea_results<-read.csv("gseaKEGG_results.csv")
  entrezid_genes<-gsea_results[4,12]
  entrezid_genes<-strsplit(entrezid_genes, split = "/")
  entrezid_genes<-unlist(entrezid_genes)
  gene<-getSYMBOL(as.character(entrezid_genes), data='org.Hs.eg')
  write.table(gene, file="hsa04360_gene_set.txt")
}


# 11. GO Gene Set Enrichment Analysis
gseaGO_analysis0<-function(foldchanges, output_dir, file_name){
  ## GSEA using gene sets 
#  gseaGO_obj <- gseGO(geneList     = foldchanges,
#                OrgDb        = org.Hs.eg.db,
#                ont          = "ALL",
#                minGSSize    = 5,
#                maxGSSize    = 500,
#                pvalueCutoff = 0.05,
#                verbose      = FALSE)
  gseaGO_obj <- gseGO(geneList= foldchanges,
                      OrgDb= org.Hs.eg.db)
  
  ## Extract the GSEA results
  gseaGO_results <- gseaGO_obj@result
  
  ## Write GSEA results to file
  file_name=(paste(output_dir, "/", file_name, sep=""))
  # View(gseaGO_results)
  write.csv(gseaGO_results, file_name)
  
  return(gseaGO_obj)
}


######### IN DEVELOPMENT #############################
# This functions aren't in the script enrichment_pathway_analysis_2.R'

heatplot1<-function(gseaGO_obj, foldchanges, type, output_dir){
## convert gene ID to Symbol
  goenx <- setReadable(gseaGO_obj, 'org.Hs.eg.db', 'ENTREZID')

  # Heatmap-like functional classification
  # look here to fix the heayplot: https://support.bioconductor.org/p/9148486/
  heatplot_file=(paste(output_dir, "/", "heatplot_", type, ".png", sep=""))
  png(file=heatplot_file, width=1000, height=700)
  p2 <- heatplot(goenx, foldChange=foldchanges, showCategory=5)
  try(print(p2), silent = TRUE)
  dev.off()
}

# 13. Upsetplot
#Upset plot

gsea_upsetPlot<-function(gseaGO_obj, output_dir){
  gseaGO_obj <- setReadable(gseaGO_obj, 'org.Hs.eg.db', 'ENTREZID')
  gseaGO_obj2 <- pairwise_termsim(gseaGO_obj)
  upsetplot_file=(paste(output_dir, "/", "gsea_upsetplot.png", sep=""))
  png(file=upsetplot_file, width=1000, height=700)
  p1<-upsetplot(gseaGO_obj2)
  try(print(p1), silent = TRUE)
  dev.off()
}

# upset plot for enrichment analysis
ora_upsetPlot<-function(go_obj, output_dir){
  
  go_obj <- pairwise_termsim(go_obj)
  upsetplot_file=(paste(output_dir, "/", "ora_upsetplot.png", sep=""))
  png(file=upsetplot_file, width=1000, height=700)
  p3<-upsetplot(go_obj)
  try(print(p3), silent = TRUE)
  dev.off()
  
}


