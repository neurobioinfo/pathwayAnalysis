#############################################################################################
## Running enrichment pathway analysis
#
# This workflow uses the gene ontology from Gene Ontology (GO) project. 
# It performs over-representation (ORA) and path enrichment (GSEA) analysis of GO terms and KEEG pathways
# in a list of significant genes using “ClusterProfiler” .
# ORA analysis performs statistical enrichment  analysis using hypergeometric testing of significant DE genes against 
# a background DE gene list.
# 
# Input: a significant DE gene list and a background DE gene list 
#
# You can run this workflow in a Rstudio session in your PC or, if you're using Compute Canada
# Clusters you can run the workflow in an interactive R session.
#
### 1. SETTING THE ENVIRONMENT
### Follow all setting environment commands (steps 1-3)  in the README.md file before running this script
#
## 2. After setting the environment , run next commands following special instructions for either case:
#     Running the analysis in your local Rstudio version or in CC clusters
#
# 3. If you're using your PC RScript version skip this step. 
#    If you're using CC clusters, load the bioconductor module,
#    and launch R:
#salloc --time=02:00:00 -N 1 --ntasks-per-node=8 --mem-per-cpu=4G --account=rrg-grouleau-ac
#module load StdEnv/2020  gcc/9.3.0 r-bundle-bioconductor/3.14
#R
#
# 4. Install required packages: tidyverse, DOSE, pathview, clusterProfiler, enrichplot, org.Hs.eg.db, ggupset:
# If you're using your PC, open a RStudio session and run this workflow beginning in this step.
# In RStudio: open the run_enrichment_pathway_analysis.R script and begin analysis running following commands in the script.
# If you're running in a CC cluster, run following commands in your interactive session.

########################################################################################

# 0. Set the <project_directory> as working directory and source the enrichment_pathway_analysis_2.R script.
setwd("/home/geo/Documents/neuro/projects/Rhalena_scRNA_DGE_analysis/NPC_cell_based_analysis/DAneurons")
set.seed(1)
source("/home/geo/Documents/neuro/projects/enrichment_pathway/pathwayAnalysis/R/enrichment_pathway_analysis_2.R")

DGE_universe_file="data/<DEG_file_name.csv>"
DGE_significant_file="data/DEG_sig_file_name.csv"
annotationType="<annotatyon type SYMBOL or ENTREZ>"
inputType="DGE package name: Deseq2, Seurat, sleuth"
go_ora_results="results/GO_ORA"
go_gsea_results="results/GO_GSEA"
kegg_ora_results="results/KEGG_ORA"
kegg_gsea_results="results/KEGG_GSEA"

## ANALYSIS BEGINS HERE

# Generate EnrichGo objects and significant GO terms files for three Go categories:
#    biological processes (BP), Molecular function (MF) and cellular components (CC)

# 1. Read universe gene list
DEG_universe<-read.csv(DGE_universe_file, header = TRUE, sep = ",", dec=".")

# Rename gene_symbol and fold change columns
DEG_universe<-rename_columns(DEG_universe, inputType)
head(DEG_universe)

# 2. Subset differentially expressed genes
DEG_list <- subset(DEG_universe, padj < 0.05)
write.csv(as.data.frame(DEG_list), file=DGE_significant_file)
# Use this comand to read an existing significant DGE list
# DEG_list<-read.csv(DGE_significant_file, header = TRUE, sep = ",", dec=".")

# 3. Extract Genes ids list
genes_universe<-extract_genesID(DEG_universe)
genes_list<-extract_genesID(DEG_list)

# 4. Run GO Overrepresentation analysis. this function uses "enrichGO" from clusterProfiler
go_BP<-go_enrichment_obj(genes_list, genes_universe, "BP", 0.05, 0.6, go_ora_results, annotationType)
go_MF<-go_enrichment_obj(genes_list, genes_universe, "MF", 0.05, 0.6, go_ora_results, annotationType)
go_CC<-go_enrichment_obj(genes_list, genes_universe, "CC", 0.05, 0.6, go_ora_results, annotationType)
go_ALL<-go_enrichment_obj(genes_list, genes_universe, "ALL", 0.05, 0.6, go_ora_results, annotationType)

# 5. GO enrichment analysis 
## extract the log2 fold changes from our results table 
dim(DEG_universe)
res_entrez<-add_entrezid(DEG_universe)
foldchanges<-name_foldchanges(res_entrez, annotationType)

gseaGo_BP<-gseaGO_analysis(foldchanges, "BP", 0.5, go_gsea_results, annotationType)
gseaGo_MF<-gseaGO_analysis(foldchanges, "MF", 0.5, go_gsea_results, annotationType)
gseaGo_CC<-gseaGO_analysis(foldchanges, "CC", 0.5, go_gsea_results, annotationType)
gseaGo_ALL<-gseaGO_analysis(foldchanges, "ALL", 0.5, go_gsea_results, annotationType)

## 6. Generate dotplots, emaplots and cnetplots
# i. Generate doptplots
go_dotPlot(go_BP, go_MF, go_CC, go_ALL, go_ora_results, "ORA")
go_dotPlot(gseaGo_BP, gseaGo_MF, gseaGo_CC, gseaGo_ALL, go_gsea_results, "GSEA")

# In RStudio: Look at the dotplots and choose the go categories to show in emaplots and cnetplots.
#     In CC clusters: download dotplots and loot  at them to chose go categories to show in emaplots and cnetplots.
#                     Run this only if you're working in CC clusters. Without closing your current R session, 
#                     open a new terminal, and from a local terminal run next command to download dotplots:
#                      scp <CC_user_ID>@<CC_cluster_name>.computecanada.ca:<output_directory>/BP_go_dotplot.png.csv <local_project_directory>
#                      scp <CC_user_ID>@<CC_cluster_name>.computecanada.ca:<output_directory>/MF_go_dotplot.png.csv <local_project_directory>
#                      scp <CC_user_ID>@<CC_cluster_name>.computecanada.ca:<output_directory>/CC_go_dotplot.png.csv <local_project_directory>
                        
# Open the dotplots and chose the go terms of each category to show in the emaplot and cnetplot.
#   You can also chose to show a number of top significant go terms
#   In your Rstudio session, or in your R interactive session, create vectors with chosen categories. Examples:                         
#   Example: CC_cats=c("presynapse", "integral component of presynaptic membrane", 
#          "intrinsic component of presynaptic membrane", "presynaptic membrane")

# ii. Generate enrichment GO plot (emaplot)
# ORA number of term
BP_cats=dim(go_BP)[1]
MF_cats=dim(go_MF)[1]
CC_cats=dim(go_CC)[1]
ALL_cats=dim(go_ALL)[1]
go_emaPlot(go_BP, go_MF, go_CC, go_ALL, go_ora_results, BP_cats, MF_cats, CC_cats, ALL_cats, "ORA")

# GSEA number of term
BP_gsea_cats=dim(gseaGo_BP)[1]
MF_gsea_cats=dim(gseaGo_MF)[1]
CC_gsea_cats=dim(gseaGo_CC)[1]
ALL_gsea_cats=dim(gseaGo_ALL)[1]
go_emaPlot(gseaGo_BP, gseaGo_MF, gseaGo_CC, gseaGo_ALL, go_gsea_results, 
           BP_gsea_cats, MF_gsea_cats, CC_gsea_cats, ALL_gsea_cats, "GSEA")


# iv. Generate Cnetplots:
go_cnetPlot(DEG_list, go_BP, go_MF, go_CC, go_ALL, go_ora_results, BP_cats, MF_cats, CC_cats, ALL_cats, "ORA")
go_cnetPlot(DEG_list, gseaGo_BP, gseaGo_MF, gseaGo_CC, gseaGo_ALL, go_gsea_results, 
            BP_gsea_cats, MF_gsea_cats, CC_gsea_cats, ALL_gsea_cats, "GSEA")


# 7. GSEA on KEEG pathways
# Functional class scoring (FCS) tools, such as GSEA, most often use the gene-level
# statistics or log2 fold changes for all genes from the differential expression results.
# Then it looks whether gene sets for particular biological pathways are enriched
# among the large positive or negative fold changes.
# 
# Its important to run this analysis with all genes from the DEG analysis.
# It is possible to run this analysis with a subset of genes, but this reduces power test.

# i Run GSEA analysis, this function uses gseKEGG from "ClusterProfiler" to find
# KEGG pathways, and it uses Pathview to generate pathway images
res_entrez<-add_entrezid(DEG_universe)
fc_entrez<-name_foldchanges(res_entrez, "ENTREZ")

gseaKEGG<-gseaKEGG_analysis(fc_entrez, kegg_gsea_results, "gseaKEGG_pathways.csv")

# ii Look at the pathways csv file and generate a list of interesting pathways
# example
pathways<-c("hsa00515", "hsa04512")

# iii Create GSEAplot and KEGG image for chosen pathways
# This function uses and it uses Pathview to generate pathway images
go_gseaKEGGplot(gseaKEGG, fc_entrez, pathways, kegg_gsea_results)

# 8. ORA analysis on KEGG pathway
# You can use this function to test if enrichment of KEEG pathways in a subset of genes. 
# For example, you can run this function for a list of significant DEG genes.
# i. Create a gene Rentrez IDs lists from the DEG subset list 
res_entrez_subset<-add_entrezid(DEG_list)
res_entrez<-add_entrezid(DEG_list)
dim(res_entrez_subset)
# ii Run the enrichKEGG analysis 
KEGGenrich<-enrichKEGG_analysis(res_entrez_subset, res_entrez, kegg_ora_results)


###################################################################
#Upset plots
###################################################

# 9. heatplot 
# Gsea heatplot
heatplot1(gseaGo_ALL, foldchanges, "gsea", results)
#GO heatplot
heatplot1(go_ALL, foldchanges, "ora", results)

# 10. Upsetplot
gsea_upsetPlot(gseaGo_ALL, results)
ora_upsetPlot(go_ALL, results)

