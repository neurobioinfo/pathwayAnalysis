## Running enrichment pathway analysis

# This workflow uses the gene ontology from Gene Ontology (GO) project. 
# It performs an over-representation analysis of GO terms in a list of significant genes using “ClusterProfiler” .
# "ClusterProfiler" performs statistical enrichment  analysis using hypergeometric testing of significant DE genes against 
# a background DE gene list.
# 
# Input: a significant DE gene list and a background DE gene list 
 
# You can run this workflow in a Rstudio session in your PC or, if you're using Compute Canada
# Clusters you can run the workflow in an interactive R session.

### 1. SETTING THE ENVIRONMENT
### Follow all setting environment commands (steps 1-3) are in the README.md file before running this script

## 2. After setting the environment , run next commands following special instructions for either case:
#     Running the analysis in your local Rstudio version or in CC clusters

# 3. If you're using your PC RScript version skip this step. 
#    If you're using CC clusters, open a Salloc session, load the bioconductor module,
#    and launch R:
salloc --time=02:00:00 -N 1 --ntasks-per-node=8 --mem-per-cpu=4G --account=rrg-grouleau-ac
module load StdEnv/2020  gcc/9.3.0 r-bundle-bioconductor/3.14
R

# 4. Install required packages: tidyverse, DOSE, pathview, clusterProfiler, enrichplot, org.Hs.eg.db:
# If you're using your PC, open a RStudio session and run this workflow beginning in this step.
# In RStudio: open the run_enrichment_pathway_analysis.R script and begin analysis running following commands in the script.
# If you're running in a CC cluster, run following commands in your interactive session.

packages<-c(`tidyverse`, `DOSE`, `pathview`, `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`)
lapply(packages, library, character.only = TRUE)

# 5. Set the <project_directory> as working directory and source the enrichment_pathway_analysis_2.R script.
setwd("<project_directory>")
source("scripts/enrichment_pathway_analysis_2.R")
#    This script contains custom functions to run clusterProfiler for a pathway enrichment analysis

## ANALYSIS BEGINS HERE

# 6. Generate EnrichGo objects and significant GO terms files for three Go categories:
#    biological processes (BP), Molecular function (MF) and cellular components (CC)

# i. Read universe gene list
DEG_universe<-read.csv("data/<all_genes_DEG_list.csv>", header = TRUE, sep = ",", dec=".")

# ii. Read significant gene list
DEG_list<-read.csv("data/<significant_genes_DEG_list.csv>", header = TRUE, sep = ",", dec=".")

# iii. Extract Genes ids list
genes_universe<-extract_genesID(DEG_universe, Deseq2)
genes_list<-extract_genesID(DEG_list, Deseq2)

# iv. Run GO enrichment analysis. this function uses "enrichGO" from clusterProfiler
annotationType="SYMBOL"
go_BP<-go_enrichment_obj(genes_list, genes_universe, "BP", 0.05, "<output_directory>/BP_go.csv", annotationType)
go_MF<-go_enrichment_obj(genes_list, genes_universe, "MF", 0.05, "<output_directory>/MF_go.csv", annotationType)
go_CC<-go_enrichment_obj(genes_list, genes_universe, "CC", 0.05, "<output_directory>/CC_go.csv", annotationType)
go_ALL<-go_enrichment_obj(genes_list, genes_universe, "ALL", 0.05, "<output_directory>/ALL_go.csv", annotationType)

## 7. Generate dotplots, emaplots and cnetplots
# i. Generate doptplots
go_dotPlot(go_BP, go_MF, go_CC, go_ALL, "output_directory")

# ii. In RStudio: Look at the dotplots and choose the go categories to show in emaplots and cnetplots.
#     In CC clusters: download dotplots and loot  at them to chose go categories to show in emaplots and cnetplots.
#                     Run this only if you're working in CC clusters. Without closing your current R session, 
#                     open a new terminal, and from a local terminal run next command to download dotplots:
                      scp <CC_user_ID>@<CC_cluster_name>.computecanada.ca:<output_directory>/BP_go_dotplot.png.csv <local_project_directory>
                      scp <CC_user_ID>@<CC_cluster_name>.computecanada.ca:<output_directory>/MF_go_dotplot.png.csv <local_project_directory>
                      scp <CC_user_ID>@<CC_cluster_name>.computecanada.ca:<output_directory>/CC_go_dotplot.png.csv <local_project_directory>
                        
# iii. Open the dotplots and chose the go terms of each category to show in the emaplot and cnetplot.
#   You can also chose to show a number of top significant go terms
#   In your Rstudio session, or in your R interactive session, create vectors with chosen categories. Examples:                         

BP_cats=c("learning or memory", "central nervous system neuron differentiation",
          "synaptic vesicle exocytosis", "chemical synaptic transmission", "regulation of neurotransmitter levels")
MF_cats=5
CC_cats=c("presynapse", "integral component of presynaptic membrane", 
          "intrinsic component of presynaptic membrane", "presynaptic membrane")
ALL_cats=20

# iv. Generate enrichment GO plot (emaplot)
go_emaPlot(go_BP, go_MF, go_CC, go_ALL, "<output_directory>", BP_cats, MF_cats, CC_cats, ALL_cats)

# v. Generate Cnetplots:
go_cnetPlot(DEG_list, go_BP, go_MF, go_CC, "<output_directory>", BP_cats, MF_cats, CC_cats)


# 8. Gsea analysis
# Functional class scoring (FCS) tools, such as GSEA, most often use the gene-level
# statistics or log2 fold changes for all genes from the differential expression results.
# Then it looks whether gene sets for particular biological pathways are enriched
# among the large positive or negative fold changes.
# 
# Its important to run this analysis with all genes from the DEG analysis.
# It is possible to run this analysis with a subset of genes, but this reduces power test.
# Below, there is another method to test a subset of genes

# 0. Rename gene_symbol and fold change columns
DEG_universe<-rename_geneSymbol_column(DEG_universe, 1, "Deseq2")
DEG_list<-rename_geneSymbol_column(DEG_list, 1,	"Deseq2")

DEG_universe<-rename_foldChange_column(DEG_universe, 3, "Deseq2")
DEG_list<-rename_foldChange_column(DEG_list, 3, "Deseq2")

# i. Create a new list with gene Entrez IDs and expression fold changes from all genes DEG list
res_entrez<-add_entrezid(DEG_universe)
foldchanges<-name_foldchanges(res_entrez)

# ii Run GSEA analysis, this function uses gseKEGG from "ClusterProfiler" to find
# KEGG pathways, and it uses Pathview to generate pathway images
gseaKEGG<-gseaKEGG_analysis(foldchanges, "<output_directory>", "gseaKEGG_pathways.csv")

# iii Look at the pathways csv file and generate a list of interesting pathways
example
pathways<-c("hsa04360", "hsa05168", ...)

# iv Create GSEAplot and KEGG image for chosen pathways
# This function uses and it uses Pathview to generate pathway images
go_gseaKEGGplot(gseaKEGG, foldchanges, pathways, "<output_directory>")

# 9. KEGG enrichment analysis
# You can use this function to test if enrichment of KEEG pathways in a subset of genes. 
# For example, you can run this function for a list of significant DEG genes.
# i. Create a gene Rentrez IDs lists from the DEG subset list 
res_entrez_subset<-add_entrezid(DEG_list)

# ii Run the enrichKEGG analysis which uses 
KEGGenrich<-enrichKEGG_analysis(res_entrez_subset, "<output_directory>")
