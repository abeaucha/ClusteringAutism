# ----------------------------------------------------------------------------
# mouse_cluster_enrichment.R
# Authors: Antoine Beauchamp, Jacob Ellegood
# Created: November 8th, 2023


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(tmod))


# Functions ------------------------------------------------------------------

# NOTE: These functions were copied from 
# /projects/jacob/ClusteringAutism_125Models_Mar2020/Code/Clustering_Functions.R

GetValidGeneClusterList<-function(Cluster_Dir="Data/Outputs/Clusters",
                                  boot="N",
                                  boot_pick=nboot){
  
  if(boot=="Y"){
    cl <- read.csv(glue("{Cluster_Dir}/Boot/Clusters_{boot_pick}.csv"))
  } else {
    cl <- read.csv(glue("{Cluster_Dir}/Clusters.csv"))
  }
  
  cl <- cl %>% dplyr::mutate(Model=X)%>% dplyr::select(-X)
  
  # Get gene names
  cl$Gene <- sapply(strsplit(as.character(cl$Model),"\\("), "[[", 1)
  cl$Gene[cl$Gene=="Andr"] <- "Ar"
  cl$Gene[cl$Gene=="Caspr2"] <- "Cntnap2"
  cl$Gene[cl$Gene=="Dat"] <- "Slc6a3"
  cl$Gene[cl$Gene=="Mor"] <- "Oprm1"
  cl$Gene[cl$Gene=="Nl1"] <- "Nlgn1"
  cl$Gene[cl$Gene=="Nl3"] <- "Nlgn3"
  cl$Gene[cl$Gene=="Nrxn1a"] <- "Nrxn1"
  cl$Gene[cl$Gene=="Sert"] <- "Slc6a4"
  cl$Gene[cl$Gene=="Pcdh"] <- "Pcdhga3"
  cl$Gene[cl$Gene=="Chd7;En1Cre"] <- "Chd7"
  cl$Gene[cl$Gene=="Ube3a.2"] <- "Ube3a"
  cl$Gene[cl$Gene=="FusDelta14"] <- "Fus"
  cl$Gene[cl$Gene=="Nr1a"] <- "Nmdar1"
  cl$Gene[cl$Model=="itsn1(+/+);itsn2(-/-)"]<- "itsn2"
  cl$Gene[cl$Model=="Snf2H(+/+);Snf2L(-/-);emxcre"]<-"Snf2l"
  cl$Gene[cl$Gene=="Snf2L"] <- "Smarca1"
  cl$Gene[cl$Gene=="Snf2H"] <- "Smarca5"
  cl$Gene[cl$Model=="Gsk3(a)"]<- "Gsk3A"
  cl$Gene[cl$Model=="Gsk3(B)"]<- "Gsk3B"
  
  # Remove CNVs, chromosomal, and behavioural models
  valid_indices <- (! (cl$Gene %in% c("15q11-13", "16p11.2", "22q11.2", "XO", "Btbr", "Balbc", "MAR","15q25","TCDD","VPA","BtbrTT")))
  cl_filt <- cl[valid_indices, ]
  return(cl_filt)
}


GetGeneNeighbourhood<-function(GeneScore="400",
                               total_clusters="10",
                               stringdb_version = "11.5",
                               output_dir="Data/Outputs/NeighbourhoodInfo/"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  all_df_list <- list()
  gx_data_list <- list()
  gscore <- GeneScore
  
  if (stringdb_version == "11.5") {
    stringdb_url <- "https://version-11-5.string-db.org/"
  } else if (stringdb_version == "12.0") {
    stringdb_url <- "https://string-db.org/"
  } else {
    stop()
  }
  
  for (num_clusters in 2:total_clusters) {
    all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", num_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", num_clusters))
    
    # Remove duplicates where multiple genes are in the same cluster
    all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
    table(all_df$cluster)
    
    # Zero based index
    all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
    
    # Get STRINGids and merge to all_df
    api_call <- paste(stringdb_url, "api/tsv/get_string_ids?identifiers=", paste(all_df$gene, collapse = "%0D"), "&species=10090&limit=1", sep="")
    string_ids <- read_tsv(api_call) %>% dplyr::select(queryIndex, stringId, preferredName)
    all_df %<>% left_join(string_ids, by="queryIndex")
    
    # Remove organism ID
    all_df$fixed_id <- all_df$stringId %>% strsplit("\\.") %>% map_chr(`[`, 2)
    write.csv(file = glue("{output_dir}/all_gene_data.csv"), x = all_df)
    
    ##########################################
    # This is where you can modify the gene score for the interactions - Normally set to 900.  Other common ones are 400 and 700
    
    # Get interaction data
    api_call <- paste(stringdb_url, "api/tsv/interaction_partners?identifiers=", paste(all_df$stringId, collapse = "%0D"), "&species=10090&required_score=",gscore, sep="")
    gx_data <- read_tsv(api_call)
    gx_data %<>% left_join(all_df %>% dplyr::select(cluster, stringId), by=c("stringId_A"="stringId")) %>% dplyr::rename(cluster_A=cluster)
    
    # Write out gene neighbourhoods for each cluster
    base_set <- unique(c(gx_data$preferredName_A, gx_data$preferredName_B))
    write.table(base_set, file=glue("{output_dir}/base_set_full_neighbourhood_{GeneScore}.txt"), quote = F, append=F, row.names = F, col.names = F)
    for (clust in 1:num_clusters) {
      target_set <- gx_data %>% filter(cluster_A==clust) %>% dplyr::select(preferredName_A, preferredName_B) %>% gather %>% .$value %>% unique 
      target_outfile <- glue("{output_dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood.txt")
      write.table(target_set, file=target_outfile, quote = F, append=F, row.names = F, col.names = F)
      
      target_set_only_models <- all_df %>% filter(cluster==clust) %>% .$gene
      target_outfile_only_models <- glue("{output_dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_only_models.txt")
      write.table(target_set_only_models, file=target_outfile_only_models, quote = F, append=F, row.names = F, col.names = F)
    }
    
    for (clust in 1:num_clusters) {
      target_set <- gx_data %>% filter(cluster_A!=clust) %>% dplyr::select(preferredName_A, preferredName_B) %>% gather %>% .$value %>% unique 
      target_outfile <- glue("{output_dir}/target_set_cluster{num_clusters}_Not_{clust}_{GeneScore}_neighbourhood.txt")
      write.table(target_set, file=target_outfile, quote = F, append=F, row.names = F, col.names = F)
      
      target_set_only_models <- all_df %>% filter(cluster!=clust) %>% .$gene
      target_outfile_only_models <- glue("{output_dir}/target_set_cluster{num_clusters}_Not_{clust}_{GeneScore}_only_models.txt")
      write.table(target_set_only_models, file=target_outfile_only_models, quote = F, append=F, row.names = F, col.names = F)
    }
    
    all_df_list[[paste0("Cluster", num_clusters)]] <- all_df
    gx_data_list[[paste0("Cluster", num_clusters)]] <- gx_data
  }
  
  listfname<-glue("{output_dir}/gene_interaction_tables_for_all_clusters_{GeneScore}.RData")
  save(list = c("all_df_list", "gx_data_list"), file = listfname)
}


get_homologs <- function(genes, species, ordered=T) {
  if (!exists("hom", envir = globalenv())) {
    hom <<- read_tsv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt")
  }
  species <- switch (species,
                     mouse={"mouse, laboratory"},
                     human={"human"}
  )
  
  if (ordered) {
    genes_human <- genes_mouse <- c()
    prog <- txtProgressBar(max=length(genes), style = 3)
    for (i in 1:length(genes)) {
      g <- genes[i]
      homologene_ids <- hom$`DB Class Key`[which((hom$Symbol==g) & hom$`Common Organism Name`==species)]
      genes_mouse <- c(genes_mouse, hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="mouse, laboratory"))])
      genes_human <- c(genes_human, hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="human"))])
      setTxtProgressBar(prog, i)
    }  
    close(prog)
  } else {
    # Faster, but unordered
    homologene_ids <- hom$`DB Class Key`[which((hom$Symbol %in% genes) & hom$`Common Organism Name`==species)]
    genes_mouse <- hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="mouse, laboratory"))]
    genes_human <- hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="human"))]
  }
  
  return(list(input_genes=genes, mouse_genes=genes_mouse, human_genes=genes_human))  
}


GetNeighbourhoodEnrichment<-function(GeneScore="900",
                                     total_clusters="10",
                                     Bader_List="Data/Raw/Human_Reactome_March_01_2020_symbol.gmt",
                                     background_set = "Data/Raw/sagittal_gene_table_normalized_filtered.csv",
                                     Neighbour_dir="Data/Outputs/NeighbourhoodInfo/",
                                     output_dir="Data/Outputs/NeighbourhoodEnrichment/"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", total_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", total_clusters))
  all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
  all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
  
  # First get different base sets
  
  base_brain <- read.csv(background_set, header = T)$msg.genes.acronym %>% as.character
  base_models_neighbourhood <- read.csv(glue("{Neighbour_dir}/{GeneScore}/base_set_full_neighbourhood_{GeneScore}.txt"), header = F)$V1  %>% as.character
  base_models_only <- all_df$gene
  base <- list(brain=base_brain, models_neighbourhood=base_models_neighbourhood, models_only=base_models_only)
  gmt_file<-Bader_List
  MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")
  
  for (num_clusters in 2:total_clusters) {
    
    target <- list(models_neighbourhood=list(), models_only=list())
    for (clust in 1:num_clusters) {
      
      if (file.size(glue("{Neighbour_dir}/{GeneScore}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood.txt")) == 0) {
        next
      } else {
        target_set_cl <- read.csv(glue("{Neighbour_dir}/{GeneScore}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood.txt"), header = F)$V1 %>% as.character
        target$models_neighbourhood[[clust]] <- target_set_cl
        
        target_set_cl <- read.csv(glue("{Neighbour_dir}/{GeneScore}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_only_models.txt"), header = F)$V1 %>% as.character
        target$models_only[[clust]] <- target_set_cl
      }
    }
    
    #####################
    # Choose gene sets
    
    for (clust in 1:num_clusters) {
      base_set <- base$brain
      target_set <- target$models_neighbourhood[[clust]]
      
      if (length(target_set) == 0){
        next
      } else {
        
        print(paste("Working on num_clusters:", num_clusters, "| cluster =", clust, "| score =", GeneScore))
        
        # Modified for the new Bader List
        
        enrichment_object <-  MappingsGMT
        fg <- get_homologs(as.character(target_set), "mouse", ordered = F)$human_genes
        bg <- get_homologs(as.character(base_set), "mouse", ordered = F)$human_genes
        
        tmodObj <- tmodHGtest(fg, bg, qval = 1.1, mset = enrichment_object, filter = F, order.by = "pval")
        tmodObj %<>% rowwise() %>% mutate(P.Value = P.Value) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
        tmodObj %<>% arrange(P.Value)
        tmodObj$rank <- 1:nrow(tmodObj)
        tmodObj %<>% mutate(NLQ = -log(adj.P.Val))
        tmodObj %<>% mutate(cluster = clust)
        tmodObj %<>% dplyr::select(cluster, rank, Title, NLQ, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
        
        write.csv(tmodObj, file = glue("{output_dir}/NewBader_enrichment_clusterneighbourhood_vs_brain_all_{num_clusters}_{clust}_{GeneScore}.csv"), row.names = F)
      }
    }
  }
}


# Main -----------------------------------------------------------------------

PROJECTPATH <- Sys.getenv("PROJECTPATH")

# Directory containing mouse cluster files
cluster_dir <- "/projects/jacob/ClusteringAutism_125Models_Mar2020/Data/Outputs/Clusters_Paper/"

# Enrichment output directory
enrichment_dir <- file.path(PROJECTPATH, "data/mouse/enrichment/")
if (!dir.exists(enrichment_dir)) {
  dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
}

# Import mouse clusters and filter for models to use for enrichment
cl_filt <- GetValidGeneClusterList(Cluster_Dir = cluster_dir,
                                   boot = "N", boot_pick = "100")

# Parameters for enrichment pipeline

# Max number of clusters
nk_max <- 10

# Gene score for StringDB
# gene_scores <- c(850, 900, 950)
gene_scores <- seq(400, 950, by = 50)

# Version of StringDB
# stringdb_versions <- c("11.5", "12.0")
stringdb_versions <- "12.0"

# Version of Bader pathways
# bader_versions <- c(2020, 2023)
bader_versions <- 2023

# Background set for enrichment analysis
background_set <- file.path(enrichment_dir, "sagittal_gene_table_normalized_filtered.csv")

params <- expand_grid(stringdb = stringdb_versions,
                      bader = bader_versions,
                      score = gene_scores)
for (i in 1:nrow(params)) {
  
  bader_version <- params[[i, "bader"]]
  stringdb_version <- params[[i, "stringdb"]]
  gene_score <- params[[i, "score"]]
  
  message(
    paste(paste("StringDB version:", stringdb_version),
          paste("Bader version:", bader_version),
          paste("Gene score:", gene_score),
          sep = "\n")
  )
  
  output_dir <- paste("StringDB", stringdb_version, "Bader", 
                      bader_version, sep = "_")
  output_dir <- file.path(enrichment_dir, output_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  if (bader_version == 2020) {
    bader_list <- file.path(enrichment_dir, "Human_Reactome_March_01_2020_symbol.gmt")
  } else if (bader_version == 2023) {
    bader_list <- file.path(enrichment_dir, "Human_Reactome_October_01_2023_symbol.gmt")  
  } else {
    stop()
  }
  
  GetGeneNeighbourhood(GeneScore = gene_score, 
                       total_clusters = nk_max, 
                       stringdb_version = stringdb_version,
                       output_dir = file.path(output_dir, "NeighbourhoodInfo"))
  
  GetNeighbourhoodEnrichment(GeneScore = gene_score,
                             total_clusters = nk_max,
                             Bader_List = bader_list,
                             background_set = background_set,
                             Neighbour_dir = file.path(output_dir, "NeighbourhoodInfo"),
                             output_dir = file.path(output_dir, "NeighbourhoodEnrichment"))
  
}
