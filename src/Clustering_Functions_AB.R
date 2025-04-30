suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(RMINC))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(SNFtool))



MakeModelList<-function(scanbase_scans_file="Data/Resources/scanbase_40um - Scans_22July19.csv",
                        scanbase_studies_file="Data/Resources/scanbase_40um - Studies_Feb2023.csv",
                        scanbase_genotypes_file="Data/Resources/scanbase_40um - Genotypes_Feb2023.csv",
                        scanbase_focus="Autism",
                        scanbase_sample_size_threshold=6,
                        base_directory="Data/Outputs/ModelInfo"){
  
  dir.create(base_directory, recursive = TRUE, showWarnings = FALSE)
  
  scans <- read.csv(scanbase_scans_file)
  studies <- read.csv(scanbase_studies_file)
  genotypes_file <- read.csv(scanbase_genotypes_file)
  
  # Use only autism studies
  studies %<>% filter(Focus==scanbase_focus)
  
  # Keep only scans in autism studies
  scans %<>% filter(Study_Name %in% studies$Study_Name)
  genotypes_file %<>% filter(Study_Name %in% studies$Study_Name)
  
  scans %<>% filter(Genotype_Code %in% genotypes_file$Genotype_Code)
  scans %<>% filter(Study_Name %in% genotypes_file$Study_Name)
  
  
  # Make genotype code a character
  scans %<>% mutate(Study_Name=as.character(Study_Name),
                    Genotype_Code=as.character(Genotype_Code))
  
  
  # Remove sex tags in scan genotype code
  # Why is this even here? There is a variable associated with sex
  scans %<>% mutate(Genotype_Code=(Genotype_Code %>%
                                     gsub("_M", "", .) %>%
                                     gsub("_F", "", .)))
  
  # Fix Pten_Kim study
  scans$Study_Name[which(scans$Study_Name=="PTEN_Kim" & endsWith(scans$Genotype_Code, "_131"))] <- "PTEN_Kim_131"
  scans$Study_Name[which(scans$Study_Name=="PTEN_Kim" & endsWith(scans$Genotype_Code, "_167"))] <- "PTEN_Kim_167"
  
  # Relabel genotypes for consistency
  wt_codes <- c("B6", "Wt", "WT_131", "WT_167", "Chd7++", "Chd7ff", "CON","Kctd13Wt_LatWt","WND","WNA","Snf2L_WT","Snf2Hcko-nestincre_WT")
  scans$Genotype_Code[which(scans$Genotype_Code %in% wt_codes)] <- "WT"
  
  het_codes <- c("HT", "HETero", "HT_131", "HT_167", "Het", "Hetero")
  scans$Genotype_Code[which(scans$Genotype_Code %in% het_codes)] <- "HET"
  
  hom_codes <- c("Hom", "HOMo","Homo")
  scans$Genotype_Code[which(scans$Genotype_Code %in% hom_codes)] <- "HOM"
  
  hem_codes <- c("Hemi","Hem","HMZ","Hmz")
  scans$Genotype_Code[which(scans$Genotype_Code %in% hem_codes)] <- "HEM"
  
  mut_codes <- c("Mut")
  scans$Genotype_Code[which(scans$Genotype_Code %in% mut_codes)] <- "MUT"
  
  # Sample size data per genotype
  table(scans$Study_Name, scans$Genotype_Code)
  
  # Remove studies with fewer than 8 wildtypes
  scans %<>% filter(Study_Name %in% names(which(table(scans$Study_Name, scans$Genotype_Code)[,"WT"] >= scanbase_sample_size_threshold)))
  
  #Removes scans with bad label files based on NA values in later steps.
  #scans<-scans[-c(1674,2402,2417,3097,3357,3387,3675),]
  
  write.csv(scans, file = glue("{base_directory}/scanbase_scans_autism_filtered_Feb2023.csv"))
  
  unique_studies <- unique(scans$Study_Name)
  model_list <- list()
  for (i in seq_along(unique_studies)) {
    study <- unique_studies[i]
    study_scans <- scans %>% filter(Study_Name==study)
    unique_genotypes <- setdiff(unique(study_scans$Genotype_Code), "WT")
    genotypes <- c()
    for (j in seq_along(unique_genotypes)) {
      genotype <- unique_genotypes[j]
      N <- nrow(study_scans %>% filter(Genotype_Code==genotype))
      if (N >= scanbase_sample_size_threshold) {
        genotypes <- c(genotypes, genotype)
      }
    }
    if (length(genotypes) >= 1) {
      model_list[[study]] <- list(wildtype="WT",
                                  genotypes=genotypes)
    }
  }
  
  save(list = "model_list", file = glue("{base_directory}/model_list_Feb2023.RData"))
  return(model_list)
}


## Making the Effect Size Images

MakeEffectSizeMaps<-function(jdtype,model_list,resolution="200",
                             dir_determinants="Data/Raw/jacobian_determinants",
                             output_phenotype_dir="Data/Outputs/EffectSizeMaps",
                             base_directory="Data/Outputs/ModelInfo",
                             boot="Y",
                             num_boot="100"){
  
  scans<-read.csv(glue("{base_directory}/scanbase_scans_autism_filtered_Feb2023.csv"))
  
  if(boot=="Y"){
    output_phenotype_dir<-glue("{output_phenotype_dir}/{resolution}/Boot")
  } else {
    output_phenotype_dir<-glue("{output_phenotype_dir}/{resolution}")
  }
  
  dir.create(output_phenotype_dir,recursive = TRUE,showWarnings = FALSE)
  mvol <- mincGetVolume(glue("{dir_determinants}/scanbase_second_level-nlin-3_mask_{resolution}um.mnc"))
  nvoxels <- length(which(mvol > 0.5))
  scans<-read.csv(file = glue("{base_directory}/scanbase_scans_autism_filtered_Feb2023.csv"))
  
  if(startsWith(tolower(jdtype),"a")){
    jddir <- glue("abs_{resolution}")
    jdtype <- "Absolute"
  }else{
    jddir <- glue("rel_{resolution}")
    jdtype <- "Relative"
  }
  
  effect_size_data_matrix <- matrix(nrow=length(as.character(unlist(sapply(model_list, "[[", 2)))), ncol=nvoxels)
  effect_size_data_matrix_rownames <- c()
  
  if (boot=="Y"){
    for (nboot in 1:num_boot) {
      for (model in names(model_list)) {
        genotypes <- model_list[[model]]$genotypes
        
        study_scans_wt <- scans %>% filter(Study_Name==model, Genotype_Code %in% c("WT"))
        wt_data_matrix <- matrix(nrow=nrow(study_scans_wt), ncol = nvoxels)
        
        for (i in 1:nrow(study_scans_wt)) {
          filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_wt$Mouse_ID)[i]}_{jddir}.mnc")
          wt_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
        }
        
        rowboot<-sample(nrow(wt_data_matrix),size=nrow(wt_data_matrix),replace=TRUE)
        wt_data_matrix<-wt_data_matrix[rowboot,]
        
        wt_mean <- colMeans(wt_data_matrix)
        # wt_sd <- colSd(wt_data_matrix)
        # Equivalent to: apply(wt_data_matrix, MARGIN=2, FUN=function(x) {mean(x)})
        wt_sd<-apply(wt_data_matrix, MARGIN=2, sd)
        
        for (genotype in genotypes) {
          
          print(paste("Working on model:", model, ",", "genotype:", genotype, ",", "num_boot:", nboot))
          
          # Make a dataframe with only wildtypes and genotype of interest
          study_scans_mut <- scans %>% filter(Study_Name==model, Genotype_Code %in% c(genotype))
          mut_data_matrix <- matrix(nrow=nrow(study_scans_mut), ncol = nvoxels)
          
          rowboot<-sample(nrow(mut_data_matrix),size=nrow(mut_data_matrix),replace=TRUE)
          mut_data_matrix<-mut_data_matrix[rowboot,]
          
          for (i in 1:nrow(study_scans_mut)) {
            filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_mut$Mouse_ID)[i]}_{jddir}.mnc")
            mut_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
          }
          
          mut_mean <- colMeans(mut_data_matrix)
          effect_size <- (mut_mean - wt_mean) / wt_sd
          
          filenameES <- glue("{output_phenotype_dir}/{as.character(model)}_{as.character(genotype)}_ES_{jdtype}_{resolution}_{nboot}.mnc")
          
          ovol <- mvol
          ovol[] <- 0
          ovol[mvol > 0.5] <- effect_size
          #ovol[is.infinite()]<-0
          ovol[is.na(ovol)]<-0
          mincWriteVolume(ovol, filenameES)
          
        }
      }
    } 
  } else {
    k <- 1
    for (model in names(model_list)) {
      genotypes <- model_list[[model]]$genotypes
      
      study_scans_wt <- scans %>% filter(Study_Name==model, Genotype_Code %in% c("WT"))
      wt_data_matrix <- matrix(nrow=nrow(study_scans_wt), ncol = nvoxels)
      
      for (i in 1:nrow(study_scans_wt)) {
        filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_wt$Mouse_ID)[i]}_{jddir}.mnc")
        wt_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
      }
      
      wt_mean <- colMeans(wt_data_matrix)
      # wt_sd <- colSd(wt_data_matrix)
      # Equivalent to: apply(wt_data_matrix, MARGIN=2, FUN=function(x) {mean(x)})
      wt_sd<-apply(wt_data_matrix, MARGIN=2, sd)
      
      for (genotype in genotypes) {
        
        message(paste("Working on model:", model, ",", "genotype:", genotype))
        
        # Make a dataframe with only wildtypes and genotype of interest
        study_scans_mut <- scans %>% filter(Study_Name==model, Genotype_Code %in% c(genotype))
        mut_data_matrix <- matrix(nrow=nrow(study_scans_mut), ncol = nvoxels)
        
        for (i in 1:nrow(study_scans_mut)) {
          filename <- glue("{dir_determinants}/{jddir}/{as.character(study_scans_mut$Mouse_ID)[i]}_{jddir}.mnc")
          mut_data_matrix[i,] <- mincGetVolume(filename)[mvol > 0.5]
        }
        
        mut_mean <- colMeans(mut_data_matrix)
        effect_size <- (mut_mean - wt_mean) / wt_sd
        
        filenameES <- glue("{output_phenotype_dir}/{as.character(model)}_{as.character(genotype)}_ES_{jdtype}_{resolution}.mnc")
        
        ovol <- mvol
        ovol[] <- 0
        ovol[mvol > 0.5] <- effect_size
        #ovol[is.infinite()]<-0
        ovol[is.na(ovol)]<-0
        mincWriteVolume(ovol, filenameES, clobber = TRUE)
        
        effect_size_data_matrix[k,] <- effect_size
        effect_size_data_matrix_rownames <- c(effect_size_data_matrix_rownames, paste(model, genotype, sep="_"))
        k <- k + 1
        
        
      }
    }
    
  }
  
  if(boot!="Y"){
    rownames(effect_size_data_matrix)<-effect_size_data_matrix_rownames
    save(file=glue("{output_phenotype_dir}/ES_Matricies_{jdtype}.RData"),effect_size_data_matrix)
    return(effect_size_data_matrix)
  }
}


SNFCombine<-function(Data1=effect_size_data_matrix_abs,Data2=effect_size_data_matrix_rel,
                     K=10,alpha=0.5,T=20,distfunc="cor",
                     output_dir="Data/Outputs/SNFMatrix"){
  
  dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)
  
  if(distfunc=="cor"){
    Dist1 = (1-cor(t(Data1)))
    Dist2 = (1-cor(t(Data2)))
  } else {
    Dist1 = (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2)
    Dist2 = (dist2(as.matrix(Data2),as.matrix(Data2)))^(1/2)
  }
  
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  W = SNF(list(W1,W2), K, T)
  
  save(file=glue("{output_dir}/WMatrix.RData"),W)
  return(W)
}

###########################################
#
# Bootstrapping SNF Matricies
#

SNFCombineBoot<-function(K=10,alpha=0.5,T=20,distfunc="cor",
                         resolution="200",
                         dir_determinants="Data/Raw/new_jacobbian_determinants",
                         output_dir="Data/Outputs/SNFMatrix_Paper/",
                         boot_effect_dir="Data/Outputs/EffectSizeMaps_Paper/",
                         num_boot="100"){
  
  output_dir<-glue("{output_dir}/{resolution}/Boot")
  boot_effect_dir<-glue("{boot_effect_dir}/{resolution}/Boot")
  dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)
  dir.create(boot_effect_dir,recursive = TRUE,showWarnings = FALSE)
  
  mvol <- mincGetVolume(glue("{dir_determinants}/scanbase_second_level-nlin-3_mask_{resolution}um.mnc"))
  nvoxels <- length(which(mvol > 0.5))
  
  Data1 <- effect_size_data_matrix_abs
  Data2 <- effect_size_data_matrix_rel
  
  rm(effect_size_data_matrix_abs,effect_size_data_matrix_rel)
  gc()
  
  for(nboot in 1:num_boot){
    
    for (i in 1:length(rownames(Data1))){
      filename_abs <- glue("{boot_effect_dir}/{rownames(Data1)[i]}_ES_Absolute_{resolution}_{nboot}.mnc")
      filename_rel <- glue("{boot_effect_dir}/{rownames(Data2)[i]}_ES_Relative_{resolution}_{nboot}.mnc")
      Data1[i,] <- mincGetVolume(filename_abs)[mvol > 0.5]
      Data2[i,] <- mincGetVolume(filename_rel)[mvol > 0.5]
    }
    
    if(distfunc=="cor"){
      Dist1 = (1-cor(t(Data1)))
      Dist2 = (1-cor(t(Data2)))
    } else {
      Dist1 = (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2)
      Dist2 = (dist2(as.matrix(Data2),as.matrix(Data2)))^(1/2)
    }
    
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    
    W = SNF(list(W1,W2), K, T)
    
    save(file=glue("{output_dir}/WMatrix_{nboot}.RData"),W)
  }
}



# Creating and Outputing Cluster Information

CreateClusters<-function(num_clusters="10",cluster_method="spectral",
                         output_dir="Data/Outputs/Clusters",
                         NameFile="Data/Raw/Names.csv"){
  
  dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)
  
  if(cluster_method == "spectral"){
    for (j in 2:num_clusters) {
      C=j
      group = spectralClustering(W,C)
      groupname <- paste0("group", as.character(j))
      assign(groupname,group)
      if (j == 2) {
        AllClusters<-data.frame(rownames(W),group)
        colnames(AllClusters)<-c("Model",groupname)
      }
      else {
        group<-data.frame(group)
        colnames(group)<-groupname
        AllClusters<-cbind(AllClusters,group)
      }
    }
  }
  NewNames<-read.csv(NameFile)
  AllClusters <- as_tibble(AllClusters)
  AllClusters$Model <- NewNames$NewName
  Clusters <- column_to_rownames(AllClusters, "Model")
  # Clusters<-AllClusters[,2:num_clusters]
  # rownames(Clusters) <-NewNames$NewName
  write.csv(Clusters,file=glue("{output_dir}/Clusters.csv"))
  return(Clusters)
  
}

BootClusters<-function(num_clusters="10",cluster_method="spectral",
                       boot_matrix_dir="Data/Outputs/SNFMatrix/200/Boot/",
                       output_dir="Data/Outputs/Clusters",
                       NameFile="Data/Raw/Names.csv",
                       num_boot="100"){
  
  output_dir<-glue("{output_dir}/Boot")
  dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)
  
  for (nboot in 1:num_boot){
    load(file=glue("{boot_matrix_dir}/WMatrix_{nboot}.RData"))
    
    if(cluster_method == "spectral"){
      for (j in 2:num_clusters) {
        C=j
        group = spectralClustering(W,C)
        groupname <- paste0("group", as.character(j))
        assign(groupname,group)
        if (j == 2) {
          AllClusters<-data.frame(rownames(W),group)
          colnames(AllClusters)<-c("Model",groupname)
        }
        else {
          group<-data.frame(group)
          colnames(group)<-groupname
          AllClusters<-cbind(AllClusters,group)
        }
      }
    }
    NewNames<-read.csv(NameFile)
    Clusters<-AllClusters[,2:num_clusters]
    rownames(Clusters) <-NewNames$NewName
    write.csv(Clusters,file=glue("{output_dir}/Clusters_{nboot}.csv"))
    
  }
}

# Outputting voxel maps from clusters
# Need to have effect size matricies loaded and the SNF W matrix

CreateClusterAnatomyMaps<-function(num_clusters="10",cluster_method="spectral",
                                   average_kind="median",
                                   volume_type="absolute",resolution="50",
                                   dir_determinants="Data/Raw/jacobian_determinants",
                                   output_dir="Data/Outputs/Clusters"){
  
  if(startsWith(tolower(volume_type),"a")){
    ES_Frame <- data.frame(effect_size_data_matrix_abs)
    volume_type <- "abs"
  }else{
    ES_Frame <- data.frame(effect_size_data_matrix_rel)
    volume_type <- "rel"
  }
  
  mvol <- mincGetVolume(glue("{dir_determinants}/scanbase_second_level-nlin-3_mask_{resolution}um.mnc"))
  nvoxels <- length(which(mvol > 0.5))
  
  for (j in 2:num_clusters) {
    C=j
    group = spectralClustering(W,C)
    groupname <- paste0("group", as.character(j))
    labels = spectralClustering(W, C)
    labelsname <- paste0("labels", as.character(j))
    assign(groupname,group)
    assign(labelsname,labels)
    
    Group <- matrix(nrow=C, ncol=nvoxels)
    
    for (i in 1:C) {
      
      if(average_kind=="mean"){
        Group[i,] <- apply(subset(ES_Frame,labels==i), MARGIN=2, mean)
      } else {
        Group[i,] <- apply(subset(ES_Frame,labels==i), MARGIN=2, median)
      }
      
      filenameES <- glue("{output_dir}/Group_{i}_Clusternum_{j}_ES_{volume_type}_{resolution}_{average_kind}.mnc")
      
      ovol <- mvol
      ovol[] <- 0
      ovol[mvol > 0.5] <- Group[i,]
      mincWriteVolume(ovol, filenameES)
    }
  }
}
  
# For Loop Plotting/Saving of Neuroanatomy Groups and Data.

PlotClusterAnatomyMaps<-function(num_clusters="10",
                                 average_kind="median",
                                 volume_type="absolute",
                                 resolution="200",
                                 dir_determinants="Data/Raw/jacobian_determinants",
                                 output_dir="Data/Outputs/Clusters",
                                 slice_dir="Axial"){
  
  anatVol<-mincGetVolume(glue("{dir_determinants}/scanbase_second_level-nlin-3_{resolution}um.mnc"))
  
  if(startsWith(tolower(volume_type),"a")){
    volume_type <- "abs"
  }
  else{
    volume_type <- "rel"
  }
  
  if(startsWith(tolower(slice_dir),"a")){
    dim<-3
    slice_dir<-"axial"
  }
  else if (startsWith(tolower(slice_dir),"c")){
    dim<-2
    slice_dir<-"coronal"
  }
  else{
    dim<-1
    slice_dir<-"sagittal"
  }
  
  for (num_clust in 2:num_clusters) {    
    
    if (resolution=="200") {
      
      plt <- sliceSeries(nrow =1,begin=24,end=26,dimension=dim) %>% 
        anatomy(mincArray(anatVol), low=300, high=2000)
      
    } else {
      
      plt <- sliceSeries(nrow =1,begin=90,end=150,dimension=dim) %>% 
        anatomy(mincArray(anatVol), low=300, high=2000)
      
    }
    
    for (clust in 1:num_clust) {
      
      ovol<-mincGetVolume(glue("{output_dir}/Group_{clust}_Clusternum_{num_clust}_ES_{volume_type}_{resolution}_{average_kind}.mnc"))
      
      plt <- plt %>% 
        sliceSeries() %>%
        anatomy() %>%
        overlay(mincArray(ovol),0.5,2,symmetric=TRUE) %>% 
        addtitle(glue("Group {clust} {volume_type}")) #%>% 
    }
    
    cairo_pdf(glue("{output_dir}/SliceSeries_{num_clust}_Clusters_{resolution}_{average_kind}_{volume_type}_{slice_dir}.pdf"), height=10,width=5*as.numeric(num_clusters))
    plt %>% legend("Effect Size") %>%
      draw()
    dev.off()
  }
}

##############################
# Bader Analysis 

CreateTermsList<-function(SFARI_List="Data/Raw/SFARI-Gene_genes_02-11-2020release_03-02-2020export-1.csv",
                          Bader_List="Data/Raw/Human_Reactome_March_01_2020_symbol.gmt",
                          TermNum="10"){
  
  SFARI<-read.csv(SFARI_List)
  gmt_file<-Bader_List
  MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")

  df <- MappingsGMT$MODULES2GENES %>%
    map_dfr(function(l) {
      tibble(count=length(intersect(l, SFARI$gene.symbol)))
    }) %>%
  mutate(name=names(MappingsGMT$MODULES2GENES))
  
  # Select terms
  df_filt <- df %>% filter(count >= TermNum) 
  terms <- as.character(df_filt$name)
  return(terms)
}
  

###########################
# Get Homologs

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


#####################
# Get and Filter the Data for the Clusters

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

#####################################
# One hop across all model genes ----
#####################################

GetGeneNeighbourhood<-function(GeneScore="950",
                               total_clusters="10",
                               output_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                               boot="Y",
                               nboot="1000"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  
  if (boot=="Y"){
    output_dir <- glue("{output_dir}/Boot/")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  all_df_list <- list()
  gx_data_list <- list()
  gscore <- GeneScore
  
  if (boot=="Y"){
    
    for (bootpick in 146:nboot){
      
      for (num_clusters in 2:total_clusters) {
        all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", num_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", num_clusters))
        
        # Remove duplicates where multiple genes are in the same cluster
        all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
        table(all_df$cluster)
        
        # Zero based index
        all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
        
        # Get STRINGids and merge to all_df
        #api_call <- paste("https://version-11-5.string-db.org/api/tsv/get_string_ids?identifiers=", paste(all_df$gene, collapse = "%0D"), "&species=10090&limit=1", sep="")
        api_call <- paste("https://string-db.org/api/tsv/get_string_ids?identifiers=", paste(all_df$gene, collapse = "%0D"), "&species=10090&limit=1", sep="")
        string_ids <- read_tsv(api_call) %>% dplyr::select(queryIndex, stringId, preferredName)
        all_df %<>% left_join(string_ids, by="queryIndex")
        
        # Remove organism ID
        all_df$fixed_id <- all_df$stringId %>% strsplit("\\.") %>% map_chr(`[`, 2)
        
        for (clust in 1:num_clusters){
          all_df_cl <- all_df %>% filter(cluster==clust)
          rnum <- nrow(all_df_cl)
          all_df_cl <- sample_n(all_df_cl,rnum,replace=TRUE)
          
          if(clust==1){
            all_df_cl_new <- all_df_cl
          } else {
            all_df_cl_new <- rbind(all_df_cl_new,all_df_cl)
          }
        }
        all_df<-all_df_cl_new
        
        write.csv(file = glue("{output_dir}/all_gene_data_{bootpick}.csv"), x = all_df)
        
        ##########################################
        # This is where you can modify the gene score for the interactions - Normally set to 900.  Other common ones are 400 and 700
        
        # Get interaction data
        #api_call <- paste("https://version-11-5.string-db.org/api/tsv/interaction_partners?identifiers=", paste(all_df$stringId, collapse = "%0D"), "&species=10090&required_score=",gscore, sep="")
        api_call <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=", paste(all_df$stringId, collapse = "%0D"), "&species=10090&required_score=",gscore, sep="")
        gx_data <- read_tsv(api_call)
        gx_data %<>% left_join(all_df %>% dplyr::select(cluster, stringId), by=c("stringId_A"="stringId")) %>% dplyr::rename(cluster_A=cluster)
        
        # Write out gene neighbourhoods for each cluster
        base_set <- unique(c(gx_data$preferredName_A, gx_data$preferredName_B))
        write.table(base_set, file=glue("{output_dir}/base_set_full_neighbourhood_{GeneScore}_{bootpick}.txt"), quote = F, append=F, row.names = F, col.names = F)
        
        for (clust in 1:num_clusters) {
          target_set <- gx_data %>% filter(cluster_A==clust) %>% dplyr::select(preferredName_A, preferredName_B) %>% gather %>% .$value %>% unique 
          target_outfile <- glue("{output_dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood_{bootpick}.txt")
          write.table(target_set, file=target_outfile, quote = F, append=F, row.names = F, col.names = F)
          
          target_set_only_models <- all_df %>% filter(cluster==clust) %>% .$gene
          target_outfile_only_models <- glue("{output_dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_only_models_{bootpick}.txt")
          write.table(target_set_only_models, file=target_outfile_only_models, quote = F, append=F, row.names = F, col.names = F)
        }
        
        for (clust in 1:num_clusters) {
          target_set <- gx_data %>% filter(cluster_A!=clust) %>% dplyr::select(preferredName_A, preferredName_B) %>% gather %>% .$value %>% unique 
          target_outfile <- glue("{output_dir}/target_set_cluster{num_clusters}_Not_{clust}_{GeneScore}_neighbourhood_{bootpick}.txt")
          write.table(target_set, file=target_outfile, quote = F, append=F, row.names = F, col.names = F)
          
          target_set_only_models <- all_df %>% filter(cluster!=clust) %>% .$gene
          target_outfile_only_models <- glue("{output_dir}/target_set_cluster{num_clusters}_Not_{clust}_{GeneScore}_only_models_{bootpick}.txt")
          write.table(target_set_only_models, file=target_outfile_only_models, quote = F, append=F, row.names = F, col.names = F)
        }
        
        all_df_list[[paste0("Cluster", num_clusters)]] <- all_df
        gx_data_list[[paste0("Cluster", num_clusters)]] <- gx_data
      }
      
      listfname<-glue("{output_dir}/gene_interaction_tables_for_all_clusters_{GeneScore}_{bootpick}.RData")
      save(list = c("all_df_list", "gx_data_list"), file = listfname)
    }
    
  } else {
    
    for (num_clusters in 2:total_clusters) {
      all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", num_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", num_clusters))
      
      # Remove duplicates where multiple genes are in the same cluster
      all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
      table(all_df$cluster)
      
      # Zero based index
      all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
      
      # Get STRINGids and merge to all_df
      #api_call <- paste("https://version-11-5.string-db.org/api/tsv/get_string_ids?identifiers=", paste(all_df$gene, collapse = "%0D"), "&species=10090&limit=1", sep="")
      api_call <- paste("https://string-db.org/api/tsv/get_string_ids?identifiers=", paste(all_df$gene, collapse = "%0D"), "&species=10090&limit=1", sep="")
      string_ids <- read_tsv(api_call) %>% dplyr::select(queryIndex, stringId, preferredName)
      all_df %<>% left_join(string_ids, by="queryIndex")
      
      # Remove organism ID
      all_df$fixed_id <- all_df$stringId %>% strsplit("\\.") %>% map_chr(`[`, 2)
      
      write.csv(file = glue("{output_dir}/all_gene_data.csv"), x = all_df)
      
      ##########################################
      # This is where you can modify the gene score for the interactions - Normally set to 900.  Other common ones are 400 and 700
      
      # Get interaction data
      #api_call <- paste("https://version-11-5.string-db.org/api/tsv/interaction_partners?identifiers=", paste(all_df$stringId, collapse = "%0D"), "&species=10090&required_score=",gscore, sep="")
      api_call <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=", paste(all_df$stringId, collapse = "%0D"), "&species=10090&required_score=",gscore, sep="")
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
}

###########################
# Enrichment tests
# First get different base sets


GetNeighbourhoodEnrichment<-function(GeneScore="950",
                                     total_clusters="10",
                                     Bader_List="Data/Raw/Human_Reactome_October_01_2023_symbol.gmt",
                                     Neighbour_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                                     output_dir="Data/Outputs/NeighbourhoodEnrichment_Paper/",
                                     boot="Y",
                                     nboot="1000"){
  
  if (boot=="Y"){
    
    
    output_dir <- glue("{output_dir}/{GeneScore}/Boot/")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", total_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", total_clusters))
    all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
    all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
    
    for (bootpick in 1:nboot){
      
      # First get different base sets
      
      base_brain <- read.csv("Data/Raw/sagittal_gene_table_normalized_filtered.csv", header = T)$msg.genes.acronym %>% as.character
      base_models_neighbourhood <- read.csv(glue("{Neighbour_dir}/{GeneScore}/Boot/base_set_full_neighbourhood_{GeneScore}_{bootpick}.txt"), header = F)$V1  %>% as.character
      base_models_only <- all_df$gene
      base <- list(brain=base_brain, models_neighbourhood=base_models_neighbourhood, models_only=base_models_only)
      gmt_file<-Bader_List
      MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")
      
      
      for (num_clusters in 2:total_clusters) {
        
        target <- list(models_neighbourhood=list(), models_only=list())
        for (clust in 1:num_clusters) {
          
          if (file.size(glue("{Neighbour_dir}/{GeneScore}/Boot/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood_{bootpick}.txt")) == 0) {
            next
          } else {
            target_set_cl <- read.csv(glue("{Neighbour_dir}/{GeneScore}/Boot/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood_{bootpick}.txt"), header = F)$V1 %>% as.character
            target$models_neighbourhood[[clust]] <- target_set_cl
            
            target_set_cl <- read.csv(glue("{Neighbour_dir}/{GeneScore}/Boot/target_set_cluster{num_clusters}_{clust}_{GeneScore}_only_models_{bootpick}.txt"), header = F)$V1 %>% as.character
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
            
            write.csv(tmodObj, file = glue("{output_dir}/NewBader_enrichment_clusterneighbourhood_vs_brain_all_{num_clusters}_{clust}_{GeneScore}_{bootpick}.csv"), row.names = F)
          }
        }
      }
    }
    
  } else {
    
    output_dir <- glue("{output_dir}/{GeneScore}")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", total_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", total_clusters))
    all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
    all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
    
    # First get different base sets
    
    base_brain <- read.csv("Data/Raw/sagittal_gene_table_normalized_filtered.csv", header = T)$msg.genes.acronym %>% as.character
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
}

###########################
# Enrichment tests
# First get different base sets
GetSingleGeneNeighbourhood<-function(GeneScore="400",
                                     total_clusters="10",
                                     output_dir="Data/Outputs/SingleGeneNeighbourhood_150/"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  all_df_list <- list()
  gx_data_list <- list()
  gscore <- GeneScore
  
  all_df <- cl_filt %>% dplyr::select(Model, Gene) %>% dplyr::rename(gene=Gene)
  
  # Remove duplicates where multiple genes are in the same cluster
  all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, sep = "-")) %>% .$dup_select),]
  
  # Zero based index
  all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
  
  # Get STRINGids and merge to all_df
  api_call <- paste("https://string-db.org/api/tsv/get_string_ids?identifiers=", paste(all_df$gene, collapse = "%0D"), "&species=10090&limit=1", sep="")
  string_ids <- read_tsv(api_call) %>% dplyr::select(queryIndex, stringId, preferredName)
  all_df %<>% left_join(string_ids, by="queryIndex")
  
  # Remove organism ID
  all_df$fixed_id <- all_df$stringId %>% strsplit("\\.") %>% map_chr(`[`, 2)
  write.csv(file = glue("{output_dir}/all_gene_data.csv"), x = all_df)
  
  ##########################################
  # This is where you can modify the gene score for the interactions - Normally set to 900.  Other common ones are 400 and 700
  
  # Get interaction data
  api_call <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=", paste(all_df$stringId, collapse = "%0D"), "&species=10090&required_score=",gscore, sep="")
  gx_data <- read_tsv(api_call)
  gx_data %<>% left_join(all_df %>% dplyr::select(gene,stringId), by=c("stringId_A"="stringId")) 
  
  # Write out gene neighbourhoods for each cluster
  base_set <- unique(c(gx_data$preferredName_A, gx_data$preferredName_B))
  write.table(base_set, file=glue("{output_dir}/base_set_full_neighbourhood_{GeneScore}.txt"), quote = F, append=F, row.names = F, col.names = F)
  
  for (sgene in 1:nrow(all_df)) {
    target_set <- gx_data %>% filter(gene==all_df$gene[sgene]) %>% dplyr::select(preferredName_A, preferredName_B) %>% gather %>% .$value %>% unique 
    target_outfile <- glue("{output_dir}/target_set_cluster_{all_df$gene[sgene]}_{GeneScore}_neighbourhood.txt")
    write.table(target_set, file=target_outfile, quote = F, append=F, row.names = F, col.names = F)
  }
}


GetSingleGeneNeighbourhoodEnrichment<-function(GeneScore="400",
                                               total_clusters="10",
                                               Bader_List="Data/Raw/Human_Reactome_March_01_2020_symbol.gmt",
                                               Neighbour_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                                               output_dir="Data/Outputs/NeighbourhoodEnrichment/SingleGene/"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", total_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", total_clusters))
  all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
  all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
  
  # First get different base sets
  
  base_brain <- read.csv("Data/Raw/sagittal_gene_table_normalized_filtered.csv", header = T)$msg.genes.acronym %>% as.character
  base_models_neighbourhood <- read.csv(glue("{Neighbour_dir}/{GeneScore}/base_set_full_neighbourhood_{GeneScore}.txt"), header = F)$V1  %>% as.character
  base_models_only <- all_df$gene
  base <- list(brain=base_brain, models_neighbourhood=base_models_neighbourhood, models_only=base_models_only)
  gmt_file<-Bader_List
  MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")
  
  for (SingleGene in 1:nrow(all_df)){
    base_set <- base$brain
    sgene<-all_df$gene[SingleGene]
    
    if(file.size(file=glue("{Neighbour_dir}/{GeneScore}/target_set_cluster_{sgene}_{GeneScore}_neighbourhood.txt"))==0){
      next
    }
    
    target_set <- read.csv(glue("{Neighbour_dir}/{GeneScore}/target_set_cluster_{sgene}_{GeneScore}_neighbourhood.txt"), header = F)$V1 %>% as.character
    
    enrichment_object <-  MappingsGMT
    fg <- get_homologs(as.character(target_set), "mouse", ordered = F)$human_genes
    bg <- get_homologs(as.character(base_set), "mouse", ordered = F)$human_genes
    
    if(any( fg %in% (enrichment_object$GENES$ID %>% as.character))==FALSE){
      next
    }
    
    tmodObj <- tmodHGtest(fg, bg, qval = 1.1, mset = enrichment_object, filter = F, order.by = "pval")
    tmodObj %<>% rowwise() %>% mutate(P.Value = P.Value) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
    tmodObj %<>% arrange(P.Value)
    tmodObj$rank <- 1:nrow(tmodObj)
    tmodObj %<>% mutate(NLQ = -log(adj.P.Val))
    #tmodObj %<>% mutate(cluster = clust)
    tmodObj %<>% dplyr::select(rank, Title, NLQ, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
    
    write.csv(tmodObj, file = glue("{output_dir}/NewBader_enrichment_singlegeneneighbourhood_vs_brain_all_{sgene}_{GeneScore}.csv"), row.names = F)
    
    
  }
}


##################################################
# Get enrichment for not clusters
##################################################

GetNeighbourhoodEnrichment_NotCluster<-function(GeneScore="900",
                                     total_clusters="10",
                                     Bader_List="Data/Raw/Human_Reactome_March_01_2020_symbol.gmt",
                                     Neighbour_dir="Data/Outputs/NeighbourhoodInfo/",
                                     output_dir="Data/Outputs/NeighbourhoodEnrichment/"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", total_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", total_clusters))
  all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
  all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
  
  # First get different base sets
  
  base_brain <- read.csv("Data/Raw/sagittal_gene_table_normalized_filtered.csv", header = T)$msg.genes.acronym %>% as.character
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
      
      for (i in 1:num_clusters) {
        if (i==clust) next 
        if (i == 1){
         target_set<-target$models_neighbourhood[[i]]
        } else {
          target_set <- c(target_set,target$models_neighbourhood[[i]])
        }
      }
      
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
        
        write.csv(tmodObj, file = glue("{output_dir}/NewBader_enrichment_clusterneighbourhood_vs_brain_all_{num_clusters}_Not_{clust}_{GeneScore}.csv"), row.names = F)
      }
    }
  }
}



ClusterPairwiseComp<-function(GeneScore=gscore,
                              total_clusters="10",
                              Bader_List="Data/Raw/Human_Reactome_March_01_2020_symbol.gmt",
                              TargetSet_Dir="Data/Outputs/NeighbourhoodInfo",
                              Enrichment_Dir="Data/Outputs/NeighbourhoodEnrichment",
                              output_dir="Data/Outputs/ClusterPairwiseComp"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  Enrichment_Dir <- glue("{Enrichment_Dir}/{GeneScore}")
  TargetSet_Dir <- glue("{TargetSet_Dir}/{GeneScore}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", total_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", total_clusters))
  all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
  all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
  
  # First get different base sets
  
  base_brain <- read.csv("Data/Raw/sagittal_gene_table_normalized_filtered.csv", header = T)$msg.genes.acronym %>% as.character
  base_models_neighbourhood <- read.csv(glue("{TargetSet_Dir}/base_set_full_neighbourhood_{GeneScore}.txt"), header = F)$V1  %>% as.character
  base_models_only <- all_df$gene
  base <- list(brain=base_brain, models_neighbourhood=base_models_neighbourhood, models_only=base_models_only)
  base_set <- base$brain
  
  gmt_file<-Bader_List
  MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")
  
  
  SFARI<-read.csv("Data/Raw/SFARI-Gene_genes_02-11-2020release_03-02-2020export-1.csv")
  gmt_file <- Bader_List
  MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")
  
  terms <- TermsList
  
  for (num_clusters in 2:total_clusters) {
    print(glue("Working on cluster: {num_clusters}"))
    
    result_all_clusters <- list()
    target_set <- list()
    for (i in 1:num_clusters) {
      clust <- i
      if (file.size(glue("{TargetSet_Dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood.txt")) == 0) {
        next
      } else {
        result_all_clusters[[i]] <- read.csv(file = glue("{Enrichment_Dir}/NewBader_enrichment_clusterneighbourhood_vs_brain_all_{num_clusters}_{clust}_{GeneScore}.csv"))
        target_set[[clust]] <- read.csv(file = glue("{TargetSet_Dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood.txt"), header = F)$V1
      }
    }
    result_all_clusters <- do.call("rbind", result_all_clusters)
    ok_clust<-target_set %>% map_lgl(function(x) {!is.null(x)}) %>% which
    
    result_individual_clusters <- NULL
    for (clust in ok_clust) {
      target_set_cl <- target_set[[clust]]
      
      enrichment_object <- MappingsGMT
      fg <- get_homologs(as.character(target_set_cl), "mouse", ordered = F)$human_genes
      bg <- get_homologs(as.character(base_set), "mouse", ordered = F)$human_genes
      
      tmodObj <- tmodHGtest(fg, bg, qval = 1.1, mset = enrichment_object, filter = F, order.by = "pval")
      tmodObj %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
      tmodObj %<>% arrange(P.Value)
      tmodObj$rank <- 1:nrow(tmodObj)
      tmodObj %<>% mutate(NLQ = -log(adj.P.Val))
      tmodObj %<>% mutate(cluster = clust)
      tmodObj %<>% dplyr::select(cluster, rank, Title, NLQ, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
      
      if (exists("result_individual_clusters")) {
        result_individual_clusters <- rbind(result_individual_clusters, tmodObj)
      } else {
        result_individual_clusters <- tmodObj
      }
      
    }
    
    results_pairwise_tests <- NULL
    for (k in 1:length(terms)) {
      out_term <- result_individual_clusters %>% filter(ID==terms[k])  
      out_term_df <- NULL
      
      for (ip in 2:length(ok_clust)) {
        for (jp in 1:(i-1)) {
          i<-ok_clust[ip]
          j<-ok_clust[jp]
          bi <- out_term$b[which(out_term$cluster==i)] # Number of target genes in term
          bj <- out_term$b[which(out_term$cluster==j)] # Number of target genes in term
          Bi <- out_term$B[which(out_term$cluster==i)] # Number of genes associated with term
          Bj <- out_term$B[which(out_term$cluster==j)] # Number of genes associated with term
          ni <- out_term$n[which(out_term$cluster==i)] # Number of target genes
          nj <- out_term$n[which(out_term$cluster==j)] # Number of target genes
          Ni <- out_term$N[which(out_term$cluster==i)] # Number of base genes
          Nj <- out_term$N[which(out_term$cluster==j)] # Number of base genes
          ctable <- matrix(c(bi, bj, Bi-bi, Bj-bj), nrow=2, dimnames = list(paste("Cluster", c(i,j)), c("b", "B-b")))
          #ctable <- matrix(c(bi/ni, bj/nj, (Bi/Ni)-(bi/ni), (Bj/Nj)-(bj/nj)), nrow=2, dimnames = list(paste("Cluster", c(i,j)), c("Proportion of input genes in module", "Proportion of GO genes in base")))
          
          # Deal with zero genes in both sets
          if (any(c(any(rowSums(ctable)==0), any(colSums(ctable)==0)))) {
            out_term_df_ij <- data.frame(term=terms[k], 
                                         Title=unique(out_term$Title),
                                         cluster_i=i, 
                                         cluster_j=j, 
                                         test=paste0("Prob[proportion(Cluster ", i,") > proportion(Cluster ", j, ")]"), 
                                         estimate=0,
                                         P.Value=1,
                                         bi=bi,
                                         Bi=Bi,
                                         ni=ni,
                                         Ni=Ni,
                                         bj=bj,
                                         Bj=Bj,
                                         nj=nj,
                                         Nj=Nj,
                                         proportion_i=(bi/Bi),
                                         proportion_j=(bj/Bj))
            #proportion_i=((bi/ni)/(Bi/Ni)),
            #proportion_j=((bj/nj)/(Bj/Nj)))
          } else {
            ctest <- exact.test(data=ctable, method="z-pooled", model="Binomial", cond.row = TRUE, alternative = "two.sided", to.plot = F, ref.pvalue = TRUE) # Use two-sided tests so we can compare i>j and j>i
            ctest$statistic
            ctest$estimate
            ctest$p.value
            ctest$parameter
            out_term_df_ij <- data.frame(term=terms[k], 
                                         Title=unique(out_term$Title),
                                         cluster_i=i, 
                                         cluster_j=j, 
                                         test=paste0("Prob[proportion(Cluster ", i,") > proportion(Cluster ", j, ")]"), 
                                         estimate=ctest$estimate,
                                         P.Value=ctest$p.value,
                                         bi=bi,
                                         Bi=Bi,
                                         ni=ni,
                                         Ni=Ni,
                                         bj=bj,
                                         Bj=Bj,
                                         nj=nj,
                                         Nj=Nj,
                                         proportion_i=(bi/Bi),
                                         proportion_j=(bj/Bj))
            #proportion_i=((bi/ni)/(Bi/Ni)),
            #proportion_j=((bj/nj)/(Bj/Nj)))
          }
          
          if (exists("out_term_df")) {
            out_term_df <- rbind(out_term_df, out_term_df_ij)
          } else {
            out_term_df <- out_term_df_ij
          }
        }
      }
      out_term_df %<>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr"))
      rownames(out_term_df) <- NULL
      if (exists("results_pairwise_tests")) {
        results_pairwise_tests <- rbind(results_pairwise_tests, out_term_df)
      } else {
        results_pairwise_tests <- out_term_df
      }
      
    }
    
    # Readjust p value
    
    # Subset for terms with differences
    #terms_with_differences <- results_pairwise_tests %>% filter(adj.P.Val <= 0.05) %>% .$term %>% as.character %>% unique
    #results_pairwise_tests_differences <- results_pairwise_tests %>% filter(term %in% terms_with_differences)
    #results_pairwise_tests_differences %<>% arrange(adj.P.Val)
    results_pairwise_tests_differences <- results_pairwise_tests %>% mutate(adj.P.Val=p.adjust(P.Value, method = "fdr"))
    results_pairwise_tests_differences <- results_pairwise_tests_differences %>% filter(adj.P.Val <= 0.05)
    results_pairwise_tests_differences %<>% arrange(adj.P.Val)
    
    # Number of times a cluster shows up -- validate that one cluster isn't overrepresented, potentially because of larger cluster size
    table(c(results_pairwise_tests_differences$cluster_i, results_pairwise_tests_differences$cluster_j))
    
    # Venn
    results_pairwise_tests_clustered <- results_pairwise_tests %>% filter(adj.P.Val <= 0.05)
    results_pairwise_tests_clustered %<>% mutate(A=ifelse(estimate > 0, cluster_i, cluster_j), B=ifelse(estimate < 0, cluster_i, cluster_j), estimate=abs(estimate)) %>% dplyr::select(term, Title, A, B, estimate, P.Value, adj.P.Val)
    write.csv(x=results_pairwise_tests_clustered, file=glue("{output_dir}/NewBader_pairwise_tests_significant_ordered_{num_clusters}.csv"))
    
    results_pairwise_tests_clustered_counting <- results_pairwise_tests_clustered
    
    A_count_table <- table(results_pairwise_tests_clustered_counting$Title, results_pairwise_tests_clustered_counting$A)
    ok_title_A <- which(A_count_table == (length(ok_clust) - 1), arr.ind = T)    
    A_counts <- paste(rownames(A_count_table)[ok_title_A[,"row"]], colnames(A_count_table)[ok_title_A[,"col"]], sep="_")
    results_pairwise_tests_clustered_majordiff <- results_pairwise_tests_clustered %>% mutate(A_count=paste(Title, A, sep="_")) %>% filter(A_count %in% A_counts) %>% dplyr::select(-A_count)
    write.csv(x=results_pairwise_tests_clustered_majordiff, file=glue("{output_dir}/NewBader_pairwise_tests_significant_ordered_{num_clusters}_majordiff.csv"))
  }
}

#########################################
#Pairwise Comparisons against the whole rest of the group

ClusterPairwiseCompvsALL<-function(GeneScore=gscore,
                              total_clusters="10",
                              Bader_List="Data/Raw/Human_Reactome_March_01_2020_symbol.gmt",
                              TargetSet_Dir="Data/Outputs/NeighbourhoodInfo",
                              Enrichment_Dir="Data/Outputs/NeighbourhoodEnrichment",
                              output_dir="Data/Outputs/ClusterPairwiseCompvsALL"){
  
  output_dir <- glue("{output_dir}/{GeneScore}")
  Enrichment_Dir <- glue("{Enrichment_Dir}/{GeneScore}")
  TargetSet_Dir <- glue("{TargetSet_Dir}/{GeneScore}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_df <- cl_filt %>% dplyr::select(Model, Gene, paste0("group", total_clusters)) %>% dplyr::rename(gene=Gene, cluster=paste0("group", total_clusters))
  all_df <- all_df[!duplicated(all_df %>% mutate(dup_select=paste(gene, cluster, sep = "-")) %>% .$dup_select),]
  all_df$queryIndex <- seq(from=0, to=(dim(all_df)[1]-1), by=1)
  
  # First get different base sets
  
  base_brain <- read.csv("Data/Raw/sagittal_gene_table_normalized_filtered.csv", header = T)$msg.genes.acronym %>% as.character
  base_models_neighbourhood <- read.csv(glue("{TargetSet_Dir}/base_set_full_neighbourhood_{GeneScore}.txt"), header = F)$V1  %>% as.character
  base_models_only <- all_df$gene
  base <- list(brain=base_brain, models_neighbourhood=base_models_neighbourhood, models_only=base_models_only)
  base_set <- base$brain
  
  gmt_file<-Bader_List
  MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")
  
  
  SFARI<-read.csv("Data/Raw/SFARI-Gene_genes_02-11-2020release_03-02-2020export-1.csv")
  gmt_file <- Bader_List
  MappingsGMT <- tmodImportMSigDB(gmt_file, format = "gmt")
  
  terms <- TermsList
  
  for (num_clusters in 2:total_clusters) {
    print(glue("Working on cluster: {num_clusters}"))
    
    result_all_clusters <- list()
    target_set <- list()
    target_set_not <- list()
    result_not_clusters <- list()
    for (i in 1:num_clusters) {
      clust <- i
      if (file.size(glue("{TargetSet_Dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood.txt")) == 0) {
        next
      } else {
        result_all_clusters[[i]] <- read.csv(file = glue("{Enrichment_Dir}/NewBader_enrichment_clusterneighbourhood_vs_brain_all_{num_clusters}_{clust}_{GeneScore}.csv"))
        result_not_clusters[[i]] <- read.csv(file = glue("{Enrichment_Dir}/NewBader_enrichment_clusterneighbourhood_vs_brain_all_{num_clusters}_Not_{clust}_{GeneScore}.csv"))
        target_set[[clust]] <- read.csv(file = glue("{TargetSet_Dir}/target_set_cluster{num_clusters}_{clust}_{GeneScore}_neighbourhood.txt"), header = F)$V1
        target_set_not[[clust]] <- read.csv(file = glue("{TargetSet_Dir}/target_set_cluster{num_clusters}_Not_{clust}_{GeneScore}_neighbourhood.txt"), header = F)$V1
      }
    }
    result_all_clusters <- do.call("rbind", result_all_clusters)
    result_not_clusters <- do.call("rbind", result_not_clusters)
    ok_clust<-target_set %>% map_lgl(function(x) {!is.null(x)}) %>% which
    
    result_individual_clusters <- NULL
    result_individual_not_clusters <- NULL

    for (clust in ok_clust) {
      target_set_cl <- target_set[[clust]]
      
      enrichment_object <- MappingsGMT
      fg <- get_homologs(as.character(target_set_cl), "mouse", ordered = F)$human_genes
      bg <- get_homologs(as.character(base_set), "mouse", ordered = F)$human_genes
      
      tmodObj <- tmodHGtest(fg, bg, qval = 1.1, mset = enrichment_object, filter = F, order.by = "pval")
      tmodObj %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
      tmodObj %<>% arrange(P.Value)
      tmodObj$rank <- 1:nrow(tmodObj)
      tmodObj %<>% mutate(NLQ = -log(adj.P.Val))
      tmodObj %<>% mutate(cluster = clust)
      tmodObj %<>% dplyr::select(cluster, rank, Title, NLQ, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
      
      if (exists("result_individual_clusters")) {
        result_individual_clusters <- rbind(result_individual_clusters, tmodObj)
      } else {
        result_individual_clusters <- tmodObj
      }
      
      target_set_not_cl <- target_set_not[[clust]]
      
      enrichment_object <- MappingsGMT
      fg <- get_homologs(as.character(target_set_not_cl), "mouse", ordered = F)$human_genes
      bg <- get_homologs(as.character(base_set), "mouse", ordered = F)$human_genes
      
      tmodObj <- tmodHGtest(fg, bg, qval = 1.1, mset = enrichment_object, filter = F, order.by = "pval")
      tmodObj %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) #tmod runs one-sided tests
      tmodObj %<>% arrange(P.Value)
      tmodObj$rank <- 1:nrow(tmodObj)
      tmodObj %<>% mutate(NLQ = -log(adj.P.Val))
      tmodObj %<>% mutate(cluster = clust)
      tmodObj %<>% dplyr::select(cluster, rank, Title, NLQ, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
      
      if (exists("result_individual_not_clusters")) {
        result_individual_not_clusters <- rbind(result_individual_not_clusters, tmodObj)
      } else {
        result_individual_not_clusters <- tmodObj
      }
      
  }
    
    results_pairwise_tests <- NULL
    for (k in 1:length(terms)) {
      out_term <- result_individual_clusters %>% filter(ID==terms[k])
      out_term_not <- result_individual_not_clusters %>% filter(ID==terms[k])
      out_term_df <- NULL
      
      for (ip in 1:length(ok_clust)) {
          i<-ok_clust[ip]
          bi <- out_term$b[which(out_term$cluster==i)] # Number of target genes in term
          bj <- out_term_not$b[which(out_term_not$cluster==i)] # Number of target genes in term
          Bi <- out_term$B[which(out_term$cluster==i)] # Number of genes associated with term
          Bj <- out_term_not$B[which(out_term_not$cluster==i)] # Number of genes associated with term
          ni <- out_term$n[which(out_term$cluster==i)] # Number of target genes
          nj <- out_term_not$n[which(out_term_not$cluster==i)] # Number of target genes
          Ni <- out_term$N[which(out_term$cluster==i)] # Number of base genes
          Nj <- out_term_not$N[which(out_term_not$cluster==i)] # Number of base genes
          ctable <- matrix(c(bi, bj, Bi-bi, Bj-bj), nrow=2, dimnames = list(c(paste("Cluster", i), paste("Not Cluster",i)), c("b", "B-b")))
          #ctable <- matrix(c(bi/ni, bj/nj, (Bi/Ni)-(bi/ni), (Bj/Nj)-(bj/nj)), nrow=2, dimnames = list(paste("Cluster", c(i,j)), c("Proportion of input genes in module", "Proportion of GO genes in base")))
          
          # Deal with zero genes in both sets
          if (any(c(any(rowSums(ctable)==0), any(colSums(ctable)==0)))) {
            out_term_df_ij <- data.frame(term=terms[k], 
                                         Title=unique(out_term$Title),
                                         cluster_i=i,  
                                         test=paste0("Prob[proportion(Cluster ", i,") > proportion(Not Cluster ", i, ")]"), 
                                         estimate=0,
                                         P.Value=1,
                                         bi=bi,
                                         Bi=Bi,
                                         ni=ni,
                                         Ni=Ni,
                                         bj=bj,
                                         Bj=Bj,
                                         nj=nj,
                                         Nj=Nj,
                                         proportion_i=(bi/Bi),
                                         proportion_j=(bj/Bj))
            #proportion_i=((bi/ni)/(Bi/Ni)),
            #proportion_j=((bj/nj)/(Bj/Nj)))
          } else {
            ctest <- exact.test(data=ctable, method="z-pooled", model="Binomial", cond.row = TRUE, alternative = "two.sided", to.plot = F, ref.pvalue = TRUE) # Use two-sided tests so we can compare i>j and j>i
            ctest$statistic
            ctest$estimate
            ctest$p.value
            ctest$parameter
            out_term_df_ij <- data.frame(term=terms[k], 
                                         Title=unique(out_term$Title),
                                         cluster_i=i,  
                                         test=paste0("Prob[proportion(Cluster ", i,") > proportion(Not Cluster ", i, ")]"), 
                                         estimate=ctest$estimate,
                                         P.Value=ctest$p.value,
                                         bi=bi,
                                         Bi=Bi,
                                         ni=ni,
                                         Ni=Ni,
                                         bj=bj,
                                         Bj=Bj,
                                         nj=nj,
                                         Nj=Nj,
                                         proportion_i=(bi/Bi),
                                         proportion_j=(bj/Bj))
            #proportion_i=((bi/ni)/(Bi/Ni)),
            #proportion_j=((bj/nj)/(Bj/Nj)))
          }
          
          if (exists("out_term_df")) {
            out_term_df <- rbind(out_term_df, out_term_df_ij)
          } else {
            out_term_df <- out_term_df_ij
          }
        
      }
      out_term_df %<>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr"))
      rownames(out_term_df) <- NULL
      if (exists("results_pairwise_tests")) {
        results_pairwise_tests <- rbind(results_pairwise_tests, out_term_df)
      } else {
        results_pairwise_tests <- out_term_df
      }
      
    }
    
    # Readjust p value
    
    # Subset for terms with differences
    #terms_with_differences <- results_pairwise_tests %>% filter(adj.P.Val <= 0.05) %>% .$term %>% as.character %>% unique
    #results_pairwise_tests_differences <- results_pairwise_tests %>% filter(term %in% terms_with_differences)
    #results_pairwise_tests_differences %<>% arrange(adj.P.Val)
    results_pairwise_tests_differences <- results_pairwise_tests %>% mutate(adj.P.Val=p.adjust(P.Value, method = "fdr"))
    results_pairwise_tests_differences <- results_pairwise_tests_differences %>% filter(adj.P.Val <= 0.05) %>% filter(estimate > 0)
    results_pairwise_tests_differences %<>% arrange(adj.P.Val)
    

    write.csv(x=results_pairwise_tests_differences, file=glue("{output_dir}/NewBader_pairwise_tests_significant_ordered_{num_clusters}.csv"))
   
    #results_pairwise_tests_clustered <- results_pairwise_tests %>% filter(adj.P.Val <= 0.05)
    #results_pairwise_tests_clustered %<>% mutate(A=ifelse(estimate > 0, cluster_i, cluster_j), B=ifelse(estimate < 0, cluster_i, cluster_j), estimate=abs(estimate)) %>% dplyr::select(term, Title, A, B, estimate, P.Value, adj.P.Val)
    #write.csv(x=results_pairwise_tests_clustered, file=glue("{output_dir}/NewBader_pairwise_tests_significant_ordered_{num_clusters}.csv"))
    
    }
}



#####################################
# Comparing the Different Scores

SummarizeAcrossScores<-function(Pairwise_Dir="Data/Outputs/ClusterPairwiseComp",
                                total_clusters="10"){
  
  score_list<-list()
  for (gscore in c(400,700,900)){
    for (num_clusters in 2:total_clusters){
      dat<-read.csv(glue("{Pairwise_Dir}/{gscore}/NewBader_pairwise_tests_significant_ordered_{num_clusters}_majordiff.csv"))
      dat<-dat %>% dplyr::select(Title,A)%>% unique()
      
      if (nrow(dat) > 0){
        dat$B <- num_clusters
      } else {
        dat <- NULL
      }
      
      score_list[[num_clusters]]<-dat
    }
    
    score_list_frame<-do.call(rbind.data.frame,score_list)
    score_name<-glue("Score{gscore}_list_frame")
    assign(score_name,score_list_frame)
  }
  
  summary_list <- list()
  for (num_clusters in 2:10){
    cluster_list<-list()
    for (clust in 1:num_clusters){
      
      dat1 <- Score400_list_frame %>% filter(B==num_clusters & A==clust)
      dat2 <- Score700_list_frame %>% filter(B==num_clusters & A==clust)
      dat3 <- Score900_list_frame %>% filter(B==num_clusters & A==clust)
      
      dat_list<-list()
      for (x in 1:nrow(dat1)){
        term <- dat1 %>% dplyr::select(Title) %>% dplyr::slice(x)
        if (nrow(term) > 0){
          blah1 <- term[1,1] %in% dat1$Title 
          blah2 <- term[1,1] %in% dat2$Title
          blah3 <- term[1,1] %in% dat3$Title
          if (blah1 == TRUE && blah2 == TRUE && blah3 == TRUE){
            dat_list[[x]] <- data.frame(term,clust,num_clusters)
          } else {
            dat_list[[x]] <- NULL
          }
        }
      }
      dat_frame <- do.call(rbind.data.frame,dat_list)
      cluster_list[[clust]] <- dat_frame
    }
    cluster_frame<-do.call(rbind.data.frame,cluster_list)
    summary_list[[num_clusters]] <- cluster_frame
  }
  
  summary_frame<-do.call(rbind.data.frame,summary_list)
  write.csv(summary_frame,glue("{Pairwise_Dir}/Summary_Across_Scores.csv"))
  return(summary_frame)
}


  