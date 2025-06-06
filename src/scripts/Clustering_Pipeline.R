library(tidyverse)
library(magrittr)
library(RMINC)
library(glue)
library(MRIcrotome)
library(grid)
library(SNFtool)
library(tmod)
library(Exact)

setwd("/projects/jacob/ClusteringAutism_125Models_Mar2020/")
source("Code/Clustering_Functions.R")

model_list <- MakeModelList(scanbase_scans_file="Data/Resources/scanbase_40um - Scans_31Jan22.csv",
                            scanbase_studies_file="Data/Resources/scanbase_40um - Studies_Feb2023.csv",
                            scanbase_genotypes_file="Data/Resources/scanbase_40um - Genotypes_Feb2023.csv",
                            scanbase_focus="Autism",
                            scanbase_sample_size_threshold=6,
                            base_directory="Data/Outputs/ModelInfo_Paper")


effect_size_data_matrix_abs <- MakeEffectSizeMaps(jdtype="Absolute",model_list=model_list,resolution="50",
                                              dir_determinants="Data/Raw/new_jacobbian_determinants",
                                              output_phenotype_dir="Data/Outputs/EffectSizeMaps_Paper",
                                              base_directory="Data/Outputs/ModelInfo_Paper",
                                              boot="N")

effect_size_data_matrix_rel <- MakeEffectSizeMaps(jdtype="Relative",model_list=model_list,resolution="50",
                                                  dir_determinants="Data/Raw/new_jacobbian_determinants",
                                                  output_phenotype_dir="Data/Outputs/EffectSizeMaps_Paper",
                                                  base_directory="Data/Outputs/ModelInfo_Paper",
                                                  boot="N")

#######
#
# Run These if you want to skip the above two steps

load(file="Data/Outputs/EffectSizeMaps_Paper/200/ES_Matricies_Absolute.RData")
effect_size_data_matrix_abs <- effect_size_data_matrix
load(file="Data/Outputs/EffectSizeMaps_Paper/200/ES_Matricies_Relative.RData")
effect_size_data_matrix_rel <- effect_size_data_matrix
rm(effect_size_data_matrix)
gc()

effect_size_data_matrix_abs[is.nan(effect_size_data_matrix_abs)] <- 0
effect_size_data_matrix_rel[is.nan(effect_size_data_matrix_rel)] <- 0

#
#######

# Bootstrapping the Effect Size Data Matrix

MakeEffectSizeMaps(jdtype="Relative",model_list=model_list,resolution="200",
                   dir_determinants="Data/Raw/new_jacobbian_determinants",
                   output_phenotype_dir="Data/Outputs/EffectSizeMaps_Paper",
                   base_directory="Data/Outputs/ModelInfo_Paper",
                   boot="Y",
                   num_boot="100")

W<- SNFCombine(Data1=effect_size_data_matrix_abs,Data2=effect_size_data_matrix_rel,
               K=10,alpha=0.5,T=20,distfunc="cor",
               output_dir="Data/Outputs/SNFMatrix_Paper")

load(file="Data/Outputs/SNFMatrix_Paper/WMatrix.RData")

# Bootstrapping the SNF W Matrix

SNFCombineBoot(K=10,alpha=0.5,T=20,distfunc="cor",
               resolution="200",
               dir_determinants="Data/Raw/new_jacobbian_determinants",
               output_dir="Data/Outputs/SNFMatrix_Paper/",
               boot_effect_dir="Data/Outputs/EffectSizeMaps_Paper/",
               num_boot="100")

Clusters<-CreateClusters(num_clusters="10",cluster_method="spectral",
                         output_dir="Data/Outputs/Clusters_Paper",
                         NameFile="Data/Raw/Names_Paper.csv")

BootClusters(num_clusters="10",cluster_method="spectral",
             boot_matrix_dir="Data/Outputs/SNFMatrix/200/Boot/",
             output_dir="Data/Outputs/Clusters",
             NameFile="Data/Raw/Names.csv",
             num_boot="100")


CreateClusterAnatomyMaps(num_clusters="10",cluster_method="spectral",
                         average_kind="median",
                         volume_type="absolute",resolution="50",
                         dir_determinants="Data/Raw/new_jacobbian_determinants/",
                         output_dir="Data/Outputs/Clusters_Paper/")

CreateClusterAnatomyMaps(num_clusters="10",cluster_method="spectral",
                         average_kind="median",
                         volume_type="relative",resolution="50",
                         dir_determinants="Data/Raw/new_jacobbian_determinants/",
                         output_dir="Data/Outputs/Clusters_Paper/")

PlotClusterAnatomyMaps(num_clusters="10",
                       average_kind="median",
                       volume_type="relative",
                       resolution="50",
                       dir_determinants="Data/Raw/new_jacobbian_determinants",
                       output_dir="Data/Outputs/Clusters_Paper",
                       slice_dir="Axial")

###########################################
# Bader Analysis

TermsList<-CreateTermsList(SFARI_List="Data/Raw/SFARI-Gene_genes_02-11-2020release_03-02-2020export-1.csv",
                          Bader_List="Data/Raw/Human_Reactome_March_01_2020_symbol.gmt",
                          TermNum="10")

#TermsList<-CreateTermsList(SFARI_List="Data/Raw/SFARI-Gene_genes_02-11-2020release_03-02-2020export-1.csv",
#                           Bader_List="Data/Raw/Human_GO_AllPathways_withPFOCR_with_GO_iea_October_01_2024_symbol.gmt",
#                           TermNum="10")

TermsList<-c("GAP JUNCTION TRAFFICKING AND REGULATION%REACTOME%R-HSA-157858.1",
             "TRANSMISSION ACROSS CHEMICAL SYNAPSES%REACTOME%R-HSA-112315.5",
             "PROTEIN-PROTEIN INTERACTIONS AT SYNAPSES%REACTOME DATABASE ID RELEASE 71%6794362",
             "ADHERENS JUNCTIONS INTERACTIONS%REACTOME%R-HSA-418990.2",
             "TIGHT JUNCTION INTERACTIONS%REACTOME DATABASE ID RELEASE 71%420029",
             "MAPK FAMILY SIGNALING CASCADES%REACTOME DATABASE ID RELEASE 71%5683057",
             "SIGNALING BY ERBB2%REACTOME DATABASE ID RELEASE 71%1227986",
             "SIGNALING BY ERBB4%REACTOME DATABASE ID RELEASE 71%1236394",
             "SIGNALING BY WNT%REACTOME%R-HSA-195721.5",
             "CA2+ PATHWAY%REACTOME DATABASE ID RELEASE 71%4086398",
             "SIGNALING BY NOTCH%REACTOME DATABASE ID RELEASE 71%157118",
             "MTOR SIGNALLING%REACTOME%R-HSA-165159.5",
             "SIGNALING BY VEGF%REACTOME DATABASE ID RELEASE 71%194138",
             "SIGNALING BY HEDGEHOG%REACTOME DATABASE ID RELEASE 71%5358351",
             "AXON GUIDANCE%REACTOME DATABASE ID RELEASE 71%422475",
             "CHROMATIN ORGANIZATION%REACTOME DATABASE ID RELEASE 71%4839726",
             "GENERIC TRANSCRIPTION PATHWAY%REACTOME%R-HSA-212436.9",
             "GENE EXPRESSION (TRANSCRIPTION)%REACTOME DATABASE ID RELEASE 71%74160",
             "LONG-TERM POTENTIATION%REACTOME DATABASE ID RELEASE 71%9620244",
             "SIGNALING BY GPCR%REACTOME%R-HSA-372790.4")

TermsList<-c("SYNAPTIC SIGNALING%GOBP%GO:0099536",
             "CHROMATIN ORGANIZATION%GOBP%GO:0006325",
             "NEGATIVE REGULATION OF TRANSLATION%GOMF%GO:0017148",
             "TRANSCRIPTION COREGULATOR ACTIVITY%GOMF%GO:0003712",
             "NEGATIVE REGULATION OF CELL-SUBSTRATE ADHESION%GOBP%GO:0010812",
             "PROTEASOME ACCESSORY COMPLEX%GOCC%GO:0022624",
             "POSTIVE REGULATION OF MYOTIC CELL CYCLE%GOBP%GO:0045931")


cl_filt<-GetValidGeneClusterList(Cluster_Dir="Data/Outputs/Clusters_Paper",
                                  boot="N",
                                  boot_pick="100")

for (gscore in c(400,700,900,950)){
  
  GetGeneNeighbourhood(GeneScore=gscore, total_clusters="10",
                       output_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                       boot="N",
                       nboot="1000")
  
  GetNeighbourhoodEnrichment(GeneScore=gscore,
                             total_clusters="10",
                             Bader_List="Data/Raw/Human_Reactome_October_01_2023_symbol.gmt",
                             Neighbour_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                             output_dir="Data/Outputs/NeighbourhoodEnrichment_Paper/",
                             boot="N",
                             nboot="1000")
  
  ClusterPairwiseComp(GeneScore=gscore,
                      total_clusters="10",
                      Bader_List="Data/Raw/Human_Reactome_October_01_2023_symbol.gmt",
                      TargetSet_Dir="Data/Outputs/NeighbourhoodInfo_Paper",
                      Enrichment_Dir="Data/Outputs/NeighbourhoodEnrichment_Paper",
                      output_dir="Data/Outputs/ClusterPairwiseComp_Paper")
}


GetGeneNeighbourhood(GeneScore=gscore,
                     total_clusters="10",
                     output_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                     boot="Y",
                     nboot="1000")

GetNeighbourhoodEnrichment(GeneScore=gscore,
                           total_clusters="10",
                           Bader_List="Data/Raw/Human_GO_AllPathways_withPFOCR_with_GO_iea_October_01_2024_symbol.gmt",
                           Neighbour_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                           output_dir="Data/Outputs/NeighbourhoodEnrichment_Paper_GO/",
                           boot="N",
                           nboot="1000")

ClusterPairwiseComp(GeneScore=gscore,
                    total_clusters="10",
                    Bader_List="Data/Raw/Human_GO_AllPathways_withPFOCR_with_GO_iea_October_01_2024_symbol.gmt",
                    TargetSet_Dir="Data/Outputs/NeighbourhoodInfo_Paper",
                    Enrichment_Dir="Data/Outputs/NeighbourhoodEnrichment_Paper_GO/",
                    output_dir="Data/Outputs/ClusterPairwiseComp_Paper_GO")

#DRUG STUFF

GetNeighbourhoodEnrichment(GeneScore=950,
                           total_clusters="10",
                           Bader_List="Data/Raw/Human_DrugBank_all_symbol.gmt",
                           Neighbour_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                           output_dir="Data/Outputs/NeighbourhoodEnrichmentDrug/",
                           boot="Y",
                           nboot="1000")



for (gscore in c(400,700,900,950)){
  
  GetSingleGeneNeighbourhood(GeneScore=gscore,
                             total_clusters="10",
                             output_dir="Data/Outputs/SingleGeneNeighbourhood_Paper/")
  
  GetSingleGeneNeighbourhoodEnrichment(GeneScore=gscore,
                                       total_clusters="10",
                                       Bader_List="Data/Raw/Human_Reactome_October_01_2023_symbol.gmt",
                                       Neighbour_dir="Data/Outputs/SingleGeneNeighbourhood_Paper/",
                                       output_dir="Data/Outputs/NeighbourhoodEnrichment_Paper/SingleGene/")
  
  GetNeighbourhoodEnrichment_NotCluster(GeneScore=gscore,
                                        total_clusters="10",
                                        Bader_List="Data/Raw/Human_Reactome_October_01_2023_symbol.gmt",
                                        Neighbour_dir="Data/Outputs/NeighbourhoodInfo_Paper/",
                                        output_dir="Data/Outputs/NeighbourhoodEnrichment_Paper/")
  
  ClusterPairwiseCompvsALL(GeneScore=gscore,
                           total_clusters="10",
                           Bader_List="Data/Raw/Human_Reactome_October_01_2023_symbol.gmt",
                           TargetSet_Dir="Data/Outputs/NeighbourhoodInfo_Paper",
                           Enrichment_Dir="Data/Outputs/NeighbourhoodEnrichment_Paper",
                           output_dir="Data/Outputs/ClusterPairwiseCompvsALL_Paper")
}



Score_Summary<-SummarizeAcrossScores(Pairwise_Dir="Data/Outputs/ClusterPairwiseComp_Paper",
                                total_clusters="10")
  

