library(tidyverse)

hom <- read_tsv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt") %>% 
  select(DBClassKey = `DB Class Key`,
         Species = `Common Organism Name`,
         Gene = Symbol, 
         EntrezID = `EntrezGene ID`) %>% 
  mutate(Species = ifelse(Species == "human", "human", "mouse"))

hom_human <- hom %>% 
  filter(Species == "human") %>% 
  select(-Species, -EntrezID) %>% 
  rename(Human = Gene)

hom_mouse <- hom %>% 
  filter(Species == "mouse") %>% 
  select(-Species, -EntrezID) %>% 
  rename(Mouse = Gene)

hom_joined <- inner_join(hom_human, hom_mouse, by = "DBClassKey")



sum(hom_human$`DB Class Key` %in% hom_mouse$`DB Class Key`)
sum(hom_mouse$`DB Class Key` %in% hom_human$`DB Class Key`)

hom

hom %>% 
  select(`DB Class Key`, `Common Organism Name`) 
  
nrow(hom_human)
  
  homologene_ids <- hom$`DB Class Key`[which((hom$Symbol %in% genes) & hom$`Common Organism Name`==species)]
genes_mouse <- hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="mouse, laboratory"))]
genes_human <- hom$Symbol[which((hom$`DB Class Key` %in% homologene_ids) & (hom$`Common Organism Name`=="human"))]


