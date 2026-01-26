library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)
library(data.table)
library(stringr)

##### NOTES #####
#This script seeks to visualize structurally complex HDRs. This is not a fully automated process and requires attention.
#Some genomes will display unique structural differences that lack reliable anchors of homologous sequences that map adjacent to a selected reference genomic region.
#In result, each region can require manual tweaks to the selected reference boudaries in order to find the homologous anchors.
#Additionally, to ensure we are not artificially assembling a haplotype using contigs mapped from different chromosomes, we can only confidently visualize regions that 
#are encompassed by a single contig mapping between the reference and wild strain genomes.
#Since strain de-novo genome assemblies vary in contiguity and completeness, some genomes must be manually dropped to avoid visualizing an incomplete or fragmented haplotype.
#Lines of code preceded with "MODIFY THIS" allow you explore regions adjacent to the initially selected reference genomic region 
#or select the strains that will be displayed when genome mapping issues are identified in the initial diagnostic plots.
#################

#read all pairwise genome coordinate comparisons
transformed_coords <- readr::read_tsv("../working_data/CBCN_nucmer_db_20250603.tsv",col_names = F) 
colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")
transformed_coords <- transformed_coords %>% dplyr::filter(!STRAIN=="AF16" & !STRAIN=="JU1422")

#read concatentated gene models of every genome (L3 features are removed)
gffCat <- readr::read_tsv("../working_data/CBCN_L1L2_master.tsv", col_names = F)
colnames(gffCat) <- c("seqid","source","type","start","end","score","strand","phase","attributes","STRAIN")
gffCat <- gffCat %>% dplyr::filter(!STRAIN=="AF16.WBPS19" & !STRAIN=="JU1422.WBPS19") %>% dplyr::mutate(STRAIN=ifelse(STRAIN=="QX1410.curated","QX1410",STRAIN))

#read ortholog relationships among gene models
orthos <- readr::read_tsv("../working_data/CBCN_orthogroups.tsv") %>% dplyr::rename(QX1410=QX1410.curated.longest.protein)
strainCol <- colnames(orthos)
strainCol_c1 <- gsub(".braker.longest.protein","",strainCol)
strainCol_c2 <- gsub(".longest.protein","",strainCol_c1)
colnames(orthos) <- strainCol_c2
orthos <- orthos %>% dplyr::select(-AF16.WBPS19,-JU1422.WBPS19) #other reference genomes (AF16, C. briggsae; JU1422, C. nigoni) can be included by dropping this line
#orthos <- orthos %>%dplyr::select(Orthogroup,CB4856,N2)
#strainCol_c2 <-colnames(orthos)

#list of nigoni strains included in orthofinder
#mainly included for exploratory purposes, they are later omitted from the plots due to lack of homology in chromosomal arms between both C.b. and C.n.
#This will change in the future, once t
nigonis <- c("JU1422","ECA2852","ECA2857","EG5268","JU1418","JU2617","JU1419","JU2484","JU4356","NIC2143","NIC2150","NIC2152","VSL2202","VX153","YR106","ZF1220")
#refs <- c("MY681","JU3207","JU2536")

#instead of targeting a region at the alignment stage, we will construct a map of anchors across a defined set of samples. 
#These anchors can then be queried for visualization.




