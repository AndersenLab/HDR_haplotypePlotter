library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)
library(data.table)
library(stringr)
library(scales)
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
#This will change in the future, once the mapping and anchoring logic and sample selection is refined
nigonis <- c("JU1422","ECA2852","ECA2857","EG5268","JU1418","JU2617","JU1419","JU2484","JU4356","NIC2143","NIC2150","NIC2152","VSL2202","VX153","YR106","ZF1220")
#refs <- c("MY681","JU3207","JU2536")

refrence_bins <- readr::read_tsv("../working_data/QX1410_genomic_windows.1kb.bed",col_names = c("CHROM","START","END"))

#instead of targeting a region at the alignment stage, we will construct a map of anchors across a defined set of samples. 
#These anchors can then be queried for visualization.
#first, some exploratory work

#high confidence intervals
#we drop the nigonis for now
transformed_coords_hidy <- transformed_coords %>% 
  dplyr::filter(IDY>=95) %>%
  dplyr::filter(L1 >5e3 & L2 > 5e3) %>%
  dplyr::filter(!(STRAIN %in% nigonis))

#intervals to data table
dt <- as.data.table(transformed_coords_hidy)
dt[, CHROM := REF]

#we don't assume orientation of reference contig
dt[, `:=`(
  ref_start = pmin(S1, E1),
  ref_end   = pmax(S1, E1)
)]

# keep S2/E2 and HIFI so they survive into ov
dt <- dt[, .(CHROM, STRAIN, IDY, ref_start, ref_end, HIFI, S2, E2)]

#genomic bins to data table
bins <- as.data.table(refrence_bins)
setnames(bins, c("START", "END"), c("bin_start", "bin_end"))

#set keys for overlap 
setkey(bins, CHROM, bin_start, bin_end)
setkey(dt, CHROM, ref_start, ref_end)

#find overlaps
ov <- foverlaps(
  x = dt,
  y = bins,
  by.x = c("CHROM", "ref_start", "ref_end"),
  by.y = c("CHROM", "bin_start", "bin_end"),
  type = "any",
  nomatch = 0L
)

#set midpoint
ov[, bin_mid := (bin_start + bin_end) / 2]

# alignment length on the WI genome 
ov[, aln_len := abs(E2 - S2)]

# sort so the best candidate is first within each bin
setorder(ov, CHROM, STRAIN, bin_start, -aln_len,-IDY)

# keep one row per bin (drops remaining ties as duplicates)
#we then explore data at the alignment level (rather than bin level)
om <- ov[, .SD[1], by = .(CHROM, STRAIN, bin_start, bin_end, bin_mid)] 

nstr=length(unique(om$STRAIN))

om_freq <- as.data.frame(om) %>%
  dplyr::group_by(CHROM,bin_start,bin_end) %>%
  dplyr::mutate(freq=(n()/nstr)*100) %>%
  dplyr::ungroup()

hm<- om %>% 
  dplyr::select(CHROM,STRAIN,IDY,ref_start,ref_end,HIFI,S2,E2,aln_len) %>%
  dplyr::group_by(STRAIN,CHROM)%>%
  dplyr::distinct(HIFI,S2,.keep_all = T)


# order strains (optional but usually nicer)
strain_levels <- om_freq %>%
  dplyr::group_by(STRAIN) %>%
  summarise(mean_idy = mean(IDY, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_idy)) %>%
  pull(STRAIN)

df <- om_freq %>%
  dplyr::mutate(
    STRAIN = factor(STRAIN, levels = strain_levels),
    y = as.numeric(STRAIN)
  )


ggplot(df) +
  geom_rect(aes(
    xmin = bin_start,
    xmax = bin_end,
    ymin = y - 0.45,
    ymax = y + 0.45,
    fill = freq
  )) +
  facet_grid(CHROM ~ ., scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("red", "yellow", "green"),
    limits = c(0, 100),
    oob = scales::squish,
    name = "Frequency"
  ) +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(df$STRAIN)),
    labels = levels(df$STRAIN),
    expand = c(0.01, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "STRAIN") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank())


ggplot(df%>% dplyr::filter(CHROM=="V")) +
  geom_rect(aes(
    xmin = bin_start,
    xmax = bin_end,
    ymin = y - 0.45,
    ymax = y + 0.45,
    fill = freq
  )) +
  facet_grid(CHROM ~ ., scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("red", "yellow", "green"),
    limits = c(0, 100),
    oob = scales::squish,
    name = "Frequency"
  ) +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0),
    limits=c(15e6,19e6)
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(df$STRAIN)),
    labels = levels(df$STRAIN),
    expand = c(0.01, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "STRAIN") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank())


#df for ggplot
df <- as.data.frame(hm) 

# order strains (optional but usually nicer)
strain_levels <- df %>%
  group_by(STRAIN) %>%
  summarise(mean_idy = mean(IDY, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_idy)) %>%
  pull(STRAIN)

df <- df %>%
  mutate(
    STRAIN = factor(STRAIN, levels = strain_levels),
    y = as.numeric(STRAIN)
  )

#genome accessibility map
ggplot(df) +
  geom_rect(aes(
    xmin = ref_start,
    xmax = ref_end,
    ymin = y - 0.45,
    ymax = y + 0.45,
    fill = HIFI
  )) +
  facet_grid(CHROM ~ ., scales = "free_x", space = "free_x") +
  # scale_fill_gradientn(
  #   colours = c("red", "yellow", "green"),
  #   limits = c(96, 100),
  #   oob = scales::squish,
  #   name = "IDY"
  # ) +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(df$STRAIN)),
    labels = levels(df$STRAIN),
    expand = c(0.01, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "STRAIN") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank())

#other diagnostic plots
cowplot::plot_grid(
ggplot() + geom_segment(data=hm %>%dplyr::filter(CHROM=="I" & STRAIN=="QG4097"),aes(y=S2,yend=E2,x=ref_start,xend=ref_end,color=HIFI))+ylim(0,5e6),
ggplot() + geom_segment(data=transformed_coords %>%dplyr::filter(REF=="I" & STRAIN=="QG4097"),aes(y=S2,yend=E2,x=S1,xend=E1,color=IDY))+ylim(0,5e6),nrow=2)


ggplot() + geom_rect(data=transformed_coords %>%dplyr::filter(REF=="I" & STRAIN=="QG4097"),aes(xmin=S1/1e6,xmax=E1/1e6,ymin=IDY+0.1,ymax=IDY-0.1,fill=HIFI))
