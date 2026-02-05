#I attempt to use namespaces at all times. This library import list might be dropped in the future once appropriate namespaces across the entire codebase are ensured.
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(data.table)
library(stringr)
library(scales)

##### NOTES #####
#This script seeks to visualize structurally complex HDRs.
#Some genomes will display unique structural differences that lack reliable anchors of homologous sequences that map adjacent to a selected reference genomic region.
#In result, each region can require manual tweaks to the selected reference boudaries in order to find the homologous anchors.
#This version attempts to automatically select boundaries given a target region
#Since strain de-novo genome assemblies vary in contiguity and completeness, some genomes may be dropped to avoid visualizing an incomplete or fragmented haplotype.
#Anchor-related parameters can help tune the stringency of the alignment selection and the distance from the target region, allowing for the optimal haplotypePlot to be generated.
#currently this script is under development
#test data is currently too large for github, and will temporarily only exist locally
#appropriately sized test data will be made available once working prod version 1 is finished
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
#this is a quirk of my input data
#mainly included for exploratory purposes, they are later omitted from the plots due to lack of homology in chromosomal arms between both C.b. and C.n.
#This will change in the future, once the mapping and anchoring logic and sample selection is refined
nigonis <- c("JU1422","ECA2852","ECA2857","EG5268","JU1418","JU2617","JU1419","JU2484","JU4356","NIC2143","NIC2150","NIC2152","VSL2202","VX153","YR106","ZF1220")
#refs <- c("MY681","JU3207","JU2536")

reference_bins <- readr::read_tsv("../working_data/QX1410_genomic_windows.1kb.bed",col_names = c("CHROM","START","END"))

#instead of targeting a region at the alignment stage, we will construct a map of anchors across a defined set of samples. 
#These anchors can then be queried for visualization.
#first, some exploratory work

#high confidence intervals
#we drop the nigonis for now
transformed_coords_hidy <- transformed_coords %>% 
  dplyr::filter(IDY>=95) %>%
  dplyr::filter(L1 >5e3 & L2 > 5e3) %>%
  dplyr::filter(!(STRAIN %in% nigonis))

transformed_coords_anyidy <- transformed_coords %>% 
  dplyr::filter(L1 >1e3 & L2 > 1e3) %>%
  dplyr::filter(!(STRAIN %in% nigonis))

# transformed_coords_hidy %>%
#   distinct(HIFI, LENQ) %>%
#   ggplot(aes(x = LENQ/1e3)) +
#   geom_histogram(binwidth = 1) +
#   theme_classic() +
#   xlim(0,100)

hifi_lengths <- transformed_coords_hidy %>%
  distinct(STRAIN, HIFI, LENQ)

# thresholds in kb
N_vals <- seq(10, 50, by = 10)

strain_stats <- dplyr::distinct(transformed_coords_hidy, STRAIN, HIFI, LENQ) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::arrange(dplyr::desc(LENQ), .by_group = TRUE) %>%
  dplyr::mutate(total_bp = sum(LENQ), half_bp = total_bp * 0.5, cum_bp = cumsum(LENQ)) %>%
  dplyr::summarise(n_HIFI = dplyr::n(),
                   total_bp = dplyr::first(total_bp),
                   N50_bp = LENQ[which(cum_bp >= half_bp)[1]],
                   .groups = "drop")

plot_df <- dplyr::distinct(transformed_coords_hidy, STRAIN, HIFI, LENQ) %>%
  dplyr::left_join(strain_stats, by = "STRAIN") %>%
  tidyr::crossing(N_kb = N_vals) %>%
  dplyr::mutate(under_N_flag = dplyr::if_else(LENQ <= N_kb * 1000, 1L, 0L),
                under_N_bp = LENQ * under_N_flag) %>%
  dplyr::group_by(STRAIN, N_kb) %>%
  dplyr::summarise(under_N_bp = sum(under_N_bp),
                   total_bp = dplyr::first(total_bp),
                   n_HIFI = dplyr::first(n_HIFI),
                   N50_bp = dplyr::first(N50_bp),
                   proportion = (under_N_bp / total_bp) * 100,
                   .groups = "drop")

strain_order <- plot_df %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(N_range = max(proportion) - min(proportion), .groups = "drop") %>%
  dplyr::arrange(desc(N_range)) %>%
  dplyr::pull(STRAIN)

plot_df %>%
  dplyr::mutate(STRAIN = factor(STRAIN, levels = strain_order)) %>%
  ggplot() +
  geom_point(aes(x = STRAIN, y = proportion, color = as.character(N_kb)), size = 3) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45,hjust = 1))+
  labs(x = "Strain",
       y = "Genome percent with contigs under N kb",
       color = "N (kb)")
    
plot_df %>%
  mutate(STRAIN = factor(STRAIN, levels = strain_order)) %>%
  ggplot() +
  geom_vline(xintercept =unique(plot_df$N50_bp)/1e6, linetype="dashed",color="grey50")+
  geom_point(aes(x = N50_bp/1e6, y = proportion, color = as.character(N_kb)), size = 3) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45,hjust = 1))+
  labs(x = "N50",
       y = "Genome percent with contigs under N kb",
       color = "N (kb)")


#intervals to data table
dt <- as.data.table(transformed_coords_hidy)
dt[, CHROM := REF]

#we don't assume orientation of reference contig
dt[, `:=`(
  ref_start = pmin(S1, E1),
  ref_end   = pmax(S1, E1)
)]

# keep S2/E2 and HIFI so they survive into ov
dt <- dt[, .(CHROM, STRAIN, IDY, ref_start, ref_end, HIFI, S2, E2,LENQ,L2)]

#genomic bins to data table
bins <- as.data.table(reference_bins)
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


#we repeat steps for anyidy
#intervals to data table
dt_anyidy<- as.data.table(transformed_coords_anyidy)
dt_anyidy[, CHROM := REF]

#we don't assume orientation of reference contig
dt_anyidy[, `:=`(
  ref_start = pmin(S1, E1),
  ref_end   = pmax(S1, E1)
)]

# keep S2/E2 and HIFI so they survive into ov
dt_anyidy <- dt_anyidy[, .(CHROM, STRAIN, IDY, ref_start, ref_end, HIFI, S2, E2,LENQ,L2)]
setkey(dt_anyidy, CHROM, ref_start, ref_end)

#find overlaps
ov_anyidy <- foverlaps(
  x = dt_anyidy,
  y = bins,
  by.x = c("CHROM", "ref_start", "ref_end"),
  by.y = c("CHROM", "bin_start", "bin_end"),
  type = "any",
  nomatch = 0L
)

#set midpoint
ov_anyidy[, bin_mid := (bin_start + bin_end) / 2]

# alignment length on the WI genome 
ov_anyidy[, aln_len := abs(E2 - S2)]

# sort so the best candidate is first within each bin
setorder(ov_anyidy, CHROM, STRAIN, bin_start, -aln_len,-IDY)

# keep one row per bin (drops remaining ties as duplicates)
#we then explore data at the alignment level (rather than bin level)
om_anyidy <- ov_anyidy[, .SD[1], by = .(CHROM, STRAIN, bin_start, bin_end, bin_mid)] 
om_anyidy_freq <- as.data.frame(om_anyidy) %>%
  dplyr::group_by(CHROM,bin_start,bin_end) %>%
  dplyr::mutate(freq=(n()/nstr)*100) %>%
  dplyr::ungroup()


hm<- om %>% 
  dplyr::select(CHROM,STRAIN,IDY,ref_start,ref_end,HIFI,S2,E2,aln_len,LENQ) %>%
  dplyr::group_by(STRAIN,CHROM)%>%
  dplyr::distinct(HIFI,S2,.keep_all = T)


# order strains (optional but usually nicer)
strain_levels <- om_freq %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(mean_idy = mean(IDY, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(mean_idy)) %>%
  dplyr::pull(STRAIN)

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

#target region position
#input args
hap_start = 2e6
hap_end = 3e6
hap_chrom = "I"

#offset parameter (how far away from your boundaries are you willing to expand to search for mappable anchors?)
#will be in args
offset = 5e4

#offset mode defines the search strategy when an offset is allowed
#if "bidirectional", offset can find mappable anchors both outside or inside of the region
#if "expand", offset can find mappable anchors only outside of the region
#will be in args
offset_mode = "bidirectional"

off_start = hap_start - offset
off_end = hap_end + offset

ref_genes <- gffCat %>%
  dplyr::filter(type=="gene" & STRAIN=="QX1410") %>%
  dplyr::mutate(refStart=hap_start,refEnd=hap_end) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("Name","post"),sep=';biotype=') %>%
  #tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>% #some of these commented lines can be included depending on the GFF format
  #dplyr::mutate(seqname=paste0(seqname)) %>%
  #tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("ID=gene:","",Name)) %>% 
  dplyr::mutate(seqname=Name) %>%
  dplyr::select(seqid,start,end,strand,Name,seqname,STRAIN,refStart,refEnd,seqname) 

ref_tran <- gffCat %>%
  dplyr::filter(type=="mRNA" & STRAIN=="QX1410" & !seqid=="MtDNA") %>%
  tidyr::separate(attributes, into=c("ID","Parent","biotype","alias","tag"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=transcript:","",ID)) %>%
  dplyr::mutate(alias=gsub("locus=","",alias)) %>%
  dplyr::mutate(Parent=gsub("Parent=gene:","",Parent)) %>%
  dplyr::mutate(alias=ifelse(alias=="hypothetical_protein",Parent,alias)) %>%
  dplyr::filter(Parent %in% ref_genes$Name) %>%
  dplyr::select(ID,Parent,alias) %>%
  dplyr::rename(tranname=ID) %>%
  dplyr::mutate(tranname=paste0("transcript_",tranname)) %>%
  dplyr::left_join(ref_genes,by=c('Parent'='Name'))

#get alt gene names/aliases
aliases <- ref_tran %>% dplyr::select(seqname,tranname,alias)

#minor diagnostic plot to visualize the REF loci captured by the HDR
#this is your REF haplotype
gene_dirs <- ggplot() + 
  geom_hline(yintercept = 1.72) +
  geom_hline(yintercept = 1.22) +
  geom_text(data=ref_tran %>% dplyr::filter(strand=="+" & seqid==hap_chrom),aes(x = start + ((end - start) / 2), y = 1.75,
                label = ifelse(strand == "+", "\u25B6", "\u25C0")),
            size = 4) +
  geom_text(data=ref_tran %>% dplyr::filter(strand=="-" & seqid==hap_chrom),aes(x = start + ((end - start) / 2), y = 1.25,
                             label = ifelse(strand == "+", "\u25B6", "\u25C0")),
            size = 4) +
  geom_rect(data=ref_tran %>% dplyr::filter(seqid==hap_chrom),aes(xmin=start,xmax=end,ymin=0.75,ymax=0.5))+
  geom_vline(xintercept = c(hap_start,hap_end),linetype='dashed')+
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank()) +
  ylab("GENES")+
  scale_x_continuous(expand = c(0.01, 0)) +
  coord_cartesian(xlim = c(hap_start, hap_end))+
  scale_y_continuous(expand = c(0.01, 0), limits = c(0.5,2))


#df for ggplot
df_hm <- as.data.frame(hm) 

# order strains (optional but usually nicer)
strain_levels <- df_hm %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(mean_idy = mean(IDY, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(mean_idy)) %>%
  dplyr::pull(STRAIN)

df_hm <- df_hm %>%
  dplyr::mutate(
    STRAIN = factor(STRAIN, levels = strain_levels),
    y = as.numeric(STRAIN)
  )

#genome accessibility map
target_aln_plt <- ggplot(df_hm %>% dplyr::filter(CHROM==hap_chrom)) +
  geom_rect(aes(
    xmin = ref_start,
    xmax = ref_end,
    ymin = y - 0.45,
    ymax = y + 0.45,
    fill = HIFI
  )) +
  facet_grid(CHROM ~ ., scales = "free_x", space = "free_x") +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)) +
  scale_y_continuous(
    breaks = seq_along(levels(df$STRAIN)),
    labels = levels(df$STRAIN),
    expand = c(0.01, 0)) +
  labs(x = "Physical position (Mb)", y = "STRAIN") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(hap_start,hap_end))

cowplot::plot_grid(gene_dirs,target_aln_plt,nrow=2,align = "v",axis = "lr",rel_heights = c(0.1,1))


# om_noPCG <- om_freq %>%
#   dplyr::filter((bin_start >= hap_start & bin_start < hap_end) | (bin_end > hap_start & bin_end <= hap_end)) %>%
#   dplyr::filter(CHROM==hap_chrom)

ggplot() + 
  annotate("rect",xmin = hap_start,xmax = hap_end,ymin = Inf,ymax = -Inf,fill = "grey50",alpha = 0.4)+
  geom_segment(data=om_freq %>%
                 dplyr::filter(STRAIN %in% c("VX34") & CHROM==hap_chrom)%>%
                 dplyr::group_by(STRAIN,bin_start,bin_end) %>%
                 dplyr::ungroup(),aes(x=ref_start,xend=ref_end,y=S2,yend=E2,color=IDY)) +
  scale_color_gradientn(colours = c("red", "yellow", "green"),
                        limits = c(95, 100),
                        oob = scales::squish,
                        name = "% Identity") +
  theme_classic() +
  facet_wrap(~STRAIN,scales = "free_y",nrow=4) +
  coord_cartesian(xlim=c(hap_start,hap_end))+
  scale_x_continuous(expand = c(0.01,0),
                     labels = label_number(scale = 1e-6, suffix = " Mb")) +
  scale_y_continuous(expand = c(0.01,0),
                     labels = label_number(scale = 1e-6, suffix = " Mb"))+
  xlab("")+
  ylab("Contig position")

om_noPCG_toReorient <- om_freq %>%
  dplyr::select(STRAIN, HIFI, S2, E2, LENQ) %>% 
  dplyr::distinct(STRAIN, HIFI, S2, E2, .keep_all = TRUE) %>%
  dplyr::mutate(INV = if_else(S2 > E2, "INV", "NOINV")) %>%
  dplyr::group_by(STRAIN, HIFI) %>%
  dplyr::mutate(inv_len   = sum(LENQ[INV == "INV"], na.rm = TRUE),
                total_len = sum(LENQ, na.rm = TRUE),
                REOR      = (inv_len / total_len) > 0.60) %>%
  dplyr::ungroup() 

#this returns true inversions, if any
#could be written out for other purposes
inversions <- om_noPCG_toReorient %>%
  dplyr::filter((REOR==F&INV=="INV")|(REOR==T&INV=="NOINV")) %>%
  dplyr::distinct(STRAIN,HIFI,S2,E2)

# om_noPCG_reoriented <- om_noPCG %>%
#   dplyr::left_join(om_noPCG_toReorient %>% 
#                      dplyr::distinct(STRAIN,HIFI,REOR),by=c("STRAIN","HIFI")) ##### THIS MIGHT NOT BE NEEDED AT THIS STAGE

toFilter <- as.data.table(om_freq)
regs <- as.data.table(ref_tran)

regs2 <- regs[, .(CHROM = seqid, start, end)]
bins <- bins[, .(CHROM,START = bin_start, END = bin_end)]
setkey(regs2, CHROM, start, end)
setkey(bins, CHROM, START, END)

hits <- foverlaps(
  x = bins,                      # intervals to keep/drop
  y = regs2,                     # intervals to overlap against
  by.x = c("CHROM","START","END"),
  by.y = c("CHROM","start","end"),
  type = "any",                  # any overlap
  nomatch = 0L                   # only return overlapping rows
)

drop_bins <- unique(hits[, .(CHROM, bin_start = START, bin_end = END)])

om_noPCG_noOverlap <- toFilter[!drop_bins, on = .(CHROM, bin_start, bin_end)] %>%
  dplyr::mutate(
    STRAIN = factor(STRAIN, levels = strain_levels),
    y = as.numeric(STRAIN))

target_aln_filt_plt <- ggplot() +
  geom_rect(data=ref_tran %>% dplyr::filter(seqid==hap_chrom),aes(xmin=start,xmax=end,ymin=Inf,ymax=-Inf),alpha=0.7)+
  geom_rect(data=om_noPCG_noOverlap%>% dplyr::filter(CHROM==hap_chrom),
            aes(xmin = bin_start,
                xmax = bin_end,
                ymin = y - 0.45,
                ymax = y + 0.45,
                fill = HIFI)) +
  facet_grid(CHROM ~ ., scales = "free_x", space = "free_x") +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)) +
  scale_y_continuous(
    breaks = seq_along(levels(df$STRAIN)),
    labels = levels(df$STRAIN),
    expand = c(0.01, 0)) +
  labs(x = "Physical position (Mb)", y = "STRAIN") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(hap_start,hap_end))

cowplot::plot_grid(gene_dirs,target_aln_filt_plt,nrow=2,align = "v",axis = "lr",rel_heights = c(0.1,1))

om_limits <- as.data.frame(om_noPCG_noOverlap) %>%
  dplyr::filter((bin_start >= off_start & bin_start < off_end) | (bin_end > off_start & bin_end <= off_end)) %>%
  dplyr::filter(CHROM==hap_chrom) %>%
  dplyr::distinct(CHROM,bin_start,bin_end,freq) %>%
  dplyr::filter(freq==100) %>%
  dplyr::mutate(start_dist=bin_start-hap_start,
                end_dist=bin_end-hap_end) %>%
  dplyr::mutate(start_dir=ifelse(start_dist<0,"NEG","POS"),
                end_dir=ifelse(end_dist<0,"NEG","POS")) %>%
  dplyr::mutate(anchor_status=ifelse(start_dir=="NEG" & start_dist > -offset,"START_ANCHOR",
                                    ifelse(start_dir=="POS" & start_dist < offset,"START_ANCHOR",
                                           ifelse(end_dir=="NEG" & end_dist > -offset,"END_ANCHOR",
                                                  ifelse(end_dir=="POS" & end_dist < offset,"END_ANCHOR","NO_ANCHOR"))))) %>%
  dplyr::filter(anchor_status!="NO_ANCHOR")

if (offset_mode=="bidirectional") {
  anchors <- om_limits %>%
    dplyr::group_by(anchor_status) %>%
    dplyr::mutate(dist_key = abs(dplyr::if_else(anchor_status == "START_ANCHOR", start_dist, end_dist)),
                  pref_key = dplyr::if_else(anchor_status == "START_ANCHOR", start_dir == "NEG", end_dir == "POS")) %>%
    dplyr::arrange(dist_key,
                   dplyr::desc(pref_key),
                   .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dist_key, -pref_key) 
  
  anchor_bound <- as.data.frame(om_anyidy_freq) %>%
    dplyr::filter(CHROM==hap_chrom) %>%
    dplyr::filter(bin_start>=min(anchors$bin_start) & bin_end <= max(anchors$bin_end))
  
  test<- anchor_bound
    
  } else if (offset_mode=="expand") {
    anchors <- om_limits %>%
      dplyr::filter((start_dir=="NEG" & anchor_status == "START_ANCHOR") | (end_dir=="POS" & anchor_status == "END_ANCHOR")) %>%
      dplyr::group_by(anchor_status) %>%
      dplyr::mutate(dist_key = abs(dplyr::if_else(anchor_status == "START_ANCHOR", start_dist, end_dist)),
                    pref_key = dplyr::if_else(anchor_status == "START_ANCHOR", start_dir == "NEG", end_dir == "POS")) %>%
      dplyr::arrange(dist_key,
                     dplyr::desc(pref_key),
                     .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(-dist_key, -pref_key)
    
    anchor_bound <- as.data.frame(om_anyidy_freq) %>%
      dplyr::filter(CHROM==hap_chrom) %>%
      dplyr::filter(bin_start>=min(anchors$bin_start) & bin_end <= max(anchors$bin_end))

  } else {
    cat("!!! OFFSET_MODE NOT SET !!!")
}

bits_to_mat <- function(bitstrings) {
  # returns a matrix: nBins x S (0/1)
  do.call(rbind, strsplit(bitstrings, "", fixed = TRUE)) |>
    apply(2, as.integer)
}

if (nrow(anchors)<2) {
  cat("A PAIR OF ANCHORS NEAR YOUR BOUNDARIES USING THE DEFINED OFFSET MODE AND LENGTH WAS NOT FOUND.\nCONSIDER CHANGING YOUR OFFSET MODE TO BIDIRECTIONAL OR INCREASING THE OFFSET LENGTH.\nSWITCHING ANCHOR SEARCH STRATEGY TO DROP SAMPLES UNTIL WE FIND A PAIR OF ANCHORS.")
  strain_levels <- unique(transformed_coords_hidy$STRAIN)
  S <- length(strain_levels)
  
  # choose a drop budget
  N <- S/2
  min_keep <- S - N
  
  strain_index <- tibble(
    STRAIN = strain_levels,
    bit_pos = seq_len(S)  # 1..S
  )
  
  om_wBinBits <- tibble::as_tibble(om_noPCG_noOverlap) %>%
    dplyr::filter((bin_start >= off_start & bin_start < off_end) | (bin_end > off_start & bin_end <= off_end)) %>%
    dplyr::filter(CHROM==hap_chrom) %>%
    dplyr::distinct(CHROM, bin_start, bin_end, STRAIN) %>%  # presence, de-dup
    dplyr::left_join(strain_index, by = "STRAIN") %>%
    dplyr::group_by(CHROM, bin_start, bin_end) %>%
    dplyr::summarise(
      bitstring = {bits <- rep.int("0", S)
      bits[bit_pos] <- "1"
      paste0(bits, collapse = "")},
      .groups = "drop") %>%
    dplyr::filter((bin_start >= off_start & bin_start <= hap_start+offset) | (bin_end <= off_end & bin_end >= hap_end-offset))
  
  # test case to force missing anchors 
  # all_ones <- paste0(rep("1", S), collapse = "")
  # om_test <- om_wBinBits %>%
  #   dplyr::filter(bitstring != all_ones)
  
  start_cand <- om_wBinBits %>%
    dplyr::filter(between(bin_start, hap_start - offset, hap_start + offset)) %>%
    dplyr::mutate(dist_start = abs(bin_start - hap_start)) %>%
    dplyr::arrange(dist_start)
  
  end_cand <- om_wBinBits %>%
    dplyr::filter(between(bin_end, hap_end - offset, hap_end + offset)) %>%
    dplyr::mutate(dist_end = abs(bin_end - hap_end)) %>%
    dplyr::arrange(dist_end)
  
  #sample by site matrices of candidate start anchors 
  start_mat <- bits_to_mat(start_cand$bitstring)
  #sample by site matrices  of candidate end anchors
  end_mat <- bits_to_mat(end_cand$bitstring)
  
  #matrix product yields sample counts per pair of start and end anchors
  keep_mat <- start_mat %*% base::t(end_mat) 
  
  #we search for lowest drop budget in matrix product
  #we don't search under the drop budget threshold S - N
  #map the start and end position of each candidate anchor, as well as distances relative to the target boundaries
  pairs <- tidyr::expand_grid(i = base::seq_len(base::nrow(start_cand)),
                              j = base::seq_len(base::nrow(end_cand))) %>%
    dplyr::mutate(keep = base::as.integer(keep_mat[base::cbind(i,j)]),
                  dropped = S - keep) %>%
    dplyr::filter(keep >= min_keep) %>%
    dplyr::mutate(start_bin_start = start_cand$bin_start[i],
                  start_bin_end = start_cand$bin_end[i],
                  end_bin_start = end_cand$bin_start[j],
                  end_bin_end = end_cand$bin_end[j],
                  dist_start = start_cand$dist_start[i],
                  dist_end = end_cand$dist_end[j],
                  total_dist = dist_start + dist_end) %>%
    dplyr::arrange(dropped,total_dist,dist_start,dist_end)
  
  #select best pair (lowest drop rate, highest proximity)
  #there could be ties given distances are +/- from target region
  #might revisit to improve search strategy for further criteria
  #for now, ties are probably equivalent in value to the user, so we slice the first
  best_pair <- pairs %>%
    dplyr::slice(1)
  
  #rewrite anchors, might not be needed
  anchors <- as.data.frame(om_noPCG_noOverlap) %>%
    dplyr::filter(CHROM==hap_chrom) %>%
    dplyr::distinct(CHROM,bin_start,bin_end,freq) %>%
    dplyr::filter(bin_start==best_pair$start_bin_start[[1]] | bin_start==best_pair$end_bin_start[[1]]) %>%
    dplyr::mutate(anchor_status=ifelse(bin_start==min(bin_start),"START_ANCHOR","END_ANCHOR"))
  
  i_best <- best_pair$i[[1]]
  j_best <- best_pair$j[[1]]
  
  common_vec <- start_mat[i_best, ] * end_mat[j_best, ]  # 1 if present in BOTH
  dropped_pos <- base::which(common_vec == 0L)
  
  #find the strains that are dropped with the selected anchor pair
  dropped_strains <- strain_index$STRAIN[dropped_pos]
  #pull all information within selected anchors
  #drop strains
  #overwrite anchor_bound
  anchor_bound <- as.data.frame(om_anyidy_freq) %>%
    dplyr::filter(CHROM==hap_chrom & !(STRAIN %in% dropped_strains)) %>%
    dplyr::filter(bin_start>=best_pair$start_bin_start[[1]] & bin_end <= best_pair$end_bin_end[[1]])
} else if (nrow(anchors)>2) {
  cat("FOUND A LIL' BUG IN PRIMARY ANCHOR SELECTION. FIX IT!")
} 


anchor_positions <- anchors %>%
  dplyr::summarise(
    xmin = min(bin_start),
    xmax = max(bin_end)
  ) %>%
  tidyr::pivot_longer(cols = everything(), values_to = "x")

gene_dirs2 <- ggplot() + 
  geom_hline(yintercept = 1.72) +
  geom_hline(yintercept = 1.22) +
  geom_text(data = ref_tran %>% dplyr::filter(strand == "+" & seqid == hap_chrom),
            aes(x = start + ((end - start) / 2),
                y = 1.75,
                label = "\u25B6"),size = 4) +
  geom_text(data = ref_tran %>% dplyr::filter(strand == "-" & seqid == hap_chrom),
            aes(x = start + ((end - start) / 2),
                y = 1.25,
                label = "\u25C0"),size = 4) +
  geom_rect(data = ref_tran %>% dplyr::filter(seqid == hap_chrom),
            aes(xmin = start, xmax = end, ymin = 0.75, ymax = 0.5)) +
  geom_text(data = anchor_positions,
            aes(x = x, y = 3),
            label = "\u25BC", 
            size = 7,
            color="red") +
  geom_vline(xintercept = c(hap_start, hap_end), linetype = "dashed") +
  geom_vline(xintercept = anchor_positions$x, linetype = "dashed",color="red") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank()
  ) +
  ylab("GENES") +
  scale_x_continuous(expand = c(0.01, 0)) +
  coord_cartesian(xlim = c(off_start, off_end)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0.5, 3))

cowplot::plot_grid(gene_dirs2+coord_cartesian(xlim = c(1.95e6, 2.05e6)),
                   target_aln_filt_plt+coord_cartesian(xlim = c(1.95e6, 2.05e6))+
                     geom_vline(xintercept = hap_start, linetype = "dashed") +
                     geom_vline(xintercept = anchor_positions$x, linetype = "dashed",color="red") ,
                   nrow=2,align = "v",axis = "lr",rel_heights = c(0.1,1))

cowplot::plot_grid(gene_dirs2+coord_cartesian(xlim = c(2.95e6, 3.1e6)),
                   target_aln_filt_plt+coord_cartesian(xlim = c(2.95e6, 3.1e6))+
                     geom_vline(xintercept = hap_end, linetype = "dashed")+  
                     geom_vline(xintercept = anchor_positions$x, linetype = "dashed",color="red") ,
                   nrow=2,align = "v",axis = "lr",rel_heights = c(0.1,1))

cowplot::plot_grid(gene_dirs2+coord_cartesian(xlim = c(min(anchor_positions$x), max(anchor_positions$x))),
                   target_aln_plt+coord_cartesian(xlim = c(min(anchor_positions$x), max(anchor_positions$x)))+
                     geom_vline(xintercept = hap_start, linetype = "dashed") +
                     geom_vline(xintercept = anchor_positions$x, linetype = "dashed",color="red"),
                   nrow=2,align = "v",axis = "lr",rel_heights = c(0.1,1))

  

anchor_bound_rle <- anchor_bound %>%
  dplyr::arrange(STRAIN, bin_start) %>%
  dplyr::group_by(STRAIN) %>%                       # or group_by(CHROM, STRAIN)
  dplyr::mutate(hifi_run = rleid(HIFI)) %>%         # 1,2,3... whenever HIFI changes
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,hifi_run) %>%
  dplyr::summarise(ref_start=min(bin_start),
                   ref_end=max(bin_end),
                   CHROM=first(CHROM),
                   HIFI=first(HIFI),
                   hifi_start=min(S2,E2),
                   hifi_end=max(S2,E2),
                   LENQ=first(LENQ),
                   L2=sum(unique(L2)),
                   STRAIN=first(STRAIN),
                   .groups = "drop") %>%
  dplyr::ungroup() %>%
  dplyr::arrange(STRAIN,ref_start)

# anchor_bound_rle <- anchor_bound %>%
#   dplyr::group_by(STRAIN,HIFI) %>%
#   dplyr::mutate(meta_start=min(ref_start)) %>%
#   dplyr::ungroup() %>%
#   dplyr::arrange(STRAIN,meta_start,bin_start) %>%
#   dplyr::group_by(STRAIN) %>%                       # or group_by(CHROM, STRAIN)
#   dplyr::mutate(hifi_run = rleid(HIFI)) %>%         # 1,2,3... whenever HIFI changes
#   dplyr::ungroup() %>%
#   dplyr::group_by(STRAIN,hifi_run) %>%
#   dplyr::summarise(ref_start=min(bin_start),
#                    ref_end=max(bin_end),
#                    CHROM=first(CHROM),
#                    HIFI=first(HIFI),
#                    hifi_start=min(S2,E2),
#                    hifi_end=max(S2,E2),
#                    LENQ=first(LENQ),
#                    L2=sum(L2),
#                    STRAIN=first(STRAIN),
#                    .groups = "drop") %>%
#   dplyr::ungroup() %>%
#   dplyr::arrange(STRAIN,ref_start)



rle_df <- anchor_bound_rle %>%
  dplyr::arrange(STRAIN, CHROM, ref_start) %>%
  dplyr::group_by(STRAIN, HIFI) %>%
  dplyr::mutate(hifi_first = min(ref_start, na.rm = TRUE)) %>%  # first appearance (within strain)
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(
    y = hifi_run   # 1..K unique HIFI per strain
  ) %>%
  dplyr::ungroup()

strain_levels <- rle_df %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(n_lanes = max(y), .groups = "drop") %>%
  dplyr::arrange(desc(n_lanes)) %>%
  dplyr::pull(STRAIN) 

rle_df <- rle_df %>%
  dplyr::mutate(STRAIN = factor(STRAIN, levels = strain_levels))

rle_plt <- ggplot(rle_df) +
  geom_rect(aes(
    xmin = ref_start,
    xmax = ref_end,
    ymin = y - 0.45,
    ymax = y + 0.45,
    fill = HIFI
  )) +
  facet_grid(STRAIN ~ CHROM, scales = "free_y", space = "free_y") +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = function(lims) seq(ceiling(lims[1]), floor(lims[2])),
    expand = c(0.01, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "HIFI lane") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(hap_start, hap_end))
rle_plt


# ctg_alignable <- transformed_coords %>%
#   dplyr::group_by(STRAIN, HIFI) %>%
#   dplyr::summarise(alignable_bases=sum(L2),
#                    total_bases=first(LENQ),
#                    HIFI=first(HIFI),
#                    STRAIN=first(STRAIN),.groups = "drop") %>%
#   dplyr::mutate(alignable_pct=alignable_bases/total_bases)


# 
# rle_metrics <- rle_df %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(
#     y_max = max(y, na.rm = TRUE),
#     
#     # special case: fully colinear (single lane)
#     is_colinear = (y_max == 1L),
#     
#     # orientation-safe contig coords
#     hifi_lo = pmin(hifi_start, hifi_end, na.rm = TRUE),
#     hifi_hi = pmax(hifi_start, hifi_end, na.rm = TRUE),
#     
#     is_terminal = !is_colinear & (y == 1L | y == y_max),
#     is_internal = !is_colinear & (y > 1L & y < y_max),
#     
#     # distances to contig ends
#     dist_to_left  = hifi_lo - 1,
#     dist_to_right = LENQ - hifi_hi,
#     
#     end_min_dist = if_else(
#       is_terminal,
#       pmin(dist_to_left, dist_to_right, na.rm = TRUE),
#       NA_real_
#     ),
#     
#     aligned_span = abs(hifi_end - hifi_start),
#     
#     internal_prop = if_else(
#       is_internal & (LENQ > 1),
#       aligned_span / (LENQ - 1),
#       NA_real_
#     )
#   ) %>%
#   dplyr::ungroup() #%>%
#   # dplyr::left_join(ctg_alignable, by=c("STRAIN","HIFI")) %>%
#   # dplyr::mutate(align_ratio = aligned_span/alignable_bases)
rle_term_groups <- rle_df %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(
    y_max = max(y, na.rm = TRUE),
    is_terminal_row = (y == 1L | y == y_max)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN, HIFI) %>%
  dplyr::mutate(
    is_terminal_contig = any(is_terminal_row, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

nonterminal_summary <- rle_term_groups %>%
  dplyr::filter(!is_terminal_contig) %>%
  dplyr::group_by(STRAIN, HIFI) %>%
  dplyr::summarise(
    CHROM=first(CHROM),
    L2_sum = sum(L2, na.rm = TRUE),
    LENQ   = first(LENQ),
    L2_per_LENQ = L2_sum / LENQ,
    n_rows = n(),
    .groups = "drop"
  ) 

ctg_alignable <- om_anyidy %>%
  dplyr::group_by(STRAIN, HIFI) %>%
  dplyr::summarise(alignable_bases=sum(unique(L2)),
                   total_bases=first(LENQ),
                   HIFI=first(HIFI),
                   STRAIN=first(STRAIN),.groups = "drop") %>%
  dplyr::mutate(alignable_pct=alignable_bases/total_bases)

rle_term_groups_wintprop <- rle_term_groups %>%
  dplyr::left_join(nonterminal_summary %>% dplyr::select(STRAIN,HIFI,L2_per_LENQ),by=c("STRAIN","HIFI"))%>%
  dplyr::left_join(ctg_alignable %>% dplyr::select(STRAIN,HIFI,alignable_bases,alignable_pct),by=c("STRAIN","HIFI")) 


rle_scatter <-  ggplot() + 
  geom_point(data=rle_term_groups_wintprop %>% dplyr::filter(is_terminal_contig==F),aes(y=L2_per_LENQ/alignable_pct,x=L2,color=STRAIN)) +
  ylab("Percent contig aligned to genome / Percent contig aligned to region")
  
bad_groups <- rle_term_groups_wintprop %>%
  dplyr::filter(L2_per_LENQ / alignable_pct < 0.75) %>%
  dplyr::distinct(STRAIN, HIFI)

rle_df_filtered <- rle_df %>%
  dplyr::anti_join(bad_groups, by = c("STRAIN", "HIFI")) %>%
  dplyr::mutate(
    hifi_lo = pmin(hifi_start, hifi_end, na.rm = TRUE),
    hifi_hi = pmax(hifi_start, hifi_end, na.rm = TRUE)
  ) %>%
  dplyr::arrange(STRAIN, ref_start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(hifi_run = data.table::rleid(HIFI)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN, hifi_run) %>%
  dplyr::summarise(
    ref_start = min(ref_start, na.rm = TRUE),
    ref_end   = max(ref_end,   na.rm = TRUE),
    CHROM     = dplyr::first(CHROM),
    HIFI      = dplyr::first(HIFI),
    hifi_start = min(hifi_lo, na.rm = TRUE),
    hifi_end   = max(hifi_hi, na.rm = TRUE),
    LENQ      = dplyr::first(LENQ),
    L2        = sum(unique(L2), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(STRAIN, ref_start)

strain_levels <- rle_df_filtered %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(n_lanes = max(hifi_run), .groups = "drop") %>%
  dplyr::arrange(desc(n_lanes)) %>%
  dplyr::pull(STRAIN) 

rle_df_filtered <- rle_df_filtered %>%
  dplyr::mutate(STRAIN = factor(STRAIN, levels = strain_levels))

rle_plt <- ggplot(rle_df_filtered) +
  geom_rect(aes(
    xmin = ref_start,
    xmax = ref_end,
    ymin = hifi_run - 0.45,
    ymax = hifi_run + 0.45,
    fill = HIFI
  )) +
  facet_grid(STRAIN ~ CHROM, scales = "free_y", space = "free_y") +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = function(lims) seq(ceiling(lims[1]), floor(lims[2])),
    expand = c(0.01, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "HIFI lane") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(hap_start, hap_end))
rle_plt


rle_plt2 <- ggplot(rle_term_groups_wintprop %>% dplyr::filter(is_terminal_contig==F,STRAIN=="BRC20075")) +
  geom_rect(aes(
    xmin = ref_start,
    xmax = ref_end,
    ymin = L2_per_LENQ - 0.1,
    ymax = L2_per_LENQ + 0.1,
    fill = HIFI
  )) +
  facet_grid(STRAIN ~ CHROM) +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "Percent contig aligned to genome / Percent contig aligned to region") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )+
  coord_cartesian(xlim = c(off_start, off_end))




rle_metrics <- rle_df %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(
    y_max = max(y, na.rm = TRUE),
    is_colinear = (y_max == 1L)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    row_span = abs(hifi_end - hifi_start),
    hifi_at_ref_start = dplyr::if_else(hifi_start <= hifi_end, hifi_start, hifi_end),
    hifi_at_ref_end   = dplyr::if_else(hifi_start <= hifi_end, hifi_end,   hifi_start)
  ) %>%
  dplyr::group_by(STRAIN, HIFI) %>%
  dplyr::mutate(
    aligned_span = sum(row_span, na.rm = TRUE),
    
    # internal-most alignments across the whole contig
    most_internal_left_pos  = hifi_at_ref_start[which.min(ref_start)],
    most_internal_right_pos = hifi_at_ref_end[which.max(ref_end)]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(
    is_terminal = !is_colinear & (y == 1L | y == y_max),
    is_internal = !is_colinear & (y > 1L & y < y_max),
    
    terminal_boundary_pos = dplyr::case_when(
      is_terminal & y == 1L    ~ most_internal_right_pos,  # interior-facing boundary for left terminal
      is_terminal & y == y_max ~ most_internal_left_pos,   # interior-facing boundary for right terminal
      TRUE ~ NA_real_
    ),
    
    end_min_dist = dplyr::if_else(
      is_terminal,
      pmin(terminal_boundary_pos - 1, LENQ - terminal_boundary_pos, na.rm = TRUE),
      NA_real_
    ),
    
    internal_prop = dplyr::if_else(
      is_internal & (LENQ > 1),
      aligned_span / (LENQ - 1),
      NA_real_
    )
  ) %>%
  dplyr::ungroup()
internal_rects <- rle_metrics %>%
  dplyr::filter(!is_colinear, !is.na(internal_prop)) %>%
  dplyr::mutate(metric = internal_prop)

p_internal <- ggplot(internal_rects %>% dplyr::filter(STRAIN=="BRC20075")) +
  geom_rect(aes(
    xmin = ref_start,
    xmax = ref_end,
    ymin = metric-0.01,
    ymax = metric+0.01,
    fill=HIFI
  ), alpha = 0.9) +
  facet_wrap(~STRAIN) +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "Difference in proportion aligned") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(off_start, off_end)) +
  geom_vline(xintercept = anchor_positions$x, linetype = "dashed",color="red")
p_internal

endpoint_rects <- rle_metrics %>%
  dplyr::filter(!is_colinear, is_terminal, !is.na(end_min_dist)) %>%
  dplyr::mutate(metric_kb = end_min_dist / 1e3)

p_endpoint <- ggplot() +
  geom_point(data=endpoint_rects %>% dplyr::filter(hifi_run==1)%>% dplyr::mutate(side="L"),aes(
    x = ref_end,
    y= metric_kb,
    color=side
  ), alpha = 0.9) +
  geom_point(data=endpoint_rects %>% dplyr::filter(hifi_run!=1) %>% dplyr::mutate(side="R"),aes(
    x = ref_end,
    y= metric_kb,
    color=side
  ), alpha = 0.9) +
  facet_grid(STRAIN ~ CHROM) +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    labels = label_number(suffix = " kb"),
    expand = c(0.1, 0)
  ) +
  labs(x = "Physical position (Mb)", y = "min endpoint dist") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(off_start, off_end))

p_endpoint





ggplot() + 
  annotate("rect",xmin = hap_start,xmax = hap_end,ymin = Inf,ymax = -Inf,fill = "grey50",alpha = 0.4)+
  geom_segment(data=transformed_coords_anyidy %>%
                 dplyr::filter(STRAIN %in% c("JU3202") & REF==hap_chrom)%>%
                 #dplyr::group_by(STRAIN,bin_start,bin_end) %>%
                 dplyr::ungroup(),aes(x=S1,xend=E1,y=S2,yend=E2,color=HIFI)) +
  # scale_color_gradientn(colours = c("red", "yellow", "green"),
  #                       limits = c(95, 100),
  #                       oob = scales::squish,
  #                       name = "% Identity") +
  theme_classic() +
  facet_wrap(~STRAIN,scales = "free_y",nrow=4) +
  coord_cartesian(xlim=c(hap_start,hap_end))+
  scale_x_continuous(expand = c(0.01,0),
                     labels = label_number(scale = 1e-6, suffix = " Mb")) +
  scale_y_continuous(expand = c(0.01,0),
                     labels = label_number(scale = 1e-6, suffix = " Mb"),
                     limits= c(3.6e6,4.7e6))+
  xlab("")+
  ylab("Contig position")

ref_tran_reg <- ref_tran %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | (end >= hap_start & end <= hap_end))  %>%
  dplyr::filter(seqid==hap_chrom)





####
#may need a conditional to ensure the boundary gene is fully accessible
####
HV_genelist <- ref_tran_reg$tranname 

#other diagnostic plots
cowplot::plot_grid(
ggplot() + geom_segment(data=hm %>%dplyr::filter(CHROM=="I" & STRAIN=="QG4097"),aes(y=S2,yend=E2,x=ref_start,xend=ref_end,color=HIFI))+ylim(0,5e6),
ggplot() + geom_segment(data=transformed_coords %>%dplyr::filter(REF=="I" & STRAIN=="QG4097"),aes(y=S2,yend=E2,x=S1,xend=E1,color=IDY))+ylim(0,5e6),nrow=2)


ggplot() + geom_rect(data=transformed_coords %>%dplyr::filter(REF=="I" & STRAIN=="QG4097"),aes(xmin=S1/1e6,xmax=E1/1e6,ymin=IDY+0.1,ymax=IDY-0.1,fill=HIFI))


target <- as.data.frame(om_freq) %>%
  dplyr::filter(CHROM==hap_chrom & bin_start > hap_start & bin_end < hap_end)

exmp1 <- ggplot() + 
  annotate("rect",
           xmin = hap_start,
           xmax = hap_end,
           ymin = Inf,
           ymax = -Inf,
           fill = "grey70",
           alpha = 0.4)+
  geom_rect(data=target %>%
              dplyr::group_by(STRAIN,bin_start,bin_end) %>%
              dplyr::ungroup(),aes(xmin=bin_start,xmax=bin_end,ymin=freq+0.2,ymax=freq-0.2,fill=freq)) +
  theme_classic()+
  scale_fill_gradientn(
    colours = c("red", "#D8C36A", "forestgreen"),
    limits = c(40, 100),
    oob = scales::squish,
    name = "Frequency"
  ) +
  scale_x_continuous(expand = c(0.01,0))+
  xlab("Physical position")+
  ylab("Frequency")
  



exmp2 <- ggplot() + 
  annotate("rect",xmin = hap_start,xmax = hap_end,ymin = Inf,ymax = -Inf,fill = "grey50",alpha = 0.4)+
  geom_segment(data=target %>%
                dplyr::filter(STRAIN %in% c("VX34"))%>%
                dplyr::group_by(STRAIN,bin_start,bin_end) %>%
                dplyr::ungroup(),aes(x=ref_start,xend=ref_end,y=S2,yend=E2,color=IDY)) +
  scale_color_gradientn(colours = c("red", "yellow", "green"),
                        limits = c(95, 100),
                        oob = scales::squish,
                        name = "% Identity") +
  theme_classic() +
  facet_wrap(~STRAIN,scales = "free_y",nrow=4) +
  coord_cartesian(xlim=c(hap_start,hap_end))+
  scale_x_continuous(expand = c(0.01,0),
                     labels = label_number(scale = 1e-6, suffix = " Mb")) +
  scale_y_continuous(expand = c(0.01,0),
                     labels = label_number(scale = 1e-6, suffix = " Mb"))+
  xlab("")+
  ylab("Contig position")

cowplot::plot_grid(exmp2,exmp1,nrow=2,align="v",axis = "lr",rel_heights = c(1,0.6))


target %>%
  #dplyr::group_by(STRAIN,HIFI) %>%
  dplyr::filter(STRAIN=="BRC20075") %>%
  ggplot() +
  geom_rect(aes(xmin=ref_start,xmax=ref_end,ymin=0,ymax=1,fill=HIFI)) +
  xlim(hap_start,hap_end)



