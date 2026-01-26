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
transformed_coords <- readr::read_tsv("../../processed_data/gene_diversity/CBCN_nucmer_db_20250603.tsv",col_names = F) 
colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")
transformed_coords <- transformed_coords %>% dplyr::filter(!STRAIN=="AF16" & !STRAIN=="JU1422")

#read concatentated gene models of every genome (L3 features are removed)
gffCat <- readr::read_tsv("../../processed_data/gene_diversity/CBCN_L1L2_master.tsv", col_names = F)
colnames(gffCat) <- c("seqid","source","type","start","end","score","strand","phase","attributes","STRAIN")
gffCat <- gffCat %>% dplyr::filter(!STRAIN=="AF16.WBPS19" & !STRAIN=="JU1422.WBPS19") %>% dplyr::mutate(STRAIN=ifelse(STRAIN=="QX1410.curated","QX1410",STRAIN))

#read ortholog relationships among gene models
#orthos <- readr::read_tsv("./input/Orthogroups.tsv")
orthos <- readr::read_tsv("../../processed_data/gene_diversity/CBCN_orthogroups.tsv") %>% dplyr::rename(QX1410=QX1410.curated.longest.protein)
strainCol <- colnames(orthos)
strainCol_c1 <- gsub(".braker.longest.protein","",strainCol)
strainCol_c2 <- gsub(".longest.protein","",strainCol_c1)
colnames(orthos) <- strainCol_c2
orthos <- orthos %>% dplyr::select(-AF16.WBPS19,-JU1422.WBPS19) #other reference genomes (AF16, C. briggsae; JU1422, C. nigoni) can be included by dropping this line
#orthos <- orthos %>%dplyr::select(Orthogroup,CB4856,N2)
#strainCol_c2 <-colnames(orthos)

#list of nigoni strains included in orthofinder
#mainly included for exploratory purposes, they are later omitted from the plots due to lack of homology in chromosomal arms between both C.b. and C.n.
nigonis <- c("JU1422","ECA2852","ECA2857","EG5268","JU1418","JU2617","JU1419","JU2484","JU4356","NIC2143","NIC2150","NIC2152","VSL2202","VX153","YR106","ZF1220")
#refs <- c("MY681","JU3207","JU2536")

######### MODIFY THIS ##########
#set your target HDR coordinates
# Figure S17a - 70 kb
hdr_chrom = "V"
hdr_start_pos = 838000
hdr_end_pos = 894000

#Figure S17b - 110kb
# hdr_chrom = "II"
# hdr_start_pos = 12500000
# hdr_end_pos = 12610000

#Figure S17c - 175kb
# hdr_chrom = "I"
# hdr_start_pos = 12469000
# hdr_end_pos = 12644000

#Figure S16 - 120 kb
# hdr_chrom = "I"
# hdr_start_pos = 11880000
# hdr_end_pos = 12000000

#offset lets you explore adjacent regions
offset = 0
hap_chrom = hdr_chrom
hap_start = hdr_start_pos - offset
hap_end = hdr_end_pos + offset 
############################### 

#use reference coordinates from g2g alginments to pull the contigs that contain the alt haplotypes for the HDR
hap_coords <- transformed_coords %>%
  dplyr::filter((REF == hap_chrom & hap_start >= S1 & hap_start <= E1 ) | 
                  (REF == hap_chrom & hap_end >= S1 & hap_end <= E1) | 
                  (REF == hap_chrom & S1 >= hap_start & E1 <= hap_end)) %>%
  dplyr::mutate(inv=ifelse(S2>E2,T,F))  %>%
  dplyr::mutate(St2=ifelse(inv==T,E2,S2),Et2=ifelse(inv==T,S2,E2))
  # dplyr::select(-S2,-E2) %>%
  # dplyr::rename(S2=newS2,E2=newE2)

# naive visualization of g2g alignments for the target region
# multiple contigs can map to the REF region, we need to filter those secondary alignments!
ggplot(hap_coords) + 
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("REF genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

#keep only the contig with the largest extent of alignment with the REF HDR
tigFilt <- hap_coords %>% 
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,HIFI) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L1)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::mutate(rangeDiff=(max(S1,E1)-min(S1,E1))-(hap_end-hap_start)) %>%
  dplyr::ungroup()

#important diagnostic plot
#after the optimal contig is selected, we can see if the HDR boundaries are within the alignment
#for proper visualization, the HDR (grey box in plot) needs to be encompassed by the selected contig
#otherwise some haplotypes will be truncated 
#a step to automatically drop genomes with incomplete coverage of the HDR could be added in the future
#have in mind that the alignments look fragmented because of sequence divergence, but the contig is linear in the original genome assembly
ggplot(tigFilt) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  #geom_rect(xmin=16219634/1e6,xmax=16221917/1e6,ymin=-Inf,ymax=Inf,aes(fill="glc-1")) +
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("QX1410 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') #+
  #scale_fill_manual(values=c("lightgrey"="lightgrey","glc-1"="pink"))

#keep the set of alignments with the largest span (i.e. removes small distant alignments)
tigFilt2 <- tigFilt %>%
  dplyr::arrange(St2) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDiff=lead(St2)-Et2) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 5e4,1,0)) %>% #CAN MODIFY THIS - this is the maximum physical distance allowed between wild genome alignments
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup()

#visualization after small distant alignments are removed
ggplot(tigFilt2) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("QX1410 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

trim_spacer = 1e3
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the WI genome)
tigTrim <- tigFilt2 %>%
  dplyr::arrange(STRAIN,S1) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rboundDist=max(S1,E1)-hap_end) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & inv==F,(E2-(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & inv==T,(E2+(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(E1=ifelse(rboundDist>trim_spacer,(E1-(rboundDist-trim_spacer)),E1)) %>%
  dplyr::mutate(lboundDist=hap_start-min(S1,E1)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & inv==F,(S2+(lboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & inv==T,(S2-(lboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(S1=ifelse(lboundDist>trim_spacer,(S1+(lboundDist-trim_spacer)),S1))

ggplot(tigTrim) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  #geom_rect(xmin=16219634/1e6,xmax=16221917/1e6,ymin=-Inf,ymax=Inf,aes(fill="glc-1")) +
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("QX1410 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

#get the minimum and maximum boundary of the WILD genome alignments that contain the HDR
HV_boundary <- tigTrim %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::select(-leadDiff,-inv) %>%
  dplyr::mutate(inv=ifelse(sum(E2)-sum(S2) >0,F,T)) %>%
  dplyr::mutate(refStart=min(S1,E1),refEnd=max(S1,E1)) %>%
  dplyr::mutate(boundStart=min(S2,E2), boundEnd=max(S2,E2)) %>%
  dplyr::distinct(STRAIN, .keep_all = T) %>%
  dplyr::select(HIFI,boundStart,boundEnd,STRAIN,REF,refStart,refEnd,inv) %>%
  dplyr::ungroup() %>%
  dplyr::rename(boundChrom=HIFI)

#filter the concatenated GFF to extract the gene models of each WILD genome contig boundary
wild_genes <- gffCat %>%
  dplyr::filter(type=="gene" & !(STRAIN=="QX1410")) %>%
  dplyr::mutate(attributes=gsub(";","",attributes)) %>% 
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::select(attributes,seqid,start,end,strand,STRAIN) %>%
  dplyr::rename(Name=attributes)
  
wild_tran <-  gffCat  %>%
  dplyr::filter(type=="mRNA" & !(STRAIN=="QX1410")) %>%
  tidyr::separate(attributes,into=c("tranname","Parent"),sep=";Parent=") %>%
  dplyr::mutate(Parent=gsub(";","",Parent)) %>%
  dplyr::mutate(tranname=gsub("ID=","",tranname)) %>%
  dplyr::select(tranname,Parent,STRAIN) %>%
  dplyr::left_join(wild_genes,by=c("Parent"="Name","STRAIN")) %>%
  dplyr::left_join(HV_boundary,by="STRAIN")
  
#extract the REF genes 
QXStart = min(HV_boundary$refStart)
QXEnd = max(HV_boundary$refEnd)
QX_genes <- gffCat %>%
  dplyr::filter(type=="gene" & STRAIN=="QX1410") %>%
  dplyr::mutate(refStart=QXStart,refEnd=QXEnd) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("Name","post"),sep=';biotype=') %>%
  #tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>% #some of these commented lines can be included depending on the GFF format
  #dplyr::mutate(seqname=paste0(seqname)) %>%
  #tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("ID=gene:","",Name)) %>% 
  dplyr::mutate(seqname=Name) %>%
  dplyr::select(seqid,start,end,strand,Name,seqname,STRAIN,refStart,refEnd,seqname) 

#extract the REF protein-coding transcripts
QX_tran <- gffCat %>%
  dplyr::filter(type=="mRNA" & STRAIN=="QX1410" & !seqid=="MtDNA") %>%
  tidyr::separate(attributes, into=c("ID","Parent","biotype","alias","tag"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=transcript:","",ID)) %>%
  dplyr::mutate(alias=gsub("locus=","",alias)) %>%
  dplyr::mutate(Parent=gsub("Parent=gene:","",Parent)) %>%
  dplyr::mutate(alias=ifelse(alias=="hypothetical_protein",Parent,alias)) %>%
  dplyr::filter(Parent %in% QX_genes$Name) %>%
  dplyr::select(ID,Parent,alias) %>%
  dplyr::rename(tranname=ID) %>%
  dplyr::mutate(tranname=paste0("transcript_",tranname)) %>%
  dplyr::left_join(QX_genes,by=c('Parent'='Name'))

QX_tran_reg <- QX_tran %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | (end >= hap_start & end <= hap_end))  %>%
  dplyr::filter(seqid==hap_chrom)

# #get gene list
HV_genelist <- QX_tran_reg$tranname 
#get alt gene names/aliases
aliases <- QX_tran %>% dplyr::select(seqname,tranname,alias)

#minor diagnostic plot to visualize the REF loci captured by the HDR
#this is your REF haplotype
ggplot(QX_tran_reg) + geom_rect(aes(xmin=start,xmax=end,ymin=1,ymax=2))


#filter orthologous groups using REF genes
#this will establish your orthology relationships between REF and WILD haplotypes
all_orthos_unnest <- orthos %>%
  dplyr::mutate(QX1410 = strsplit(as.character(QX1410), ",")) %>%
  tidyr::unnest(QX1410) %>%
  dplyr::mutate(QX1410=trimws(QX1410)) %>%
  dplyr::mutate(na_count = rowSums(is.na(.))) %>%
  dplyr::filter(na_count < length(strainCol_c2) - 2) %>%
  dplyr::left_join(QX_tran %>% dplyr::select(tranname,seqid,seqname,start,end),by=c("QX1410"="tranname"))  %>%
  dplyr::select(-na_count)

filtOrthos <- orthos %>%
  dplyr::filter(grepl(paste(HV_genelist,collapse="|"),QX1410)) %>%
  dplyr::mutate(QX1410 = strsplit(as.character(QX1410), ",")) %>%
  tidyr::unnest(QX1410) %>%
  dplyr::mutate(QX1410=trimws(QX1410)) %>%
  dplyr::left_join(QX_tran_reg %>% dplyr::select(tranname,seqid,seqname,start,end) %>% dplyr::mutate(og_loc="in_region"),by=c("QX1410"="tranname")) 

inreg_orthos <- filtOrthos %>% dplyr::filter(!is.na(seqid)) 
outreg_orthos <- filtOrthos %>% dplyr::filter(is.na(seqid)) %>% 
  dplyr::select(-seqid,-seqname,-start,-end,-og_loc) %>%
  dplyr::left_join(QX_tran %>% dplyr::select(tranname,seqid,seqname,start,end,refStart,refEnd) %>% 
                     dplyr::mutate(refChrom=hap_chrom) %>%
                     dplyr::mutate(og_loc="out_region"),by=c("QX1410"="tranname")) %>%
  dplyr::filter(seqid==refChrom) %>%
  dplyr::mutate(start_dist=abs(refStart-end),end_dist=abs(start-refEnd)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(min_dist_bases=min(start_dist,end_dist)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(updown=ifelse(end < refStart,"upstream","downstream")) %>%
  dplyr::mutate(status=ifelse(min_dist_bases <2e4,"out_expand","outside"))

if (nrow(outreg_orthos  %>% dplyr::filter(min_dist_bases < 2e4)) > 0) {
  print("WARNING: There is at least one paralog that is within 10 kb of a gene within your defined boundary in QX1410. Your boundary will be automatically expanded to include:")
  print(outreg_orthos %>% dplyr::select(seqid,seqname,start,end,QX1410,min_dist_bases) %>% dplyr::filter(min_dist_bases < 2e4))
}

#generate a lookup table (all_ortho_pairs) which contains all pairwise gene orthologs between REF and WILD
orthoList <- list()
orthoList_bound <- list()
orthoList_raw <- list()
strainCol_iter <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup","QX1410","JU1422.WBPS19","AF16.WBPS19","JU1422")]

for (i in 1:length(strainCol_iter)) {
  
  id=strainCol_iter[i]
  raw_tmp <- orthos %>%
    dplyr::select(Orthogroup,strainCol_iter[i],QX1410) %>%
    dplyr::mutate(str=!!sym(strainCol_iter[i])) %>%
    dplyr::mutate(str = strsplit(as.character(str), ",")) %>%
    tidyr::unnest(str) %>%
    dplyr::mutate(str=trimws(str)) %>%
    dplyr::filter(!is.na(QX1410)) %>%
    dplyr::select(-strainCol_iter[i]) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    dplyr::mutate(has_any_ortho=T) %>%
    dplyr::left_join(wild_tran,by=c("STRAIN","str"="tranname"))
    #dplyr::select(-QX1410)
  
  orthoList_raw[[i]] <- raw_tmp
    
  print(paste0("Mapped orthologs for ",i,"/",length(strainCol_iter)," strains."))
  tmp <- rbind(inreg_orthos %>% dplyr::mutate(status="within") %>% dplyr::select(Orthogroup,strainCol_iter[i],QX1410,seqid,seqname,start,end,og_loc,status),outreg_orthos %>% 
                 dplyr::select(Orthogroup,strainCol_iter[i],QX1410,seqid,seqname,start,end,og_loc,status)) %>%
    dplyr::select(Orthogroup,QX1410,strainCol_iter[i],og_loc,status) %>%
    dplyr::rename(tmpSel=strainCol_iter[i]) %>%
    dplyr::mutate(newSel = strsplit(as.character(tmpSel), ",")) %>%
    tidyr::unnest(newSel) %>%
    dplyr::mutate(newSel=trimws(newSel)) %>%
    dplyr::select(-tmpSel) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    tidyr::separate(newSel,into=c("Name","tnum"),sep="\\.",remove = F) %>%
    dplyr::select(Orthogroup,newSel,Name,STRAIN,QX1410,-tnum,og_loc,status) %>%
    dplyr::rename(tranname=newSel,Parent=Name) %>%
    dplyr::left_join(wild_tran,by=c("tranname","Parent","STRAIN")) 

  orthoList[[i]] <- tmp
  boundg <- tmp %>% 
                            dplyr::filter(og_loc=="in_region" | status=="out_expand") %>%
                            dplyr::select(-og_loc,-status) %>%
                            dplyr::filter(seqid==boundChrom) %>%
                            dplyr::mutate(og_loc=ifelse(((start >= boundStart & start <= boundEnd) | (end >= boundStart & end <= boundEnd)),"in_region","out_region")) %>%
                            dplyr::mutate(start_dist=ifelse(og_loc=="out_region",abs(boundStart-end),NA),end_dist=ifelse(og_loc=="out_region",abs(start-boundEnd),NA)) %>%
                            dplyr::rowwise() %>%
                            dplyr::mutate(min_dist_bases=min(start_dist,end_dist)) %>%
                            dplyr::ungroup() %>%
                            dplyr::mutate(status=ifelse(og_loc=="in_region","within",ifelse(min_dist_bases < 10000 & !is.na(min_dist_bases),"out_expand","outside")))
  check <- boundg %>% dplyr::filter(status=="out_expand")
  
  if (nrow(check) > 0) {
    print(paste0("WARNING: There is at least one paralog that is within 10 kb of a gene within your derived boundary in ",strainCol_iter[[i]],". Your boundary will be automatically expanded to include:"))
    print(check %>% dplyr::select(seqid,tranname,start,end,strand,STRAIN,min_dist_bases) %>% dplyr::distinct(tranname,.keep_all = T))
    
    sorter <- check %>% dplyr::mutate(updown=ifelse(min_dist_bases==start_dist,"upstream","downstream"))
    upstream <- sorter %>% dplyr::filter(updown=="upstream")
    downstream <- sorter %>% dplyr::filter(updown=="downstream")
    
    if (nrow(upstream) > 0) {
      outer_lim <- max(upstream$end)
      inner <- boundg %>% dplyr::filter(status=="within") %>% dplyr::arrange(start) %>% dplyr::filter(start==min(start))
      inner_lim <- min(inner$start)
      seqid_match <- as.character(unique(inner$seqid))
      extension <- raw_tmp %>%
        dplyr::filter(seqid==seqid_match & start > outer_lim & end <inner_lim) %>%
        dplyr::rename(tranname=str) %>%
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,QX1410,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(QX1410 = strsplit(as.character(QX1410), ",")) %>%
        tidyr::unnest(QX1410) %>%
        dplyr::mutate(QX1410=trimws(QX1410))
        
      boundg_inc <- rbind(extension, boundg %>% dplyr::select(-start_dist,-end_dist,-min_dist_bases))
      orthoList_bound[[i]]  <- boundg_inc %>% dplyr::arrange(start)
    } 
    
    if(nrow(downstream) > 0) {
      outer_lim <- min(downstream$end)
      inner <- boundg %>% dplyr::filter(status=="within") %>% dplyr::arrange(start) %>% dplyr::filter(end==max(end))
      inner_lim <- max(inner$end)
      seqid_match <- as.character(unique(inner$seqid))
      extension <- raw_tmp %>%
        dplyr::filter(seqid==seqid_match & start > inner_lim & end < outer_lim) %>%
        dplyr::rename(tranname=str) %>%
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,QX1410,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(QX1410 = strsplit(as.character(QX1410), ",")) %>%
        tidyr::unnest(QX1410) %>%
        dplyr::mutate(QX1410=trimws(QX1410))
      
      boundg_inc <- rbind(extension, boundg %>% dplyr::select(-start_dist,-end_dist,-min_dist_bases))
      orthoList_bound[[i]]  <- boundg_inc %>% dplyr::arrange(start)
    } 
    
    
  } else {
   orthoList_bound[[i]] <- boundg %>% 
     dplyr::select(-start_dist,-end_dist,-min_dist_bases) %>% dplyr::arrange(start)
 }
}

all_ortho_pairs  <- ldply(orthoList,data.frame) 
all_ortho_pairs_bound_pre <-ldply(orthoList_bound,data.frame) %>% 
  dplyr::filter(!status=="outside") 

corr_jumps <- all_ortho_pairs_bound_pre %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::arrange(STRAIN,start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDist=lead(start)-start) %>%
  dplyr::mutate(leadDist=ifelse(is.na(leadDist),0,leadDist)) %>%
  dplyr::mutate(jump=ifelse(lag(leadDist)>5e4 & lag(status)=="within","JUMP","NOJUMP")) %>%
  dplyr::mutate(jump=ifelse(is.na(jump),"NOJUMP",jump)) %>%
  dplyr::mutate(jumpID=rleid(jump)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,jumpID) %>%
  dplyr::mutate(jgroup_size=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(jgroup_size==max(jgroup_size)) %>%
  dplyr::mutate(keep=T)


all_ortho_pairs_bound <- all_ortho_pairs_bound_pre %>%
  dplyr::arrange(STRAIN,start)

all_ortho_pairs_raw <- ldply(orthoList_raw,data.frame) %>% dplyr::select(STRAIN,str,has_any_ortho) %>% dplyr::rename(tranname=str) 

new_boundaries_WI <-  all_ortho_pairs_bound %>%
  dplyr::select(seqid,start,end,STRAIN,tranname,Parent) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(minStart=min(start), maxEnd=max(end)) %>%
  dplyr::distinct(tranname,.keep_all = T) %>%
  dplyr::filter(start==minStart | end==maxEnd) %>%
  dplyr::mutate(gene2gene=paste(Parent,collapse="-")) %>%
  dplyr::distinct(minStart,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::select(seqid,minStart,maxEnd,STRAIN,gene2gene)

QX_expand <- rbind(inreg_orthos %>% 
                       dplyr::mutate(status="within"),outreg_orthos %>% 
                       dplyr::select(-refStart,-refEnd,-refChrom,-start_dist,-end_dist,-min_dist_bases,-updown)) %>%
                       #dplyr::select(Orthogroup,strainCol_iter[i],QX,seqid,seqname,start,end,og_loc,status)) %>%
  dplyr::filter(!status=="outside")

new_boundaries_QX <- QX_expand %>%
  dplyr::mutate(minStart=min(start), maxEnd=max(end)) %>%
  dplyr::distinct(QX1410,.keep_all = T) %>% 
  dplyr::filter(start==minStart | end==maxEnd) %>%
  dplyr::mutate(gene2gene=paste(seqname,collapse="-")) %>%
  dplyr::mutate(STRAIN="QX1410") %>%
  dplyr::distinct(minStart,.keep_all = T) %>%
  dplyr::select(seqid,minStart,maxEnd,STRAIN,gene2gene)

new_boundaries <- rbind(new_boundaries_WI,new_boundaries_QX) %>%
  dplyr::rename(boundStart=minStart,boundEnd=maxEnd)

#find the bound genes for each strain that are not orthologous
boundGenes <- rbind(wild_tran %>% 
                      dplyr::select(-boundChrom,-boundStart,-boundEnd,-REF,-refStart,-refEnd,-inv) %>% 
                      dplyr::mutate(alias=NA),
                    QX_tran %>% dplyr::select(tranname,seqname,STRAIN,seqid,start,end,strand,alias) %>%
                      dplyr::rename(Parent=seqname)) %>%
  dplyr::left_join(new_boundaries,by=c("STRAIN","seqid")) %>%
  dplyr::filter(!is.na(boundStart)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(start >= boundStart & start <= boundEnd) %>%
  dplyr::ungroup()


QX_ad <- boundGenes %>% 
  dplyr::filter(STRAIN=="QX1410") %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(tranname %in% all_orthos_unnest$QX1410,T,F)) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(tranname %in% all_ortho_pairs_bound$QX1410,T,F)) %>%
  dplyr::arrange(start) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup()

reassess_distal <- QX_ad %>% dplyr::filter(has_any_ortho==T & has_bound_ortho==F) 

g_count <- length(unique(QX_ad$Parent))

allstrain_list <- unique(boundGenes$STRAIN)[!unique(boundGenes$STRAIN) %in% nigonis]
#lineages <- readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byLineage_GeoLocAdmCol_20250909.tsv") %>% dplyr::filter(isotype %in% allstrain_list)
#lineage_order <- c("QX1410",lineages$isotype[!lineages$isotype %in% c("QX1410")])

######### MODIFY THIS ##########
#select and organize how genomes are displayed
desired_strains <- c("QX1410","QG1005","QG2964","ED3102","QG2902","BRC20492","BRC20530","ECA2670","NIC1667", "ECA2666","JU1348")

#the ref genome (QX1410 is omitted from WI_ad as it is handled separately)
#MY681 has been flagged for removal due to a potential mixture during library prepped
#Additional strains can be removed from display (e.g.: lack of genome coverage on adjacent regions to the HDR - anchoring regions to haplotype of interest is not properly mapped)
always_omit_WI_ad <- c("QX1410","MY681")
###############################

WI_ad <- boundGenes %>% 
  dplyr::filter(STRAIN %in% desired_strains & !STRAIN %in% always_omit_WI_ad) %>% #handle selection
  #dplyr::filter(!STRAIN %in% nigonis & !STRAIN %in% always_omit_WI_ad) #if no selection/order is desired, then simply display all briggsae genomes (Figure S16)
  dplyr::left_join(all_ortho_pairs_raw,by=c("STRAIN","tranname")) %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(is.na(has_any_ortho),F,has_any_ortho)) %>%
  dplyr::left_join(all_ortho_pairs_bound %>% dplyr::select(tranname,STRAIN,QX1410,status) %>% dplyr::filter(QX1410 %in% QX_ad$tranname),by=c("STRAIN","tranname")) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(!is.na(status),T,F)) %>%
  dplyr::select(-status,-has_any_ortho) %>% 
  dplyr::rename(QX_name=QX1410) %>%
  dplyr::left_join(aliases %>% dplyr::select(-seqname),by=c("QX_name"="tranname")) %>%
  dplyr::mutate(alias.x=alias.y) %>%
  dplyr::select(-alias.y,-QX_name) %>%
  dplyr::rename(alias=alias.x) %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(HV_boundary %>% dplyr::select(STRAIN,inv),by="STRAIN") %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(n_gene=n_distinct(tranname)) %>%
  dplyr::mutate(start_sort = ifelse(rep(dplyr::first(inv), dplyr::n()), -start, start)) %>%
  dplyr::arrange(start_sort, .by_group = TRUE) %>%
  dplyr::mutate(first_gene=dplyr::first(alias)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::arrange(first_gene, desc(n_gene), .by_group = TRUE) %>% 
  dplyr::mutate(order_gene = dplyr::first(first_gene),
         order_num = dplyr::first(n_gene)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(STRAIN = factor(STRAIN, levels = rev(desired_strains))) %>% #can be commented out if no order is desired
  dplyr::arrange(STRAIN, order_gene, order_num) %>%
  dplyr::select(-order_gene, -order_num) %>%
  dplyr::mutate(g_diff = abs(n_gene-g_count)) %>%
  dplyr::mutate(y_pos=rleid(STRAIN)) %>%
  dplyr::select(-g_diff,-start_sort,-first_gene,-inv,-n_gene)

QX_ad_corr <- QX_ad %>%
  dplyr::mutate(y_pos=max(WI_ad$y_pos)+1)


all_ad <- rbind(QX_ad_corr,WI_ad) %>% 
  dplyr::arrange(STRAIN,start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(gen_pos=rleid(Parent)) %>%
  dplyr::mutate(shift=min(start)) %>%
  dplyr::mutate(end=end-min(start),start=start-min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(has_any_ortho==T & has_bound_ortho ==T,2,ifelse(has_any_ortho==T,1,0))) %>%
  dplyr::left_join(HV_boundary %>% dplyr::select(STRAIN,inv),by="STRAIN") %>%
  dplyr::mutate(inv=ifelse(is.na(inv),F,inv)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(bound_corr=max(end)) %>%
  dplyr::mutate(start=ifelse(inv==T,abs(start-bound_corr),start)) %>%
  dplyr::mutate(end=ifelse(inv==T,abs(end-bound_corr),end))

hlines <- new_boundaries %>% 
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,shift) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN")

segments <- all_ortho_pairs_bound %>%
  dplyr::select(STRAIN,Parent,start,end,QX1410,strand) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::distinct(Parent,QX1410,.keep_all = T) %>%
  dplyr::left_join(QX_tran %>% 
                     dplyr::rename(start_QX=start,end_QX=end,strand_QX=strand,chrom_QX=seqid,QXid=STRAIN) %>% 
                     dplyr::select(tranname,chrom_QX,start_QX,end_QX,strand_QX,alias,seqname,QXid),
                   by=c("QX1410"="tranname")) %>% 
  dplyr::filter(QX1410 %in% QX_ad$tranname) %>%
  #dplyr::filter(chrom_QX==seqid & start_QX > (boundStart-5e4) & end_QX < (boundEnd+5e4)) #%>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::rename(WI_y_pos=y_pos) %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by=c("QXid"="STRAIN")) %>%
  dplyr::rename(QX_y_pos=y_pos) %>%
  dplyr::mutate(QX_shift=min(start_QX),WI_shift=min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(WI_x_pos=(start+((end-start)/2))-WI_shift, QX_x_pos=(start_QX+((end_QX-start_QX)/2)-QX_shift)) %>%
  dplyr::mutate(WI_y_pos=WI_y_pos+0.2,QX_y_pos=QX_y_pos-0.2) %>%
  dplyr::distinct(STRAIN,Parent,seqname,.keep_all = T) %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::mutate(n1=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,seqname) %>%
  dplyr::mutate(n2=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(n1>1 | n2>1,"multi_copy","single_copy")) #%>%
  #dplyr::mutate(QX_y_pos=ifelse(strand_QX=="+",QX_y_pos,QX_y_pos-0.4),WI_y_pos=ifelse(strand=="+",WI_y_pos,WI_y_pos-0.4))

plot_ad <- all_ad %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::filter(col==max(col)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::mutate(class=ifelse(col==0,"no_known_ortho",ifelse(col==1,"has_distal_ortho","has_local_ortho"))) 

all_hap <- ggplot() +
  geom_segment(data=hlines,aes(x=boundStart-shift,xend=boundEnd-shift,y=y_pos,yend=y_pos))+
  #geom_segment(data=segments,aes(x=WI_x_pos,xend=N2_x_pos,y=WI_y_pos,yend=N2_y_pos,col=col)) +
  geom_rect(data=plot_ad, aes(xmin=start,xmax=end,ymin=y_pos+0.2,ymax=y_pos-0.2,fill=class),color="black") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(),
        axis.ticks = element_blank()) +
  #scale_color_manual(values=c("multi_copy"="black","single_copy"="grey")) +
  scale_fill_manual(values=c("has_local_ortho"="grey","has_distal_ortho"="black","no_known_ortho"="red","no_known_allelic_CB"="blue")) +
  scale_x_continuous(expand = c(0.01,0)) +
  ylab("")#+
#geom_segment(data=tigs,aes(x=S2_adj,xend=E2_adj,y=0.7,yend=0.7,col=INV))
all_hap

plot_ad <- all_ad %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::filter(col==max(col)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::mutate(class=ifelse(col==0,"no_known_ortho",ifelse(col==1,"has_distal_ortho","has_local_ortho"))) 

# Exclude "non-ortho" from the trapezium joining
plot_ad_filtered <- plot_ad %>% 
  dplyr::mutate(alias=ifelse(is.na(alias),"non-ortho",alias)) %>%
  filter(alias != "non-ortho")

# Join filtered data frames for many-to-many connections
trapeziums <- inner_join(
  plot_ad_filtered, plot_ad_filtered,
  by = "alias",
  suffix = c("_upper", "_lower"),
  relationship = "many-to-many"
) %>% 
  filter(y_pos_upper - y_pos_lower == 1)

# Create trapezium polygons using min/max for x-coordinates so that start/end orientation is corrected.
trapezium_polys <- trapeziums %>% 
  dplyr::rowwise() %>%
  do({
    # Calculate corrected x coordinates for the upper rectangle
    x_left_upper <- min(.$start_upper, .$end_upper)
    x_right_upper <- max(.$start_upper, .$end_upper)
    
    # Calculate corrected x coordinates for the lower rectangle
    x_left_lower <- min(.$start_lower, .$end_lower)
    x_right_lower <- max(.$start_lower, .$end_lower)
    
    data.frame(
      alias = .$alias,
      group = paste(.$alias, .$y_pos_upper, sep = "_"),
      x = c(x_left_upper, x_right_upper, x_right_lower, x_left_lower),
      y = c(.$y_pos_upper - 0.2,  # bottom edge of the upper rectangle
            .$y_pos_upper - 0.2,
            .$y_pos_lower + 0.2,  # top edge of the lower rectangle
            .$y_pos_lower + 0.2)
    )
  }) %>%
  dplyr::ungroup()

# Extract unique aliases at y_pos max+1 in order of increasing start position
ordered_aliases <- plot_ad %>%
  filter(y_pos == max(plot_ad$y_pos)) %>%
  arrange(start) %>%
  pull(alias) %>%
  unique()

# Reorder the factor levels so that the legend follows the ordered aliases
plot_ad <- plot_ad %>%
  mutate(alias = factor(alias, levels = ordered_aliases))

# Also update any other data frames with alias info, e.g. trapezium_polys:
trapezium_polys <- trapezium_polys %>%
  mutate(alias = factor(alias, levels = ordered_aliases))

# Shuffle the assignment of colors to the ordered aliases
set.seed(123)  # for reproducibility
shuffled_aliases <- sample(ordered_aliases)

# Generate colors using hcl.colors() for the shuffled aliases
default_colors <- setNames(hcl.colors(length(shuffled_aliases), "Dark 3"), shuffled_aliases)

# But to keep the legend order as ordered_aliases, we re-map these colors back:
final_colors <- default_colors[ordered_aliases]

# Optionally, if you have the "non-ortho" alias (or any other), add it explicitly:
final_colors <- c(final_colors, "Unknown gene" = "darkgrey")


plot_ad_segments <- plot_ad %>%
  dplyr::mutate(
    # Adjust strand logic if inverted
    strand_logic = case_when(
      strand == "+" & !inv ~ "+",
      strand == "-" & !inv ~ "-",
      strand == "+" & inv  ~ "-",
      strand == "-" & inv  ~ "+"
    ),
    seg_color = ifelse(strand_logic == "+", "black", "red"),
    
    x_start = start,
    x_end   = end,
    y_seg   = y_pos - 0.25  # just under the geom_rect (geom_rect is y_pos ± 0.2)
  )
#apply final_colors in your ggplot scale:
all_hap_bg <- ggplot() +
  geom_segment(data = hlines, 
               aes(x = boundStart - shift, xend = boundEnd - shift, y = y_pos, yend = y_pos)) +
  geom_polygon(data = trapezium_polys, 
               aes(x = x, y = y, group = group, fill = alias)) +
  geom_rect(data = plot_ad %>% dplyr::mutate(alias=ifelse(is.na(alias),"Unknown gene",as.character(alias))),
            aes(xmin = start, xmax = end, ymin = y_pos + 0.2, ymax = y_pos - 0.2, fill = alias),color = "black") +
  #geom_segment(
  #  data = plot_ad_segments,
  #  aes(x = x_start, xend = x_end, y = y_seg, yend = y_seg, color = seg_color),
  #  linewidth = 0.5,
  #  inherit.aes = FALSE
  #) + # orientation of the genes can be displayed with this geom_segment() call
  scale_y_continuous(
    expand = c(0.01, 0),
    breaks = hlines$y_pos,
    labels = hlines$STRAIN
  ) +
  scale_x_continuous(expand = c(0.01, 0),labels = function(x) x / 1000) +
  scale_fill_manual(values = final_colors,
                    breaks = names(final_colors)) +
  scale_color_identity()  +
  #scale_color_manual(values = c("+"="black","-"="red")) +
  labs(fill="Reference\ngene")+
  xlab("Physical distance (kb)") +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8),  # adjust size as needed
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(),
    axis.title.x = element_text(),
    legend.position='none'
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    # legend.key.size = unit(0.4, "lines"),
    # legend.text = element_text(size = 6),
    # legend.title = element_text(size = 7) # legend can be enabled with these theme() paramters
  ) +
  guides(fill = guide_legend(nrow = 6, byrow = TRUE))

#save png
#used for Figure S16
#outfile_png <- paste0("../../processed_data/gene_diversity/HDR_",hdr_chrom,"_",hdr_start_pos,"_",hdr_end_pos,"_",length(desired_strains),"rg",".png")
#ggsave(plot = all_hap_bg,filename = outfile_png,width = 7.5,height = 8.5,device = "png",units = "in",dpi = 600,bg = "white")

#store ggplot in Rds for concatenating with other plots
outfile_rd <- paste0("../../processed_data/gene_diversity/HDR_",hdr_chrom,"_",hdr_start_pos,"_",hdr_end_pos,"_",length(desired_strains),"rg",".Rds")
saveRDS(all_hap_bg, file=outfile_rd)
