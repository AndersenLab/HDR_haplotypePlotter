library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)

setwd("/vast/eande106/projects/Nicolas/github/HDR_haplotypePlotter/")

#read collapsed reference HDRs
collapsed_ff <- readr::read_tsv("./input/HDR_5kbclust_collapsed_wFreq.tsv")

#read strain-specific HDRs
all_SR_calls <- readr::read_tsv("./input/HDR_allStrain_5kbclust_1IBfilt.tsv")
#all_LR_calls <- readr::read_tsv("/projects/b1059/projects/Nicolas/hyperdivergent_regions/elegans/HDR_LRcalls_95idy.tsv")

#read all pairwise genome coordinate comparisons
transformed_coords <- readr::read_tsv("./input/all_hifi_nucmer_CE.coords",col_names = F) 
colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")

#read concatentated gene models of every genome
gffCat <- ape::read.gff("./input/all_LRiso.gff")

#read ortholog relationships among gene models
orthos <- readr::read_tsv("./input/Orthogroups.tsv")
strainCol <- colnames(orthos)

#set your target HDR coordinates
#SEA-1 from Lee et al. 2021 is currently displayed
#these could be arguments in the future
hap_chrom = "II"
hap_start = 3667179
hap_end = 3701405

#some other examples that were tried below
# 6 HAPLOT
# hap_chrom = "II"
# hap_start = 3271510
# hap_end = 3415738

# 7 HAPLOT
# hap_chrom = "V"      
# hap_start = 20193463
# hap_end = 20267244 

#use reference coordinates from g2g alginments to pull the contigs that contain the alt haplotypes for the HDR
hap_coords <- transformed_coords %>%
  dplyr::filter((REF == hap_chrom & hap_start >= S1 & hap_start <= E1 ) | 
                  (REF == hap_chrom & hap_end >= S1 & hap_end <= E1) | 
                  (REF == hap_chrom & S1 >= hap_start & E1 <= hap_end)) %>%
  dplyr::mutate(inv=ifelse(S2>E2,T,F)) 
  # dplyr::select(-S2,-E2) %>%
  # dplyr::rename(S2=newS2,E2=newE2)

# naive visualization of g2g alignments for the target region
# multiple contigs may map to the REF region, we need to filter those!
ggplot(hap_coords) + geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=inv)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("REF genome position (Mb)") +
  ylab("WILD contig position (Mb)")

#keep only the contig with the largest extent of alignment with the REF HDR
tigFilt <- hap_coords %>% 
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HIFI) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L2)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::ungroup()

#important diagnostic plot
#after the optimal contig is selected, we can see if the HDR boundaries are within the alignment
#for proper visualization, the HDR (grey box in plot) needs to be encompassed by the selected contig
#otherwise some haplotypes will be truncated 
#a step to drop genomes with incomplete coverage of the HDR could be added
#SEA-1 locus behaves well
#have in mind that the alignments look fragmented because of sequence divergence, but the contig is linear in the genome file
ggplot(tigFilt) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("REF genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA))

#get the minimum and maximum boundary of the WILD genome alignments that contain the HDR
HV_boundary <- tigFilt %>%
  dplyr::mutate(refStart=min(S1),refEnd=max(E1)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(boundStart=min(S2), boundEnd=max(E2)) %>%
  dplyr::distinct(STRAIN, .keep_all = T) %>%
  dplyr::select(HIFI,boundStart,boundEnd,STRAIN,REF,refStart,refEnd) %>%
  dplyr::ungroup()

#filter the concatenated GFF to extract the gene models of each WILD genome contig boundary
filtGff <- gffCat %>%
  dplyr::filter(type=="gene") %>%
  dplyr::filter(seqid %in% HV_boundary$HIFI) %>%
  dplyr::left_join(HV_boundary,by=c('seqid'='HIFI')) %>%
  dplyr::filter((start >= boundStart & start <= boundEnd) | 
                  (end >= boundStart & end <= boundEnd)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ngene=n()) %>%
  dplyr::arrange(ngene) %>%
  dplyr::mutate(gid=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(boundSize=abs(boundStart-boundEnd))

#extract the REF genes 
N2Start = filtGff$refStart[1]
N2End = filtGff$refEnd[1]
N2_genes <- gffCat %>%
  dplyr::filter(seqid==hap_chrom & type=="gene") %>%
  dplyr::mutate(refStart=N2Start,refEnd=N2End) %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | 
                  (end >= hap_start & end <= hap_end)) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=';sequence_name=') %>%
  tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>%
  dplyr::mutate(seqname=paste0("N2.",seqname)) %>%
  tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("Name=","",Name))  

#extract the REF protein-coding transcripts
N2_tran <- gffCat %>%
  dplyr::filter(seqid==hap_chrom & type=="mRNA") %>%
  tidyr::separate(attributes, into=c("ID","Parent","Name","wormpep","locus"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=Transcript:","",ID)) %>%
  tidyr::separate(ID,into = c("fosmid","tseqID",'tranum'),sep="\\.",remove = F) %>%
  dplyr::mutate(tseqname=paste0("N2.",fosmid,".",tseqID)) %>%
  dplyr::mutate(Parent=gsub("Parent=Gene:","",Parent)) %>%
  dplyr::filter(Parent %in% N2_genes$Name) %>%
  dplyr::select(tseqname,Parent) %>% 
  dplyr::rename(tranname=tseqname) %>%
  dplyr::left_join(N2_genes,by=c('Parent'='Name'))

#get gene list
HV_genelist <- N2_genes$seqname

#get alt gene names/aliases
aliases <- N2_genes %>%
  dplyr::select(seqname,post2) %>%
  tidyr::separate(post2,into=c("rem","aliases"),sep=";Alias=") %>%
  dplyr::select(-rem) %>%
  tidyr::separate(aliases,into=c('locus_name',"alias2"),sep=',') %>%
  dplyr::select(-alias2)

#minor diagnostic plot to visualize the REF loci captured by the HDR
#this is your REF haplotype
ggplot(N2_genes) + geom_rect(aes(xmin=start,xmax=end,ymin=1,ymax=2))


#filter orthologous groups using REF genes
#this will establish your orthology relationships between REF and WILD haplotypes
filtOrthos <- orthos %>%
  dplyr::filter(grepl(paste(HV_genelist,collapse="|"),N2.c_elegans.PRJNA13758.WS270.protein.fa_longest_isoforms)) %>%
  dplyr::mutate(N2 = strsplit(as.character(N2.c_elegans.PRJNA13758.WS270.protein.fa_longest_isoforms), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end),by=c("N2"="tranname")) %>%
  dplyr::filter(!is.na(seqid))

#generate a lookup table (all_ortho_pairs) which contains all pairwise gene orthologs between REF and WILD
orthoList <- list()
for (i in 2:length(strainCol)) {
  tmp <- filtOrthos %>%
    dplyr::select(N2,strainCol[i]) %>%
    dplyr::rename(tmpSel=strainCol[i]) %>%
    dplyr::mutate(newSel = strsplit(as.character(tmpSel), ",")) %>%
    tidyr::unnest(newSel) %>%
    dplyr::mutate(newSel=trimws(newSel)) %>%
    dplyr::select(-tmpSel) %>%
    tidyr::separate(newSel,into = c("STRAIN","geneid",'tranid'),sep="\\.")
  
  orthoList[[i-1]] <- tmp
}
all_ortho_pairs  <- ldply(orthoList,data.frame) %>%
  dplyr::mutate(geneid=ifelse(STRAIN=='N2',paste0(geneid,".",tranid),geneid))

#prep GFF fields for join with ortholog table
gffFilt_ortho <- gffCat %>%
  dplyr::filter((!source=="WormBase") & type=="gene") %>%
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::rename(geneid=attributes) %>%
  tidyr::separate(seqid,into = c("STRAIN","tigID1",'tigID2'),sep = '\\.',remove = F)

# #join ortholog table to GFF coordinates
#seqid delineates the contig (thus the wild strain), which sets the Y position in the haplotype plot
#keep the longest transcript per gene
orthoPairs_wcoord <- all_ortho_pairs %>%
  dplyr::left_join(gffFilt_ortho,by = c("STRAIN","geneid")) %>%
  dplyr::arrange(seqid) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(ypos=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(gsize==max(gsize)) %>%
  dplyr::ungroup()

#similar to DF above, but
#some contigs are inverted, and gene models need to be adjusted for orientation
#some contigs have partial alignments that contain genes further away from the local HDR
#using lead() we can remove genes that have large (>5e4) jumps in position
#this trimming process may be the cause of misbehavior in plots of other HDRs
orthoPairs_wcoord_trimmer <- all_ortho_pairs %>%
  dplyr::left_join(gffFilt_ortho,by = c("STRAIN","geneid")) %>%
  dplyr::arrange(seqid) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(ypos=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(gsize==max(gsize)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(hap_coords %>%
                     dplyr::select(HIFI,inv) %>%
                     dplyr::distinct(HIFI,.keep_all=T),by=c("seqid"="HIFI")) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(newStart=ifelse(inv==T,abs(end-(max(end))),start),newEnd=ifelse(inv==T,abs(start-(max(end))),end)) %>%
  #dplyr::select(-start,-end) %>%
  #dplyr::rename(start=newStart,end=newEnd) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(seqid,newStart) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDist=abs(newEnd-lead(newStart))) %>%
  dplyr::mutate(trm=ifelse(leadDist>5e4 | lead(leadDist) >5e4,"R","NR")) %>% #this is the problematic line
  dplyr::filter(is.na(trm) | trm=="NR")

#pull the boundaries of the ortholgous genes
pullBound <- orthoPairs_wcoord_trimmer %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(boundStart = min(start), boundEnd = max(end)) %>%
  dplyr::distinct(seqid,.keep_all = T) %>%
  dplyr::select(STRAIN,seqid,boundStart,boundEnd) %>%
  dplyr::ungroup()

#find the bound genes for each strain that are not orthologous
boundGenes <- gffCat %>%
  dplyr::filter(type=="gene" & seqid %in% pullBound$seqid) %>%
  dplyr::left_join(pullBound,by="seqid") %>%
  dplyr::group_by(seqid) %>%
  dplyr::filter(start >= boundStart & start <= boundEnd) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(geneid=gsub("ID=","",attributes)) %>%
  dplyr::left_join(orthoPairs_wcoord,by=c('geneid','STRAIN')) %>%
  dplyr::filter(is.na(N2)) %>%
  dplyr::select(seqid.x, start.x, end.x, geneid) %>%
  dplyr::rename(seqid=seqid.x,start=start.x,end=end.x) %>%
  dplyr::mutate(N2="non-ortho")


#bind ortholgous and non orthologous genes for each strain, and transform the coordinates of genes to a midle axis center
plotCoords <- rbind(boundGenes,orthoPairs_wcoord %>% dplyr::select(seqid,start,end,geneid,N2)) %>%
  dplyr::mutate(ortho_status=ifelse(N2=='non-ortho',F,T)) %>%
  dplyr::left_join(hap_coords %>%
                     dplyr::select(HIFI,inv) %>%
                     dplyr::distinct(HIFI,.keep_all=T),by=c("seqid"="HIFI")) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(newStart=ifelse(inv==T,abs(end-(max(end))),start),newEnd=ifelse(inv==T,abs(start-(max(end))),end)) %>%
  dplyr::select(-start,-end) %>%
  dplyr::rename(start=newStart,end=newEnd) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(ortho_status,seqid,start) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(leadDist=abs(end-lead(start))) %>%
  dplyr::mutate(trm=ifelse(leadDist>5e4 | lag(leadDist) >5e4,"R","NR")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.na(trm) | trm=="NR") %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(newstart=start-min(start),newend=end-min(start)) %>%
  dplyr::ungroup() %>%
  tidyr::separate(seqid,into = c("STRAIN","tig1","tig2"),sep = "\\.",remove = F) %>%
  dplyr::select(-tig1,-tig2,-geneid,-inv,-start,-end,-trm,-leadDist)

#get N2 genes
N2ad <- N2_genes %>%
  dplyr::mutate(STRAIN="N2",ortho_status=T) %>%
  dplyr::rename(N2=seqname) %>%
  dplyr::mutate(shift=min(start)) %>%
  dplyr::mutate(newstart=start-shift,newend=end-shift) %>%
  dplyr::select(seqid,STRAIN,N2,ortho_status,newstart,newend,shift) 

#store the N2 coordinate shift
N2_shift <- c("N2",unique(N2ad$shift))
#remove it from DF
N2ad <- N2ad %>% dplyr::select(-shift)

# get SR HVR calls
hd_reg <- all_SR_calls %>%
  dplyr::filter(CHROM==hap_chrom) %>%
  dplyr::filter((minStart >= hap_start & minStart<= hap_end) | (maxEnd<= hap_end & maxEnd >= hap_start)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(divSize==max(divSize))

#get reference-collapsed calls
hd_collapse <- collapsed_ff %>%
  dplyr::filter(CHROM==hap_chrom) %>%
  dplyr::filter((minStart >= hap_start & minStart<= hap_end) | (maxEnd<= hap_end & maxEnd >= hap_start)) %>%
  dplyr::mutate(minStart=minStart-as.numeric(N2_shift[2]),maxEnd=maxEnd-as.numeric(N2_shift[2])) %>%
  dplyr::mutate(STRAIN="N2")



########################### PLOT ALL POSSIBLE HAP ##############################

#gene positions
plotCoords2 <- rbind(N2ad,plotCoords) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(minStart=min(newstart),maxEnd=max(newend)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(midpoint=(minStart+maxEnd)/2) %>%
  dplyr::mutate(centeredStart=newstart-midpoint,centeredEnd=newend-midpoint) %>%
  dplyr::filter(!is.na(seqid))

#center line positions
hlines <- plotCoords2 %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(hlineStart=min(centeredStart)-1000,
                hlineEnd=max(centeredEnd)+1000) %>%
  dplyr::select(STRAIN,hlineStart,hlineEnd)

#lab positions
labs <- plotCoords2 %>%
  dplyr::mutate(labStart=min(centeredStart)+0.2*(min(centeredStart))) %>%
  dplyr::select(STRAIN,labStart,seqid) %>%
  dplyr::distinct(seqid,.keep_all = T) 

#plot
allhap<- ggplot() +
  geom_segment(data=hlines,aes(x=hlineStart,xend=hlineEnd,y=1,yend=1)) +
  geom_rect(data=plotCoords2 %>% dplyr::filter(ortho_status==T),aes(xmin=centeredStart,xmax=centeredEnd,ymin=1-0.5,ymax=1+0.5,fill=N2,color="black")) +
  geom_rect(data=plotCoords2 %>% dplyr::filter(ortho_status==F),aes(xmin=centeredStart,xmax=centeredEnd,ymin=1-0.5,ymax=1+0.5,color='non-ortho')) +
  geom_text(data=labs,aes(x=labStart,y=1,label=STRAIN)) +
  #facet_wrap(~factor(STRAIN, levels=c('N2','NIC2','ECA36','ECA396', 'MY2693', 'EG4725','JU2600','JU310','MY2147','JU1400','NIC526','JU2526','QX1794','CB4856','XZ1516','DL238')),ncol=1) +
  facet_wrap(~factor(STRAIN, levels=c('N2','NIC2','ECA36' ,'DL238',"CB4856",'EG4725','MY2147','NIC526','JU2600','JU310','JU1400','XZ1516','JU2526',"QX1794","MY2693","ECA396")),ncol=1) +
  theme(strip.text = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank()) +
  scale_color_manual(values = c('non-ortho'='darkgrey')) 
allhap
########################### NEW 6HAP LOCUS ##############################
# 
#   
# plotCoords3 <- rbind(N2ad %>% dplyr::filter(!(N2=="N2.F02E11.3")),plotCoords) %>%
#     dplyr::group_by(STRAIN) %>%
#     dplyr::mutate(minStart=min(newstart),maxEnd=max(newend)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(midpoint=(minStart+maxEnd)/2) %>%
#     dplyr::mutate(centeredStart=newstart-midpoint,centeredEnd=newend-midpoint) %>%
#     dplyr::filter(!is.na(seqid)) %>%
#     dplyr::filter(STRAIN =="N2" | STRAIN == "CB4856" | STRAIN=="DL238"| STRAIN=="XZ1516"| STRAIN=="QX1794"| STRAIN=="MY2693")
# 
# midpoints <- plotCoords3 %>% 
#   dplyr::select(STRAIN,midpoint) %>%
#   dplyr::distinct(STRAIN,.keep_all = T)
# 
# hlines <- plotCoords3 %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(hlineStart=min(centeredStart)-1000,
#                 hlineEnd=max(centeredEnd)+1000) %>%
#   dplyr::select(STRAIN,hlineStart,hlineEnd)
# 
# labs <- plotCoords3 %>%
#   dplyr::mutate(labStart=min(centeredStart)+0.2*(min(centeredStart))) %>%
#   dplyr::select(STRAIN,labStart,seqid) %>%
#   dplyr::distinct(seqid,.keep_all = T) 
# 
# ggplot() +
#   geom_rect(data=hd_collapse,aes(xmin=minStart,xmax=maxEnd,ymin=1-0.5,ymax=1+0.5,color='lightgrey',alpha=0.2)) +
#   geom_segment(data=hlines,aes(x=hlineStart,xend=hlineEnd,y=1,yend=1)) +
#   geom_rect(data=plotCoords3 %>% dplyr::filter(ortho_status==T),aes(xmin=centeredStart,xmax=centeredEnd,ymin=1-0.5,ymax=1+0.5,fill=N2,color="black")) +
#   geom_rect(data=plotCoords3 %>% dplyr::filter(ortho_status==F),aes(xmin=centeredStart,xmax=centeredEnd,ymin=1-0.5,ymax=1+0.5,color='non-ortho')) +
#   geom_text(data=labs,aes(x=labStart,y=1,label=STRAIN)) +
#   facet_wrap(~factor(STRAIN, levels=c('N2','CB4856','XZ1516','DL238','QX1794','MY2693')),ncol=1) +
#   theme(strip.text = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         panel.background = element_blank()) +
#   scale_color_manual(values = c('non-ortho'='darkgrey'))

############ PLOT SEA-1 3HAP LOCUS ##############
#here I have previous knowledge of the haplotype strycture of this region
#I've collapsed the plot to the coordinates of the three core haplotypes

plotCoords4 <- rbind(N2ad,plotCoords) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(minStart=min(newstart),maxEnd=max(newend)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(midpoint=(minStart+maxEnd)/2) %>%
  dplyr::mutate(centeredStart=newstart-midpoint,centeredEnd=newend-midpoint) %>%
  dplyr::filter(!is.na(seqid)) %>%
  dplyr::filter(STRAIN =="N2" | STRAIN == "CB4856" | STRAIN=="DL238") %>%
  dplyr::left_join(aliases,by=c("N2"="seqname"))

midpoints <- plotCoords4 %>% 
  dplyr::select(STRAIN,midpoint) %>%
  dplyr::distinct(STRAIN,.keep_all = T)

hlines <- plotCoords4 %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(hlineStart=min(centeredStart)-1000,
                hlineEnd=max(centeredEnd)+1000) %>%
  dplyr::select(STRAIN,hlineStart,hlineEnd)

labs <- plotCoords4 %>%
  dplyr::mutate(labStart=min(centeredStart)+0.2*(min(centeredStart))) %>%
  dplyr::select(STRAIN,labStart,seqid) %>%
  dplyr::distinct(seqid,.keep_all = T) 


# test<- data.frame(STRAIN=c("N2","CB4856","DL238"),newLab=c("N2\nNIC2\nECA36\nECA396",
#                                                            "QX1794\tJU1400\nMY2693\tJU2526\nJU2600\tJU310\nMY2147\tNIC526\nCB4856\tEG4725\nXZ1516\t            \n",
#                                                            "DL238"),newLab2=c(NA,"MY2147\nNIC526\nCB4856\nEG4725\nXZ1516\n",NA))
#test<- data.frame(STRAIN=c("N2","CB4856","DL238"),newLab=c("N2\nNIC2\nECA36\nECA396",
#                                                           "QX1794\nJU1400\nMY2693\nJU2526\nJU2600\nJU310",
#                                                           "DL238"),newLab2=c(NA,"MY2147\nNIC526\nCB4856\nEG4725\nXZ1516\n",NA))

# labs <- labs %>%
#   dplyr::left_join(test,by="STRAIN")

regDef_WI <- plotCoords4 %>%
  dplyr::filter((!STRAIN=="N2") & (N2=="N2.F19B10.9" | N2=="N2.F40H7.5")) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(regStart=min(centeredStart),regEnd=max(centeredEnd)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,.keep_all = T) %>%
  dplyr::select(STRAIN,regStart,regEnd)

hap3 <- ggplot() +
  geom_rect(data=hd_collapse,aes(xmin=minStart-as.numeric(midpoints[1,2]),xmax=maxEnd-as.numeric(midpoints[1,2]),ymin=1-0.5,ymax=1+0.5),fill='lightgrey')+
  geom_rect(data=regDef_WI,aes(xmin=regStart,xmax=regEnd,ymin=1-0.5,ymax=1+0.5),fill='lightgrey')+
  geom_segment(data=hlines,aes(x=hlineStart,xend=hlineEnd,y=1,yend=1)) +
  geom_rect(data=plotCoords4 %>% dplyr::filter(ortho_status==T),aes(xmin=centeredStart,xmax=centeredEnd,ymin=1-0.5,ymax=1+0.5,fill=locus_name),color="black") +
  geom_rect(data=plotCoords4 %>% dplyr::filter(ortho_status==F),aes(xmin=centeredStart,xmax=centeredEnd,ymin=1-0.5,ymax=1+0.5),fill='darkgrey',color="black") +
  #geom_text(data=labs %>% dplyr::filter(!STRAIN=="CB4856"),aes(x=labStart-2000,y=1,label=newLab),size=2.2,hjust=0.5,fontface='bold') +
  #geom_text(data=labs %>% dplyr::filter(STRAIN=="CB4856"),aes(x=labStart+3000,y=1,label=newLab),size=2.2,hjust=1,fontface='bold') +
  #geom_text(data=labs %>% dplyr::filter(STRAIN=="CB4856"),aes(x=labStart-2500,y=1,label=newLab2),size=2.2,hjust=1,fontface='bold') +
  facet_wrap(~factor(STRAIN, levels=c('N2','CB4856','XZ1516','DL238','QX1794','MY2693')),ncol=1) +
  theme(strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7,)) +
  #scale_color_manual(values = c('non-ortho'='darkgrey')) +
  #xlim(min(labs$labStart)-5000,max(hlines$hlineEnd)) +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  labs(fill='Locus name') 
hap3
       
reg <- paste0(hap_chrom,".",hap_start,"-",hap_end)
dir <-getwd()
f1 <- paste0(dir,"/figs/","allhap","_",reg,".png")
f2 <- paste0(dir,"/figs/","hap3","_",reg,".png")
ggsave(f1,allhap,height = 4.5,width = 8.5,units = 'in',device = 'png',dpi = 900) 
ggsave(f2,hap3,height = 3.5,width = 6.5,units = 'in',device = 'png',dpi = 900) 



# ggplot(orthoPairs_wcoord_trimmer) +
#   geom_rect(aes(xmin=start,xmax=end,ymin=1-0.5,ymax=1+0.5,fill=N2)) +
#   facet_wrap(~seqid,scales = 'free_x',ncol=1) +
#   theme(strip.text = element_blank())

