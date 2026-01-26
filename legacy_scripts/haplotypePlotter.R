library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)
library(data.table)
library(stringr)

#setwd("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/glc1_variation/HDR_haplotypePlotter")

#read collapsed reference HDRs
collapsed_ff <- readr::read_tsv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/elegans/HDR_5kbclust_collapsed.tsv")

#read strain-specific HDRs
all_SR_calls <- readr::read_tsv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/elegans/HDR_allStrain_5kbclust_1IBfilt.tsv")
#all_LR_calls <- readr::read_tsv("/projects/b1059/projects/Nicolas/hyperdivergent_regions/elegans/HDR_LRcalls_95idy.tsv")

#read all pairwise genome coordinate comparisons
transformed_coords <- readr::read_tsv("/vast/eande106/projects/Nicolas/c.elegans/reference_genealn/N2vCB/N2_hifi_transformed2.tsv",col_names = F) %>% dplyr::mutate(STRAIN="CB4856")
#transformed_coords <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/asssemblies/nucmer_runs/all_WI_transformed.tsv",col_names = F) 
colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")

#read concatentated gene models of every genome
gffCat1 <- ape::read.gff("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/braker_runs/gff/CB4856.braker.gff3") %>% dplyr::mutate(STRAIN="CB4856")
#gffCat1 <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/braker_runs/merged_gff/all_WI_braker.clean.gff", col_names = F)
#colnames(gffCat1) <- c("seqid","source","type","start","end","score","strand","phase","attributes","STRAIN")
gffCat2 <- ape::read.gff("/vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3") %>% dplyr::mutate(STRAIN="N2")
gffCat <- rbind(gffCat1,gffCat2)

#read ortholog relationships among gene models
#orthos <- readr::read_tsv("./input/Orthogroups.tsv")
orthos <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/N2xCB/OrthoFinder/Results_Mar17/Orthogroups/Orthogroups.tsv")
#orthos <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/prot/OrthoFinder/Results_Mar20/Orthogroups/Orthogroups.tsv")
strainCol <- colnames(orthos)
strainCol_c1 <- gsub(".braker.protein","",strainCol)
strainCol_c2 <- gsub("_WS283.protein","",strainCol_c1)
colnames(orthos) <- strainCol_c2
#orthos <- orthos %>%dplyr::select(Orthogroup,CB4856,N2)
#strainCol_c2 <-colnames(orthos)

#set your target HDR coordinates
#SEA-1 from Lee et al. 2021 is currently displayed
# #these could be arguments in the future
# hap_chrom = "II"
# hap_start = 3667179
# hap_end = 3701405

#GLC-1
hdr_chrom = "II"
hdr_start_pos = 2365967
hdr_end_pos = 2669007

#offset lets you explore adjacent (non-HDR sequences) - set to 0 if not needed
offset = 0
hap_chrom = hdr_chrom
hap_start = hdr_start_pos - offset
hap_end = hdr_end_pos + offset


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
  #geom_rect(xmin=16219634/1e6,xmax=16221917/1e6,ymin=-Inf,ymax=Inf,aes(fill="glc-1")) +
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') #+
  #scale_fill_manual(values=c("lightgrey"="lightgrey","glc-1"="pink"))


tigFilt2 <- tigFilt %>%
  dplyr::arrange(S2) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDiff=lead(S2)-E2) %>%
  dplyr::mutate(lagDiff=E2-lag(S2)) %>%
  dplyr::mutate(lagDiff=ifelse(is.na(lagDiff),0,lagDiff)) %>%
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(drop=ifelse(lagDiff >1.5e5 & lag(leadDiff) > 1.5e5,T,F)) %>% 
  dplyr::mutate(droprle=rleid(drop)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,droprle) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(E2-S2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize)

ggplot(tigFilt2) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  #geom_rect(xmin=16219634/1e6,xmax=16221917/1e6,ymin=-Inf,ymax=Inf,aes(fill="glc-1")) +
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

#get the minimum and maximum boundary of the WILD genome alignments that contain the HDR
HV_boundary <- tigFilt2 %>%
  dplyr::select(-lagDiff,-leadDiff,-drop) %>%
  dplyr::mutate(refStart=min(S1,E1),refEnd=max(S1,E1)) %>%
  dplyr::mutate(boundStart=min(S2,E2), boundEnd=max(S2,E2)) %>%
  dplyr::distinct(STRAIN, .keep_all = T) %>%
  dplyr::select(HIFI,boundStart,boundEnd,STRAIN,REF,refStart,refEnd) %>%
  dplyr::ungroup() %>%
  dplyr::rename(boundChrom=HIFI)




#filter the concatenated GFF to extract the gene models of each WILD genome contig boundary
wild_genes <- gffCat %>%
  dplyr::filter(type=="gene" & !(STRAIN=="N2")) %>%
  dplyr::mutate(attributes=gsub(";","",attributes)) %>% 
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::select(attributes,seqid,start,end,strand,STRAIN) %>%
  dplyr::rename(Name=attributes)
  
wild_tran <-  gffCat  %>%
  dplyr::filter(type=="mRNA" & !(STRAIN=="N2")) %>%
  tidyr::separate(attributes,into=c("tranname","Parent"),sep=";Parent=") %>%
  dplyr::mutate(Parent=gsub(";","",Parent)) %>%
  dplyr::mutate(tranname=gsub("ID=","",tranname)) %>%
  dplyr::select(tranname,Parent,STRAIN) %>%
  dplyr::left_join(wild_genes,by=c("Parent"="Name","STRAIN")) %>%
  dplyr::left_join(HV_boundary,by="STRAIN")
  
#extract the REF genes 
N2Start = HV_boundary$refStart[1]
N2End = HV_boundary$refEnd[1]
N2_genes <- gffCat %>%
  dplyr::filter(type=="gene" & STRAIN=="N2") %>%
  dplyr::mutate(refStart=N2Start,refEnd=N2End) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=';sequence_name=') %>%
  tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>%
  #dplyr::mutate(seqname=paste0(seqname)) %>%
  tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("Name=","",Name)) %>% 
  dplyr::select(seqid,start,end,strand,Name,rest2,seqname,STRAIN,refStart,refEnd,seqname) %>%
  dplyr::mutate(rest2=ifelse(grepl("locus",rest2),gsub("locus=","",rest2),seqname)) %>%
  dplyr::rename(alias=rest2)

#extract the REF protein-coding transcripts
N2_tran <- gffCat %>%
  dplyr::filter(type=="mRNA" & STRAIN=="N2") %>%
  tidyr::separate(attributes, into=c("ID","Parent","Name","wormpep","locus"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=Transcript:","",ID)) %>%
  tidyr::separate(ID,into = c("fosmid","tseqID",'tranum'),sep="\\.",remove = F) %>%
  dplyr::mutate(tseqname=paste0(fosmid,".",tseqID,".",tranum)) %>%
  dplyr::mutate(Parent=gsub("Parent=Gene:","",Parent)) %>%
  dplyr::filter(Parent %in% N2_genes$Name) %>%
  dplyr::select(tseqname,Parent) %>%
  dplyr::rename(tranname=tseqname) %>%
  dplyr::mutate(tranname=paste0("Transcript_",tranname)) %>%
  dplyr::left_join(N2_genes,by=c('Parent'='Name'))


##
# dplyr::filter((start >= hap_start & start <= hap_end) | 
#                 (end >= hap_start & end <= hap_end)) %>%
##
N2_tran_reg <- N2_tran %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | (end >= hap_start & end <= hap_end))  %>%
  dplyr::filter(seqid==hap_chrom)


# #get gene list
HV_genelist <- N2_tran_reg$tranname 
#get alt gene names/aliases
aliases <- N2_tran %>% dplyr::select(seqname,tranname,alias)

#minor diagnostic plot to visualize the REF loci captured by the HDR
#this is your REF haplotype
ggplot(N2_tran_reg) + geom_rect(aes(xmin=start,xmax=end,ymin=1,ymax=2))


#filter orthologous groups using REF genes
#this will establish your orthology relationships between REF and WILD haplotypes

all_orthos_unnest <- orthos %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::mutate(na_count = rowSums(is.na(.))) %>%
  dplyr::filter(na_count < length(strainCol_c2) - 2) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end),by=c("N2"="tranname"))  %>%
  dplyr::select(-na_count)

filtOrthos <- orthos %>%
  dplyr::filter(grepl(paste(HV_genelist,collapse="|"),N2)) %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::left_join(N2_tran_reg %>% dplyr::select(tranname,seqid,seqname,start,end) %>% dplyr::mutate(og_loc="in_region"),by=c("N2"="tranname")) 

inreg_orthos <- filtOrthos %>% dplyr::filter(!is.na(seqid)) 
outreg_orthos <- filtOrthos %>% dplyr::filter(is.na(seqid)) %>% 
  dplyr::select(-seqid,-seqname,-start,-end,-og_loc) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end,refStart,refEnd) %>% 
                     dplyr::mutate(refChrom=hap_chrom) %>%
                     dplyr::mutate(og_loc="out_region"),by=c("N2"="tranname")) %>%
  dplyr::filter(seqid==refChrom) %>%
  dplyr::mutate(start_dist=abs(refStart-end),end_dist=abs(start-refEnd)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(min_dist_bases=min(start_dist,end_dist)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(updown=ifelse(end < refStart,"upstream","downstream")) %>%
  dplyr::mutate(status=ifelse(min_dist_bases <10000,"out_expand","outside"))

if (nrow(outreg_orthos  %>% dplyr::filter(min_dist_bases < 10000)) > 0) {
  print("WARNING: There is at least one paralog that is within 10 kb of a gene within your defined boundary in N2. Your boundary will be automatically expanded to include:")
  print(outreg_orthos %>% dplyr::select(seqid,seqname,start,end,N2,min_dist_bases) %>% dplyr::filter(min_dist_bases < 10000))
}

#generate a lookup table (all_ortho_pairs) which contains all pairwise gene orthologs between REF and WILD
orthoList <- list()
orthoList_bound <- list()
orthoList_raw <- list()
strainCol_iter <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup","N2")]

for (i in 1:length(strainCol_iter)) {
  
  id=strainCol_iter[i]
  raw_tmp <- orthos %>%
    dplyr::select(Orthogroup,strainCol_iter[i],N2) %>%
    dplyr::mutate(str=!!sym(strainCol_iter[i])) %>%
    dplyr::mutate(str = strsplit(as.character(str), ",")) %>%
    tidyr::unnest(str) %>%
    dplyr::mutate(str=trimws(str)) %>%
    dplyr::filter(!is.na(N2)) %>%
    dplyr::select(-strainCol_iter[i]) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    dplyr::mutate(has_any_ortho=T) %>%
    dplyr::left_join(wild_tran,by=c("STRAIN","str"="tranname"))
    #dplyr::select(-N2)
  
  orthoList_raw[[i]] <- raw_tmp
    
  print(paste0("Mapped orthologs for ",i,"/",length(strainCol_iter)," strains."))
  tmp <- rbind(inreg_orthos %>% dplyr::mutate(status="within") %>% dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status),outreg_orthos %>% 
                 dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status)) %>%
    dplyr::select(Orthogroup,N2,strainCol_iter[i],og_loc,status) %>%
    dplyr::rename(tmpSel=strainCol_iter[i]) %>%
    dplyr::mutate(newSel = strsplit(as.character(tmpSel), ",")) %>%
    tidyr::unnest(newSel) %>%
    dplyr::mutate(newSel=trimws(newSel)) %>%
    dplyr::select(-tmpSel) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    tidyr::separate(newSel,into=c("Name","tnum"),sep="\\.",remove = F) %>%
    dplyr::select(Orthogroup,newSel,Name,STRAIN,N2,-tnum,og_loc,status) %>%
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
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,N2,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
        tidyr::unnest(N2) %>%
        dplyr::mutate(N2=trimws(N2))
        
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
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,N2,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
        tidyr::unnest(N2) %>%
        dplyr::mutate(N2=trimws(N2))
      
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
  dplyr::arrange(STRAIN,start)# %>%
  #dplyr::left_join(corr_jumps %>% dplyr::select(STRAIN,Parent,keep),by=c("STRAIN","Parent")) %>%
  #dplyr::filter(!is.na(keep))

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

N2_expand <- rbind(inreg_orthos %>% 
                       dplyr::mutate(status="within"),outreg_orthos %>% 
                       dplyr::select(-refStart,-refEnd,-refChrom,-start_dist,-end_dist,-min_dist_bases,-updown)) %>%
                       #dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status)) %>%
  dplyr::filter(!status=="outside")

new_boundaries_N2 <- N2_expand %>%
  dplyr::mutate(minStart=min(start), maxEnd=max(end)) %>%
  dplyr::distinct(N2,.keep_all = T) %>% 
  dplyr::filter(start==minStart | end==maxEnd) %>%
  dplyr::mutate(gene2gene=paste(seqname,collapse="-")) %>%
  dplyr::mutate(STRAIN="N2") %>%
  dplyr::distinct(minStart,.keep_all = T) %>%
  dplyr::select(seqid,minStart,maxEnd,STRAIN,gene2gene)

new_boundaries <- rbind(new_boundaries_WI,new_boundaries_N2) %>%
  dplyr::rename(boundStart=minStart,boundEnd=maxEnd)

#find the bound genes for each strain that are not orthologous
boundGenes <- rbind(wild_tran %>% 
                      dplyr::select(-boundChrom,-boundStart,-boundEnd,-REF,-refStart,-refEnd) %>% 
                      dplyr::mutate(alias=NA),
                    N2_tran %>% dplyr::select(tranname,seqname,STRAIN,seqid,start,end,strand,alias) %>%
                      dplyr::rename(Parent=seqname)) %>%
  dplyr::left_join(new_boundaries,by=c("STRAIN","seqid")) %>%
  dplyr::filter(!is.na(boundStart)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(start >= boundStart & start <= boundEnd) %>%
  dplyr::ungroup()


N2_ad <- boundGenes %>% 
  dplyr::filter(STRAIN=="N2") %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(tranname %in% all_orthos_unnest$N2,T,F)) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(tranname %in% all_ortho_pairs_bound$N2,T,F)) %>%
  dplyr::arrange(start) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup()

reassess_distal <- N2_ad %>% dplyr::filter(has_any_ortho==T & has_bound_ortho==F) 
distal_ortho <- inreg_orthos %>% 
  dplyr::filter(N2 %in% reassess_distal$tranname) %>% 
  dplyr::mutate(comma_count = stringr::str_count(CB4856, ",")+1) %>%
  dplyr::group_by(CB4856) %>%
  dplyr::mutate(comma_count=sum(comma_count)) %>%
  dplyr::filter(comma_count > 1)

N2_ad_corr <- N2_ad %>%
  dplyr::mutate(has_any_ortho=ifelse(tranname %in% distal_ortho$N2,F,has_any_ortho))

WI_ad <- boundGenes %>% 
  dplyr::filter(!STRAIN=="N2") %>%
  dplyr::left_join(all_ortho_pairs_raw,by=c("STRAIN","tranname")) %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(is.na(has_any_ortho),F,has_any_ortho)) %>%
  dplyr::left_join(all_ortho_pairs_bound %>% dplyr::select(tranname,STRAIN,N2,status) %>% dplyr::filter(N2 %in% N2_ad$tranname),by=c("tranname","STRAIN")) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(!is.na(status),T,F)) %>%
  dplyr::select(-status,-has_any_ortho,-N2) %>% 
  dplyr::group_by(Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(n_gene=n_distinct(Parent)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(n_gene) %>%
  dplyr::mutate(y_pos=rleid(STRAIN)) %>%
  dplyr::select(-n_gene) #%>%
  #dplyr::mutate(y_pos=ifelse(strand=="+",y_pos+0.2,y_pos-0.2))

N2_ad_corr <- N2_ad_corr %>%
  dplyr::mutate(y_pos=max(WI_ad$y_pos)+1)

reassess_distal_WI <- WI_ad %>% dplyr::filter(has_any_ortho==T & has_bound_ortho==F) 
raw_WI <- ldply(orthoList_raw,data.frame) %>% dplyr::rename(tranname=str) 
distal_ortho_WI <- all_orthos_unnest %>% 
  dplyr::filter(grepl(paste(reassess_distal_WI$tranname,collapse="|"),CB4856)) %>% 
  dplyr::mutate(comma_count_N2 = stringr::str_count(N2, ",")+1) %>%
  dplyr::mutate(comma_count_CB = stringr::str_count(CB4856, ",")+1) %>% 
  dplyr::filter(comma_count_N2 > 1 | comma_count_CB > 1 ) %>%
  dplyr::mutate(CB4856 = strsplit(as.character(CB4856), ",")) %>%
  tidyr::unnest(CB4856) %>%
  dplyr::mutate(CB4856=trimws(CB4856)) 

WI_ad_corr <- WI_ad %>%
  dplyr::mutate(has_any_ortho=ifelse(tranname %in% distal_ortho_WI$CB4856,F,has_any_ortho))

all_ad <- rbind(N2_ad_corr,WI_ad_corr) %>% 
  dplyr::arrange(STRAIN,start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(gen_pos=rleid(Parent)) %>%
  dplyr::mutate(shift=min(start)) %>%
  dplyr::mutate(end=end-min(start),start=start-min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(has_any_ortho==T & has_bound_ortho ==T,2,ifelse(has_any_ortho==T,1,0))) 


hlines <- new_boundaries %>% 
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,shift) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN")

segments <- all_ortho_pairs_bound %>%
  dplyr::select(STRAIN,Parent,start,end,N2,strand) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::distinct(Parent,N2,.keep_all = T) %>%
  dplyr::left_join(N2_tran %>% 
                     dplyr::rename(start_N2=start,end_N2=end,strand_N2=strand,chrom_N2=seqid,N2id=STRAIN) %>% 
                     dplyr::select(tranname,chrom_N2,start_N2,end_N2,strand_N2,alias,seqname,N2id),
                   by=c("N2"="tranname")) %>% 
  dplyr::filter(N2 %in% N2_ad$tranname) %>%
  #dplyr::filter(chrom_N2==seqid & start_N2 > (boundStart-5e4) & end_N2 < (boundEnd+5e4)) #%>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::rename(WI_y_pos=y_pos) %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by=c("N2id"="STRAIN")) %>%
  dplyr::rename(N2_y_pos=y_pos) %>%
  dplyr::mutate(N2_shift=min(start_N2),WI_shift=min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(WI_x_pos=(start+((end-start)/2))-WI_shift, N2_x_pos=(start_N2+((end_N2-start_N2)/2)-N2_shift)) %>%
  dplyr::mutate(WI_y_pos=WI_y_pos+0.2,N2_y_pos=N2_y_pos-0.2) %>%
  dplyr::distinct(STRAIN,Parent,seqname,.keep_all = T) %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::mutate(n1=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,seqname) %>%
  dplyr::mutate(n2=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(n1>1 | n2>1,"multi_copy","single_copy")) #%>%
  #dplyr::mutate(N2_y_pos=ifelse(strand_N2=="+",N2_y_pos,N2_y_pos-0.4),WI_y_pos=ifelse(strand=="+",WI_y_pos,WI_y_pos-0.4))

plot_ad <- all_ad %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::filter(col==max(col)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::mutate(class=ifelse(col==0,"no_known_ortho",ifelse(col==1,"has_distal_ortho","has_local_ortho"))) %>%
  dplyr::mutate(class=ifelse(class=="no_known_ortho" & STRAIN=="N2","no_known_allelic_N2",ifelse(class=="no_known_ortho" & STRAIN=="CB4856","no_known_allelic_CB",class)))

plot_ad2 <- all_ad %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::filter(col==max(col)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::mutate(class=ifelse(col==0,"no_known_ortho",ifelse(col==1,"has_distal_ortho","has_local_ortho"))) %>%
  dplyr::mutate(class=ifelse(class=="no_known_ortho" & STRAIN=="N2","no_known_allelic_N2",ifelse(class=="no_known_ortho" & STRAIN=="CB4856","no_known_allelic_CB",class))) %>%
  dplyr::mutate(class=ifelse(tranname=="Transcript_F11A5.10.1"|tranname=="g6321.t1","glc-1",class))

#CB_shift <- all_ad %>% dplyr::filter(STRAIN=="CB4856")
# tigs <- tigFilt %>%
#   dplyr::arrange(S2) %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(leadDiff=lead(S2)-S2) %>%
#   dplyr::mutate(lagDiff=S2-lag(S2)) %>%
#   dplyr::mutate(lagDiff=ifelse(is.na(lagDiff),0,lagDiff)) %>%
#   dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
#   dplyr::mutate(drop=ifelse(lagDiff >2e5 & lag(leadDiff) > 2e5,T,F)) %>%
#   dplyr::filter(drop==F)%>%
#   dplyr::ungroup() %>%
#  dplyr::select(-lagDiff,-leadDiff,-drop) %>% dplyr::filter(STRAIN=="CB4856") %>% dplyr::select(S2,E2) %>% dplyr::mutate(INV=ifelse(S2>E2,"INV","NORM")) %>% dplyr::mutate(S2_adj=S2-unique(CB_shift$shift),E2_adj=E2-unique(CB_shift$shift))

all_hap <- ggplot() +
  geom_segment(data=hlines,aes(x=boundStart-shift,xend=boundEnd-shift,y=y_pos,yend=y_pos))+
  geom_segment(data=segments,aes(x=WI_x_pos,xend=N2_x_pos,y=WI_y_pos,yend=N2_y_pos,col=col)) +
  geom_rect(data=plot_ad, aes(xmin=start,xmax=end,ymin=y_pos+0.2,ymax=y_pos-0.2,fill=class),color="black") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text("77 isotype strains"),
        axis.ticks = element_blank()) +
  scale_color_manual(values=c("multi_copy"="black","single_copy"="grey")) +
  scale_fill_manual(values=c("has_local_ortho"="grey","has_distal_ortho"="black","no_known_allelic_N2"="#DB6333","no_known_allelic_CB"="blue")) +
  scale_x_continuous(expand = c(0.01,0)) #+
  #ylab("77 isotype strains")#+
  #geom_segment(data=tigs,aes(x=S2_adj,xend=E2_adj,y=0.7,yend=0.7,col=INV))

all_hap2 <- ggplot() +
  geom_segment(data=hlines,aes(x=boundStart-shift,xend=boundEnd-shift,y=y_pos,yend=y_pos))+
  geom_segment(data=segments,aes(x=WI_x_pos,xend=N2_x_pos,y=WI_y_pos,yend=N2_y_pos,col=col)) +
  geom_rect(data=plot_ad, aes(xmin=start,xmax=end,ymin=y_pos+0.2,ymax=y_pos-0.2,fill=class),color="black") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text("77 isotype strains"),
        axis.ticks = element_blank()) +
  scale_color_manual(values=c("multi_copy"="grey","single_copy"="black")) +
  scale_fill_manual(values=c("has_local_ortho"="grey","has_distal_ortho"="black","no_known_allelic_N2"="#DB6333","no_known_allelic_CB"="blue")) +
  scale_x_continuous(expand = c(0.01,0)) #+
#ylab("77 isotype strains")#+
#geom_segment(data=tigs,aes(x=S2_adj,xend=E2_adj,y=0.7,yend=0.7,col=INV))
all_hap
all_hap2


# all_hap3 <- ggplot() +
#   geom_segment(data=hlines,aes(x=boundStart-shift,xend=boundEnd-shift,y=y_pos,yend=y_pos))+
#   geom_segment(data=segments,aes(x=WI_x_pos,xend=N2_x_pos,y=WI_y_pos,yend=N2_y_pos,col=col)) +
#   geom_rect(data=plot_ad2, aes(xmin=start,xmax=end,ymin=y_pos+0.2,ymax=y_pos-0.2,fill=class),color="black") +
#   theme(legend.title = element_blank(),
#         panel.background = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.title.x = element_blank(),
#         #axis.title.y = element_text("77 isotype strains"),
#         axis.ticks = element_blank()) +
#   scale_color_manual(values=c("multi_copy"="black","single_copy"="grey")) +
#   scale_fill_manual(values=c("has_local_ortho"="grey","has_distal_ortho"="black","no_known_allelic_N2"="#DB6333","no_known_allelic_CB"="blue","glc-1"="green")) +
#   scale_x_continuous(expand = c(0.01,0)) #+
# #ylab("77 isotype strains")#+
# #geom_segment(data=tigs,aes(x=S2_adj,xend=E2_adj,y=0.7,yend=0.7,col=INV))
# 
# all_hap4 <- ggplot() +
#   geom_segment(data=hlines,aes(x=boundStart-shift,xend=boundEnd-shift,y=y_pos,yend=y_pos))+
#   geom_segment(data=segments,aes(x=WI_x_pos,xend=N2_x_pos,y=WI_y_pos,yend=N2_y_pos,col=col)) +
#   geom_rect(data=plot_ad2, aes(xmin=start,xmax=end,ymin=y_pos+0.2,ymax=y_pos-0.2,fill=class),color="black") +
#   theme(legend.title = element_blank(),
#         panel.background = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.title.x = element_blank(),
#         #axis.title.y = element_text("77 isotype strains"),
#         axis.ticks = element_blank()) +
#   scale_color_manual(values=c("multi_copy"="grey","single_copy"="black")) +
#   scale_fill_manual(values=c("has_local_ortho"="grey","has_distal_ortho"="black","no_known_allelic_N2"="#DB6333","no_known_allelic_CB"="blue","glc-1"="green")) +
#   scale_x_continuous(expand = c(0.01,0))

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc1_QTL_CBxN2_CNVPAV.png",all_hap, device = 'png',dpi=900,width = 15,height = 2.7,units = 'in')
# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc1_QTL_CBxN2_CNVPAV_rev.png",all_hap2, device = 'png',dpi=900,width = 15,height = 2.7,units = 'in')
# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc1_QTL_CBxN2_CNVPAV_glc1.png",all_hap3, device = 'png',dpi=900,width = 15,height = 2.7,units = 'in')
# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc1_QTL_CBxN2_CNVPAV_glc1_rev.png",all_hap4, device = 'png',dpi=900,width = 15,height = 2.7,units = 'in')

ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/IIL_2.36-2.66_CBxN2_CNVPAV.png",all_hap, device = 'png',dpi=900,width = 15,height = 2.7,units = 'in')
ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/IIL_2.36-2.66_CBxN2_CNVPAV_rev.png",all_hap2, device = 'png',dpi=900,width = 15,height = 2.7,units = 'in')
  
