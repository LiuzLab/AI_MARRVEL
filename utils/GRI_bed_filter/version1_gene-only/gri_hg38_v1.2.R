#GRI.V1=protein_coding gene + HGMD mutations (50 flanking bp) + Cinvar mutations (50 flanking bp)

#library(GenomicRanges)
library(dplyr)
library(VariantAnnotation) # https://bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.html#variant-call-format-vcf-files
library(data.table)
# #check file md5
# file_path <- "clinvar.vcf.gz"
# md5_checksum <- digest::digest(file_path, algo = "md5", file = TRUE)
# md5file<- readLines('clinvar.vcf.gz.md5')
# md5_expected<-strsplit(md5file[1], " ")[[1]][1]
# print(md5_checksum==md5_expected)
# 
# # uncompress .gz files
# folder_path <- ""
# gz_files <- list.files(path = folder_path, pattern = "\\.gz$", full.names = TRUE)
# for (file in gz_files) {
#   R.utils::gunzip(file, remove = FALSE) # Set remove = TRUE if you want to delete the .gz file after uncompressing
# }

# define some functions
coverage_bed <- function(bed_file, sequence_report='sequence_report.tsv'){
  hg.seq<-readr::read_tsv(sequence_report)
  # bed<-readr::read_tsv('exome_filter.bed', col_names = FALSE)
  ##calculate hg seq size
  seq.24<-hg.seq[c(1:24),]
  seq.size<-sum(seq.24$`Seq length`)
  ##calculate bed size
  bed.size<-sum(bed_file[,4])
  ##get percentage
  size.per<-bed.size/seq.size
  print(paste(paste0(size.per*100,"%"),bed.size,seq.size))
}

# 1. Use gene info from GENCODE (hg38) for GRI ----------------------------
#read gtf file from gencode: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz
gtf_file <- "gencode.v46.basic.annotation.gtf"
gtf_data <- rtracklayer::import(gtf_file)
gtf_dt <- as.data.table(as.data.frame(gtf_data))
##include the chromosome ending positions: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
##downloaded tsv file contains sequence information
hg.seq<-readr::read_tsv('sequence_report.tsv')
seq.24<-hg.seq[c(1:24),]
chromosome_sizes<-seq.24 |> 
  dplyr::select(`Chromosome name`,`Seq length`) |> 
  dplyr::mutate(seqnames=paste0('chr',`Chromosome name`))
colnames(chromosome_sizes)[2]<-'chr_end'
gtf_dt_size <- merge(gtf_dt, chromosome_sizes, by="seqnames") #removed MT
#select protein coding genes
genes <- gtf_dt_size[type=="gene" & gene_type=="protein_coding"]
#extract gene location info and save into bed format
##first add promoter 1kb regions
##add source column
bed_coding <- genes |> 
  dplyr::select(
    chr = seqnames,
    start = start,  # Use the upstream 1kb start for BED
    end = end,
    name = gene_name,                 # Use gene name as the BED name field      
    strand = strand,
  ) |> 
  mutate(
    source = "GencodeV46.protein_coding"
  )
##merge overlapped region
gr <- GRanges(
  seqnames = bed_coding$chr,
  ranges = IRanges(start = bed_coding$start, end = bed_coding$end)
)
merged_gr <- reduce(gr)
bed_coding.dt<-as.data.table(merged_gr)
coverage_bed(bed_coding.dt)

# 2. Use HGMD 2022 (hg38) for GRI -----------------------------------------
#read hgmd vcf file
vcf_hg38<-'HGMD_Pro_2022.2_hg38.vcf'
vcf_data <- readVcf(vcf_hg38, genome = "hg38")
## Extract fixed fields and convert to data.table
fixed_dt <- as.data.table(vcf_data@fixed)
info_dt <- as.data.table(vcf_data@info)
range_dt <- as.data.table(vcf_data@rowRanges)
## merged useful columns
hgmd<-cbind(range_dt,info_dt)
hgmd_select<-hgmd |> 
  dplyr::select(seqnames,start,end,GENE,STRAND,CLASS,MUT) |> 
  mutate(source="HGMD_hg38_2022")
## Obtain row numebrs of "Rejected" mutations in HGMD and exclude them in range_dt
## add chr size
hgmd_noR <- hgmd_select[-which(hgmd_select$CLASS=='R'),]
hgmd_noR<-hgmd_noR |>  
  dplyr::mutate(seqnames=paste0('chr',seqnames))
hgmd_noR_size<-merge(hgmd_noR, chromosome_sizes[,c(2:3)], by='seqnames')
##add 50bp flanking region
hgmd_noR_size_flank<-hgmd_noR_size |> 
  dplyr::mutate(final_start=ifelse(start > 50, start - 50, 0),
         final_end=ifelse((end+50)>chr_end, chr_end, end+50),
         .after = seqnames )
##generate bed 
bed_hgmd<- hgmd_noR_size_flank |> 
  dplyr::select(
         chr=seqnames,
         start=final_start,
         end=final_end,
         name = GENE,                 # Use gene name as the BED name field      
         strand = STRAND,
         source
         )
##merge bed files from step1 and step2, and reduce them into non-overlapping regions
bed_gencode_hgmd<-rbind(bed_coding, bed_hgmd, fill=TRUE) |>
  arrange(chr, start, end)
##use genomicrange reduce function to create non-overlapping regions:
##https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html#basic-interval-operations-for-granges-objects
gr <- GRanges(
  seqnames = bed_gencode_hgmd$chr,
  ranges = IRanges(start = bed_gencode_hgmd$start, end = bed_gencode_hgmd$end)
)
merged_gr <- reduce(gr)
bed_gencode_hgmd.dt<-as.data.table(merged_gr)
coverage_bed(bed_gencode_hgmd.dt)

# 3. Include ClinVar (hg38) for GRI ----------------------------------------
#https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
#process vcf data
vcf_hg38<-'clinvar.vcf'
vcf_data <- readVcf(vcf_hg38, genome = "hg38")
##obtain pathogenic variants
info_df<-as.data.frame(info(vcf_data))
range_df <-as.data.frame(vcf_data@rowRanges)
colnames(info_df)
unique(info_df$CLNSIG)
clinvar<-cbind(range_df,info_df)
clinvar_patho<- clinvar |> 
  filter(grepl('pathogenic', CLNSIG))
##remove MT
clinvar_patho<-clinvar_patho |> 
  filter(seqnames!='MT')
clinvar_patho<-clinvar_patho |>  
  mutate(seqnames=paste0('chr',seqnames))
##add chromosome size column
clinvar_patho_size<-merge(clinvar_patho,chromosome_sizes[,c(2,3)],by='seqnames')
##add flanking region
clinvar_patho_size_flank<-clinvar_patho_size |> 
  mutate(final_start=ifelse(start > 50, start - 50, 0),
         final_end=ifelse((end+50)>chr_end, chr_end, end+50),
         .after = seqnames )
##generate bed and merge step3 with previously generated bed file from step1 and step2
##generate bed 
bed_clinvar <- clinvar_patho_size_flank %>% 
  dplyr::select(chr=seqnames,
         start=final_start,
         end=final_end,
         name = GENEINFO,             
         strand = strand
         ) |> 
  mutate(source="Clinvar_hg38_080624")
##merge and reduce
bed_gencode_hgmd_clinvar <- rbind(bed_gencode_hgmd, bed_clinvar) |> 
  arrange(chr,start,end)
gr <- GRanges(
  seqnames = bed_gencode_hgmd_clinvar$chr,
  ranges = IRanges(start = bed_gencode_hgmd_clinvar$start, end = bed_gencode_hgmd_clinvar$end)
)
merged_gr <- reduce(gr)
bed_gencode_hgmd_clinvar.dt<-as.data.table(merged_gr)
coverage_bed(bed_gencode_hgmd_clinvar.dt)
#write the bed file
write.table(bed_gencode_hgmd_clinvar.dt[,c(1:3)], file='genomic_region_filter.hg38.v1.bed', sep="\t", quote=FALSE, col.names = FALSE, row.names=FALSE)
