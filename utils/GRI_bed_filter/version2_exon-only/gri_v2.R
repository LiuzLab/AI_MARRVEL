# BED region =  all exons
#             + HGMD mutations (50 flanking bp)
#             + Cinvar mutations (50 flanking bp)
#             + spliceAI data (50 flanking bp)
#
# 1. gencode data:
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.basic.annotation.gtf.gz
#
# 2. HGMD data:
# HGMD_Pro_2024.3_hg38.vcf HGMD_Pro_2024.3_hg19.vcf
#
# 3. Clinvar data:
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/
#
# 4. spliceAI data v1.3 is on skyriver server:
# /home/zhijiany/data-transfer/AIM/spliceAI-2019/exome_filter

### Notes: 1. Please prepare/download all the data files before running the script
###        2. Decompress the gz files to avoid error when loading the data

library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(httr)
library(dplyr)
library(readr)
library(data.table)

generate_gri_v2 <- function(genome_version = "hg38") {
    # data preparation --------------------------------------------------------

    ## construct splice data file name
    snv <- paste0("snv.", genome_version, ".vcf")
    indel <- paste0("indel.", genome_version, ".vcf")
    columnnames <- unlist(strsplit("CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO", "\t"))

    ## load spliceAI data
    splice_snv <- read_table(snv,
        comment = "#",
        col_names = columnnames,
        col_types = "cncccccc"
    )

    splice_indel <- read_table(indel,
        comment = "#",
        col_names = columnnames,
        col_types = "cncccccc"
    )


    ## hgmd local file name
    hgmd_file <- paste0("HGMD_Pro_2024.3_", genome_version, ".vcf")

    ## use different data based on genome version
    if (genome_version == "hg38") {
        ## load chr size
        hg38 <- BSgenome.Hsapiens.UCSC.hg38
        size <- seqlengths(hg38)
        chr_size <-
            data.frame(seqnames = names(size), chr_end = as.integer(size))[c(1:24), ]

        ## load gencode data
        gtf_data <- rtracklayer::import("gencode.v46.basic.annotation.gtf")

        ## load hgmd data
        hgmd_data <- readVcf(hgmd_file, genome = "hg38")

        ## load clinvar data
        clinvar_data <- readVcf("clinvar.vcf.gz", genome = "hg38")
    } else if (genome_version == "hg19") {
        ## load chr size
        hg19 <- BSgenome.Hsapiens.UCSC.hg19
        size <- seqlengths(hg19)
        chr_size <-
            data.frame(seqnames = names(size), chr_end = as.integer(size))[c(1:24), ]

        ## load gencode data
        gtf_data <- rtracklayer::import("gencode.v46lift37.basic.annotation.gtf")

        ## load hgmd data
        hgmd_data <- readVcf(hgmd_file, genome = "hg19")

        ## load clinvar data
        clinvar_data <- readVcf("clinvar.vcf.gz", genome = "hg19")
    }


    ## function to check bed file coverage
    coverage_bed <- function(bed_file, chromosome_size = chr_size) {
        ## calculate hg chr bp length
        seq.size <- sum(chromosome_size[, 2])
        ## calculate bed size
        bed.size <- sum(bed_file[, 4])
        ## get percentage
        size.per <- bed.size / seq.size
        print(paste(
            paste0("bed file coverage is: ", size.per * 100, "% ;"),
            bed.size,
            "bp out of",
            seq.size,
            "bp"
        ))
    }

    # 1. include gencode v46 annotation ----------------------------

    ## convert gtf data to data.table and add chr size
    gtf_dt <- as.data.table(gtf_data)
    gtf_dt_size <-
        merge(gtf_dt, chr_size, by = "seqnames") # no MT chr

    ## select exon only
    exons <-
        gtf_dt_size[type == "exon"]

    ## extract gene location info and save into bed format
    ## add source column
    bed_coding <- exons |>
        dplyr::select(
            chr = seqnames,
            start = start,
            end = end,
            name = gene_name,
            strand = strand,
        ) |>
        mutate(source = "GencodeV46.exon")

    ## merge overlapped region
    gr <- GRanges(
        seqnames = bed_coding$chr,
        ranges = IRanges(start = bed_coding$start, end = bed_coding$end)
    )

    merged_gr <- reduce(gr)
    bed_coding.dt <- as.data.table(merged_gr)
    coverage_bed(bed_coding.dt) # about 4% coverage


    # 2. include hgmd 2022 variants -----------------------------------------

    ## extract fields and convert to data.table
    info_dt <- as.data.table(info(hgmd_data))
    range_dt <- as.data.table(rowRanges(hgmd_data))

    ## merged useful columns
    hgmd <- cbind(range_dt, info_dt) |>
        dplyr::select(seqnames, start, end, GENE, STRAND, CLASS, MUT) |>
        mutate(source = "HGMD_2024")

    ## exclude "rejected variants" in HGMD database
    hgmd_noR <- hgmd[-which(hgmd$CLASS == "R"), ]
    hgmd_noR <- hgmd_noR |>
        dplyr::mutate(seqnames = paste0("chr", seqnames))

    ## add 50bp flanking region
    hgmd_noR <- merge(hgmd_noR, chr_size, by = "seqnames")
    hgmd_noR_flank <- hgmd_noR |>
        dplyr::mutate(
            final_start = ifelse(start > 50, start - 50, 1),
            final_end = ifelse((end + 50) > chr_end, chr_end, end + 50),
            .after = seqnames
        )

    ## generate bed
    bed_hgmd <- hgmd_noR_flank |>
        dplyr::select(
            chr = seqnames,
            start = final_start,
            end = final_end,
            name = GENE,
            strand = STRAND,
            source
        )

    ## merge bed files from step1 and step2, and reduce them into non-overlapping regions
    bed_gencode_hgmd <- rbind(bed_coding, bed_hgmd, fill = TRUE) |>
        arrange(chr, start, end)

    gr <- GRanges(
        seqnames = bed_gencode_hgmd$chr,
        ranges = IRanges(start = bed_gencode_hgmd$start, end = bed_gencode_hgmd$end)
    )

    merged_gr <- reduce(gr)
    bed_gencode_hgmd.dt <- as.data.table(merged_gr)
    coverage_bed(bed_gencode_hgmd.dt) # about 4.16% coverage

    # 3. Include clinvar variants----------------------------------------

    ## obtain pathogenic variants
    info_df <-
        as.data.frame(clinvar_data@info) # faster to load as dataframe than datatable
    range_df <- as.data.frame(clinvar_data@rowRanges)
    clinvar <- cbind(range_df, info_df)
    clinvar_patho <- clinvar |>
        filter(grepl("pathogenic", CLNSIG))

    ## remove MT chromosome
    ## rename seqnames
    clinvar_patho <- clinvar_patho |>
        filter(seqnames != "MT")
    clinvar_patho <- clinvar_patho |>
        mutate(seqnames = paste0("chr", seqnames))

    ## add chromosome size column
    clinvar_patho <- merge(clinvar_patho, chr_size, by = "seqnames")

    ## add flanking region
    clinvar_patho_flank <- clinvar_patho |>
        mutate(
            final_start = ifelse(start > 50, start - 50, 1),
            final_end = ifelse((end + 50) > chr_end, chr_end, end + 50),
            .after = seqnames
        )

    ## generate bed
    bed_clinvar <- clinvar_patho_flank |>
        dplyr::select(
            chr = seqnames,
            start = final_start,
            end = final_end,
            name = GENEINFO,
            strand = strand
        ) |>
        mutate(source = "Clinvar_latest") # the Clinvar link provides latest vcf file

    ## merge and reduce
    bed_gencode_hgmd_clinvar <-
        rbind(bed_gencode_hgmd, bed_clinvar) |>
        arrange(chr, start, end)

    gr <- GRanges(
        seqnames = bed_gencode_hgmd_clinvar$chr,
        ranges = IRanges(start = bed_gencode_hgmd_clinvar$start, end = bed_gencode_hgmd_clinvar$end)
    )

    merged_gr <- reduce(gr)
    bed_gencode_hgmd_clinvar.dt <- as.data.table(merged_gr)

    ## check bed coverage
    coverage_bed(bed_gencode_hgmd_clinvar.dt) # coverage about 4.22%



    # 4. Include spliceAI data ----------------------------------------------------
    # Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">

    # merge snv and indel and add 50 flanking bp
    splice <- rbind(splice_snv, splice_indel) |>
        mutate(seqnames = paste0("chr", CHROM))
    splice <- merge(splice, chr_size, by = "seqnames")

    splice_flank <- splice |>
        mutate(
            chr = seqnames,
            start = ifelse(POS - 50 < 1, 1, POS - 50),
            end = ifelse(POS + 50 > chr_end, chr_end, POS + 50),
            source = "spliceAI_v1.3"
        ) |>
        dplyr::select(chr, start, end, source)

    ## merge and reduce
    bed_gencode_hgmd_clinvar_splice <-
        rbind(bed_gencode_hgmd_clinvar, splice_flank, fill = TRUE) |>
        arrange(chr, start, end)

    gr <- GRanges(
        seqnames = bed_gencode_hgmd_clinvar_splice$chr,
        ranges = IRanges(start = bed_gencode_hgmd_clinvar_splice$start, end = bed_gencode_hgmd_clinvar_splice$end)
    )

    merged_gr <- reduce(gr)
    splice.dt <- as.data.table(merged_gr)
    print("Final coverage is below:")
    coverage_bed(splice.dt) # coverage is about (4.36% without 50bp flanks)  5.75%

    # convert start position from 1-based to 0-based of bed format ------------
    splice.dt <- splice.dt |>
        mutate(start = pmax(0, start - 1))

    bed_gencode_hgmd_clinvar_splice <- bed_gencode_hgmd_clinvar_splice |>
        mutate(start = pmax(0, start - 1))

    # save the file
    # write the bed file
    file_name <- paste0("exon_only.", genome_version, ".bed")
    write_tsv(splice.dt[, c(1:3)], file = file_name, col_names = FALSE)

    # save the initial bed files which are not reduced, which contains the source column to indicate where that region is from
    file_name <- paste0("exon_only_source.", genome_version, ".bed")
    write_tsv(bed_gencode_hgmd_clinvar_splice,
        file = file_name,
        col_names = FALSE
    )
}

# Apply generate_gri_geneonly function to generate bed files --------------
# generate_gri_v2("hg19")
generate_gri_v2("hg38")