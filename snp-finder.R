library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
# library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# url <- "http://www.affymetrix.com/analysis/downloads/na32/genotyping/GenomeWideSNP_6.na32.annot.csv.zip"
# download.file(url, destfile = "data/affy.zip")
# unzip("data/affy.zip")
# z <- read_csv("data/NIHMS958804-supplement-Supplementary_Table4", skip = 7)


sz_snps <-
  read_csv("data/NIHMS958804-supplement-Supplementary_Table4", skip = 7) %>% 
  janitor::clean_names() %>%  
  select(chr_id = chromosome, 
         rs_name = index_snp_db_snp_b141, 
         pos_int = start_bp, 
         pos_end_int = end_bp, 
         allele_a = a1, 
         allele_b = a2, 
         number_of_sn_ps, length_kb) %>% 
  mutate(rs_id = str_remove(rs_name, "rs"))

all_snps <- read_tsv("data/snps.txt") %>% 
  filter(nitems > 100)

# z <- z1 %>% 
#   filter(rs_id %in% snps$snp_id)

sz_snps_granges <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, 
              ids = sz_snps$rs_name, 
              ifnotfound = "drop")

all_snps_granges <- GRanges(seqnames = all_snps$chr_id, 
               ranges = IRanges(start = all_snps$pos_int-1000, 
                                end = all_snps$pos_end_int+1000), 
               strand = all_snps$strand, 
               snp = all_snps$rs_id)

sz_overlaps <- findOverlaps(sz_snps_granges, 
                            all_snps_granges, 
                            type = "any")
sz <- ubsetByOverlaps(sz_snps_granges, all_snps_granges)
as.data.frame(z)$RefSNP_id
