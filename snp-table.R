library(tidyverse)
library(rsnps)
library(gt)

z <- tribble(~affy, ~snp,
             8212413, 2456973,
             10343904, 10148671,
             12069402, 783540,
             12084869, 950169,
             12176372, 4702,
             13716002, 1620977,
             13716502, 34305371,
             14930173, 301800,
             15652810, 2905426,
             18819114, 6704641,
             20324519, 10496091,
             22398323, 35761247,
             22442521, 1080500,
             22568310, 1452075,
             23047667, 13107325,
             23965547, 1106568,
             24061823, 851,
             24502599, 3101246,
             25886669, 12522290,
             27133315, 16867576,
             30229962, 140,
             31197888, 11993663,
             31756553, 4129585,
             35945093, 112635299,
             37123094, 117074560) %>% 
  mutate(snp = str_c("rs", snp))

asd <- read_csv("data/tabula-NIHMS1015648-supplement-Supplementary_Data (3).csv") %>% 
  janitor::clean_names() %>% 
  mutate(snp = str_remove_all(snp, "'"))

my_sz <- genes::pardinas_snp()

z <- z %>% mutate(sz = snp %in% my_sz$snps,
                  asd = snp %in% asd$snp)

z1 <- rsnps::ncbi_snp_query(snps = z$snp)
z1 %>% left_join(z, by = c("query" = "snp")) %>% 
  select(-class, -rsid, -assembly, -ref_seq, -minor, -hgvs) %>% 
  rename(snp = query,
         old = ancestral_allele,
         new = variation_allele) %>% 
  gt()  %>% 
  cols_hide(c(sz, asd)) %>% 
  tab_style(
    style = cell_fill(color = "#E0FFE0"),
    locations = cells_body(
      rows = sz == TRUE)
  ) %>% 
  tab_style(
    style = cell_fill(color = "#E0E0FF"),
    locations = cells_body(
      rows = asd == TRUE)
  ) %>% 
  tab_options(table.font.size = 10)
