library(tidyverse)
library(gridExtra)
library(grid)
library(googlesheets4)
opt <- list(
  run="SF07",
  nmol_conc=150
)




probe_ov_url <- "https://docs.google.com/spreadsheets/d/1vNJssytzfJEsmYrIdLDZi3v_V6m-sokpKpTKtSpLXbk/edit?gid=1469504260#gid=1469504260"

sf_plan_url <- "https://docs.google.com/spreadsheets/d/1irN9c8recUKXVmJwj8HFg_S2YxxrPnLqldEXZ9ddcTw/edit?gid=1206615644#gid=1206615644"

probe_ov <- probe_ov_url %>% read_sheet(sheet="ext_probes", col_types="c") 

sf_plan <- read_sheet(sf_plan_url, sheet=opt$run)


volume_plan <-
  sf_plan %>%
  ## remove unused wells
  filter(!is.na(vol)) %>%
  group_by(probes) %>%
  summarize(
    wells=paste(well, collapse=", "),
    tot_vol=sum(vol)
  )


all_probes <- unlist(str_split(sf_plan$probes, ", ")) %>% 
  .[which(!is.na(.))] %>%
  unique()

primary_ov_pre <- probe_ov %>%
  filter(ext_probe_id %in% all_probes) %>%
  mutate(molconc=as.numeric(molar_conc_nmol_L)) %>%
  select(ext_probe_id, molconc) %>%
  ## join indfo about total volume per probe
  left_join(
    volume_plan %>%
      group_by(wells) %>%
      reframe(
        ext_probe_id=str_split(probes, ", ") %>% unlist(),
        tot_vol=tot_vol
      ),
    by = join_by(ext_probe_id)
  ) %>%
  mutate(
    vol_probe=(opt$nmol_conc*tot_vol)/molconc
  ) %>%
  select(-tot_vol, -molconc) %>%
  pivot_wider(id_cols="wells", names_from = "ext_probe_id", values_from = "vol_probe") %>%
  
  
  bind_cols(
    .,
    tibble(probe_vol=rowSums(select(.,-wells), na.rm=T))
  ) %>%
  left_join(
    volume_plan
    ,.,
    by=c("wells")
  ) %>%
  mutate(
    probe_vol=ifelse(is.na(probe_vol), 0, probe_vol),
    formamide=0.5*tot_vol,
    fish_mix_4x=0.25*tot_vol,
    nfw=tot_vol-(formamide+fish_mix_4x+probe_vol)
  )%>%
  mutate_at(.vars=c(all_probes, "nfw"), .funs=round, digits=1) %>%
  mutate_all(~ ifelse(is.na(.), "-", .)) %>%
  select(-probes, -probe_vol)


names_converter <- tibble(old=names(primary_ov_pre)) %>%
  left_join(probe_ov %>% select(old=ext_probe_id, label)) %>%
  mutate(
    new=ifelse(is.na(label), old, label)
  )


primary_ov <- primary_ov_pre %>%
  set_names(nm=names_converter$new)


secondary_ov <- 
sf_plan %>%
  rowwise() %>%
  mutate(primers=paste(unlist(str_extract_all(probes, paste0("p",seq(10,60),collapse="|"))),collapse=", "),
         cond_pairs_id=str_remove(condition, "_T$|_C$")) %>%
  group_by(cond_pairs_id) %>%
  summarize(
    pair_vol=sum(vol),
    primers=primers %>% .[which(.!="NA")],
    wells=paste(well, collapse=", ")
  ) %>%
  group_by(primers) %>%
  summarize(
    tot_vol=sum(pair_vol),
    wells=paste(wells, collapse=", ")
  ) %>%
  group_by(wells) %>%
  reframe(
    tot_vol=tot_vol, 
    primer=str_split(primers, ", ") %>% unlist()
  ) %>%
  mutate(
    primer_vol=0.01*tot_vol,
    pbs_10x=0.1*tot_vol
  ) %>%
  pivot_wider(id_cols=c(wells, tot_vol, pbs_10x),
              names_from = primer,
              values_from=primer_vol) %>%
  bind_cols(
    .,
    tibble(probe_vol=rowSums(select(.,-wells, -tot_vol), na.rm=T))
  ) %>%
  mutate(nfw=tot_vol-probe_vol) %>%
  mutate_all(~ ifelse(is.na(.), "-", .)) %>%
  select(-probe_vol)






title1 <- textGrob(
  paste(opt$run, "primary hybridisation setup"),
  gp = gpar(fontsize = 16, fontface = "bold")
)


tbl1_grob <- tableGrob(
  primary_ov,
  rows = NULL,
  theme = ttheme_default(base_size = 9)
)
  

title2 <- textGrob(
  paste(opt$run, "Secondary fluorostaining setup"),
  gp = gpar(fontsize = 16, fontface = "bold")
)


tbl2_grob <- tableGrob(
  secondary_ov,
  rows = NULL,
  theme = ttheme_default(base_size = 9)
)
  
pdf("/g/schwab/marco/table.pdf", width = 10, height = 4)


  grid.arrange(
    title1,
    tbl1_grob,
    title2,
    tbl2_grob,
    heights=c(0.1,1,0.1,1)
  )
dev.off()






