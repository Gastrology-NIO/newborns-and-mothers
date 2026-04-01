ps_genus_B<-subset_samples(ps_genus, Type =="B")
sample_data(ps_genus_B)$mergedLocStad<-mapply(paste0, sample_data(ps_genus_B)$Location,";", sample_data(ps_genus_B)$Stadium)
unique(sample_data(ps_genus_B)$mergedLocStad)-> locations
for (location in locations){
for (location2 in locations){
if (location != location2) {
if (!(location == "Cheek;Mother" && location2 =="Cervix;Mother")) {
if (!(location2 == "Cheek;Mother" && location =="Cervix;Mother")){
if (!(location == "Stool;Mother" && location2 =="Cervix;Mother")){
if (!(location2 == "Stool;Mother" && location =="Cervix;Mother")){
if (!(location == "Stomach;Newborn" && location2 =="Cervix;Mother")){
if (!(location2 == "Stomach;Newborn" && location =="Cervix;Mother")){
ps_tmp<-subset_samples(ps_genus_B, mergedLocStad %in% c(location, location2))
output_ancomb2_ark5= ancombc2(data = ps_tmp, tax_level = "genus",
p_adj_method = "BH",
fix_formula = "mergedLocStad",
prv_cut=0.2,
pseudo_sens =FALSE)
result<-output_ancomb2_ark5$res
tmp_tax=data.frame(tax_table(ps_tmp))
otu_tab<-t(as.data.frame(as.matrix(otu_table(ps_tmp))))
tmp_res<-merge(result,tmp_tax, by.x="taxon",by.y=0)
tmp_res<-merge(tmp_res, otu_tab, by.x="taxon",by.y=0)
write.csv2(tmp_res, paste0(folder_out,"B/differential_ancom_deseq2_",location,"_", location2, "_B.csv"))
}}}}}}}}}


ps_genus_K<-subset_samples(ps_genus, Type =="K")
sample_data(ps_genus_K)$mergedLocStad<-mapply(paste0, sample_data(ps_genus_K)$Location,";", sample_data(ps_genus_K)$Stadium)
unique(sample_data(ps_genus_K)$mergedLocStad)-> locations
for (location in locations){
for (location2 in locations){
if (location != location2) {
if (!(location == "Cheek;Mother" && location2 =="Cervix;Mother")) {
if (!(location2 == "Cheek;Mother" && location =="Cervix;Mother")){
if (!(location == "Stool;Mother" && location2 =="Cervix;Mother")){
if (!(location2 == "Stool;Mother" && location =="Cervix;Mother")){
if (!(location == "Stomach;Newborn" && location2 =="Cervix;Mother")){
if (!(location2 == "Stomach;Newborn" && location =="Cervix;Mother")){
ps_tmp<-subset_samples(ps_genus_K, mergedLocStad %in% c(location, location2))
output_ancomb2_ark5= ancombc2(data = ps_tmp, tax_level = "genus",
p_adj_method = "BH",
fix_formula = "mergedLocStad",
prv_cut=0.2,
pseudo_sens =FALSE)
result<-output_ancomb2_ark5$res
tmp_tax=data.frame(tax_table(ps_tmp))
otu_tab<-t(as.data.frame(as.matrix(otu_table(ps_tmp))))
tmp_res<-merge(result,tmp_tax, by.x="taxon",by.y=0)
tmp_res<-merge(tmp_res, otu_tab, by.x="taxon",by.y=0)
write.csv2(tmp_res, paste0(folder_out,"differential_ancom_deseq2_",location,"_", location2, "_K.csv"))


for (location in locations){
  for (location2 in locations){
  read.csv2()
}
}

  

  sample_data(ps_genus)$mergedLocStad<-mapply(paste0, sample_data(ps_genus)$Location,";", sample_data(ps_genus)$Stadium)
unique(sample_data(ps_genus)$mergedLocStad)-> locations
for (location in locations){
ps_tmp<-subset_samples(ps_genus, mergedLocStad %in% c(location))
output_ancomb2_ark5= ancombc2(data = ps_tmp, tax_level = "genus",
p_adj_method = "BH",
fix_formula = "Type",
prv_cut=0.2,
pseudo_sens =FALSE)
result<-output_ancomb2_ark5$res
tmp_tax=data.frame(tax_table(ps_tmp))
otu_tab<-t(as.data.frame(as.matrix(otu_table(ps_tmp))))
tmp_res<-merge(result,tmp_tax, by.x="taxon",by.y=0)
tmp_res<-merge(tmp_res, otu_tab, by.x="taxon",by.y=0)
write.csv2(tmp_res, paste0(folder_out,"B vs K/differential_ancom_deseq2_",location, ".csv"))
}
}}}}}}}}}






ps_genus<-tax_glom(ps_all, "genus")
ps_genus_K<-subset_samples(ps_genus, Type =="K")
unique(sample_data(ps_genus_K)$mergedLocStad)-> locations
df_K_result <- data.frame(matrix(0, nrow = length(locations), ncol = length(locations)))
colnames(df_K_result) <- locations
rownames(df_K_result) <- locations
folder_out<-"/home/chwast/nio/noworodki/new_result/noworodki i matki/v3 - matki i noworodki razem/differential_bacteria/K/"
for (location in locations){
  for (location2 in locations){
    if (location!=location2){
      if (!(location == "Cheek;Mother" && location2 =="Cervix;Mother")) {
        if (!(location2 == "Cheek;Mother" && location =="Cervix;Mother")){
          if (!(location == "Stool;Mother" && location2 =="Cervix;Mother")){
            if (!(location2 == "Stool;Mother" && location =="Cervix;Mother")){
              if (!(location == "Stomach;Newborn" && location2 =="Cervix;Mother")){
                if (!(location2 == "Stomach;Newborn" && location =="Cervix;Mother")){
    file<-read.csv2(paste0(folder_out,"differential_ancom_deseq2_",location,"_", location2, "_K.csv"))
    lfc_name<-colnames(file)[grepl("lfc_mergedLocStad", colnames(file))]
    q_name<-colnames(file)[grepl("p_mergedLocStad", colnames(file))]
    tmp_1<-file[file[,q_name] < 0.1,]
    tmp_1<-tmp_1[tmp_1[,lfc_name] < 1,]
    tmp_1<-tmp_1[tmp_1[,lfc_name] > -1,]
    df_K_result[location, location2]<-nrow(tmp_1)
    
    
                }}}}}}}}}

write.csv2(df_K_result, paste0(folder_out,"pval_0.1_K.csv"))




ps_genus_K<-subset_samples(ps_genus, Type =="B")
unique(sample_data(ps_genus_K)$mergedLocStad)-> locations
df_K_result <- data.frame(matrix(0, nrow = length(locations), ncol = length(locations)))
colnames(df_K_result) <- locations
rownames(df_K_result) <- locations
folder_out<-"/home/chwast/nio/noworodki/new_result/noworodki i matki/v3 - matki i noworodki razem/differential_bacteria/B/"
for (location in locations){
  for (location2 in locations){
    if (location!=location2){
      if (!(location == "Cheek;Mother" && location2 =="Cervix;Mother")) {
        if (!(location2 == "Cheek;Mother" && location =="Cervix;Mother")){
          if (!(location == "Stool;Mother" && location2 =="Cervix;Mother")){
            if (!(location2 == "Stool;Mother" && location =="Cervix;Mother")){
              if (!(location == "Stomach;Newborn" && location2 =="Cervix;Mother")){
                if (!(location2 == "Stomach;Newborn" && location =="Cervix;Mother")){
                  file<-read.csv2(paste0(folder_out,"differential_ancom_deseq2_",location,"_", location2, "_B.csv"))
                  lfc_name<-colnames(file)[grepl("lfc_mergedLocStad", colnames(file))]
                  q_name<-colnames(file)[grepl("Q_mergedLocStad", colnames(file))]
                  tmp_1<-file[file[,q_name] < 0.1,]
                  tmp_1<-tmp_1[tmp_1[,lfc_name] < 1,]
                  tmp_1<-tmp_1[tmp_1[,lfc_name] > -1,]
                  df_K_result[location, location2]<-nrow(tmp_1)
                  
                  
                }}}}}}}}}

write.csv2(df_K_result, paste0(folder_out,"Qval_0.1_B.csv"))

library(patchwork)
plots <- list()

# mother stool vs neonat rectum LP 
df_long=data.frame(group=c("Barnesiella", "Butyricimonas", "Family XIII AD3011 group", "Alistipes", "[Eubacterium] eligens group", "Monoglobus", "Bacteroides", "Parabacteroides", "Sutterella_2", "Ruminococcus", "Incertae Sedis_218", "[Eubacterium] siraeum group", 
                          "Incertae Sedis_497", "Pseudomonas_2", "Rothia_2", "Porphyromonas", "Negativibacillus", "Faecalimonas", "Acutalibacter", "Intestinimonas", "Howardella", "Schaalia"), 
                   value=c(3.68825825, 3.802639199, 3.43642963, 2.532370646, 2.884823725, 3.080801731, 2.973559629, 3.238879598, 2.80509465, 2.138947076, 1.537903123, 2.263655423, 
                           -4.095654736, -5.494754708, -3.290578553, -5.313825699, -3.981574484, -3.045024343,-2.131867362,-3.078874751,-2.464869761,-2.877913715

                          ))
df_long$side<-NaN
df_long$side[df_long$value>0]<-"over-represented"
df_long$side[df_long$value<=0]<-"aunder-represented"
df_long <- df_long[df_long$group %in% c(
  "Barnesiella",
  "Butyricimonas",
  "Alistipes",
  "Bacteroides",
  "Parabacteroides",
  "Rothia_2",
  "Porphyromonas",
  "Odoribacter",
  "Butyricimonas",
  "Anaerococcus",
  "Rothia_2",
  "Porphyromonas",
  "Prevotella",
  "Gemella"
),]
df_long <- df_long %>%
  arrange(value) %>%
  mutate(group = factor(group, levels = group))
title<-"Mother stool vs neonat rectum LP"
p<-plot_enriched_bacteria(df_long, title)
plots[[1]] <- p 

# mother stool vs neonat rectum TP (Kontrola)
df_long=data.frame(group=c("Bacteroides", "Alistipes", "Barnesiella", "Parabacteroides", "Incertae Sedis_813", "Odoribacter", "Bilophila_2", "Butyricimonas", "Phascolarctobacterium", "Megasphaera", 
                          "Hoylesella", "Incertae Sedis_497", "Anaerococcus", "Rothia_2", 
  "Howardella", "Porphyromonas", "Campylobacter", "Leuconostoc", 
  "S5-A14a", "Sphingobium_2", "Peptoniphilus", "Prevotella", 
  "Streptococcus", "Schaalia", "Incertae Sedis_823", 
  "Incertae Sedis_700", "Gemella", "Incertae Sedis_310"), 
                   value=c(5.6317635, 4.204874179, 5.724571252, 3.002695564, 4.752656243, 3.29411916, 3.065767166, 3.752706889, 4.254214737, 3.854237555, 
                          -8.481330028, -5.208367693, -6.3103908, -4.173473615,
  -3.398012531, -4.769111267, -5.139803434, -4.173474956,
  -3.929030436, -4.878533337, -5.872864194, -6.214844003,
  -2.888941051, -3.139913279, -2.912656747, -2.531498773,
  -2.406810023, -1.809682004
                          ))

df_long$side<-NaN
df_long$side[df_long$value>0]<-"over-represented"
df_long$side[df_long$value<=0]<-"aunder-represented"
df_long <- df_long[df_long$group %in% c(
  "Barnesiella",
  "Butyricimonas",
  "Alistipes",
  "Bacteroides",
  "Parabacteroides",
  "Rothia_2",
  "Porphyromonas",
  "Odoribacter",
  "Butyricimonas",
  "Anaerococcus",
  "Rothia_2",
  "Porphyromonas",
  "Prevotella",
  "Gemella"
),]
df_long <- df_long %>%
  arrange(value) %>%
  mutate(group = factor(group, levels = group))

title2<-"Mother stool vs neonat rectum TP"
p2<-plot_enriched_bacteria(df_long, title2)
plots[[2]] <- p2

# Neonat placenta vs mother cervix TP
df_long=data.frame(group=c("Incertae Sedis_823", "Corynebacterium_2", "Incertae Sedis_497",
  "Rothia_2", "Peptostreptococcus", "Incertae Sedis_375",
  "Campylobacter", "Rheinheimera_2", "Bilophila_2",
  "Lactobacillus", "Clostridium", "Enhydrobacter_2",
  "Howardella", "Litchfieldia"), 
                   value=c(3.192005971, 2.348829513, 3.063650033, 3.225880615,
  2.830260626, 2.517052822, 2.461695686, 3.334372413,
  2.916675213, -8.649570477, -6.325495873, -4.211071864,
  -6.366283444, -4.614924484))

df_long$side<-NaN
df_long$side[df_long$value>0]<-"over-represented"
df_long$side[df_long$value<=0]<-"aunder-represented"

df_long <- df_long[df_long$group %in% c(
  "Rheinheimera_2",
  "Lactobacillus",
  "Clostridium",
  "Enhydrobacter_2",
  "Streptococcus",
  "Cupriavidus_2",
  "Staphylococcus",
  "Lactobacillus",
  "Enhydrobacter_2",
  "Gardnerella_2",
  "Clostridium",
  "Veillonella"
)
,]

df_long <- df_long %>%
  arrange(value) %>%
  mutate(group = factor(group, levels = group))

title2<-"Neonat placenta vs mother cervix LP"
p3<-plot_enriched_bacteria(df_long, title2)
plots[[3]] <- p3


# Neonat placenta vs mother cervix LP
df_long=data.frame(group=c(
  "Streptococcus",
  "Cupriavidus_2",
  "Staphylococcus",
  "Corynebacterium_2",
  "Prevotella",
  "Incertae Sedis_497",
  "Lactobacillus",
  "Enhydrobacter_2",
  "Gardnerella_2",
  "Clostridium",
  "Finegoldia",
  "Veillonella"
), 
                   value=c(
  4.220038,
  3.267905,
  3.451095,
  2.8449,
  2.816292,
  2.339428,
  -7.56392,
  -5.3046,
  -6.06571,
  -6.53919,
  -2.51187,
  -1.66644
))

df_long$side<-NaN
df_long$side[df_long$value>0]<-"over-represented"
df_long$side[df_long$value<=0]<-"aunder-represented"


df_long <- df_long[df_long$group %in% c(
  "Rheinheimera_2",
  "Lactobacillus",
  "Clostridium",
  "Enhydrobacter_2",
  "Streptococcus",
  "Cupriavidus_2",
  "Staphylococcus",
  "Lactobacillus",
  "Enhydrobacter_2",
  "Gardnerella_2",
  "Clostridium",
  "Veillonella"
)
,]
df_long <- df_long %>%
  arrange(value) %>%
  mutate(group = factor(group, levels = group))

title2<-"Neonat placenta vs mother cervix TP"
p4<-plot_enriched_bacteria(df_long, title2)
plots[[4]] <- p4

# Neonat stomach vs mother cheek TP
df_long=data.frame(group= c(
  "Rheinheimera_2",
  "Micrococcus_2",
  "Staphylococcus",
  "Pseudomonas_2",
  "Kocuria_2",
  "Brevundimonas_2",
  "Acinetobacter_2",
  "Enhydrobacter_2",
  "Chryseobacterium",
  "Tepidimonas_2",
  "Pseudoxanthomonas_2",
  "Methylobacterium",
  "Anaerococcus",
  "Facklamia",
  "Nocardioides_2",
  "Peptoniphilus",
  "Deinococcus",
  "Incertae Sedis_310",
  "Corynebacterium_2",
  "Enterococcus",
  "Incertae Sedis_648",
  "Bacteroides",
  "Cupriavidus_2",
  "Pedobacter",
  "Lactobacillus",
  "Comamonas_2",
  "Exiguobacterium",
  "Acidibacter_2",
  "Lactococcus",
  "Incertae Sedis_813",
  "Leyella",
  "Sphingomonas_2",
  "Moraxella_2",
  "Prevotella",
  "Rothia_2",
  "Neisseria_2",
  "Parvimonas",
  "Abiotrophia",
  "Lautropia_2",
  "Veillonella",
  "Haemophilus_2",
  "Hoylesella",
  "Campylobacter",
  "Segatella",
  "Filifactor",
  "Actinomyces_2",
  "Selenomonas",
  "Lancefieldella",
  "Aggregatibacter_2",
  "Catonella",
  "Incertae Sedis_702",
  "Porphyromonas",
  "Anaeroglobus",
  "Alloprevotella",
  "F0058",
  "Anoxybacillus",
  "Scardovia_2",
  "Schaalia",
  "Shuttleworthia",
  "Peptostreptococcus",
  "Cardiobacterium_2",
  "Gemella",
  "Streptococcus"
)
, 
                   value=c(
  5.040970751,
  4.325281842,
  4.0676956,
  4.45852139,
  4.473808821,
  5.003485204,
  4.025365869,
  5.05237412,
  5.344487794,
  4.373954428,
  4.741492242,
  4.21238448,
  3.589169658,
  4.866240663,
  4.021275846,
  3.572203984,
  4.855118749,
  4.27738857,
  2.780170802,
  2.842740257,
  3.530762568,
  3.44995154,
  2.335999058,
  3.367193554,
  2.984748409,
  1.684517173,
  2.504504383,
  3.611245812,
  2.531682003,
  2.453172967,
  3.67293393,
  3.577580851,
  1.480211936,
  -3.60947264,
  -3.65462063,
  -3.875552285,
  -3.873406355,
  -3.497734238,
  -3.881929807,
  -2.715905996,
  -3.491323086,
  -2.940959037,
  -2.638444179,
  -1.995165818,
  -3.495602514,
  -1.922295441,
  -2.377612566,
  -1.753534983,
  -2.206405678,
  -2.490990555,
  -1.920986764,
  -2.206692359,
  -2.258815627,
  -2.042906558,
  -2.616732681,
  -2.444576693,
  -2.758593865,
  -1.705684059,
  -1.843664139,
  -2.334495436,
  -1.708046344,
  -2.648721822,
  -1.888019343
)
)

df_long$side<-NaN
df_long$side[df_long$value>0]<-"over-represented"
df_long$side[df_long$value<=0]<-"aunder-represented"

gatunki <- c(
  "Micrococcus_2",
  "Staphylococcus",
  "Pseudomonas_2",
  "Kocuria_2",
  "Acinetobacter_2",
  "Methylobacterium",
  "Deinococcus",
  "Corynebacterium_2",
  "Bacteroides",
  "Prevotella",
  "Neisseria_2",
  "Parvimonas",
  "Veillonella",
  "Haemophilus_2",
  "Actinomyces_2",
  "Porphyromonas",
  "Gemella",
  "Barnesiella",
  "Butyricimonas"
)

df_long <- df_long[df_long$group %in% gatunki,]

df_long <- df_long %>%
  arrange(value) %>%
  mutate(group = factor(group, levels = group))

title2<-"Neonat stomach vs mother cheek TP"
p5<-plot_enriched_bacteria(df_long, title2)
plots[[5]] <- p5


# Neonat stomach vs mother cheek LP
df_long=data.frame(group=c(
  "Micrococcus_2",
  "Acinetobacter_2",
  "Bilophila_2",
  "Nocardioides_2",
  "Staphylococcus",
  "Pseudomonas_2",
  "Methylobacterium",
  "Deinococcus",
  "Bacteroides",
  "Lactobacillus",
  "Kocuria_2",
  "Incertae Sedis_497",
  "Enhydrobacter_2",
  "Sutterella_2",
  "Scardovia_2",
  "Corynebacterium_2",
  "Barnesiella",
  "Butyricimonas",
  "Peptoniphilus",
  "Comamonas_2",
  "Brevibacterium_2",
  "Rothia_2",
  "Neisseria_2",
  "Porphyromonas",
  "Hoylesella",
  "Prevotella",
  "Actinomyces_2",
  "Gemella",
  "Veillonella",
  "Treponema",
  "Campylobacter",
  "Kingella_2",
  "Lautropia_2"
), 
  value= c(
  4.672029609,
  6.433279942,
  4.209939866,
  5.361072212,
  3.91134557,
  3.437474433,
  4.09791136,
  4.75240512,
  2.976468637,
  3.763053505,
  3.816924981,
  3.437522801,
  3.130345079,
  2.670597547,
  3.237130067,
  2.194698684,
  3.373850836,
  3.837504838,
  2.524073495,
  2.381918758,
  4.39367517,
  -3.52754045,
  -3.910799306,
  -3.041976486,
  -2.985485278,
  -3.063042592,
  -2.753877557,
  -3.301226832,
  -2.720580739,
  -2.868067311,
  -1.93239652,
  -2.305286137,
  -2.829893653
)
)

df_long$side<-NaN
df_long$side[df_long$value>0]<-"over-represented"
df_long$side[df_long$value<=0]<-"aunder-represented"
df_long <- df_long[df_long$group %in% gatunki,]

df_long <- df_long %>%
  arrange(value) %>%
  mutate(group = factor(group, levels = group))

title2<-"Neonat stomach vs mother cheek LP"
p6<-plot_enriched_bacteria(df_long, title2)
plots[[6]] <- p6




  library(patchwork)

  p<-wrap_plots(plots, ncol = 2) +
    plot_annotation(tag_levels = "A")

ggsave(
     "different.svg",
     plot = p,
     width = 18, height = 18,     
)
plot_enriched_bacteria<-function(df_long, title){
p <- ggplot(df_long, aes(x = value, y = group)) +
  geom_col(aes(fill = value < 0)) +
  scale_fill_manual(
    values = c("TRUE" = "#EF5350", "FALSE" = "lightblue"),
    guide = "none"
  ) +
  
  facet_grid(. ~ side, scales = "free_x", space = "free_x") +
  
  coord_cartesian(clip = "off") +
  
  
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(20, 160, 20, 140),
    strip.background = element_blank(),
  ) +
  geom_text(
  aes(
    x = value + ifelse(value > 0, 0.2, -0.2),
    label = gsub("g__", "Genus: ",
            gsub("f__", "Family: ",
            gsub(";", "\n", group))),
    hjust = ifelse(value > 0, 0, 1)
  ),
  size = 3,
  lineheight = 0.7
)+ggtitle(title)
  return(p)}

# #2db62bff - zielony tp był #add8e6ff
# #e73785ff - różowy LP był #ef5350ff
