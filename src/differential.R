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

