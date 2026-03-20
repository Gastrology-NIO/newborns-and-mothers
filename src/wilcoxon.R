richness<-estimate_richness(ps_genus)
rownames(richness)<-sub("X", "", rownames(richness))

sample_data(ps_genus)$mergedLocStad<-mapply(paste0, sample_data(ps_genus)$Location,";", sample_data(ps_genus)$Stadium)
sample_data(ps_genus)$location_type<-paste(sample_data(ps_genus)$mergedLocStad, sample_data(ps_genus)$Type)
tmp<-merge(richness, sample_data(ps_genus), by=0)
locations_all<-unique(tmp$location_type)
locations_all
res_wil<-data.frame(location=c("a"),location2=c("a"),  chao1=c("a"), shannon=c("a"), observed=c("a"), fisher=c("a"), simpson=c("a"))
for (location in locations_all){
for (location2 in locations_all){
if (location!=location2){

  location_tmp<- tmp[tmp$location_type %in% c(location,location2),]
  res_chao<-pairwise.wilcox.test(location_tmp$Chao1, location_tmp$location_type, 
                                 p.adjust.method = "BH")$p.value[1]
  res_shan<-pairwise.wilcox.test(location_tmp$Shannon, location_tmp$location_type, 
                                 p.adjust.method = "BH")$p.value[1]
  # pairwise.wilcox.test(Shannon ~ Type, data = location_tmp)$p.value
  res_simpson<-pairwise.wilcox.test(location_tmp$Simpson, location_tmp$location_type, 
                                    p.adjust.method = "BH")$p.value[1]
  # pairwise.wilcox.test(Simpson ~ Type, data = location_tmp)$p.value
  res_fisher<-pairwise.wilcox.test(location_tmp$Fisher, location_tmp$location_type, 
                                   p.adjust.method = "BH")$p.value[1]
  # pairwise.wilcox.test(Fisher ~ Type, data = location_tmp)$p.value
  res_observed<-pairwise.wilcox.test(location_tmp$Observed, location_tmp$location_type, 
                                     p.adjust.method = "BH")$p.value[1]
  # pairwise.wilcox.test(Observed ~ Type, data = location_tmp)$p.value
  res_wil[nrow(res_wil) + 1,]<-c(location, location2, res_chao, res_shan, res_simpson, res_fisher, res_observed)
}}}
res_wil<-res_wil[c(2:nrow(res_wil)),]
write.csv2(res_wil, paste0(folder_out, "tests_wilcoxon.csv"), row.names = F)


