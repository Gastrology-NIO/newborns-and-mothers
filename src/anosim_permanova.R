
sample_data(ps_genus)$location_type<-paste(sample_data(ps_genus)$Location, sample_data(ps_genus)$Type, sample_data(ps_genus)$Stadium)
locations_all<-unique(sample_data(ps_genus)$location_type)
res_anosim<-data.frame(location=c("a"),location2=c("a"),statistic=c("a"), p_val=c("a"))
for (location in locations_all){
for (location2 in locations_all){
if (location!=location2){

  ps_tmp<-subset_samples(ps_genus, location_type %in% c(location, location2))
  ps_tmp <- microbiome::transform(ps_tmp, "clr")
  metadata_sub <- data.frame(sample_data(ps_tmp))
  
  anosim_test<-anosim(phyloseq::distance(ps_tmp, method="euclidean"), 
                             metadata_sub$location_type, permutations = 9999)
  res_anosim[nrow(res_anosim) + 1,]<-c(location,location2, anosim_test$statistic, anosim_test$signif)

}}}
p.adj <- p.adjust(res_anosim$p_val, method = "BH")
res_anosim <- cbind.data.frame(res_anosim, p.adj=p.adj)
res_anosim<-res_anosim[c(2:nrow(res_anosim)),]
write.csv2(res_anosim, paste0(folder_out, "anosim_bray_dist_bh.csv"))

ps_genus<-tax_glom(ps_all, "genus")
sample_data(ps_genus)$location_type<-paste(sample_data(ps_genus)$Location, sample_data(ps_genus)$Type)
physeq_clr <- microbiome::transform(ps_genus, "clr")
ps_dist_matrix <- phyloseq::distance(physeq_clr, method ="euclidean")
control_adonis<-pairwise.adonis(ps_dist_matrix, phyloseq::sample_data(ps_genus)$location_type, p.adjust.m="BH")
write.csv2(control_adonis,  paste0(folder_out, "adonis_euclidean_bh.csv"), row.names = F)
