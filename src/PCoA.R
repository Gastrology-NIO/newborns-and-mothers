


PCoA <- function(ps, folder_out){
  ps_genus<-tax_glom(ps, "genus")
  sample_data(ps_genus)$mergedLocStad<-mapply(paste0, sample_data(ps_genus)$Location,";", sample_data(ps_genus)$Stadium)
  pseq.compositional <- microbiome::transform(ps_genus, "clr")
  ps.pcoa.a<-ordinate(pseq.compositional, method="PCoA", distance="euclidean")
  p<-plot_ordination(pseq.compositional, ps.pcoa.a, color = "mergedLocStad", shape="Type", label="Organism_id") + 
    geom_point(size = 3) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  ggsave(
    paste(folder_out, "pca_all_mothers_clr_aitchison.png"),
    plot = p,
  )

  unique(sample_data(ps_genus)$mergedLocStad) -> locations
  for (location in locations){
    for (location2 in locations){
      if (location!=location2){
        ps_tmp <- subset_samples(ps_genus, mergedLocStad %in% c(location, location2))
        pseq.compositional <- microbiome::transform(ps_tmp, "clr")
        ps.pcoa.a <- ordinate(pseq.compositional, method="PCoA", distance="euclidean")
        p <- plot_ordination(pseq.compositional, ps.pcoa.a, color = "mergedLocStad", shape="Type", label="Organism_id") + 
        geom_point(size = 3)
        ggsave(
          paste0(folder_out, "PCOA_", location,"_",location2, "_clr_aitchison.png"),
          plot = p,
        )
      }
    }
  }
}
