library(ggpubr)
PCoA_B_ks_K <- function(ps_all, folder_out){
  ps_genus<-tax_glom(ps_all, "genus")
  sample_data(ps_genus)$mergedLocStad<-mapply(paste0, sample_data(ps_genus)$Location,";", sample_data(ps_genus)$Stadium)
  pseq.compositional <- microbiome::transform(ps_genus, "clr")
  ps.pcoa.a<-ordinate(pseq.compositional, method="PCoA", distance="euclidean")
  p<-plot_ordination(pseq.compositional, ps.pcoa.a, color = "mergedLocStad", shape="Type", label="Organism_id") + 
    geom_point(size = 3) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  ggsave(
    paste(folder_out, "pca_all_mothers_clr_aitchison.svg"),
    plot = p,
  )

  unique(sample_data(ps_genus)$mergedLocStad) -> locations
  for (location in locations){
        ps_tmp <- subset_samples(ps_genus, mergedLocStad %in% c(location))
        pseq.compositional <- microbiome::transform(ps_tmp, "clr")
        ps.pcoa.a <- ordinate(pseq.compositional, method="PCoA", distance="euclidean")
        p <- plot_ordination(pseq.compositional, ps.pcoa.a, color="Type", label="Organism_id") + 
        geom_point(size = 3) + ggtitle(location)
        ggsave(
          paste0(folder_out, "PcOA/PCOA_", location, "_clr_aitchison.svg"),
          plot = p,
        )

  }
}


PCoA <- function(ps_all, folder_out){
  ps_genus<-tax_glom(ps_all, "genus")
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

PCoA_selected<- function(ps_all, save_path, type){
  pairs<- data.frame(loc1=c('Stool;Mother', 'Cervix;Mother', 'Cheek;Mother', 'Mother placenta;Mother'), loc2=c('Rectum;Newborn','Placenta;Newborn','Stomach;Newborn','Placenta;Newborn'))
  ps_genus<-tax_glom(ps_all, "genus")
  
  sample_data(ps_genus)$location_type<-paste0(sample_data(ps_genus)$Location, ";", sample_data(ps_genus)$Stadium,";", sample_data(ps_genus)$Type)
  physeq_clr <- microbiome::transform(ps_genus, "clr")
  ps_dist_matrix <- phyloseq::distance(physeq_clr, method ="euclidean")
  control_adonis<-pairwise.adonis(ps_dist_matrix, phyloseq::sample_data(ps_genus)$location_type, p.adjust.m="BH")
  
  ps_genus_type<-subset_samples(ps_genus, Type %in% c("TP"))
  type<-"TP"
  #   ps_genus_type<-subset_samples(ps_genus, Type %in% c("LP"))
  # type<-"LP"
plots <- list()
  for (pair_no in 1:nrow(pairs)) {

   sample_data(ps_genus_type)$mergedLocStad<-mapply(paste0, sample_data(ps_genus_type)$Location,";", sample_data(ps_genus_type)$Stadium) 
    ps_tmp<-subset_samples(ps_genus_type, mergedLocStad %in% pairs[pair_no,])
    pseq.compositional <- microbiome::transform(ps_tmp, "clr")
    ad_res <- control_adonis[control_adonis$pair  %in% c(paste(paste0(pairs$loc1[pair_no], ";",type),"vs" ,paste0(pairs$loc2[pair_no], ";",type)), paste(paste0(pairs$loc2[pair_no], ";",type),"vs" ,paste0(pairs$loc1[pair_no], ";",type))), ]
    ad_label <- paste0(
        "R² = ", round(ad_res$R2, 3), "\n",
        "p = ", signif(ad_res$p.adjusted, 3)
    )
    ps.pcoa.a <- ordinate(pseq.compositional, method="PCoA", distance="euclidean")
    p <- plot_ordination(pseq.compositional, ps.pcoa.a, color = "mergedLocStad", label="Organism_id") + 
    geom_point(size = 3)+
  stat_ellipse(type = "norm")+
  scale_color_manual(values = c("#e73785ff", "#2db62bff"))

    if ((paste0(pairs$loc1[pair_no], ";",type) == "Cheek;Mother;TP") & (paste0(pairs$loc1[pair_no], ";",type) == "Stomach;Newborn;TP")){
        p <- p +
  annotation_custom(
    grob = textGrob(
      ad_label,
      x = unit(0.98, "npc"),   # prawy margines
      y = unit(0.98, "npc"),   # górny margines
      hjust = 1, vjust = 1,    # wyrównanie do prawej i góry
      gp = gpar(col = "black", fontsize = 10)
    )
  ) +
  coord_cartesian(clip = "off")
    } else {
          p<-p+annotation_custom(
    grob = textGrob(
      ad_label,
      x = unit(0.02, "npc"),   # lewy margines
      y = unit(0.98, "npc"),   # górny margines
      hjust = 0, vjust = 1,
      gp = gpar(col = "black", fontsize = 10)
    )
  ) +
  coord_cartesian(clip = "off")

    }

    
   plots[[pair_no]] <- p 
  }
  library(patchwork)

  p<-wrap_plots(plots, ncol = 2) +
    plot_annotation(tag_levels = "A")
  
  ggsave(
      paste0(folder_out,type, "_PCOA_clr_aitchison.svg"),
      plot = p,
        width = 12, height = 6,

  )
}
