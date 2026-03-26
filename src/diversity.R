
library(tidyverse)
library(ggsignif)
library(patchwork)


plots <- list()

# shannon <- plot_richness(ps_genus, "Location", measures = c("Shannon"))
# shannon + theme_bw() + geom_boxplot(aes(fill = Type)) +facet_wrap(~Stadium, scales="free_x")+ theme(axis.text.x = element_text(angle=90, hjust=1)) +
#   geom_jitter(width = 0.15, alpha = 0.6) 

library(dplyr)


# funkcja zamieniająca p-value na gwiazdki
p_to_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.1)   return(".")
  return("ns")
}


position <- function(df, location_name, all_loc){ 
  # obliczamy pozycje x dla Skin osobno w każdym Stadium
  for (i in 1:length(all_loc)){
    if (all_loc[i]==location_name){
        pos <- df %>%
          filter(Location == location_name) %>%
          group_by(Stadium) %>%
          summarise(
            x = i,
            y = max(value) * 1.1
          )
        return(pos)
      }
    }
}


add_clamr <- function(plot, pos, label_star){
   plot<-plot+ geom_segment(
    data = pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = pos,
    aes(x = x, y = y * 1.05, label = label_star),
    inherit.aes = FALSE
  ) 
  print(label_star)
  return(plot)
}



plot_richness_with_p_val<-function(ps_genus, measure){
    df <- plot_richness(ps_genus, "Location", measures = measure)$data
    mother<-df %>%
        filter(Stadium== "Mother") 
    mother_loc<-sort(unique(mother$Location))
    
    
    newborn<-df %>%
        filter(Stadium== "Newborn") 
    newborn_loc<-sort(unique(newborn$Location))
    
    plot <- ggplot(df, aes(Location, value, fill = Type)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, size = 1) +
      facet_wrap(~Stadium, scales = "free_x") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle(measure) +
      labs(x =NULL, y = NULL)+
      scale_fill_manual(values = c("B" = "#e73785ff",
                               "K" = "#2db62bff"))
    
    for (location_name in mother_loc){
      position_loc <- position(df, location_name, mother_loc)
      p_value <- wilcox.test(value ~ Type, data = df %>% filter(Location == location_name))$p.value
      p_value_star <- p_to_stars(p_value)

      if (p_value_star!="ns"){
          print("mother")
          print(location_name)
          print(p_value_star)    
          plot <- add_clamr(plot, position_loc, p_value_star)
        }
    }
    
    for (location_name in newborn_loc){
      position_loc <- position(df, location_name, newborn_loc)
      p_value <- wilcox.test(value ~ Type, data = df %>% filter(Location == location_name))$p.value
      p_value_star <- p_to_stars(p_value)
      if (p_value_star!="ns"){
        print("newborn")
        print(location_name)
        print(p_value_star)
        
        plot <- add_clamr(plot, position_loc, p_value_star)
        }
    }
    return(plot)
}

# e737e3ff + 2db62bff
# e73785ff

# było f8766dff 00bfc4ff
shannon<-plot_richness_with_p_val(ps_genus, "Shannon")
shannon
plots[[1]] <- shannon 

chao1<-plot_richness_with_p_val(ps_genus, "Chao1")
chao1
plots[[2]] <- chao1

simpson<-plot_richness_with_p_val(ps_genus, "Simpson")
simpson
plots[[3]] <- simpson 

observed<-plot_richness_with_p_val(ps_genus, "Observed")
observed
plots[[4]] <- observed 

fisher<-plot_richness_with_p_val(ps_genus, "Fisher")
fisher
plots[[5]] <- fisher 


  p<-wrap_plots(plots, ncol = 2) +
    plot_annotation(tag_levels = "A")
  
  ggsave(
      paste0(folder_out, "diversity.svg"),
      plot = p,
    width = 12, height = 18,
    limitsize  = FALSE
  )




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

####


library(ggpubr)

df <- plot_richness(ps_genus, "Location", measures = "Fisher")$data

ggplot(df, aes(Location, value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~Stadium, scales = "free_x") +
  
  # globalny test Kruskala dla Location
  stat_compare_means(
    aes(group = 1),
    method = "kruskal.test",
    label = "p.format"
  ) +
  
  # klamry: B vs K w obrębie każdej Location
  stat_compare_means(
    aes(group = Type),
    method = "wilcox.test",
    label = "p.format"   # <--- tu zmieniamy gwiazdki na p-value
  ) +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



## between different locations
df <- plot_richness(ps_genus, "Location", measures = "Shannon")$data

# listy porównań osobno dla Mother i Newborn
comparisons_mother <- combn(unique(df$Location[df$Stadium == "Mother"]), 2, simplify = FALSE)
comparisons_newborn <- combn(unique(df$Location[df$Stadium == "Newborn"]), 2, simplify = FALSE)


ggplot(df, aes(Location, value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  facet_wrap(~Stadium, scales = "free_x") +
  
  # globalny test Kruskala dla Location
  stat_compare_means(
    aes(group = 1),   # ważne: ignoruje Type
    method = "kruskal.test",
    label = "p.format"
  ) +
  
  # klamry dla Mother
  stat_compare_means(
    comparisons = comparisons_mother,
    method = "wilcox.test",
    label = "p.signif",
    data = subset(df, Stadium == "Mother")
  ) +
  
  # klamry dla Newborn
  stat_compare_means(
    comparisons = comparisons_newborn,
    method = "wilcox.test",
    label = "p.signif",
    data = subset(df, Stadium == "Newborn")
  ) +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
