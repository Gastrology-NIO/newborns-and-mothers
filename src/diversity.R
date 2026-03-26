
library(tidyverse)
library(ggsignif)



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
            x = 3,
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



plot_richness_with_p_val<-function(ps_genus, measure)
    df <- plot_richness(ps_genus, "Location", measures = "Shannon")$data
    mother<-df %>%
        filter(Stadium== "Mother") 
    mother_loc<-sort(unique(mother$Location))
    
    
    newborn<-df %>%
        filter(Stadium== "Newborn") 
    newborn_loc<-sort(unique(newborn$Location))
    
    shannon <- ggplot(df, aes(Location, value, fill = Type)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, size = 1) +
      facet_wrap(~Stadium, scales = "free_x") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle(measure) +
      labs(x =NULL, y = NULL)
    
    for (location_name in mother_loc){
      position_loc <- position(df, location_name, mother_loc)
      p_value <- wilcox.test(value ~ Type, data = df %>% filter(Location == location_name))$p.value
      p_value_star <- p_to_stars(p_value)
      if (p_value_star!="ns"){
            print(p_value_star)
    
      shannon <- add_clamr(shannon, position_loc, p_value_star)
        }
    }
    
    for (location_name in newborn_loc){
      position_loc <- position(df, location_name, newborn_loc)
      p_value <- wilcox.test(value ~ Type, data = df %>% filter(Location == location_name))$p.value
      p_value_star <- p_to_stars(p_value)
      if (p_value_star!="ns"){
        print(location_name)
        print(p_value_star)
        
        shannon <- add_clamr(shannon, position_loc, p_value_star)
        }
    }
    return(shannon)
}

shannon<-plot_richness_with_p_val(ps_genus, "Shannon")
    plots[[1]] <- shannon 

### chao


df <- plot_richness(ps_genus, "Location", measures = "Chao1")$data

# obliczamy pozycje x dla Skin osobno w każdym Stadium



loc_levels <- levels(factor(df$Location))
mother_placenta_pos <- df %>%
  filter(Location == "Mother placenta") %>%
  group_by(Stadium) %>%
  summarise(
    x = which(loc_levels == "Mother placenta"),
    y = max(value) * 1.1
  )



rectum_pos <- df %>%
  filter(Location == "Rectum") %>%
  group_by(Stadium) %>%
  summarise(
    x = as.numeric(factor(Location, levels = unique(df$Location[df$Stadium == Stadium]))),
    y = max(value) * 1.1
  )

placenta_pos <- df %>%
  filter(Location == "Placenta") %>%
  group_by(Stadium) %>%
  summarise(
    x = as.numeric(factor(Location, levels = unique(df$Location[df$Stadium == Stadium]))),
    y = max(value) * 1.1
  )



# p-value dla Placenta, Rectum i Mother Placenta
p_rectum <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Rectum"))$p.value
p_rectum_star <- p_to_stars(p_rectum)

p_placenta <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Placenta"))$p.value
p_placenta_star <- p_to_stars(p_placenta)

p_mother_placenta <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Mother placenta"))$p.value
p_mother_placenta_star <- p_to_stars(p_mother_placenta)

chao1 <- ggplot(df, aes(Location, value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~Stadium, scales = "free_x") +

  # klamra dla Skin w każdym panelu
  geom_segment(
    data = rectum_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = rectum_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = rectum_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = rectum_pos,
    aes(x = x, y = y * 1.05, label = p_placenta_star),
    inherit.aes = FALSE
  ) +
  # klamra dla Placenta w każdym panelu

 geom_segment(
    data = placenta_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = placenta_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = placenta_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = placenta_pos,
    aes(x = x, y = y * 1.05, label = p_placenta_star),
    inherit.aes = FALSE
  ) +

# klamra dla mother Placenta w każdym panelu

 geom_segment(
    data = mother_placenta_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = mother_placenta_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = mother_placenta_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = mother_placenta_pos,
    aes(x = x, y = y * 1.05, label = p_mother_placenta_star),
    inherit.aes = FALSE
  ) +


  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle("Chao1") +
  labs(x =NULL, y = NULL)

chao1


#####


# chao <- plot_richness(ps_genus, "Location", measures = c("Chao1"))
# chao + theme_bw() + geom_boxplot(aes(fill = Type)) +facet_wrap(~Stadium, scales="free_x")+ theme(axis.text.x = element_text(angle=90, hjust=1))
plots[[2]] <- chao1


#### simpson



df <- plot_richness(ps_genus, "Location", measures = "Simpson")$data

# obliczamy pozycje x dla Skin osobno w każdym Stadium
df$Location <- factor(df$Location, levels = sort(unique(df$Location)))

skin_pos <- df %>%
  filter(Location == "Skin") %>%
  group_by(Stadium) %>%
  summarise(
    x = as.numeric(factor(Location, levels = unique(df$Location[df$Stadium == Stadium]))),
    y = max(value) * 1.1
  )

#toto: fix x, aby był automatycznie wybierany
stool_pos <- df %>%
  filter(Location == "Stool") %>%
  group_by(Stadium) %>%
  summarise(
    x = 4,
    y = max(value) * 1.1
  )



# p-value dla Placenta, Rectum i Mother Placenta
p_skin <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Skin"))$p.value
p_skin_star <- p_to_stars(p_skin)

p_stool <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Stool"))$p.value
p_stool_star <- p_to_stars(p_stool)


simpson <- ggplot(df, aes(Location, value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~Stadium, scales = "free_x") +

  # klamra dla Skin w każdym panelu
  geom_segment(
    data = skin_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = skin_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = skin_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = skin_pos,
    aes(x = x, y = y * 1.05, label = p_skin_star),
    inherit.aes = FALSE
  ) +
  # klamra dla Placenta w każdym panelu

 geom_segment(
    data = stool_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = stool_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = stool_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = stool_pos,
    aes(x = x, y = y * 1.05, label = p_stool_star),
    inherit.aes = FALSE
  ) +



  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle("Simpson") +
  labs(x =NULL, y = NULL)

simpson
plots[[3]] <- simpson 
###

# simpson <- plot_richness(ps_genus, "Location", measures = c("Simpson"))
# simpson + theme_bw() + geom_boxplot(aes(fill = Type)) +facet_wrap(~Stadium, scales="free_x")+ theme(axis.text.x = element_text(angle=90, hjust=1))
# plots[[3]] <- simpson 


##### observed

df <- plot_richness(ps_genus, "Location", measures = "Observed")$data

# obliczamy pozycje x dla Skin osobno w każdym Stadium

loc_levels <- levels(factor(df$Location))
mother_placenta_pos <- df %>%
  filter(Location == "Mother placenta") %>%
  group_by(Stadium) %>%
  summarise(
    x = which(loc_levels == "Mother placenta"),
    y = max(value) * 1.1
  )


rectum_pos <- df %>%
  filter(Location == "Rectum") %>%
  group_by(Stadium) %>%
  summarise(
    x = as.numeric(factor(Location, levels = unique(df$Location[df$Stadium == Stadium]))),
    y = max(value) * 1.1
  )

placenta_pos <- df %>%
  filter(Location == "Placenta") %>%
  group_by(Stadium) %>%
  summarise(
    x = as.numeric(factor(Location, levels = unique(df$Location[df$Stadium == Stadium]))),
    y = max(value) * 1.1
  )


# p-value dla Placenta, Rectum i Mother Placenta
p_placenta <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Placenta"))$p.value
p_placenta_star <- p_to_stars(p_placenta)

p_rectum <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Rectum"))$p.value
p_rectum_star <- p_to_stars(p_rectum)


p_mother_placenta <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Mother placenta"))$p.value
p_mother_placenta_star <- p_to_stars(p_mother_placenta)

observed <- ggplot(df, aes(Location, value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~Stadium, scales = "free_x") +

  # klamra dla Skin w każdym panelu
  geom_segment(
    data = placenta_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = placenta_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = placenta_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = placenta_pos,
    aes(x = x, y = y * 1.05, label = p_placenta_star),
    inherit.aes = FALSE
  ) +
  # klamra dla Placenta w każdym panelu

 geom_segment(
    data = rectum_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = rectum_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = rectum_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = rectum_pos,
    aes(x = x, y = y * 1.05, label = p_rectum_star),
    inherit.aes = FALSE
  ) +
 # klamra dla Mother placenta w każdym panelu

 geom_segment(
    data = mother_placenta_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = mother_placenta_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = mother_placenta_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = mother_placenta_pos,
    aes(x = x, y = y * 1.05, label = p_mother_placenta_star),
    inherit.aes = FALSE
  ) +


  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("Observed") +
  labs(x =NULL, y = NULL)

observed


# observed <- plot_richness(ps_genus, "Location", measures = c("Observed"))
# observed + theme_bw() + geom_boxplot(aes(fill = Type)) +facet_wrap(~Stadium, scales="free_x")+ theme(axis.text.x = element_text(angle=90, hjust=1))
plots[[4]] <- observed 

### Fisher


df <- plot_richness(ps_genus, "Location", measures = "Fisher")$data

# obliczamy pozycje x dla Skin osobno w każdym Stadium


mother_placenta_pos <- df %>%
  filter(Location == "Mother placenta") %>%
  group_by(Stadium) %>%
  summarise(
    x =3,
    y = max(value) * 1.1
  )
rectum_pos <- df %>%
  filter(Location == "Rectum") %>%
  group_by(Stadium) %>%
  summarise(
    x = as.numeric(factor(Location, levels = unique(df$Location[df$Stadium == Stadium]))),
    y = max(value) * 1.1
  )

placenta_pos <- df %>%
  filter(Location == "Placenta") %>%
  group_by(Stadium) %>%
  summarise(
    x = as.numeric(factor(Location, levels = unique(df$Location[df$Stadium == Stadium]))),
    y = max(value) * 1.1
  )

loc_levels <- levels(factor(df$Location))
mother_placenta_pos <- df %>%
  filter(Location == "Mother placenta") %>%
  group_by(Stadium) %>%
  summarise(
    x = which(loc_levels == "Mother placenta"),
    y = max(value) * 1.1
  )


# p-value dla Placenta, Rectum i Mother Placenta
p_placenta <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Placenta"))$p.value
p_placenta_star <- p_to_stars(p_placenta)

p_rectum <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Rectum"))$p.value
p_rectum_star <- p_to_stars(p_rectum)


p_mother_placenta <- wilcox.test(value ~ Type, data = df %>% filter(Location == "Mother placenta"))$p.value
p_mother_placenta_star <- p_to_stars(p_mother_placenta)

fisher <- ggplot(df, aes(Location, value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~Stadium, scales = "free_x") +

  # klamra dla Skin w każdym panelu
  geom_segment(
    data = placenta_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = placenta_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = placenta_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = placenta_pos,
    aes(x = x, y = y * 1.05, label = p_placenta_star),
    inherit.aes = FALSE
  ) +
  # klamra dla Placenta w każdym panelu

 geom_segment(
    data = rectum_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = rectum_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = rectum_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = rectum_pos,
    aes(x = x, y = y * 1.05, label = p_rectum_star),
    inherit.aes = FALSE
  ) +
 # klamra dla Mother placenta w każdym panelu

 geom_segment(
    data = mother_placenta_pos,
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = mother_placenta_pos,
    aes(x = x - 0.2, xend = x - 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = mother_placenta_pos,
    aes(x = x + 0.2, xend = x + 0.2, y = y * 0.98, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = mother_placenta_pos,
    aes(x = x, y = y * 1.05, label = p_mother_placenta_star),
    inherit.aes = FALSE
  ) +


  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("Fisher") +
  labs(x =NULL, y = NULL)

fisher

# fisher <- plot_richness(ps_children_family, "Location", measures = c("Fisher"))
# fisher + theme_bw() + geom_boxplot(aes(fill = Type)) +facet_wrap(~Stadium, scales="free_x")+ theme(axis.text.x = element_text(angle=90, hjust=1))
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
