sample_data(ps2)$Stadium[sample_data(ps2)$Stadium == "dziecko"] <- "Newborn"
sample_data(ps2)$Stadium[sample_data(ps2)$Stadium == "matka"] <- "Mother"
sample_data(ps2)$Location[sample_data(ps2)$sample_type  == "skora"] <- "Skin"
sample_data(ps2)$Location[sample_data(ps2)$sample_type  == "zoladek"] <- "Stomach"
sample_data(ps2)$Location[paste0(sample_data(ps2)$sample_type,sample_data(ps2)$Stadium ) == "lozysko_noworodek_1"] <- "Placenta"
sample_data(ps2)$Location[paste0(sample_data(ps2)$sample_type,sample_data(ps2)$Stadium ) == "lozysko_matka_1"] <- "Mother placenta"
sample_data(ps2)$Location[sample_data(ps2)$sample_type  == "policzek"] <- "Cheek"
sample_data(ps2)$Location[sample_data(ps2)$sample_type  == "odbyt"] <- "Rectum"
sample_data(ps2)$Location[sample_data(ps2)$sample_type  == "stolec"] <- "Stool"
sample_data(ps2)$Location[sample_data(ps2)$sample_type  == "szyjka"] <- "Cervix"

ps_all <- subset_samples(ps2, Stadium %in% c("_noworodek_1", '_noworodek', "_matka", "_matka_1"))
