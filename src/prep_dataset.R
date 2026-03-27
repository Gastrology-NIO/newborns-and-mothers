
path<-"~/data02/Seq_AZ/2025_noworodki"
folder<-"/data/Noworodki\ 281024\ ch1/"
ch1_folder<-paste0(path, folder)
# bam <- list.files(ch1_folder, pattern = "bam", full.names = TRUE)
list.files("./data_fastq/", pattern = "bam", full.names = TRUE) -> bam


out <- filterAndTrim(bam, bam_filtered) 
erc_removed<- 100*(head(out)[,1]-head(out)[,2])/head(out)[,2]

list.files(bam_filtered)
sample.names <- sapply(strsplit(basename(bam), ".Ion"), `[`, 1)

fastq_r1_filtered <- file.path(bam_filtered, list.files(bam_filtered))
# fastq_r1_filtered <- file.path(path, "filtered", paste0(sample.names, ".fastq"))

# https://github.com/benjjneb/dada2/issues/795
err_r1 <- learnErrors(fastq_r1_filtered, multithread=TRUE)
dada_r1 <- dada(fastq_r1_filtered, err=err_r1, multithread=TRUE, pool="pseudo")



seqtab <- makeSequenceTable(dada_r1)
dim(seqtab) 
table(nchar(getSequences(seqtab)))

#Remove chimeras <- to jest raczej tu nie potrzebne
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#> dim(seqtab.nochim)
#[1]    6 4889
#> sum(seqtab.nochim)/sum(seqtab)
#[1] 0.3611021
plotErrors(err_r1, nominalQ=TRUE)
ggsave(
  "/data02/Seq_AZ/2025_noworodki/error_correction_R1.png",
  plot = last_plot(),
)
plotErrors(err_r2, nominalQ=TRUE)
ggsave(
  "/data02/Seq_AZ/2025_noworodki/error_correction_R2.png",
  plot = last_plot(),
)

#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_r1, getN), sapply(dada_r2, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("./SILVA_SSU_r138_2019.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest

ranks <- c("domain", "genus") 
taxid_tmp <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa <- taxid
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
samples.out <- rownames(seqtab.nochim)
location<-sapply(strsplit(sample.names, "_"), `[`, c(3))
sample.names <- sapply(strsplit(basename(bam), ".Ion"), `[`, 1)
organism<-sapply(strsplit(sample.names, "(_zoladek|_odbyt)"), `[`, c(1))
samdf <- data.frame(Location=location, Organism_id=organism, Tries=sample.names)
rownames(samdf) <- samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)

sample_data(ps)$is_control<-F
sample_data(ps)$is_control[sample_data(ps)$Names %in% c("KN.IonXpress_071.bam.fastq", "KP_AT1_IonXpress_007.bam.fastq")]<-T
colnames(sample_data(ps))<-c("Sample_id","sample_type","Organism_id","Tries","Stadium","Names","is_control")
tmp_sample<-as.matrix(sample_data(ps))

tmp_sample<-as.matrix(sample_data(ps)[,c('is_control','sample_type')])
tmp_sample[,2][is.na(tmp_sample[,2])]<-'Control'
tmp_sample<-as.data.frame(tmp_sample)
tmp<-as.matrix(otu_table(ps))
tmp<-as.data.frame(tmp)
orig.composition_results = micRoclean(counts = tmp,
                                      meta = tmp_sample,
                                      research_goal = 'orig.composition',
                                      control_name = 'Control')


orig.composition_results$decontaminated_count
sum(orig.composition_results$decontaminated_count)/sum(otu_table(ps))
ps2<- ps
otu_table(ps2)<-otu_table(orig.composition_results$decontaminated_count, taxa_are_rows=F) 



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
sample_data(ps_all)$Stadium[sample_data(ps_all)$Stadium %in% c("_dziecko_2", "_noworodek_1", "_noworodek")]<-"Newborn"
sample_data(ps_all)$Stadium[sample_data(ps_all)$Stadium %in% c("_matka", "_matka_1", "_matka_2")]<-"Mother"
sample_data(ps_all)$Type<-sapply(strsplit(sample_data(ps_all)$Organism_id, "_"), `[`, c(2))
ps_genus<-tax_glom(ps_all, "genus")
sample_data(ps_genus)$Type[sample_data(ps_genus)$Type=="K"]<-"TP"
sample_data(ps_genus)$Type[sample_data(ps_genus)$Type=="B"]<-"LP"
sample_data(ps_genus)$Stadium[sample_data(ps_genus)$Stadium =="Newborn"]<-"Neonate"
