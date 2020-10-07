# 1. Install packages ####
devtools::install_github("hadley/devtools")
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("dada2", "ShortRead", "Biostrings", "phyloseq", "DESeq2", "microbiome", "DECIPHER", "phangorn"))

devtools::install_github("bryandmartin/corncob")

install.packages(c("vegan", "ggplot2", "lmerTest", "lme4", "tibble", "car", "rcompanion", "emmeans", "RVAideMemoire"))

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("dada2", "ShortRead", "Biostrings", "phyloseq", "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4", "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire")
ipak(packages)
# at this step, you should see 'TRUE' under each package, if you don't, you'll need to run install.packages("packagename") and troubleshoot the install

# 2. set working directory ####
setwd("~/16S Sequences")
path<- setwd("~/Downloads/TriSocieties-Amplicon-Workshop")

# check to make sure the fastq files are in the directory
list.files(path)


# 3. Sort the forward and reverse reads #####
# read in the file names and store as F.sam and R.sam
F.sam <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
#You should see in your global env, fnFs chr[1:20] "and the path directory"
R.sam <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
#this is depedent on how your samples are named

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(F.sam), "_"), `[`, 1)
sample.names
#should see a list of the sample names


# 4. Examine the quality profiles ####

## forward reads ##

plotQualityProfile(F.sam[1:4])
# gray scale heatmap is the frequencey of each quality score at each bp p
# green line represents the mean quality score 
# orange lines represent the quartiles of the quality scores
# red line shows how many reads extend to at least that bp position

## reverse reads ##

plotQualityProfile(R.sam[1:4])
# do not be alarmed if the reverse reads have a poorer qualify profile, this is expected because of how the reads are sequenced with Illumina
# we can account for this when we filter and trim the reads

# 5. Assign where filtered samples will be sent and Filter samples####

# This will place the filtered fastq files into the 'filtered' folder 
filt_path <- file.path(path, "filtered") #you should see this in your global environment
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter forward and reverse reads
# On Windows set multithread=FALSE

# the commands included in this filter and trim are dependent on your specific data- especially the truncLen- because you want to trim your seqs based on the quality profiles you generated- quality scores of 30 + are good
# so you should choose the trunclen to reflect the quality scores of your data but you can't cut too much because you need the forward and reverse reads to overlap when you merge them later
# this information is helpful if you're using 515F/806R primers in determining bp overlap https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/appnote_miseq_16S.pdf

out <- filterAndTrim(F.sam, filtFs, R.sam, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     compress=TRUE, multithread=TRUE) 
# set multithread to FALSE if using windows
# this step takes awhile to perform
# setting maxN = 0 means that 
# truncQ (default is 2), truncates reads where quality scores drop less than or equal to the number specified
# maxEE, reads with higher expected errors than what is specified are discarded
# setting compress to TRUE will gzip the output files 

head(out) # you should see an abbreviated list of your samples, how many reads were used in the filter and trim command and how many reads were removed after the command
# if you see that a lot of your reads were removed, you may need to be less conservative in the filterandTrim command


# 6. Train Dada2 to assess errors ####
# The program has to estimate the sequencing error rate, to try and distinguish errors from true Exact Sequence Variants
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

## take a look at the modeled errors

plotErrors(errF, nominalQ=TRUE)
# what you should see is that your samples track well with the predicted errors
# you should also see that as the consensus quality score goes up, the error frequencey goes down


# 7. Dereplicate the filtered fastq files ####
# we want to dereplicate the sequences to decrease computation time for the next step
derepFs <- derepFastq(filtFs, verbose=TRUE)
#you should see "encountered XXXX unique sequences from XXX total sequences read
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# 8. Infer the sequence variants in each sample #####
# This command goes through every sequence in your sample and removes sequences that have sequencing errors
# the estimated error rates from step 5 are important here

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


dadaFs[[1]]
# output says "XXX sequence variants were inferred from XXXX input unique sequences"
# the difference in these values is the number of sequences that were removed because they were likely to include sequencing errors 
dadaRs[[1]]

# 9.  Merge the denoised forward and reverse reads: ####

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#You should see an "abundance", "forward", "reverse", "nmismatch", "nindel", "prefer" and "accept" column
#nmatch tells you how many base pairs overlapped between the forward and reverse reads- you want this number to be fairly high (50, 100) 
#you should see 0's in all the rows under "nmismatch" and "nindels"


#We can now construct a sequence table of our samples, a #
#higher-resolution version of the ESV table produced by traditional methods.#

seqtab <- makeSequenceTable(mergers)
#you will get a message "The sequences being tabled vary in length" but ignore

dim(seqtab)
#there are two values here, the first value tells you how many samples you have, which should be 20 for our samples, and the second number tells you how many ESVs you have

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))
#you will see how variable the length of your sequences are and you'll have to decide how many different lengths you're willing to include


#remove sequences less than 289 and greater than 296 bp in length
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(240,255)] 
#this removes sequences shorter than 240 or longer than 255
table(nchar(getSequences(seqtab2)))

# 10.  Remove chimeric sequences ####

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
#should get an output that says "Identified XXX bimeras out of XXX input seqs"
#so it filtered out the bimeras (chimeras) from the table

dim(seqtab.nochim)
#So, this is saying that there are 151 samples and a total of 35828 ESVs after chimeric filtering

# Percent of sequences that are non-chimeric
sum(seqtab.nochim)/sum(seqtab)
# 95.7 non chimeric


# 11. Combine sequence reductions by sample from each step #####
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
# this displays how many sequences you started with in each sample, how many were filtered by the filterAndTrim command, how many were merged, and how many were non-chimeric
# Show sequence reductions throughout pipeline


# 12. Assign taxonomy ####
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr99_v138_train_set.fa", multithread=TRUE)
# check your taxonomic classifications #
taxa.print<- taxa
rownames(taxa.print)
head(taxa.print)


# 13. Make a Phyloseq object ####

meta <- read.table("Metadata.txt", header = TRUE, row.names = 1)
asv.table<- otu_table(seqtab.nochim, taxa_are_rows=FALSE)

#Now we can make the phyloseq object
ps <- phyloseq(asv.table, tax_table(taxa), sample_data(meta))

# 14. Make a phylogenetic tree ####
# to do this, I follow the pipeline outlined in this paper: https://f1000research.com/articles/5-1492
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
# the above step performed a multiple sequence alignment
# we will use this to create a neighbor-joining tree, 

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

# we use this NJ tree to make a Generalized time-reversible with Gamma rate variation) maximum likelihood tree

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

# now we can add the tree to our phyloseq object
ps <- merge_phyloseq(ps, phy_tree(fitGTR$tree))

# 15. A little more data formatting ####
# we want to duplicate the phyloseq object, so that we can have 1 original and edit the other
ps.1 <- ps

dna <- Biostrings::DNAStringSet(taxa_names(ps.1))
names(dna) <- taxa_names(ps.1)
ps.1 <- merge_phyloseq(ps.1, dna)
taxa_names(ps.1) <- paste0("ASV", seq(ntaxa(ps.1)))
ps.1
# Now, let's make sure that the correct information is included in our phyloseq object

summary(ps.1@otu_table) # Should include asv info
ps.1@tax_table # Should include taxonomic info
ps.1@sam_data # Should reflect the mapping file that we imported

ps.2 = subset_taxa(ps.1, Kingdom == "Bacteria" | Kingdom == "Archaea")
# Removes anything not assigned to Bacteria or Archaea

# let us see if there are any samples that have a small number of sequences
sample_sums(ps.2)
# while it doesn't seem to be the case, if it were we would use the following lines to remove samples with low sequence counts
set.seed(500)
ps.pruned <- prune_samples(sample_sums(ps.2)>=2000, ps.2)
# check to see how many samples were removed 

# in some cases you may want to have a relative abundance OTU table, rather than one with absolute abundances
# the following line of code transforms the absolute abundance data into relative abundance data
ps.perc <- transform_sample_counts(ps.pruned, function(x) x / sum(x)) 

# 16. Alpha-Diversity Estimates and Plotting####
# It may be appropriate for you to rarefy the dataset before computing alpha-diversity estimates
ordered(sample_sums(ps.pruned))
# now that you can look at the range of total sequences per sample you can get a better idea of whether of not a dataset should be rarefied
# since our total sequences per sample are all within an order of magnitude from eachother, we will not rarefy
# let's calcuolate Shannon diversity and Richness (Observed)
alpha.div<-estimate_richness(ps.pruned, measures=c("Shannon", "Observed"))
even <- evenness(ps.pruned, 'pielou')
# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses
meta$Shannon <- paste(alpha.div$Shannon)
meta$Observed <- paste(alpha.div$Observed)
meta$Evenness <- even$pielou
meta$Observed <- as.numeric(meta$Observed)
meta$Shannon <- as.numeric(meta$Shannon)
meta$Block <- as.factor(meta$Block)
meta$Line <- as.factor(meta$Line)
# the below commands will perform anova's based on flavonoid production and plant genotype
# Note that I have included Block as a random factor here, thus I have to make a linear mixed effects model before running the anova
# make sure lmerTest is a loaded package or you will not see p-values
shannon.model<-lmer(Shannon ~ Line * Flavs  + (1|Block), data = meta)

# check assumptions for linear modeling
# 1. linearity
# plot fitted versus residuals
plot(shannon.model)
# plot the observed values versus the model values
# should see no pattern
# 2. homogeneity of variance
leveneTest(residuals(shannon.model) ~ meta$Line)
# if p <0.05, assumption is not met

# 3. normal distribution of residuals
qqmath(shannon.model)
# values should ~ follow the line

# assuming that the assumptions are satisfied, run the anova

anova(shannon.model) 
# note that this is a Type III Analysis of variance table

ggplot(meta, aes(Shannon, Line, fill = Line)) + geom_boxplot() + theme_classic()

# Evenness #
even.model<-lmer(Evenness ~ Line * Total_Phenolics + (1|Block), data = meta)
plot(even.model)
leveneTest(residuals(even.model) ~ meta$Line)
qqmath(even.model)

anova(even.model)

ggplot(meta, aes(Evenness, Line, fill = Line)) + geom_boxplot() + theme_classic() + 
  labs(axis.text.x = element_text(angle = -90, hjust = 0))

ggplot(meta, aes(Total_Phenolics, Evenness)) + geom_line() + theme_classic() + 
  labs(axis.text.x = element_text(angle = -90, hjust = 0))
ggplot(meta, aes(Total_Phenolics, Evenness, colour = Line)) + geom_line() + theme_classic() + 
  labs(axis.text.x = element_text(angle = -90, hjust = 0))

# post-hoc test using estimated marginal means
emmeans(even.model, pairwise ~ Line, adjust = "none")
# indicates no differences between lines

# if covariate, then we can use estimated marginal means of linear trends
emtrends(even.model, pairwise ~ Line, var = "Total_Phenolics")
# indicates that Lines do not have a different trend of Total Phenolics ~ Richness 
# we can plot this like so:
emmip(even.model, Line ~ Total_Phenolics, cov.reduce = range)
# 17. Heatmaps ####
# this command, taxa_level, comes from the microbiomeSeq package and should be cited if you use this
taxa_level <- function(physeq,which_level){
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  OTU <- otu_table(physeq)
  SAM <- sample_data(physeq)
  OTU_taxonomy <- tax_table(physeq)
  new_abund_table<-NULL
  if(which_level=="Otus"){
    OTU_tree <- phy_tree(physeq)
    new_abund_table<-OTU
  } else {
    list<-na.omit(unique(OTU_taxonomy[,which_level]))
    new_abund_table<-NULL
    for(i in list){
      rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
      tmp<-data.frame(rowSums(OTU[,rt]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  OTU<-as.data.frame(as(new_abund_table,"matrix"))
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(SAM)
  #reconstruct the phyloseq object
  physeq<-NULL
  if(which_level=="Otus"){
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }
  return(physeq)
}

ps.phylum<- taxa_level(ps.pruned, "Phylum") #this calculates the abundances of each phylum across every sample into a new phyloseq object
plot_heatmap(ps.phylum, "PCoA", "bray", "Line")

# Use this to look at genera # 
ps.genera.perc<- taxa_level(ps.pruned, "Genus")

plot_heatmap(ps.genera.perc, "PCoA", "bray", "Line", na.value = "grey",  low= "#66CCFF", high = "#000033" )

# 18. Beta-Diversity Estimates####
# there are two main things that you will likely want to analyze, a permanova and homogeneity of group dispersions 
# adonis runs a permutational MANOVA. One major assumption for this analysis is that within group similarity is similar across groups
# betadisper runs PERMDISP2 to calculate group dispersions
# we will start with permdisp and then run adonis, after creating a basic ordination plot

plot_ordination(ps.pruned, ordinate(ps.pruned, "PCoA", "bray"), color = "Line") + theme_classic() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)

bc.asv <- phyloseq::distance(otu_table(ps.pruned), "bray")
bc.phyl <- phyloseq::distance(otu_table(ps.phylum), "bray")
unifrac <- UniFrac(ps.pruned, weighted=TRUE, normalized=TRUE)

# betadisper
disp.phyl <- betadisper(bc.phyl, meta$Line)
anova(disp.phyl)

disp.uni <- betadisper(unifrac, meta$Line)
anova(disp.uni)

disp.asv <- betadisper(bc.asv, meta$Line)
disp.asv
anova(disp.asv)
# if the anova had a p <0.05, perform post-hoc analysis like so:

TukeyHSD(disp.asv) 
# creates confidence intervals between means of levels of a factor 
plot(disp.asv)

# adonis
adonis(bc.phyl ~ Line*Total_Phenolics, data = meta)

adonis(bc.asv ~ Line, data = meta)
# if adonis is significant, move on to post-hoc test
pairwise.perm.manova(bc.asv, meta$Line, nperm = 999)
# what about phylogenetic distance?
adonis(unifrac ~ Line, data = meta)

pairwise.perm.manova(unifrac, meta$Line, nperm= 999, p.method = "BH")

# 19. Plotting Relative Abundance Bar Charts####
# phylum-level
ps.phyla.perc <- taxa_level(ps.perc, "Phylum")

# identify the 10 most abundant phylum
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:10])
# now we can use the list of the top phylum and subset those from the phyloseq object
ps.phylum.10 <- prune_taxa(phylum.10, ps.phyla.perc)

melt.phylum <- psmelt(ps.phylum.10)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

ggplot(melt.phylum, aes(x = Line, y = Abundance, fill = OTU)) + theme_classic() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + labs(fill = "Phylum")

# do you notice anything ...unpleasant about this plot? - relative abundances are used and should sum to 1...

phyl.means <- aggregate(Abundance~Line+OTU, melt.phylum, FUN=mean)

ggplot(phyl.means, aes(x = Line, y = Abundance, fill = OTU)) + theme_classic() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette)

# bar plot at the genus-level
ps.gen <- taxa_level(ps.perc, "Genus")

gen.10 <- names(sort(taxa_sums(ps.gen), TRUE)[1:10])
ps.gen.10 <- prune_taxa(gen.10, ps.gen)
melt.gen <- psmelt(ps.gen.10)

gen.means <- aggregate(Abundance~Line+OTU, melt.gen, FUN=mean)

ggplot(gen.means, aes(x = Line, y = Abundance, fill = OTU)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + theme_classic()

# is there a difference in abundance of Pseudomonas between the Lines?
# extract Pseudomonas abundances from melt.gen file
pseud <- melt.gen[ which(melt.gen$OTU=='Pseudomonas'), ]

# create model
pseud.lm <- lm(Abundance ~ Line, data = pseud)
anova(pseud.lm)

emmeans(pseud.lm, pairwise ~ Line)

# 20. Mantel Tests ####
asv.table <- data.frame(otu_table(ps.pruned))

xdist.prod <- vegdist(asv.table, method = "bray")
ydist.prod <- vegdist(meta$DPPH, method = "euclid")

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.prod, method = "spear")
# mantel test with weighted UniFrac distance
vegan::mantel(unifrac, ydist.prod, method = "spear")
# we can also include more than 1 environmental variable
ydist.prod <- vegdist(meta$DPPH*meta$Flavs, method = "euclid")
vegan::mantel(unifrac, ydist.prod, method = "spear")


# 21. Differential Abundance Analysis Corn Cob ####
# differential abundance analyses compare 2 groups 
# Since we have three groups, we need to split the data 
ps.N321.N334 <- subset_samples(ps.pruned, Line == "N321" | Line == "N334")
ps.N321.N336 <- subset_samples(ps.pruned, Line == "N321" | Line == "N336")
ps.N334.N336 <- subset_samples(ps.pruned, Line == "N334" | Line == "N336")

ex1 <- differentialTest(formula = ~ Line+ Block,
                        phi.formula = ~ 1,
                        formula_null = ~ Block,
                        phi.formula_null = ~ 1,
                        data = ps.N321.N334,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex1)

ex2 <- differentialTest(formula = ~ Line+ Block,
                        phi.formula = ~ 1,
                        formula_null = ~ Block,
                        phi.formula_null = ~ 1,
                        data = ps.N321.N336,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex2)

ex3 <- differentialTest(formula = ~ Line+ Block,
                        phi.formula = ~ 1,
                        formula_null = ~ Block,
                        phi.formula_null = ~ 1,
                        data = ps.N334.N336,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex3)

# how do you know which Line has more of ASV107?

plot(asv.table$ASV107 ~ meta$Line)
aggregate(asv.table$ASV107~Line, meta, FUN=mean)

