rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Edinburgh Honours Project Write-Up')

####Libraries####
library(car)
library(beeswarm)

####Loading Datasets####
load('Honours_Project- Kai/Honours_Project--main/flora.Rda')

flora <- flora[complete.cases(flora$species),]
flora <- flora[complete.cases(flora$msystem),]

genome <- flora[complete.cases(flora$DNA2C),]
range <- flora[complete.cases(flora$total_hectads_2010_19),]

####Genome Complete Analysis####
#Analysis
shapiro.test(genome$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(genome$DNA2C)) #p<0.05 therefore non-normal
leveneTest(DNA2C~msystem, data=genome) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=genome) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=genome)
boxplot(log(DNA2C)~msystem, data=genome, add=T, 
        col="#0000ff22")

####Genome Family Specific Analysis####
genome_family <- genome[complete.cases(genome$family),]
family_freq <- as.data.frame(table(genome_family$family))
large_fams <- family_freq[which(family_freq$Freq>=20),]

##Apiaceae
Apiaceae <- genome_family[which(genome_family$family=="Apiaceae"),]
table(Apiaceae$msystem)
#Analysis
shapiro.test(Apiaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Apiaceae$DNA2C)) #p<0.05 therefore non-normal
leveneTest(DNA2C~msystem, data=Apiaceae) #p<0.05 so variances aren't homogenous
leveneTest(log(DNA2C)~msystem, data=Apiaceae) #p<0.05 so variances aren't homogenous
#kruskal.test(DNA2C~msystem, data=Apiaceae)
#Can't test because of non-homogenous variances, prob. because there's only 2 selfing species
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Apiaceae)
boxplot(log(DNA2C)~msystem, data=Apiaceae, add=T, 
        col="#0000ff22")

##Asteraceae
Asteraceae <- genome_family[which(genome_family$family=="Asteraceae"),]
table(Asteraceae$msystem)
#Analysis
shapiro.test(Asteraceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Asteraceae$DNA2C)) #p>0.05 therefore normal
aster_anova <- aov(log(DNA2C)~msystem, data=Asteraceae)
summary(aster_anova) #p>0.05 so no sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(DNA2C~msystem, data=Asteraceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Asteraceae) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Asteraceae)
boxplot(log(DNA2C)~msystem, data=Asteraceae, add=T, 
        col="#0000ff22")

##Brassicaceae
Brassicaceae <- genome_family[which(genome_family$family=="Brassicaceae"),]
table(Brassicaceae$msystem)
#Analysis
shapiro.test(Brassicaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Brassicaceae$DNA2C)) #p<0.05 therefore non-normal
leveneTest(DNA2C~msystem, data=Brassicaceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Brassicaceae) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Brassicaceae)
boxplot(log(DNA2C)~msystem, data=Brassicaceae, add=T, 
        col="#0000ff22")

##Caryophyllaceae
Caryophyllaceae <- genome_family[which(genome_family$family=="Caryophyllaceae"),]
table(Caryophyllaceae$msystem)
#Analysis
shapiro.test(Caryophyllaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Caryophyllaceae$DNA2C)) #p>0.05 therefore normal
caryo_anova <- aov(log(DNA2C)~msystem, data=Caryophyllaceae)
summary(caryo_anova) #p<0.05 so sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(DNA2C~msystem, data=Caryophyllaceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Caryophyllaceae) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(log(Caryophyllaceae$DNA2C), Caryophyllaceae$msystem, p.adjust.method = "BH")
#Mean and median gen. cross higher than mixed or gen. self
#Mean and median mixed lower than gen. cross and higher than gen. self
#Mean and median gen. self lower than gen. cross
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Caryophyllaceae)
boxplot(log(DNA2C)~msystem, data=Caryophyllaceae, add=T, 
        col="#0000ff22")

##Cyperaceae
Cyperaceae <- genome_family[which(genome_family$family=="Cyperaceae"),]
table(Cyperaceae$msystem) #No gen. selfers

##Fabaceae
Fabaceae <- genome_family[which(genome_family$family=="Fabaceae"),]
table(Fabaceae$msystem)
#Analysis
shapiro.test(Fabaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Fabaceae$DNA2C)) #p<0.05 therefore non-normal
leveneTest(DNA2C~msystem, data=Fabaceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Fabaceae) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Fabaceae)
boxplot(log(DNA2C)~msystem, data=Fabaceae, add=T, 
        col="#0000ff22")

##Lamiaceae
Lamiaceae <- genome_family[which(genome_family$family=="Lamiaceae"),]
table(Lamiaceae$msystem)
#Analysis
shapiro.test(Lamiaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Lamiaceae$DNA2C)) #p>0.05 therefore normal
lami_anova <- aov(log(DNA2C)~msystem, data=Lamiaceae)
summary(lami_anova) #p>0.05 so no sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(DNA2C~msystem, data=Lamiaceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Lamiaceae) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Lamiaceae)
boxplot(log(DNA2C)~msystem, data=Lamiaceae, add=T, 
        col="#0000ff22")

##Plantaginaceae
Plantaginaceae <- genome_family[which(genome_family$family=="Plantaginaceae"),]
table(Plantaginaceae$msystem)
#Analysis
shapiro.test(Plantaginaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Plantaginaceae$DNA2C)) #p<0.05 therefore non-normal
leveneTest(DNA2C~msystem, data=Plantaginaceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Plantaginaceae) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(Plantaginaceae$DNA2C, Plantaginaceae$msystem, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Plantaginaceae)
boxplot(log(DNA2C)~msystem, data=Plantaginaceae, add=T, 
        col="#0000ff22")
#Median gen. cross higher than mixed
#Median gen. self not different from gen. cross or mixed

##Poaceae
Poaceae <- genome_family[which(genome_family$family=="Poaceae"),]
table(Poaceae$msystem)
#Analysis
shapiro.test(Poaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Poaceae$DNA2C)) #p>0.05 therefore normal
poa_anova <- aov(log(DNA2C)~msystem, data=Poaceae)
summary(poa_anova) #p>0.05 so no sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(DNA2C~msystem, data=Poaceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Poaceae) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Poaceae)
boxplot(log(DNA2C)~msystem, data=Poaceae, add=T, 
        col="#0000ff22")

##Ranunculaceae
Ranunculaceae <- genome_family[which(genome_family$family=="Ranunculaceae"),]
table(Ranunculaceae$msystem)
#Analysis
shapiro.test(Ranunculaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Ranunculaceae$DNA2C)) #p<0.05 therefore normal
leveneTest(DNA2C~msystem, data=Ranunculaceae) #p>0.05 so variances are homogenous
kruskal.test(DNA2C~msystem, data=Ranunculaceae) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Ranunculaceae)
boxplot(log(DNA2C)~msystem, data=Ranunculaceae, add=T, 
        col="#0000ff22")

##Rosaceae
Rosaceae <- genome_family[which(genome_family$family=="Rosaceae"),]
table(Rosaceae$msystem)
#Analysis
shapiro.test(Rosaceae$DNA2C) #p<0.05 therefore non-normal
shapiro.test(log(Rosaceae$DNA2C)) #p>0.05 therefore normal
rosa_anova <- aov(log(DNA2C)~msystem, data=Rosaceae)
summary(rosa_anova) #p<0.05 so sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(DNA2C~msystem, data=Rosaceae) #p<0.05 so variances aren't homogenous
leveneTest(log(DNA2C)~msystem, data=Rosaceae) #p>0.05 so variances aren't homogenous
kruskal.test(log(DNA2C)~msystem, data=Rosaceae) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(log(Caryophyllaceae$DNA2C), Caryophyllaceae$msystem, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(log(DNA2C)~msystem, data=Caryophyllaceae)
boxplot(log(DNA2C)~msystem, data=Caryophyllaceae, add=T, 
        col="#0000ff22")
#Mean and median gen. cross higher than mixed and gen. self
#Mean and median mixed lower than gen. cross and higher than gen. self
#Mean and median gen. self lower than gen. cross and mixed


####Genome Pair Analysis####
self_cross_pairs_genome <- data.frame(self_species=c("Briza minor",
                                                     "Callitriche brutia",
                                                     "Cephalanthera damasonium",
                                                     "Juncus squarrosus",
                                                     "Luzula campestris",
                                                     "Papaver dubium",
                                                     "Senecio vulgaris",
                                                     "Trifolium dubium",
                                                     "Trifolium fragiferum",
                                                     "Trifolium glomeratum",
                                                     "Veronica polita"),
                                      cross_species=c("Briza media",
                                                      "Callitriche stagnalis",
                                                      "Cephalanthera longifolia",
                                                      "Juncus bufonius",
                                                      "Luzula forsteri",
                                                      "Papaver rhoeas",
                                                      "Senecio squalidus",
                                                      "Trifolium campestre",
                                                      "Trifolium hybridum",
                                                      "Trifolium repens",
                                                      "Veronia chamaedrys"),
                                      genus=c("Briza",
                                              "Callitriche",
                                              "Cephalanthera",
                                              "Juncus",
                                              "Luzula",
                                              "Papaver",
                                              "Senecio",
                                              "Trifolium",
                                              "Trifolium",
                                              "Trifolium",
                                              "Veronica"),
                                      self_2c=c(mean(flora$DNA2C[which(flora$species=="Briza minor")]),
                                                mean(flora$DNA2C[which(flora$species=="Callitriche brutia")]),
                                                mean(flora$DNA2C[which(flora$species=="Cephalanthera damasonium")]),
                                                mean(flora$DNA2C[which(flora$species=="Juncus squarrosus")]),
                                                mean(flora$DNA2C[which(flora$species=="Luzula campestris")]),
                                                mean(flora$DNA2C[which(flora$species=="Papaver dubium")]),
                                                mean(flora$DNA2C[which(flora$species=="Senecio vulgaris")]),
                                                mean(flora$DNA2C[which(flora$species=="Trifolium dubium")]),
                                                mean(flora$DNA2C[which(flora$species=="Trifolium fragiferum")]),
                                                mean(flora$DNA2C[which(flora$species=="Trifolium glomeratum")]),
                                                mean(flora$DNA2C[which(flora$species=="Veronica polita")])),
                                      cross_2c=c(mean(flora$DNA2C[which(flora$species=="Briza media")]),
                                                 mean(flora$DNA2C[which(flora$species=="Callitriche stagnalis")]),
                                                 mean(flora$DNA2C[which(flora$species=="Cephalanthera longifolia")]),
                                                 mean(flora$DNA2C[which(flora$species=="Juncus bufonius")]),
                                                 mean(flora$DNA2C[which(flora$species=="Luzula forsteri")]),
                                                 mean(flora$DNA2C[which(flora$species=="Papaver rhoeas")]),
                                                 mean(flora$DNA2C[which(flora$species=="Senecio squalidus")]),
                                                 mean(flora$DNA2C[which(flora$species=="Trifolium campestre")]),
                                                 mean(flora$DNA2C[which(flora$species=="Trifolium hybridum")]),
                                                 mean(flora$DNA2C[which(flora$species=="Trifolium repens")]),
                                                 mean(flora$DNA2C[which(flora$species=="Veronica chamaedrys")])),
                                      pair=c("B. minor - B. media",
                                             "C. brutia - C. stagnalis",
                                             "C. damasonium - C. longifolia",
                                             "J. squarrosus - J. bufonius",
                                             "L. campestris - L. forsteri",
                                             "P. dubium - P. rhoeas",
                                             "S. vulgaris - S. squalidus",
                                             "T. dubium - T. campestre",
                                             "T. fragigerum - T. hybridum",
                                             "T. glomeratum - T. repens",
                                             "V. polita - V. chamaedrys"))
self_cross_pairs_genome$difference_2c <- self_cross_pairs_genome$self_2c-self_cross_pairs_genome$cross_2c

library(ggplot2)
ggplot(data= self_cross_pairs_genome, aes(x=pair, y=difference_2c, fill=genus))+
  geom_bar(stat="identity")+
  xlab("")+
  ylab("Difference in 2C DNA Content (pg)")+
  labs(fill="Genus")+
  theme(axis.text.x=element_text(angle=90, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.text=element_text(face="italic"))

t.test(self_cross_pairs_genome$difference_2c, mu=0)
#p>0.05, so mean difference between pairs doesn't differ sig. from 0

####Genome Family Difference Plot####
fam_nums <- as.data.frame(table(genome$family))

mean_fam_2c_cs <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(mean_fam_2c_cs) <- c("family", "cross_av_2c", "self_av_2c")
for(i in 1:nrow(fam_nums)){
  if(fam_nums$Freq[i]>1){
    fam <- fam_nums$Var1[i]
    temp <- as.data.frame(table(genome$msystem[which(genome$family==fam_nums$Var1[i])]))
    if(nrow(as.data.frame(temp$Freq[which(temp$Var1=="generally cross")]))!=0&
       nrow(as.data.frame(temp$Freq[which(temp$Var1=="generally self")]))!=0){
      if(temp$Freq[which(temp$Var1=="generally cross")]>=1&
       temp$Freq[which(temp$Var1=="generally self")]>=1){
      family <- unique(genome$family[which(genome$family==fam_nums$Var1[i])])
      cross_av_2c <- mean(genome$DNA2C[which(genome$family==fam&genome$msystem=="generally cross")])
      self_av_2c <- mean(genome$DNA2C[which(genome$family==fam&genome$msystem=="generally self")])
      temp2 <- data.frame(t(as.data.frame(c(family, cross_av_2c, self_av_2c))))
      colnames(temp2) <- c("family", "cross_av_2c", "self_av_2c")
      mean_fam_2c_cs <- rbind(mean_fam_2c_cs, temp2)
    }}
  }
}

mean_fam_2c_cs$self_av_2c <- as.numeric(mean_fam_2c_cs$self_av_2c)
mean_fam_2c_cs$cross_av_2c <- as.numeric(mean_fam_2c_cs$cross_av_2c)
mean_fam_2c_cs$diff <- NA
for(i in 1:nrow(mean_fam_2c_cs)){
  mean_fam_2c_cs$diff <- mean_fam_2c_cs$self_av_2c - mean_fam_2c_cs$cross_av_2c
}

library(ggplot2)
dev.off()
ggplot(data=mean_fam_2c_cs, aes(x=family, y=diff, 
                             fill=family)) +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.text=element_text(face="italic")) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Difference in Genome Size")+
  scale_fill_hue(name="Family")

####Range Outliers####

####Range Complete Analysis####
#Analysis
shapiro.test(range$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(range$total_hectads_2010_19)) #p<0.05 therefore non-normal
leveneTest(total_hectads_2010_19~msystem, data=range) #p<0.05 so variances aren't homogenous
leveneTest(log(total_hectads_2010_19)~msystem, data=range) #p>0.05 so vairances are homogenous
kruskal.test(log(total_hectads_2010_19)~msystem, data=range) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(log(range$total_hectads_2010_19), range$msystem, p.adjust.method = "BH")
#Mean and median gen cross larger than gen self
#Mean and median mixed smaller than gen cross
#mean and median gen self and mixed not sig different
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=range)
boxplot(log(total_hectads_2010_19)~msystem, data=range, add=T, 
        col="#0000ff22")

####Range Family Specific Analysis####
range_family <- range[complete.cases(range$family),]
family_freq_2 <- as.data.frame(table(range_family$family))
large_fams_2 <- family_freq_2[which(family_freq_2$Freq>=20),]

#Apiaceae
Apiaceae_r <- range_family[which(range_family$family=="Apiaceae"),]
table(Apiaceae_r$msystem)
#Analysis
shapiro.test(Apiaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Apiaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
leveneTest(total_hectads_2010_19~msystem, data=Apiaceae_r) #p>0.05 so variances are homogenous
#leveneTest(log(total_hectads_2010_19)~msystem, data=Apiaceae_r) #p<0.05 so variances aren't homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Apiaceae_r) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Apiaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Apiaceae_r, add=T, 
        col="#0000ff22")

#Asteraceae
Asteraceae_r <- range_family[which(range_family$family=="Asteraceae"),]
table(Asteraceae_r$msystem)
#Analysis
shapiro.test(Asteraceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Asteraceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Asteraceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Asteraceae_r) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Asteraceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Asteraceae_r, add=T, 
        col="#0000ff22")

#Brassicaceae
Brassicaceae_r <- range_family[which(range_family$family=="Brassicaceae"),]
table(Brassicaceae_r$msystem)
#Analysis
shapiro.test(Brassicaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Brassicaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Brassicaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Brassicaceae_r) #p>0.05 therefore median genome size not sig. different
pairwise.wilcox.test(log(Brassicaceae_r$total_hectads_2010_19), Brassicaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed or gen. self
#Mean and median mixed not sig diff to gen. cross or gen. self
#Mean and median gen. self not sig diff to gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Asteraceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Asteraceae_r, add=T, 
        col="#0000ff22")

#Caryophyllaceae
Caryophyllaceae_r <- range_family[which(range_family$family=="Caryophyllaceae"),]
table(Caryophyllaceae_r$msystem)
#Analysis
shapiro.test(Caryophyllaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Caryophyllaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#caryo_anova <- aov(log(total_hectads_2010_19)~msystem, data=Caryophyllaceae_r)
#summary(caryo_anova) #p<0.05 so sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Caryophyllaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Caryophyllaceae_r) #p>0.05 therefore median genome size not sig. different
pairwise.wilcox.test(log(Caryophyllaceae_r$total_hectads_2010_19), Caryophyllaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed or gen. self
#Mean and median mixed not sig diff to gen. cross or gen. self
#Mean and median gen. self not sig diff to gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Caryophyllaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Caryophyllaceae_r, add=T, 
        col="#0000ff22")

#Cyperaceae
Cyperaceae_r <- range_family[which(range_family$family=="Cyperaceae"),]
table(Cyperaceae_r$msystem)

#Fabaceae
Fabaceae_r <- range_family[which(range_family$family=="Fabaceae"),]
table(Fabaceae_r$msystem)
#Analysis
shapiro.test(Fabaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Fabaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#fab_anova <- aov(log(total_hectads_2010_19)~msystem, data=Fabaceae_r)
#summary(fab_anova) #p<0.05 so sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Fabaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Fabaceae_r) #p>0.05 therefore median genome size not sig. different
pairwise.wilcox.test(log(Fabaceae_r$total_hectads_2010_19), Fabaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed or gen. self
#Mean and median mixed not sig diff to gen. cross or gen. self
#Mean and median gen. self not sig diff to gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Fabaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Fabaceae_r, add=T, 
        col="#0000ff22")

#Lamiaceae
Lamiaceae_r <- range_family[which(range_family$family=="Lamiaceae"),]
table(Lamiaceae_r$msystem)
#Analysis
shapiro.test(Lamiaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Lamiaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#lami_anova <- aov(log(total_hectads_2010_19)~msystem, data=Lamiaceae_r)
#summary(lami_anova) #p<0.05 so sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Lamiaceae_r) #p<0.05 so variances are not homogenous
leveneTest(log(total_hectads_2010_19)~msystem, data=Lamiaceae_r) #p<0.05 so variances are not homogenous
#kruskal.test(total_hectads_2010_19~msystem, data=Lamiaceae_r)
#pairwise.wilcox.test(log(Lamiaceae_r$total_hectads_2010_19), Lamiaceae_r$msystem, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Lamiaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Lamiaceae_r, add=T, 
        col="#0000ff22")

#Orchidaceae
Orchidaceae_r <- range_family[which(range_family$family=="Orchidaceae"),]
table(Orchidaceae_r$msystem)
#Analysis
shapiro.test(Orchidaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Orchidaceae_r$total_hectads_2010_19)) #p>0.05 therefore normal
orchid_anova <- aov(log(total_hectads_2010_19)~msystem, data=Orchidaceae_r)
summary(orchid_anova) #p>0.05 so no sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Orchidaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Orchidaceae_r) #p>0.05 therefore median genome size not sig. different
pairwise.wilcox.test(log(Orchidaceae_r$total_hectads_2010_19), Orchidaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed or gen. self
#Mean and median mixed not sig diff to gen. cross or gen. self
#Mean and median gen. self not sig diff to gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Orchidaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Orchidaceae_r, add=T, 
        col="#0000ff22")

#Orobanchaceae
Orobanchaceae_r <- range_family[which(range_family$family=="Orobanchaceae"),]
table(Orobanchaceae_r$msystem)
#Analysis
shapiro.test(Orobanchaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Orobanchaceae_r$total_hectads_2010_19)) #p>0.05 therefore normal
orobanch_anova <- aov(log(total_hectads_2010_19)~msystem, data=Orobanchaceae_r)
summary(orobanch_anova) #p>0.05 so no sig. diff between means
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Orobanchaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Orobanchaceae_r) #p>0.05 therefore median genome size not sig. different
pairwise.wilcox.test(log(Orobanchaceae_r$total_hectads_2010_19), Orobanchaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed or gen. self
#Mean and median mixed not sig diff to gen. cross or gen. self
#Mean and median gen. self not sig diff to gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Orobanchaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Orobanchaceae_r, add=T, 
        col="#0000ff22")

#Plantaginaceae
Plantaginaceae_r <- range_family[which(range_family$family=="Plantaginaceae"),]
table(Plantaginaceae_r$msystem)
#Analysis
shapiro.test(Plantaginaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Plantaginaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#plantagin_anova <- aov(log(total_hectads_2010_19)~msystem, data=Plantaginaceae_r)
#summary(plantagin_anova)
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Plantaginaceae_r) #p<0.05 so variances aren't homogenous
leveneTest(log(total_hectads_2010_19)~msystem, data=Plantaginaceae_r) #p>0.05 so variances are homogenous
kruskal.test(log(total_hectads_2010_19)~msystem, data=Plantaginaceae_r) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(log(Plantaginaceae_r$total_hectads_2010_19), Plantaginaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed or gen. self
#Mean and median mixed sig diff to gen. cross or gen. self
#Mean and median gen. self sig diff to gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Plantaginaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Plantaginaceae_r, add=T, 
        col="#0000ff22")

#Poaceae
Poaceae_r <- range_family[which(range_family$family=="Poaceae"),]
table(Poaceae_r$msystem)
#Analysis
shapiro.test(Poaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Poaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#poa_anova <- aov(log(total_hectads_2010_19)~msystem, data=Poaceae_r)
#summary(poa_anova)
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Poaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Poaceae_r) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(Poaceae_r$total_hectads_2010_19, Poaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed and higher than gen. self
#Mean and median mixed sig diff to gen. cross or gen. self
#Mean and median gen. self lower than gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Poaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Poaceae_r, add=T, 
        col="#0000ff22")

#Polygonaceae
Polygonaceae_r <- range_family[which(range_family$family=="Polygonaceae"),]
table(Polygonaceae_r$msystem)
#Analysis
shapiro.test(Polygonaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Polygonaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#polygon_anova <- aov(log(total_hectads_2010_19)~msystem, data=Polygonaceae_r)
#summary(polygon_anova)
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Polygonaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Polygonaceae_r) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(Polygonaceae_r$total_hectads_2010_19, Polygonaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed and higher than gen. self
#Mean and median mixed sig diff to gen. cross or gen. self
#Mean and median gen. self lower than gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Polygonaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Polygonaceae_r, add=T, 
        col="#0000ff22")

#Ranunculaceae
Ranunculaceae_r <- range_family[which(range_family$family=="Ranunculaceae"),]
table(Ranunculaceae_r$msystem)
#Analysis
shapiro.test(Ranunculaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Ranunculaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#ranancul_anova <- aov(log(total_hectads_2010_19)~msystem, data=Ranunculaceae_r)
#summary(ranancul_anova)
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Ranunculaceae_r) #p<0.05 so variances aren't homogenous
leveneTest(log(total_hectads_2010_19)~msystem, data=Ranunculaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Ranunculaceae_r) #p>0.05 therefore median genome size not sig. different
pairwise.wilcox.test(Ranunculaceae_r$total_hectads_2010_19, Ranunculaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed or gen. self
#Mean and median mixed not sig diff to gen. cross or gen. self
#Mean and median gen. self not sid diff to gen. cross
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Ranunculaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Ranunculaceae_r, add=T, 
        col="#0000ff22")

#Rosaceae
Rosaceae_r <- range_family[which(range_family$family=="Rosaceae"),]
table(Rosaceae_r$msystem)
#Analysis
shapiro.test(Rosaceae_r$total_hectads_2010_19) #p<0.05 therefore non-normal
shapiro.test(log(Rosaceae_r$total_hectads_2010_19)) #p<0.05 therefore non-normal
#rosa_anova <- aov(log(total_hectads_2010_19)~msystem, data=Rosaceae_r)
#summary(rosa_anova)
#KW as alternative if all test need to be the same?
leveneTest(total_hectads_2010_19~msystem, data=Rosaceae_r) #p<0.05 so variances aren't homogenous
leveneTest(log(total_hectads_2010_19)~msystem, data=Rosaceae_r) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_19~msystem, data=Rosaceae_r) #p<0.05 therefore median genome size sig. different
pairwise.wilcox.test(Rosaceae_r$total_hectads_2010_19, Rosaceae_r$msystem, p.adjust.method = "BH")
#Mean and median gen. cross not sig diff to mixed and higher than gen. self
#Mean and median mixed not sig diff to gen. cross and higher than gen. self
#Mean and median gen. self sig diff to gen. cross and mixed
#Beeswarm Plot
beeswarm(log(total_hectads_2010_19)~msystem, data=Rosaceae_r)
boxplot(log(total_hectads_2010_19)~msystem, data=Rosaceae_r, add=T, 
        col="#0000ff22")

####Range Subset Analysis####

####Range PGLS Analysis####

###Range Change Complete Analysis####
#Analysis
shapiro.test(range$total_hectad_change) #p<0.05 therefore non-normal
leveneTest(total_hectad_change~msystem, data=range) #p>0.05 so variances are homogenous
kruskal.test(total_hectad_change~msystem, data=range) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(total_hectad_change~msystem, data=range)
boxplot(total_hectad_change~msystem, data=range, add=T, 
        col="#0000ff22")
