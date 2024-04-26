library(vegan)
library (pairwiseAdonis)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-3000-23")

#TECHNICAL FACTORS

tech_map_uw <- read.table('CH_metadata_for_tech_uw.txt', header=T)
tech_map_w <- read.table('CH_metadata_for_tech_w.txt', header=T)

tech_uw_dm <- as.dist(read.table('uw-distance-matrix-tech.tsv', header=T))
tech_w_dm <- as.dist(read.table('w-distance-matrix-tech.tsv', header=T))

tech2_uw_dm <- as.dist(read.table('uw-distance-matrix-tech-time.txt', header=T))
tech2_w_dm <- as.dist(read.table('w-distance-matrix-tech-time.txt', header=T))
tech2_map_uw <- read.table('CH_metadata_for_tech_uw_time.txt', header=T)
tech2_map_w <- read.table('CH_metadata_for_tech_w_time.txt', header=T)

##beta for for technical

adonis2(tech_uw_dm~Individual+Location, data=tech_map_uw, permutations=5000)
adonis2(tech_w_dm~Individual+Location, data=tech_map_w, permutations=5000)

adonis2(tech_uw_dm~Individual+Substrate2, data=tech_map_uw, permutations=5000)
adonis2(tech_w_dm~Individual+Substrate2, data=tech_map_w, permutations=5000)

adonis2(tech_uw_dm~Individual+Weather, data=tech_map_uw, permutations=5000)
adonis2(tech_w_dm~Individual+Weather, data=tech_map_w, permutations=5000)

adonis2(tech_uw_dm~Individual+Temp, data=tech_map_uw, permutations=5000)
adonis2(tech_w_dm~Individual+Temp, data=tech_map_w, permutations=5000)

adonis2(tech_uw_dm~Individual+Time_to_collection, data=tech_map_uw, permutations=5000)
adonis2(tech_w_dm~Individual+Time_to_collection, data=tech_map_w, permutations=5000)

###only samples with time to collection included
adonis2(tech2_uw_dm~Individual+Time_to_collection, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Individual+Time_to_collection, data=tech2_map_w, permutations=5000)

adonis2(tech2_uw_dm~Individual+Temp*Time_to_collection, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Individual+Temp*Time_to_collection, data=tech2_map_w, permutations=5000)

##alpha for temperature and time

library(nlme)
library(car)
library(ggplot2)

alpha_tech<-read.table('alpha_summary_reduced.txt', header=T)

shan_temp<-lme(fixed=shannon_entropy~Temp, data=alpha_tech, 
         random= ~1|Individual)
Anova(shan_temp)

faith_temp<-lme(fixed=faith_pd~Temp, data=alpha_tech, 
           random= ~1|Individual)
Anova(faith_temp)

obs_temp<-lme(fixed=observed_features~Temp, data=alpha_tech, 
           random= ~1|Individual)
Anova(obs_temp)


shan_time<-lme(fixed=shannon_entropy~Time_to_collection, data=alpha_tech, 
               random= ~1|Individual)
Anova(shan_time)

faith_time<-lme(fixed=faith_pd~Time_to_collection, data=alpha_tech, 
                random= ~1|Individual)
Anova(faith_time)

obs_time<-lme(fixed=observed_features~Time_to_collection, data=alpha_tech, 
              random= ~1|Individual)
Anova(obs_time)

ggplot(alpha_tech, aes(x=Temp, y=faith_pd)) + geom_point()
ggplot(alpha_tech, aes(x=Time_to_collection, y=shannon_entropy)) + geom_point()


#CHIMP STATIC FACTORS -- GOING TO HAVE TO DO WITH AVERAGES
##beta diversity
adonis2(tech2_uw_dm~Temp+Time_to_collection+Group, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Temp+Time_to_collection+Group, data=tech2_map_w, permutations=5000)

adonis2(tech2_uw_dm~Temp+Time_to_collection+Health, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Temp+Time_to_collection+Health, data=tech2_map_w, permutations=5000)

adonis2(tech2_uw_dm~Temp+Time_to_collection+Age, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Temp+Time_to_collection+Age, data=tech2_map_w, permutations=5000)

adonis2(tech2_uw_dm~Temp+Time_to_collection+Sex, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Temp+Time_to_collection+Sex, data=tech2_map_w, permutations=5000)

##beta diversity with avg

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-avg")

avg_uw_dm<- as.dist(read.table('uw-avg-distance-matrix.tsv', header=T))
avg_w_dm<- as.dist(read.table('w-avg-distance-matrix.tsv', header=T))
avg_map<-read.table('CH_metadata_for_avg.txt', header=T)

avg_map$Group<-as.factor(avg_map$Group)

adonis2(avg_uw_dm~Age, data=avg_map, permutations=5000)
adonis2(avg_w_dm~Age, data=avg_map, permutations=5000)

adonis2(avg_uw_dm~Sex, data=avg_map, permutations=5000)
adonis2(avg_w_dm~Sex, data=avg_map, permutations=5000)

adonis2(avg_uw_dm~Group, data=avg_map, permutations=5000)
adonis2(avg_w_dm~Group, data=avg_map, permutations=5000)

adonis2(avg_uw_dm~Health, data=avg_map, permutations=5000)
adonis2(avg_w_dm~Health, data=avg_map, permutations=5000)

##beta diversity with t1
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-t1")

t1_uw_dm<- as.dist(read.table('uw-distance-matrix-t1.tsv', header=T))
t1_w_dm<- as.dist(read.table('w-distance-matrix-t1.tsv', header=T))
avg_map<-read.table('CH_metadata_for_avg.txt', header=T)

adonis2(t1_uw_dm~Age, data=avg_map, permutations=5000)
adonis2(t1_w_dm~Age, data=avg_map, permutations=5000)

adonis2(t1_uw_dm~Sex, data=avg_map, permutations=5000)
adonis2(t1_w_dm~Sex, data=avg_map, permutations=5000)

adonis2(t1_uw_dm~Group, data=avg_map, permutations=5000)
adonis2(t1_w_dm~Group, data=avg_map, permutations=5000)

adonis2(t1_uw_dm~Health, data=avg_map, permutations=5000)
adonis2(t1_w_dm~Health, data=avg_map, permutations=5000)

##beta diversity with tlast
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-tlast")

tlast_uw_dm<- as.dist(read.table('uw-distance-matrix-tlast.tsv', header=T))
tlast_w_dm<- as.dist(read.table('w-distance-matrix-tlast.tsv', header=T))
avg_map<-read.table('CH_metadata_for_avg.txt', header=T)

adonis2(tlast_uw_dm~Age, data=avg_map, permutations=5000)
adonis2(tlast_w_dm~Age, data=avg_map, permutations=5000)

adonis2(tlast_uw_dm~Sex, data=avg_map, permutations=5000)
adonis2(tlast_w_dm~Sex, data=avg_map, permutations=5000)

adonis2(tlast_uw_dm~Group, data=avg_map, permutations=5000)
adonis2(tlast_w_dm~Group, data=avg_map, permutations=5000)

adonis2(tlast_uw_dm~Health, data=avg_map, permutations=5000)
adonis2(tlast_w_dm~Health, data=avg_map, permutations=5000)


##alpha diversity with average
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-avg")

alpha_avg<-read.table('alpha_summary.txt', header=T)

shan_age<-lm(shannon_entropy~Age, data=alpha_avg)
Anova(shan_age)

faith_age<-lm(faith_pd~Age, data=alpha_avg)
Anova(faith_age)

obs_age<-lm(observed_features~Age, data=alpha_avg)
Anova(obs_age)



shan_sex<-lm(shannon_entropy~Sex, data=alpha_avg)
Anova(shan_sex)

faith_sex<-lm(faith_pd~Sex, data=alpha_avg)
Anova(faith_sex)

obs_sex<-lm(observed_features~Sex, data=alpha_avg)
Anova(obs_sex)



shan_group<-lm(shannon_entropy~Group, data=alpha_avg)
Anova(shan_group)

faith_group<-lm(faith_pd~Group, data=alpha_avg)
Anova(faith_group)

obs_group<-lm(observed_features~Group, data=alpha_avg)
Anova(obs_group)



shan_health<-lm(shannon_entropy~Health, data=alpha_avg)
Anova(shan_health)

faith_health<-lm(faith_pd~Health, data=alpha_avg)
Anova(faith_health)

obs_health<-lm(observed_features~Health, data=alpha_avg)
Anova(obs_health)

#LONGITUDINAL FACTORS
##beta diversity

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-3000-23")

Bristolc1<-as.factor(tech2_map_uw$Bristol_cat)
Bristolc2<-as.factor(tech2_map_w$Bristol_cat)

adonis2(tech2_uw_dm~Individual+Temp+Time_to_collection+Bristolc1, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Individual+Temp+Time_to_collection+Bristolc2, data=tech2_map_w, permutations=5000)

adonis2(tech2_uw_dm~Individual+Temp+Time_to_collection+Enclosure, data=tech2_map_uw, permutations=5000)
adonis2(tech2_w_dm~Individual+Temp+Time_to_collection+Enclosure, data=tech2_map_w, permutations=5000)

##alpha diversity
Bristolc3<-as.factor(alpha_tech$Bristol_cat)

shan_bris<-lme(fixed=shannon_entropy~Bristolc3, data=alpha_tech, 
               random= list(Individual=~1, Temp=~1, Time_to_collection=~1))
Anova(shan_bris)

faith_bris<-lme(fixed=faith_pd~Bristol, data=alpha_tech, 
                random= list(Individual=~1, Temp=~1, Time_to_collection=~1))
Anova(faith_bris)

obs_bris<-lme(fixed=observed_features~Bristol, data=alpha_tech, 
              random= list(Individual=~1, Temp=~1, Time_to_collection=~1))
Anova(obs_bris)


shan_encl<-lme(fixed=shannon_entropy~Enclosure, data=alpha_tech, 
               random= list(Individual=~1, Temp=~1, Time_to_collection=~1))
Anova(shan_encl)

faith_encl<-lme(fixed=faith_pd~Enclosure, data=alpha_tech, 
                random= list(Individual=~1, Temp=~1, Time_to_collection=~1))
Anova(faith_encl)

obs_encl<-lme(fixed=observed_features~Enclosure, data=alpha_tech, 
              random= list(Individual=~1, Temp=~1, Time_to_collection=~1))
Anova(obs_encl)

ggplot(alpha_tech, aes(x=Bristol, y=shannon_entropy)) + geom_point()


#SOCIAL GROUP MERGE
##beta for merged group merge in 2020

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-group-merge")

uw_dm_merge<-as.dist(read.table('uw-distance-matrix-merge_red.txt', header=T))
w_dm_merge<-as.dist(read.table('w-distance-matrix-merge-red.txt', header=T))
map_merge<-read.table('CH_metadata_group_merge_red.txt', header=T)


adonis2(uw_dm_merge~Individual+Temp+Time_to_collection+Group2*SampleNumber, data=map_merge, permutations=5000)
adonis2(uw_dm_merge~Individual+Group2*SampleNumber, data=map_merge, permutations=5000)
adonis2(uw_dm_merge~Group2*SampleNumber, data=map_merge, permutations=5000)
adonis2(uw_dm_merge~Individual+Temp+Time_to_collection+Group2*Month, data=map_merge, permutations=5000)
adonis2(uw_dm_merge~Individual+Group2*Month, data=map_merge, permutations=5000)
adonis2(uw_dm_merge~Group2*Month, data=map_merge, permutations=5000)

adonis2(w_dm_merge~Individual+Temp+Time_to_collection+Group2*SampleNumber, data=map_merge, permutations=5000)
adonis2(w_dm_merge~Individual+Group2*SampleNumber, data=map_merge, permutations=5000)
adonis2(w_dm_merge~Group2*SampleNumber, data=map_merge, permutations=5000)
adonis2(w_dm_merge~Individual+Temp+Time_to_collection+Group2*Month, data=map_merge, permutations=5000)
adonis2(w_dm_merge~Individual+Group2*Month, data=map_merge, permutations=5000)
adonis2(w_dm_merge~Group2*Month, data=map_merge, permutations=5000)

#CORRELATION LOOP FOR TIME, TEMP, BRISTOL
#time
library(fdrtool)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis")

cor_data <- read.csv("feature_table_for_corloop.csv", header=T)
pvaluematrix = mat.or.vec(12230,2)
fvaluematrix = mat.or.vec(12230,2)
estmatrix = mat.or.vec(12230,2)
 c = for(i in 5:12234){b=cor.test(cor_data[,2], cor_data[,i], data=cor_data)
 p = b$'p.value'
 f = b$'statistic'
 e = b$'estimate'
pvaluematrix[i-4,] = p
fvaluematrix[i-4,] = f
estmatrix[i-4,]=e}
write.table(pvaluematrix, "time_cor_p.txt")
write.table(fvaluematrix, "time_cor_f.txt")
write.table(estmatrix, "time_cor_est.txt")

p_no_na = na.omit(pvaluematrix)
time_cor_corrected = fdrtool(p_no_na[,1], statistic = "pvalue", plot = FALSE)
write.table(time_cor_corrected, "time_cor_q.txt")

#temp
cor_data <- read.csv("feature_table_for_corloop.csv", header=T)
pvaluematrix = mat.or.vec(12230,2)
fvaluematrix = mat.or.vec(12230,2)
estmatrix = mat.or.vec(12230,2)
c = for(i in 5:12234){b=cor.test(cor_data[,3], cor_data[,i], data=cor_data)
p = b$'p.value'
f = b$'statistic'
e = b$'estimate'
pvaluematrix[i-4,] = p
fvaluematrix[i-4,] = f
estmatrix[i-4,]=e}
write.table(pvaluematrix, "temp_cor_p.txt")
write.table(fvaluematrix, "temp_cor_f.txt")
write.table(estmatrix, "temp_cor_est.txt")

p_no_na = na.omit(pvaluematrix)
temp_cor_corrected = fdrtool(p_no_na[,1], statistic = "pvalue", plot = FALSE)
write.table(temp_cor_corrected, "temp_cor_q.txt")

#bristol
cor_data <- read.csv("feature_table_for_corloop.csv", header=T)
pvaluematrix = mat.or.vec(12230,2)
fvaluematrix = mat.or.vec(12230,2)
estmatrix = mat.or.vec(12230,2)
c = for(i in 5:12234){b=cor.test(cor_data[,4], cor_data[,i], data=cor_data)
p = b$'p.value'
f = b$'statistic'
e = b$'estimate'
pvaluematrix[i-4,] = p
fvaluematrix[i-4,] = f
estmatrix[i-4,]=e}
write.table(pvaluematrix, "bris_cor_p.txt")
write.table(fvaluematrix, "bris_cor_f.txt")
write.table(estmatrix, "bris_cor_est.txt")

p_no_na = na.omit(pvaluematrix)
bris_cor_corrected = fdrtool(p_no_na[,1], statistic = "pvalue", plot = FALSE)
write.table(bris_cor_corrected, "bris_cor_q.txt")


#ANCOM FOR HEALTH AND SEX AND ENCLOSURE
#health
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis")
asv = read.table("average_abc.txt", header=T, check.names=FALSE)
metadata = read.table("CH_metadata_for_avg.txt", header=T)
taxonomy = read.table("taxonomy_abc.txt", header=T)

asv_matrix = asv %>% column_to_rownames("sampleid") %>% as.matrix()
tax_matrix = taxonomy %>% column_to_rownames("feature") %>% as.matrix()
meta = metadata %>% column_to_rownames("sampleid")

ASV<-otu_table(asv_matrix, taxa_are_rows = TRUE)
TAX<-tax_table(tax_matrix)
samples<-sample_data(meta)

asv_phylo = phyloseq(ASV, TAX, samples)


health_asv = ancombc2(data=asv_phylo, fix_formula="Health",
                     p_adj_method = "fdr",
                     group = "Health", global = T)

res_h_asv<-health_asv$res
res_global_h_asv<-health_asv$res_global

write_csv(res_h_asv, "Diff_abund_health_asv.csv")
write_csv(res_global_h_asv, "Diff_abund_trash_health_global.csv")

#sex
sex_asv = ancombc2(data=asv_phylo, fix_formula="Sex",
                      p_adj_method = "fdr",
                      group = "Sex", global = T)

res_s_asv<-sex_asv$res
res_global_s_asv<-sex_asv$res_global

write_csv(res_s_asv, "Diff_abund_sex_asv.csv")
write_csv(res_global_s_asv, "Diff_abund_sex_global.csv")

#enclosure
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis")
asv2 = read.csv("overtime_abc.csv", header=T, check.names=FALSE)
metadata2 = read.table("CH_metadata_for_overtime.txt", header=T)
taxonomy = read.table("taxonomy_abc.txt", header=T)

asv_matrix2 = asv2 %>% column_to_rownames("sampleid") %>% as.matrix()
tax_matrix2 = taxonomy %>% column_to_rownames("feature") %>% as.matrix()
meta2 = metadata2 %>% column_to_rownames("sampleid")

ASV2<-otu_table(asv_matrix2, taxa_are_rows = TRUE)
TAX2<-tax_table(tax_matrix2)
samples2<-sample_data(meta2)

asv_phylo2 = phyloseq(ASV2, TAX2, samples2)
enc_asv = ancombc2(data=asv_phylo2, fix_formula="Individual+Enclosure",
                   p_adj_method = "fdr",
                   group = "Enclosure", global = T)

res_enc_asv<-enc_asv$res
res_global_enc_asv<-enc_asv$res_global

write_csv(res_enc_asv, "Diff_abund_enc_asv.csv")
write_csv(res_global_enc_asv, "Diff_abund_enc_global.csv")

#Bristol cat
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis")
asv2 = read.csv("overtime_abc.csv", header=T, check.names=FALSE)
metadata2 = read.table("CH_metadata_for_overtime.txt", header=T)
taxonomy = read.table("taxonomy_abc.txt", header=T)

asv_matrix2 = asv2 %>% column_to_rownames("sampleid") %>% as.matrix()
tax_matrix2 = taxonomy %>% column_to_rownames("feature") %>% as.matrix()
meta2 = metadata2 %>% column_to_rownames("sampleid")

ASV2<-otu_table(asv_matrix2, taxa_are_rows = TRUE)
TAX2<-tax_table(tax_matrix2)
samples2<-sample_data(meta2)

asv_phylo2 = phyloseq(ASV2, TAX2, samples2)
bris_asv = ancombc2(data=asv_phylo2, fix_formula="Individual+Bristol_cat",
                   p_adj_method = "fdr",
                   group = "Bristol_cat", global = T)

res_bris_asv<-bris_asv$res
res_global_bris_asv<-bris_asv$res_global

write_csv(res_bris_asv, "Diff_abund_bris_asv.csv")
write_csv(res_global_bris_asv, "Diff_abund_bris_global.csv")

#Figures
library(ggplot2)

##beta diversity location
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-3000-23")

tech2_uw_dm <- as.dist(read.table('uw-distance-matrix-tech-time.txt', header=T))
tech2_w_dm <- as.dist(read.table('w-distance-matrix-tech-time.txt', header=T))
tech2_map_uw <- read.table('CH_metadata_for_tech_uw_time.txt', header=T)
tech2_map_w <- read.table('CH_metadata_for_tech_w_time.txt', header=T)

u_mds<-metaMDS(tech2_uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = tech2_map_uw, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Location), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds

w_mds<-metaMDS(tech2_w_dm, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = tech2_map_w, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=location), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

##beta diversity age, sex, health, social group
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-avg")

avg_uw_dm<- as.dist(read.table('uw-avg-distance-matrix.tsv', header=T))
avg_w_dm<- as.dist(read.table('w-avg-distance-matrix.tsv', header=T))
avg_map<-read.table('CH_metadata_for_avg.txt', header=T)

avg_map$Group<-as.factor(avg_map$Group)

u_mds<-metaMDS(avg_uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Group), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_ellipse(aes(group=Group, color=Group), type='t')
unmds

w_mds<-metaMDS(avg_w_dm, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Group), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_ellipse(aes(group=Group, color=Group), type='t')
wnmds

u_mds<-metaMDS(avg_uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Sex), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_ellipse(aes(group=Sex, color=Sex), type='t')
unmds

w_mds<-metaMDS(avg_w_dm, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Sex), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_ellipse(aes(group=Sex, color=Sex), type='t')
wnmds

u_mds<-metaMDS(avg_uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Age), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds

w_mds<-metaMDS(avg_w_dm, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Age), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

u_mds<-metaMDS(avg_uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Health), size=5) +
  scale_color_manual(values = c("#4daf4a", "#377eb8","#984ea3")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds

w_mds<-metaMDS(avg_w_dm, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = avg_map, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Health), size=5)+
  scale_color_manual(values = c("#4daf4a", "#377eb8","#984ea3")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

##beta enclosure and bristol
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-3000-23")

Bristolc1<-as.factor(tech2_map_uw$Bristol_cat)
Bristolc2<-as.factor(tech2_map_w$Bristol_cat)

u_mds<-metaMDS(tech2_uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = tech2_map_uw, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Enclosure), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds

w_mds<-metaMDS(tech2_w_dm, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = tech2_map_w, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Enclosure), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

u_mds<-metaMDS(tech2_uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = tech2_map_uw, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Bristolc2), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds

w_mds<-metaMDS(tech2_w_dm, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = tech2_map_w, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Bristolc2), size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds


##alpha diversity time to collection
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-3000-23")
alpha_tech<-read.table('alpha_summary_reduced.txt', header=T)

ggplot(alpha_tech, aes(x=Time_to_collection, y=shannon_entropy)) + geom_point() +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))

##alpha diversity bristol
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Chimp_Haven/final_analysis/core-metrics-results-3000-23")
alpha_tech<-read.table('alpha_summary_reduced.txt', header=T)

Bristolc<-as.factor(alpha_tech$Bristol_cat)

ggplot(alpha_tech, aes(x=Bristolc, y=shannon_entropy)) + geom_boxplot() +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
