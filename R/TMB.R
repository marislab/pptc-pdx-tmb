###### Calculate TMB and create plots for Figures 3B and S3B########
#
#     Authors: Jo Lynne Rokita
#     Updated 2019-09-16
################################################################

# working directory (created with git clone)
mainDir <- "~/pptc-pdx-tmb/"
dataDir <- paste0(mainDir,"data/")
# set path to your git cloned repo
script.folder <- paste0(mainDir, "R/") 
##create directory for output files
ifelse(!dir.exists(file.path(paste0(mainDir, "results/"))), dir.create(file.path(paste0(mainDir, "results/"))), 
       "Directory exists!")

results.folder <- paste0(mainDir, "results/") 

##set wd
setwd(mainDir)

####Install package Dependencies if needed
source(paste0(script.folder, "install.packages.R"))

####Source plot theme script
source(paste0(script.folder, "theme.R"))

####Load histology color function
source(paste0(script.folder, "demog-color-function.R"))

####load mutations calculated from maftools
load <- read.delim(paste0(dataDir,"mutations-per-model.txt"), sep = "\t", header = T, as.is = T)

####Below are the calculations from https://github.com/marislab/create-pptc-pdx-oncoprints/blob/master/R/load-maf-maftools.R

####for muts per mb, only use substitutions 
#load$subs <- load$Missense_Mutation + load$Nonsense_Mutation

####used Roche Nimblegen VCRome v. 2.1 = 45.1 Mb capture
#load$MutperMB <- round(load$subs/45.1,2)


##add clinical data
clin = read.delim(paste0(dataDir, "pptc-pdx-clinical-web.txt"), sep = "\t", header = T, as.is = T)

#merge load and clinical
load.short <- merge(load, clin[c("Model", "Histology", "Histology.Detailed","Phase")])

#calculate N per histology 
load.freq <- ddply(load.short, .(Histology.Detailed), mutate, freq = length(Histology.Detailed))

#create facet labels
load.freq$hist_label <- paste(load.freq$Histology.Detailed, " (n = ", load.freq$freq, ")", sep = "")

#reorder full dataset by median mutational load
reorder.load <- dplyr::mutate(load.freq,
                       hist_label = reorder(hist_label, MutperMB, median))

#add the median per histology
reorder.load.med <- reorder.load %>% group_by(hist_label) %>%
  dplyr::mutate(med = median(MutperMB))


#create summary stats for supplemental table - MutperMB and Ages
summ.mutpermb <- describeBy(x=reorder.load$MutperMB, group=reorder.load$hist_label, mat = T, na.rm = T)
summ.mutpermb$vars <- NULL
write.table(summ.mutpermb[,2:ncol(summ.mutpermb)], paste0(results.folder, Sys.Date(), "-summary-mutperMB.txt"),
            sep = "\t", col.names = T, row.names = F, quote = F)
ages <- describeBy(x=clin$Age, group=clin$Histology.Detailed, mat = T, na.rm = T)
ages$vars <- NULL
#remove hist in which there is no age data
ages2 <- subset(ages, n >0)
write.table(ages2[,2:ncol(ages2)], paste0(results.folder, Sys.Date(),"-summary-ages.txt"),
            sep = "\t", col.names = T, row.names = F, quote = F)


#plot Mut per Mb by histology
breaks = c(1,5,10,50,100, 150)
pdf(paste0(results.folder, Sys.Date(), "-pdx-mutload-byhistology.pdf"), width = 22, height = 8)
ggplot(data=reorder.load.med, aes(x = reorder(Model, MutperMB), y = MutperMB, color = Histology.Detailed)) +
  geom_point(stat = "identity", size = 4, shape = 1)+
  scale_color_manual(values =histcol) +
  facet_wrap("hist_label", nrow = 1, scales="free_x") +
  geom_hline(aes(yintercept = med, group = hist_label), colour = 'black')+
  theme_Publication()+
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle=90, face = "bold", hjust = 0, vjust = 0),
        panel.border = element_rect(fill = NA, colour = "grey"),
        strip.background = element_blank(),
        legend.position = "none")+
  scale_x_discrete(breaks=NULL, expand = c(0.3, 0.3))+
  scale_y_continuous(limits=c(0,10))+
  labs(x ="", y = "Mutations per Mb")

dev.off()



###group by phase; first lump PD-PM and PD phases together as relapse (NBL patients which progress or had relapsed disease and died)
load.short$Phase2 <- ifelse((load.short$Phase == "PD-PM" | load.short$Phase == "Progressing Disease" | load.short$Phase == "Relapse"), 
                       "Relapse", load.short$Phase)
###add frequencies
load.short.freq <- ddply(load.short, .(Phase2), dplyr::mutate, freq = length(Phase2))

###select only diagnosis and relapse
load.short.freq <- subset(load.short.freq, Phase2 == "Diagnosis" | Phase2 == "Relapse")

###reorder by median
load.med <- dplyr::mutate(load.short.freq,
                    phase_label = reorder(Phase2, MutperMB, median))
load.med2 <- load.med %>% group_by(phase_label) %>%
  dplyr::mutate(med = median(MutperMB))

###subset dx and relapse only
dx.relapse <- subset(load.med2, phase_label == "Diagnosis" | phase_label == "Relapse")

###add log2 for plotting
dx.relapse$logmut <- log(dx.relapse$MutperMB,2)

dx <- subset(dx.relapse, Phase2 == "Diagnosis")
rel <- subset(dx.relapse, Phase2 == "Relapse")
median(dx$MutperMB) 
median(rel$MutperMB)

###produce figure 3B, first panel
pdf(paste0(results.folder, Sys.Date(), "-mut-load-dx-relapse-allhist.pdf"), width = 6, height = 6)
dx.relapse$facet <- "All Histologies"
print(ggviolin(dx.relapse, x = "phase_label", y = "logmut", fill = "phase_label",
               palette = dxrelcol, alpha = 0.8,#add = "boxplot",
              add.params = list(fill = "white"))+
        stat_compare_means(comparisons = list(c("Diagnosis", "Relapse")), method = "wilcox.test", label = "p.format")+ # Add significance levels
        #stat_compare_means(label.y = 5) + ##global p
        facet_grid(~facet)+
        theme_Publication() + 
        xlab("") + ylab('log2[Mutations per Mb]')+
        scale_y_continuous(limits=c(0,8)))

print(ggplot(dx.relapse, aes(x = phase_label, y = logmut, colour= phase_label))+
        #geom_violin(alpha = 0.1, colour = "black")+
        geom_beeswarm(cex=1.5, alpha = 0.8, shape=21, size =4, aes(fill=phase_label), colour = "black")+
        #geom_beeswarm(cex=1.5, alpha = 0.8, size =3)+
        geom_violin(alpha = 0.1, colour = "black")+
        scale_fill_manual(values = dxrelcol, name = "Phase of Therapy")+
        facet_grid(~facet)+
        theme_Publication() + 
        xlab("") + ylab('log2[Mutations per Mb]')+
        scale_y_continuous(limits=c(0,8)))

        
dev.off()

##compare within histology 
##remove histologies which do not have both phases
df <- as.data.frame(table(dx.relapse$Histology.Detailed, dx.relapse$phase_label))
df <- subset(df, Freq > 0)

new.df <- df %>% 
  group_by(Var1) %>% 
  summarise_all(funs(trimws(paste(., collapse = '-'))))

both.phases <- subset(new.df, Var2 == "Diagnosis-Relapse")

new <- subset(dx.relapse, Histology.Detailed %in% both.phases$Var1)
table(new$Histology.Detailed, new$phase_label)



##remove if N == 3 in any phase and must have two phases
greaterthan3 <- subset(as.data.frame(table(new$Histology.Detailed, new$phase_label)), Freq>3)
#phase.table <- as.data.frame(table(rel.sub$Histology.Detailed, rel.sub$phase_label))

unite.table <- as.data.frame(greaterthan3 %>% 
                               group_by(Var1) %>% 
                               summarise_all(funs(trimws(paste(., collapse = ', ')))))

unite.table<- subset(unite.table, Var2 == "Diagnosis, Relapse")
rel.sub <- subset(new, Histology.Detailed %in% unite.table$Var1)


rel.coll <- rel.sub %>% 
  group_by(Histology.Detailed, phase_label) %>% dplyr::mutate(n = n()) %>% 
  mutate(label = paste0('\nN = ', n))
table(rel.coll$Histology.Detailed, rel.coll$phase_label)

###produce figure 3B, histology-specific panels
pdf(paste0(results.folder, Sys.Date(), "-mut-load-dx-relapse-byhist.pdf"), width = 12, height = 6)
print(ggviolin(rel.coll, x = "phase_label", y = "logmut", fill = "phase_label",
               palette = dxrelcol, alpha = 0.8,#add = "boxplot",
               add.params = list(fill = "white"))+
        facet_wrap(~Histology.Detailed, nrow = 1) +
        stat_compare_means(comparisons = list(c("Diagnosis", "Relapse")), method = "wilcox.test", label.y = c(7,4,4,6), label = "p.format")+ # Add significance levels        theme_Publication() + 
        stat_n_text() +
        xlab("") + ylab('log2[Mutations per Mb]')+
        scale_y_continuous(limits=c(0,8)))


print(ggplot(rel.coll, aes(x = phase_label, y = logmut, colour= phase_label))+
        geom_violin(alpha = 1, colour = "black")+
        geom_beeswarm(cex=2, alpha = 0.8, shape=21, size =4, aes(fill=phase_label), colour = "black")+
        scale_fill_manual(values = dxrelcol, name = "Phase of Therapy")+
        facet_wrap(~Histology.Detailed, nrow = 1) +
        theme_Publication() + 
        xlab("") + ylab('log2[Mutations per Mb]')+
        scale_y_continuous(limits=c(0,8)))

dev.off()


###mutations from dx/relapse paired models
muts <- read.delim(paste0(dataDir, "total_Mutations-All-Dx-Relapse-Models.csv"), sep = ",")
muts$log2muts <- log(muts$Total.Mutations, 2)
muts$facet <- "All Histologies"
dx.muts <- subset(muts, Phase_1 == "Diagnosis")
rel.muts <- subset(muts, Phase_1 == "Relapse")
median(dx.muts$Total.Mutations)
median(rel.muts$Total.Mutations)

###Figure 3B last panel
pdf(paste0(results.folder, Sys.Date(), "-mut-load-dx-relapse-matchedpairs.pdf"), width = 9, height = 6)
print(ggviolin(muts, x = "Phase_1", y = "log2muts", fill = "Phase_1",
                    palette = dxrelcol, alpha = 0.8,#add = "boxplot",
                    add.params = list(fill = "white"))+
             stat_compare_means(comparisons = list(c("Diagnosis", "Relapse")), method = "wilcox.test", label.y = c(20), label = "p.format")+ # Add significance levels
             facet_grid(~facet) +##global p
             theme_Publication() + 
             xlab("") + ylab('log2[Total Mutations]'))
dev.off()


