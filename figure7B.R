library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
options(bitmapType='cairo')
options(scipen=10000)

##########################################################################################################################################
#####	ID Signature Count and Proportion Bar Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-HPRT-SIG-2020-1028/part10/output/sigProfilerExtractor.id3/output/signature/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/COSMIC_ID83_Activities_refit.txt')

output.dir = '/BiO/Research/UNIST-HPRT-SIG-2020-1028/part10/output/sigProfilerExtractor.id3/plot'
if(!dir.exists(output.dir)){dir.create(output.dir)}

output.dir = '/BiO/Research/UNIST-HPRT-SIG-2020-1028/part10/output/sigProfilerExtractor.id3/plot/indel.sig'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)
meta = fread('/BiO/Research/UNIST-HPRT-SIG-2020-1028/20251025.ju.sample.list.txt')

#no.samples = setdiff(meta$Sample, i$Samples)
#no.df = data.frame(no.samples, 0, 0)
#colnames(no.df) = colnames(i)
#i = rbind(i, no.df)

df = melt(i, id.vars='Samples')
colnames(df) = c('Sample', 'Signature', 'Count')

df = merge(df, meta, by='Sample')

prop.df = cbind(Samples=i[,1], i[,-1]/rowSums(i[,-1]))
rowSums(prop.df[,-1])

prop.df = melt(prop.df, id.vars='Samples')
colnames(prop.df) = c('Sample', 'Signature', 'Proportion')

df2 = merge(df, prop.df, by=c('Sample', 'Signature'))
df2$Proportion[is.na(df2$Proportion)] = 0

df2$Group = factor(df2$Group, levels=c('Human_0Gy', 'Human_50Gy', 'Irradiated_human_cancer'))
df2$Sample = factor(df2$Sample, levels=meta$Sample)

df2 = df2[df2$Group=='Irradiated_human_cancer',]

length(unique(df2$Signature))
#[1] 6

pdf(file=file.path(output.dir, 'sig.indel.bar.count.canonical.v4.pdf'), width=80, height=100)
ggplot(df2, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=Signature), position='stack', stat='identity') +
	scale_fill_manual(values = c("ID1" = "#d02127", "ID5" = "#f26622", "ID6" = "#efeb40", 
							"ID8" = "#97c53e", "ID9" = "#516db5", "ID12" = "#7e4e9f",
							"ID2" = "#e91e63", "ID3" = "#1b2a6b")) +
	facet_grid(~Group, scales='free_x', space='free_x') +
	theme_light(base_size=300, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', size=200, face='bold')) +
	theme(legend.key.size=unit(3, 'cm'),
		axis.text.x=element_text(colour='black', size=140, angle=90, vjust=0.5, hjust=1),
		strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'))
dev.off()

pdf(file=file.path(output.dir, 'sig.indel.bar.proportion.canonical.v4.pdf'), width=80, height=100)
ggplot(df2, aes(x=Sample, y=Proportion)) +
	geom_bar(aes(fill=Signature), position='stack', stat='identity') +
	scale_fill_manual(values = c("ID1" = "#d02127", "ID5" = "#f26622", "ID6" = "#efeb40", 
							"ID8" = "#97c53e", "ID9" = "#516db5", "ID12" = "#7e4e9f",
							"ID2" = "#e91e63", "ID3" = "#1b2a6b")) +
	facet_grid(~Group, scales='free_x', space='free_x') +
	theme_light(base_size=300, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', size=200, face='bold')) +
	theme(legend.key.size=unit(3, 'cm'),
		axis.text.x=element_text(colour='black', size=140, angle=90, vjust=0.5, hjust=1),
		strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'))
dev.off()

