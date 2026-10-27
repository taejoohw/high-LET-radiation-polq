library(data.table)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(janitor)
library(ggbreak)
library(scales)
library(ggpubr)
library(bayestestR)
library(patchwork)
library(ggh4x)

options(bitmapType='cairo')
options(scipen=10000)

##########################################################################################################################################
#####	ID Signature Ratio Bar Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-HPRT-SIG-2020-1028/part8/output/sigProfilerExtractor.ref.only.sv.id2/output/signature/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/COSMIC_ID83_Activities_refit.txt')

output.dir = '/BiO/Research/UNIST-HPRT-SIG-2020-1028/part8/output/sigProfilerExtractor.ref.only.sv.id2/plot'
if(!dir.exists(output.dir)){dir.create(output.dir)}

output.dir = '/BiO/Research/UNIST-HPRT-SIG-2020-1028/part8/output/sigProfilerExtractor.ref.only.sv.id2/plot/indel.sig'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)
df = melt(i, id.vars='Samples')
colnames(df) = c('Sample', 'Signature', 'Count')

df$Signature = gsub('ID83A', 'ID8', df$Signature)

meta = fread('/BiO/Research/UNIST-HPRT-SIG-2020-1028/20250112.sample.list.txt')
df = merge(df, meta, by='Sample')

prop.df = cbind(Samples=i[,1], i[,-1]/rowSums(i[,-1]))
rowSums(prop.df[,-1])

prop.df = melt(prop.df, id.vars='Samples')
colnames(prop.df) = c('Sample', 'Signature', 'Proportion')

prop.df$Signature = gsub('ID83A', 'ID8', prop.df$Signature)

df2 = merge(df, prop.df, by=c('Sample', 'Signature'))

length(unique(df2$Signature))
#[1] 5

group.level = c('X-ray_WT', 'X-ray_LIG4-MT', 'X-ray_POLQ-MT', 'Carbon-ion_WT', 'Carbon-ion_LIG4-MT', 'Carbon-ion_POLQ-MT')
df3 = df2[df2$Group %in% group.level,]
df3$Group = factor(df3$Group, levels=group.level)

df3$IR = gsub('\\_\\S+', '', df3$Group)
df3$Gene = gsub('X-ray_', '', df3$Group)
df3$Gene = gsub('Carbon-ion_', '', df3$Gene)
df3$Gene = gsub('Fe-ion_', '', df3$Gene)

gene.level = c('WT', 'LIG4-MT', 'POLQ-MT')
ir.level = c('X-ray', 'Carbon-ion', 'Fe-ion')

df3$Gene = factor(df3$Gene, levels=gene.level)
df3$Group = factor(df3$Group, levels=group.level)
df3$IR = factor(df3$IR, levels=ir.level)

df3 = df3[order(factor(df3$Group, levels=group.level)),]
sample.level = unique(df3$Sample)
df3$Sample = factor(df3$Sample, levels=sample.level)

facet.label = gsub('X-ray_', '', group.level)
facet.label = gsub('Carbon-ion_', '', facet.label)
facet.label = gsub('Fe-ion_', '', facet.label)
facet.label = gsub('-MT', '-/-', facet.label)

facet.name = c(ir.level, group.level)
temp.label = c(ir.level, facet.label)
temp.label = gsub('_', '\n', temp.label)

facet.label = temp.label
names(facet.label) = facet.name

colors = c('#efeb40', '#97c53e')

legend = ggplot(df3, aes(x=Sample, y=Proportion)) +
	geom_bar(aes(fill=Signature), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_nested(~IR+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=100, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(axis.text.x=element_blank(), axis.title=element_text(size=150, face='bold')) +
	theme(strip.text=element_text(colour='black', size=50)) +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60'), legend.key.size=unit(3, 'cm')) +
	theme(legend.position='bottom', legend.box='horizontal') + guides(fill=guide_legend(nrow=1))

p = ggplot(df3, aes(x=Sample, y=Proportion)) +
	geom_bar(aes(fill=Signature), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_nested(~IR+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip_nested(size='variable')) +
	theme_light(base_size=150, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(axis.text.x=element_blank(), axis.title=element_text(size=150, face='bold')) +
	theme(strip.text=element_text(colour='black', face='bold')) +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60')) +
	guides(x='none') + theme(legend.position='none')

strip.colors = c('#ef3c41', '#3762ae', '#34ae5b')

pdf(file=file.path(output.dir, 'sig.indel.sv.ratio.strip.wt.lig4.polq.bar.plot.pdf'), width=200, height=60)
g = ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k = which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- strip.colors[i]
}

pdf(file=file.path(output.dir, 'sig.indel.sv.ratio.strip.wt.lig4.polq.bar.plot.pdf'), width=200, height=60)
plot(g)
dev.off()

pdf(file=file.path(output.dir, 'sig.indel.sv.ratio.strip.wt.lig4.polq.bar.plot.legend.pdf'), width=200, height=60)
legend
dev.off()


