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

rm(list=ls())

############################################################################################################################################
#####	Indel bar count for all samples
i = fread('/BiO/Research/UNIST-HPRT-SIG-2020-1028/part8/output/sigProfilerExtractor.ref/output/SBS/sigProfilerExtractor.ref.SBS96.all')
meta = fread('/BiO/Research/UNIST-HPRT-SIG-2020-1028/20250112.sample.list.txt')

output.dir = '/BiO/Research/UNIST-HPRT-SIG-2020-1028/part8/output/sigProfilerExtractor.ref/plot/sbs/figure4'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)
i$MutationType = substr(i$MutationType, 3, 5)
i = aggregate(.~MutationType, i, sum)

df = melt(i, id.vars = 'MutationType')
colnames(df) = c('MutationType', 'Sample', 'Count')

df = merge(df, meta, by='Sample')

df$Count = as.numeric(df$Count)

df[df$Group=='WT_Ref',]
df = df[df$Group != 'WT_Ref',]
df = df[df$Sample != 'LIG4#2_C#1',]

df$IR = gsub('\\_\\S+', '', df$Group)

df = df[df$IR != 'Fe-ion',]

group.level = c('X-ray_WT', 'X-ray_LIG4-MT', 'X-ray_POLQ-MT', 'X-ray_WT_siCtrl', 'X-ray_WT_siBRCA2',
'Carbon-ion_WT', 'Carbon-ion_LIG4-MT', 'Carbon-ion_POLQ-MT', 'Carbon-ion_WT_siCtrl', 'Carbon-ion_WT_siBRCA2')

df = df[df$Group %in% group.level,]

df$Gene = gsub('X-ray_', '', df$Group)
df$Gene = gsub('Carbon-ion_', '', df$Gene)
df$Gene = gsub('Fe-ion_', '', df$Gene)

unique(df$IR)
unique(df$Gene)
unique(df$Group)

length(unique(df$Group))
length(unique(meta$Group))
setdiff(unique(meta$Group), unique(df$Group))

#####	Only for sample>=3
#df2 = df[df$Group %in% names(table(meta$Group)[table(meta$Group)>2]),]
df2 = df

###########################################################################################
#####	Draw plot
length(group.level)
setdiff(group.level, unique(df2$Group))
setdiff(unique(df2$Group), group.level)

#group.level = group.level[group.level %in% unique(df2$Group)]

setdiff(group.level, unique(df2$Group))

unique(df2$Gene)

gene.level = c('WT', 'LIG4-MT', 'POLQ-MT', 'WT_siCtrl', 'WT_siBRCA2')

setdiff(gene.level, unique(df2$Gene))
setdiff(unique(df2$Gene), gene.level)

ir.level = c('X-ray', 'Carbon-ion')

df2$Gene = factor(df2$Gene, levels=gene.level)
df2$Group = factor(df2$Group, levels=group.level)
df2$IR = factor(df2$IR, levels=ir.level)
df2$MutationType = factor(df2$MutationType, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))

df2 = df2[order(factor(df2$Group, levels=group.level)),]
sample.level = unique(df2$Sample)
df2$Sample = factor(df2$Sample, levels=sample.level)

sum.df = df2[,c('Sample', 'Count')]
sum.df = aggregate(.~Sample, sum.df, sum)
sum.df = merge(sum.df, meta, by='Sample')

table(sum.df$Group)

test.df = compare_means(Count~Group, method='t.test', data=sum.df)
test.df2 = as.data.frame(subset(test.df, group1=='X-ray_WT'))
test.df3 = as.data.frame(subset(test.df, group2=='X-ray_WT'))
test.df3 = test.df3[,c(1,3,2,4:8)]
colnames(test.df3) = colnames(test.df2)
test.df4 = rbind(test.df2, test.df3)

row.group = group.level[group.level!='X-ray_WT']
row.group = row.group[row.group %in% unique(sum.df$Group)]
test.df4 = test.df4[match(row.group, test.df4$group2),]

sum.df2 = sum.df[,c(3,2)]
mean.df = aggregate(.~Group, sum.df2, mean)

wt.mean.df = mean.df[mean.df$Group=='X-ray_WT',]
mean.df = mean.df[mean.df$Group!='X-ray_WT',]

mean.df = mean.df[match(row.group, mean.df$Group),]

mean.df$Count = mean.df$Count/wt.mean.df[1,2]

test.df5 = test.df4
test.df5$Ratio = mean.df$Count
test.df5$p.format = as.numeric(test.df5$p.format)
test.df5$colors = ifelse(test.df5$p.format>0.1, 'white',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.05, '#e899af',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.01, '#dd6688',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.001, '#d23260',
					ifelse(test.df5$Ratio>1 & test.df5$p.format<=0.001, '#C70039',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.05, '#DDF2FD',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.01, '#9BBEC8',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.001, '#427D9D', '#164863'))))))))

#####	Down
#####	0.1, 0.05, 0.01, 0.001
#####	#DDF2FD, #9BBEC8, #427D9D, #164863

#####	UP
#####	0.1, 0.05, 0.01, 0.001
#####	#e899af, #dd6688, #d23260, #C70039

pval.strips.colors = c('white', test.df5$colors)

pval.facet.label = ifelse(test.df5$p.format<0.05, paste0(test.df5$p.signif, '\n', test.df5$group2), 
				ifelse(test.df5$p.format<0.1, paste0('#', '\n', test.df5$group2), test.df5$group2))

pval.facet.label = c('WT', pval.facet.label)

pval.facet.label = gsub('X-ray_', '', pval.facet.label)
pval.facet.label = gsub('Carbon-ion_', '', pval.facet.label)
pval.facet.label = gsub('Fe-ion_', '', pval.facet.label)
pval.facet.label = gsub('-MT', '-/-', pval.facet.label)

facet.name = c(ir.level, group.level)
temp.label = c(ir.level, pval.facet.label)
temp.label = gsub('_', '\n', temp.label)

facet.label = temp.label
names(facet.label) = facet.name

hline.mean.df = aggregate(.~Group, sum.df2, mean)
hline.mean.df = hline.mean.df[order(factor(hline.mean.df$Group, levels=group.level)),]
hline.mean.df$IR = gsub('\\_\\S+', '', hline.mean.df$Group)
hline.mean.df$Group = factor(hline.mean.df$Group, levels=group.level)
hline.mean.df$IR = factor(hline.mean.df$IR, levels=ir.level)

colors = c('#0cbdeb', '#0e0e0e', '#e12c26', '#b8b3b4', '#a0ce62', '#ecc6c5')

legend = ggplot(df2, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(data=hline.mean.df, aes(yintercept=Count), linetype='dashed', colour='red', linewidth=5) +
	scale_fill_manual(values=colors) +
	facet_nested(~IR+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=100, base_family='sans') +
	ylab('Mutation Count') + xlab(NULL) +
	theme(axis.text.x=element_blank(), axis.title=element_text(size=150, face='bold')) +
	theme(strip.text=element_text(colour='black', size=50)) +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60'), legend.key.size=unit(3, 'cm'))

p = ggplot(df2, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(data=hline.mean.df, aes(yintercept=Count), linetype='dashed', colour='red', linewidth=5) +
	scale_fill_manual(values=colors) +
	facet_nested(~IR+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip_nested(size='variable')) +
	theme_light(base_size=150, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(axis.text.x=element_blank(), axis.title=element_text(size=150, face='bold')) +
	theme(strip.text=element_text(colour='black', face='bold')) +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60')) +
	guides(x='none') + 	theme(legend.position='none')

strip.colors = c(c('#ef3c41', '#3762ae'), pval.strips.colors)

pdf(file=file.path(output.dir, 'sbs.6class.barplot.figure4.pdf'), width=120, height=40)
g = ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k = which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- strip.colors[i]
}

pdf(file=file.path(output.dir, 'sbs.6class.barplot.figure4.pdf'), width=120, height=40)
plot(g)
dev.off()

pdf(file=file.path(output.dir, 'sbs.6class.barplot.legend.figure4.pdf'), width=120, height=40)
legend
dev.off()

