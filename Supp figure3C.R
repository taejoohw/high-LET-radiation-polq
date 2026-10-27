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
#####	SV bar count for all samples
vcfs = Sys.glob('/BiO/Research/UNIST-HPRT-SIG-2020-1028/part8/output/delly/filter/onlypass/*.pass.somatic.vcf')
vcfs = c(Sys.glob('/BiO/Research/UNIST-HPRT-SIG-2020-1028/part7/output/delly/filter/onlypass/*.pass.somatic.vcf'), vcfs)

meta = fread('/BiO/Research/UNIST-HPRT-SIG-2020-1028/20250112.sample.list.txt')

output.dir = '/BiO/Research/UNIST-HPRT-SIG-2020-1028/part8/output/delly/filter/onlypass/plot/figure4'
if(!dir.exists(output.dir)){dir.create(output.dir)}

no.samples = c()
samples = c()

df = data.frame()
for (vcf in vcfs) {
	data.dir = paste0('zgrep -v "^##" ', vcf)
	i = fread(data.dir, header=TRUE)
	colnames(i)[1] = 'CHROM'
	samples = append(colnames(i)[10], samples)

	if (nrow(i) == 0) {
		no.samples = append(colnames(i)[10], no.samples)
		next;
	}

	i2 = i[,c(1,2,3,4,5,8)]
	i2$Sample = colnames(i)[10]
	df = rbind(df, i2)
}

print(no.samples)
length(no.samples)
#[1] 163

df2 = df

df2$END = gsub('\\S+\\;END=', '', df2$INFO)
df2$END = gsub('\\;\\S+', '', df2$END)

df2$LEN = as.numeric(df2$END) - as.numeric(df2$POS)

df2$INSLEN = gsub('\\S+\\;INSLEN=', '', df2$INFO)
df2$INSLEN = gsub('\\;\\S+', '', df2$INSLEN)

df2$INFO = NULL

df2$ALT = ifelse(grepl('[', df2$ALT, fixed=TRUE), '<TRA>', df2$ALT)
df2$ALT = ifelse(grepl(']', df2$ALT, fixed=TRUE), '<TRA>', df2$ALT)
df2$ALT = gsub('<', '', df2$ALT)
df2$ALT = gsub('>', '', df2$ALT)
df2$ALT = ifelse(df2$ALT=='INV', 'Inversions', 
			ifelse(df2$ALT=='DUP', 'Duplications',
			ifelse(df2$ALT=='TRA', 'Translocations',
			ifelse(df2$ALT=='DEL', 'Deletions',
			ifelse(df2$ALT=='INS', 'Insertions', 0)))))

df2$LEN2 = ifelse(df2$ALT=='Insertions', df2$INSLEN, df2$LEN)
df2$LEN2 = as.numeric(df2$LEN2)

df2[df2$ALT=='Insertions',]

temp.df = df2[df2$LEN2>200 & df2$ALT=='Insertions',]
temp.df2 = df2[df2$LEN2>200 & df2$ALT=='Deletions',]

'%notin%' = Negate('%in%')
temp.df3 = df2[df2$ALT %notin% c('Insertions', 'Deletions'),]

df3 = rbind(temp.df, temp.df2)
df3 = rbind(df3, temp.df3)

table(df3$ALT)
#Deletions Inversions Duplications Translocations
#607       101        66                  47

df4 = as.data.frame(table(df3$Sample, df3$ALT))
colnames(df4) = c('Sample', 'MutationType', 'Count')

df4$Sample = as.character(df4$Sample)
df4$MutationType = as.character(df4$MutationType)

##### no over 200 insertions
#all.types.order = c('Insertions', 'Deletions', 'Inversions', 'Duplications', 'Translocations')
all.types.order = c('Deletions', 'Inversions', 'Duplications', 'Translocations')

df5 = df4
df5$Sample = as.character(df5$Sample)
no.samples = setdiff(samples, unique(df5$Sample))

sv.types = all.types.order

if(length(no.samples)!=0) {
	for (i in 1:length(no.samples)) {
		for (j in 1:length(sv.types)) {
			df5[nrow(df5)+1,] = c(no.samples[i], sv.types[j], 0)
		}
	}
}

df6 = merge(df5, meta, by='Sample')
df6$Count = as.numeric(df6$Count)

df6[df6$Group=='WT_Ref',]
df6 = df6[df6$Group != 'WT_Ref',]
df6 = df6[df6$Sample != 'LIG4#2_C#1',]

df6$IR = gsub('\\_\\S+', '', df6$Group)
df6 = df6[df6$IR != 'Fe-ion',]

group.level = c('X-ray_WT', 'X-ray_LIG4-MT', 'X-ray_POLQ-MT', 'X-ray_WT_siCtrl', 'X-ray_WT_siBRCA2',
'Carbon-ion_WT', 'Carbon-ion_LIG4-MT', 'Carbon-ion_POLQ-MT', 'Carbon-ion_WT_siCtrl', 'Carbon-ion_WT_siBRCA2')

df6 = df6[df6$Group %in% group.level,]

df6$Gene = gsub('X-ray_', '', df6$Group)
df6$Gene = gsub('Carbon-ion_', '', df6$Gene)
df6$Gene = gsub('Fe-ion_', '', df6$Gene)

unique(df6$IR)
unique(df6$Gene)
unique(df6$Group)

length(unique(df6$Group))
length(unique(meta$Group))
setdiff(unique(meta$Group), unique(df6$Group))

#####	Only for sample>=3
#df7 = df6[df6$Group %in% names(table(meta$Group)[table(meta$Group)>2]),]
df7 = df6

###########################################################################################
#####	Draw plot
length(group.level)
setdiff(group.level, unique(df7$Group))

#group.level = group.level[group.level %in% unique(df7$Group)]

setdiff(group.level, unique(df7$Group))

unique(df7$Gene)

gene.level = c('WT', 'LIG4-MT', 'POLQ-MT', 'WT_siCtrl', 'WT_siBRCA2')
ir.level = c('X-ray', 'Carbon-ion')

df7$Gene = factor(df7$Gene, levels=gene.level)
df7$Group = factor(df7$Group, levels=group.level)
df7$IR = factor(df7$IR, levels=ir.level)

df7 = df7[order(factor(df7$Group, levels=group.level)),]
sample.level = unique(df7$Sample)
df7$Sample = factor(df7$Sample, levels=sample.level)Inversions
df7$MutationType = factor(df7$MutationType, levels=c('Deletions', 'Duplications', 'Inversions', 'Translocations'))

sum.df = df7[,c('Sample', 'Count')]
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

colors = c('#8b4512', '#b02222', '#008a8a', '#678b6c')

legend = ggplot(df7, aes(x=Sample, y=Count)) +
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

p = ggplot(df7, aes(x=Sample, y=Count)) +
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

pdf(file=file.path(output.dir, 'no.read.filter.sv.5class.barplot.figure4.pdf'), width=120, height=40)
g = ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k = which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- strip.colors[i]
}

pdf(file=file.path(output.dir, 'no.read.filter.sv.5class.barplot.figure4.pdf'), width=120, height=40)
plot(g)
dev.off()

pdf(file=file.path(output.dir, 'no.read.filter.sv.5class.barplot.legend.figure4.pdf'), width=120, height=40)
legend
dev.off()

