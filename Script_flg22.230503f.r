##### glm.nb analysis of flg22 RNA-seq data
##### starting with read counts per gene data

#### load packages
library(MASS)
library(qvalue)
library(lsa)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(nlme)
library(minpack.lm)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(stringr)

#### import read-counts-per-gene dataset
rm(list=ls())

load('./data/flg22_raw_gene_counts.RData')
rm('raw_data_corrected_sample_labels')
dim(raw_data_corrected_sample_labels_nonzero_rows) # this is the read count data
#[1] 26275   357

genotypes = c("JEPS", "xEPS", "JxPS", "JExS", "JEPx", "xxPS", "xExS", "xEPx", "JxxS", "JxPx", "JExx", "xxxS", "xxPx", "xExx", "Jxxx", "xxxx", "fls2")
flg22_genotype = factor(rep(genotypes, times=21))
flg22_time = factor(rep(rep(c(0,1,2,3,5,9,18), each=17), times=3))
fix.ef = paste(as.character(flg22_genotype), as.character(flg22_time), sep='.')
fix.ef = factor(fix.ef, levels = fix.ef[1:(17*7)])

exp.dat = raw_data_corrected_sample_labels_nonzero_rows
rm('raw_data_corrected_sample_labels_nonzero_rows')  # change the object name

### 90 percentile normalization offset for each library,
cval.90 = apply(exp.dat, 2, function(x) {
  quantile(x, probs = 0.9)
})

### filter out the genes with many unreliable measurements
## remove genes that are 
## more than half of the data points are 0
## AND the 10th largest value is less than 25.
numb.0.v = apply(exp.dat, 1, function(x) {
  return( ( sum(x==0) > 0.5*ncol(exp.dat) ) &
    ( sort(x, decreasing = T)[10] < 25 ) )
} )
exp.dat = exp.dat[!numb.0.v,]
nrow(exp.dat)
#[1] 17423

### Pseudocounts addition
### to avoid glm fitting problem, for each genotype:time with all rep 0s, 
### Add pseudo count of 1 scaled to the 90%-tile value
min(cval.90)
#[1] 27  # 1 pseudocount for this minimum library
ps.c = round(cval.90/min(cval.90))

table(ps.c)
#ps.c
#  1   2   3   4   5   6   7   8   9  11 
# 21  75 136  92  16   8   5   1   2   1 

ps.c.mat = matrix(rep(ps.c, nrow(exp.dat)), nrow=nrow(exp.dat), byrow = T)
exp.one.d = exp.dat + ps.c.mat

sum(exp.dat==0)
#[1] 384237
sum(exp.one.d==0)
#[1] 0  # good
hist(as.matrix(log(exp.one.d)))  # looking reasonable

### adjust cval.90 with pseudo counts
cval.90 = cval.90 + ps.c
offset.cval.90 = log(cval.90/500)  # make 90%-tile read count = 500

##### fit glm.nb model
#### visualize all data, genotype/time/rep, for Fig 1a
### log2 transformed read count data, between-libraries normalized
exp.od.log2 = sapply(1:ncol(exp.one.d), function(x) {
  log2(exp.one.d[,x]) - log2(exp(offset.cval.90[x]))
})
ff.n1 = paste(as.character(fix.ef), rep(c('r1','r2','r3'), each=119), sep = '.')
dimnames(exp.od.log2) = list(rownames(exp.one.d), ff.n1)

### reorder columns
ff.n1 = as.character(matrix(ff.n1, nrow=3, byrow = T))
ff.n = strsplit(ff.n1, '\\.')
genot.n = sapply(ff.n, '[', 1)
genot.nu = unique(genot.n)
genot.nf = factor(genot.n, levels = genot.nu)
time.n = sapply(ff.n, '[', 2)
samp.or = order(as.numeric(time.n))
genot.nf = genot.nf[samp.or]
samp.or = samp.or[order(as.numeric(genot.nf))]


##### all genes.dyn for scaled min-to-max to 0-to-1 for real figure
all.est = exp.od.log2[,ff.n1[samp.or]]

####### heatmap for 18978 genes x 540 libraries treatment/genotype/time/rep
ams.1 = all.est
### Person correlation among genes
row_distances_subset = 1 - cor(t(ams.1)) # 85 sec
###
distances = as.dist(row_distances_subset, diag = FALSE, upper = FALSE) # 6 sec
hclust.genes.all = hclust( distances, method = "complete" ) # 8 sec
save(all.est, hclust.genes.all, file='./outputs/hclust.genes.all.samp.pearson.complete.RData')
ord = hclust.genes.all$order
data_reordered_full = ams.1[rev(ord),]

## parameters for the heat map
repeat_no1 = 4  # width of each column for mock

data_reordered_rescaled = data_reordered_full # no scaling
upr.limit = quantile(data_reordered_full, probs = 0.95) # 95 percentile as max color
lwr.limit = quantile(data_reordered_full, probs = 0.05) # 5 percentile as max color

fit_coef_rep = t(matrix( rep( c(t(data_reordered_rescaled)), each = repeat_no1), ncol = nrow(data_reordered_rescaled)))

gap.pos = 7*3*1:16*repeat_no1
data_reordered_rescaled_with_gaps = fit_coef_rep
separator=NA
for (coln in rev(gap.pos)) {
  data_reordered_rescaled_with_gaps = cbind(data_reordered_rescaled_with_gaps[,1:coln], 
                                            rep(separator, nrow(data_reordered_rescaled_with_gaps)),
                                            data_reordered_rescaled_with_gaps[,(coln+1):ncol(data_reordered_rescaled_with_gaps)])
} # 3sec
colnames(data_reordered_rescaled_with_gaps) = NULL
## add on an extra row: the bottom row is getting chopped off for mysterious reasons. This fixes the problem. (Note: there was no problem with the pdf, just with the console visualization.)
data_reordered_rescaled_with_gaps = rbind(data_reordered_rescaled_with_gaps, data_reordered_rescaled_with_gaps[nrow(data_reordered_rescaled_with_gaps),])
longData <- melt(data_reordered_rescaled_with_gaps)
time_labels <- factor(longData$Var1, levels = rev(unique(longData$Var1)))
sigalloc_labels <- factor(longData$Var2, levels = unique(longData$Var2))
value_labels = longData$value
value_labels = round(value_labels, digits = 2)
value_labels[value_labels==0] = ''
longData = cbind(longData, value_labels)
alloc_vis <- ggplot(longData,
                    aes(y=time_labels, x = sigalloc_labels, fill = value, label=value_labels))
alloc_vis <- alloc_vis + geom_tile() + ggtitle('') + xlab('') + ylab('') + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
alloc_vis <- alloc_vis + scale_fill_gradientn(colours=rev(c(rainbow(35, end=0.23),rainbow(26, start=0.46, end=0.62))), limits = c(lwr.limit, upr.limit), na.value="black")
#print(alloc_vis)

png(file = './outputs/flg22.all.samp.x.genes.all.not.scaled2.png', width = 15*900, height = 6.5*600, res=300)
print(alloc_vis) # 2 min
dev.off()

#####################################
#### fit glm.nb model for mean estimates
gene.names = rownames(exp.one.d)
date()
glm.coefs = apply(exp.one.d, 1, function(exp.d){
  glmnb1 = glm.nb(exp.d ~ -1 + fix.ef + offset(offset.cval.90),
                  link=log)
  glm.c = summary(glmnb1)$coef[,1:2]
  glm.c = glm.c / log(2)  # to convert the base of log to 2 from e.
  est = glm.c[,1]
  sem = glm.c[,2]
  fitted = glmnb1$fitted.values
  d.resid = residuals(glmnb1, type = 'deviance') / sqrt(fitted) /
     log(2) # to convert the base of log to 2 from e
  return(list(est=est, sem=sem, d.resid=d.resid))
})
date()
# ~13 min

save(glm.coefs, genotypes, fix.ef, gene.names, cval.90, ps.c,
     file='./outputs/flg22.glm.nb.220314.RData')

#######################################
###
#### select dynamical genes compared to t=0
#### q < 0.05, mean.diff > 1, at least two consecutive time points (> 1 hr)
#### including fls2
#### return significant genotypes

fef.Ed = levels(fix.ef)
fef.Ed
ffef.Ed = paste0('fix.ef', fef.Ed)
## organize it by genotype
ffef.Ed.g = lapply(genotypes, function(x) {
   ffef.Ed[grep(x,ffef.Ed)]
})
names(ffef.Ed.g) = genotypes

geno.p.md.by.dyn.genes = lapply(glm.coefs, function(x) {
   est = x$est
   sem = x$sem
   
   ea.geno = lapply(genotypes, function(y) {
      geno.n = ffef.Ed.g[[y]]
      est.0 = est[geno.n[1]]
      est.other = est[geno.n[-1]]
      sem.0 = sem[geno.n[1]]
      sem.other = sem[geno.n[-1]]
      mean.diff = est.other - est.0
      sedm = sqrt(sem.other^2 + sem.0^2)
      p.vals = 2* pnorm(abs(mean.diff)/sedm, lower.tail = F)
      return(list(mean.diff=mean.diff, p.vals=p.vals))
   })
   names(ea.geno) = genotypes
   return(ea.geno)
}) # 3 sec

### qval cutoff
all.pvals.vs0h = t(sapply(geno.p.md.by.dyn.genes, function(x) {
   as.numeric(sapply(x, function(y) y$p.vals))
}))
hist(all.pvals.vs0h, freq=F) # looking good
all.qvals.vs0h = qvalue(all.pvals.vs0h)$qvalues
upr.pval = min(sort(all.pvals.vs0h)[sort(all.qvals.vs0h) >= 0.05]) # find boundary pval for qval < 0.05
upr.pval
#[1] 0.06611285 # pval should be smaller than this - so many positives

sig.geno.by.dyn.gene = t(sapply(geno.p.md.by.dyn.genes, function(x) {
   sapply(x, function(y) {
      p.val.pass = y$p.vals < upr.pval
      p.val.pass = p.val.pass[-1] # after 1 hour only
      p.val.pass = p.val.pass[-1] * p.val.pass[-length(p.val.pass)] # 2 consecutive
      md.pass = abs(y$mean.diff) > 1
      md.pass = md.pass[-1] # after 1 hour only
      md.pass = md.pass[-1] * md.pass[-length(md.pass)] # 2 consecutive
      md.dir = y$mean.diff
      md.dir = md.dir[-1] # after 1 hour only
      md.dir = (md.dir[-1] * md.dir[-length(md.dir)]) > 0 # 2 consecutive, same direction 
      geno.p = sum(p.val.pass * md.pass * md.dir > 0)
   })
})) # 2 sec

### any single genotype
sum(apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 0))
#[1] 12537
### any two genotypes
sum(apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 1))
#[1] 11266
### any three genotypes
sum(apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 2))
#[1] 10407
# it doesn't decrease quickly. 

### Take any single genotype
genes.dyn = gene.names[apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 0)]

#############
##### compare each genotype to flg22 with same treatment at same time point
geno.p.md.by.genes.var.geno = lapply(glm.coefs, function(x) {
   est = x$est
   sem = x$sem
   est = matrix(est, nrow=length(genotypes))
   sem = matrix(sem, nrow=length(genotypes))
   est.dfls2 = est[1:16,] - matrix(rep(est[17,],each=16), nrow=16)
   sem.dfls2 = sqrt(sem[1:16,]^2 + matrix(rep(sem[17,],each=16), nrow=16) ^2)
   pval.dfls2 = 2*pnorm(abs(est.dfls2)/sem.dfls2, lower.tail = F)
   return(list(mean.diff=est.dfls2, p.vals=pval.dfls2))
}) 

### qval cutoff
all.pvals.vsr = t(sapply(geno.p.md.by.genes.var.geno, function(x) {
   as.numeric(x$p.vals)
}))
hist(all.pvals.vsr, freq=F) # looking good
all.qvals.vsr = qvalue(all.pvals.vsr)$qvalues
upr.pval = min(sort(all.pvals.vsr)[sort(all.qvals.vsr) >= 0.05]) # find boundary pval for qval < 0.05
upr.pval
#[1] 0.02570349 # pval should be smaller than this

sig.geno.by.gene.var.geno = t(sapply(geno.p.md.by.genes.var.geno, function(x) {
   p.val.pass = x$p.vals < upr.pval
   md.pass = abs(x$mean.diff) > 1
   geno.p = apply(p.val.pass * md.pass, 1, sum)
})) 

### any single genotype
sum(apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 0))
#[1] 13612
### any two genotypes
sum(apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 1))
#[1] 12250
### any three genotypes
sum(apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 2))
#[1] 11356
# it doesn't decrease quickly. 

### Take any single genotype
genes.var.geno = gene.names[apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 0)]

save(genes.dyn,genes.var.geno, file='./outputs/flg22.selected.genes.dynamics.and.geno.variation.RData')
#######

######################
########## heatmap visualization
##### all genes.dyn for scaled min-to-max to 0-to-1 for real figure
all.est = t(sapply(glm.coefs[genes.dyn], function(x) x$est))
### order the columns
all.samps = colnames(all.est)
all.est = all.est[,as.character(matrix(all.samps, ncol=17, byrow = T))]

####### heatmap for 10833 genes.dyn x 180 treatment/genotype/time
ams.1 = all.est
### Person correlation among genes
row_distances_subset = 1 - cor(t(ams.1)) 
###
distances = as.dist(row_distances_subset, diag = FALSE, upper = FALSE) 
hclust.genes.dyn = hclust( distances, method = "complete" ) 
ord = hclust.genes.dyn$order
save(all.est, hclust.genes.dyn, file='./outputs/flg22.hclust.genes.dyn.all.samp.pearson.complete.RData')
data_reordered_full = ams.1[rev(ord),]

## parameters for the heat map
repeat_no1 = 4  # width of each column for mock

geno_maxes = apply(abs(data_reordered_full), 1, max, na.rm=TRUE)
geno_mins = apply(abs(data_reordered_full), 1, min, na.rm=TRUE)
geno_maxes_matrix = matrix( rep(geno_maxes, times = ncol(data_reordered_full)), ncol = ncol(data_reordered_full) )
geno_mins_matrix = matrix( rep(geno_mins, times = ncol(data_reordered_full)), ncol = ncol(data_reordered_full) )
data_reordered_rescaled = (data_reordered_full-geno_mins_matrix) / (geno_maxes_matrix - geno_mins_matrix)
upr.limit = 1
lwr.limit = 0

fit_coef_rep = t(matrix( rep( c(t(data_reordered_rescaled)), each = repeat_no1), ncol = nrow(data_reordered_rescaled)))

gap.pos = 7*1:16*repeat_no1
data_reordered_rescaled_with_gaps = fit_coef_rep
separator=NA
for (coln in rev(gap.pos)) {
  data_reordered_rescaled_with_gaps = cbind(data_reordered_rescaled_with_gaps[,1:coln], 
                                            rep(separator, nrow(data_reordered_rescaled_with_gaps)),
                                            data_reordered_rescaled_with_gaps[,(coln+1):ncol(data_reordered_rescaled_with_gaps)])
} # 3sec
colnames(data_reordered_rescaled_with_gaps) = NULL
## add on an extra row: the bottom row is getting chopped off for mysterious reasons. This fixes the problem. (Note: there was no problem with the pdf, just with the console visualization.)
data_reordered_rescaled_with_gaps = rbind(data_reordered_rescaled_with_gaps, data_reordered_rescaled_with_gaps[nrow(data_reordered_rescaled_with_gaps),])

longData <- melt(data_reordered_rescaled_with_gaps)
time_labels <- factor(longData$Var1, levels = rev(unique(longData$Var1)))
sigalloc_labels <- factor(longData$Var2, levels = unique(longData$Var2))
value_labels = longData$value
value_labels = round(value_labels, digits = 2)
value_labels[value_labels==0] = ''
longData = cbind(longData, value_labels)
alloc_vis <- ggplot(longData,
                    aes(y=time_labels, x = sigalloc_labels, fill = value, label=value_labels))
alloc_vis <- alloc_vis + geom_tile() + ggtitle('') + xlab('') + ylab('') + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
alloc_vis <- alloc_vis + scale_fill_gradientn(colours=rev(c(rainbow(20, end=0.25),rainbow(30, start=0.4, end=0.65))), limits = c(lwr.limit, upr.limit), na.value="black")
#print(alloc_vis)

png(file = './outputs/flg22.all.samp.x.genes.dyn.scaled2.png', width = 15*300, height = 6.5*300, res=300)
print(alloc_vis) 
dev.off()


##########################################################################
##########################################################################
###################
#### distribution of all expression level values by boxplot

#### load workspace
rm(list=ls())
load('./outputs/flg22.glm.nb.220314.RData')

all.g.est = t(sapply(glm.coefs, function(x) x$est))
dim(all.g.est)
#[1] 17423   119
boxplot(all.g.est, cex=0)
#saved as "flg22.all.log2.exp.val.distr.pdf"
# no outliers shown, very consistent across samples
boxplot(as.numeric(all.g.est), cex=0)

### expression levels of FLS2, BAK1, SERK4, SERK5, BIK1, PBL1, PBL9, PBL11
### "AT5G46330", "AT4G33430", "AT2G13790", "AT2G13800", "AT2G39660", "AT3G55450", "AT1G07570", "AT5G02290"
r.g.n = c("AT5G46330", "AT4G33430", "AT2G13790", "AT2G13800", "AT2G39660", "AT3G55450", "AT1G07570", "AT5G02290")
name.rgn = c("FLS2", "BAK1", "SERK4", "SERK5", "BIK1", "PBL1", "PBL9", "PBL11")
names(name.rgn) = r.g.n
time.p=c(0,1,2,3,5,9,18)

pdf('./outputs/eight.FLS2.signaling.genes.pdf', height = 10, width = 7.5)
opar=par(mfrow=c(3,2), mar=c(4.4,2,0.5,0.5))
for (gene in r.g.n) {
   g.exp = all.g.est[gene,]
   g.exp = t(matrix(g.exp, nrow=17))
   colnames(g.exp) = genotypes
   y.min = min(g.exp)
   y.max = max(g.exp)
   y.max2 = y.min + 1.7 *(y.max-y.min)
   y.min2 = y.min - 0.3 *(y.max-y.min)
   matplot(sqrt(time.p), g.exp, type='l', col=rainbow(18, end=0.65), lty=rep(1:2, 9),
           xlab=NA, ylab=NA,
           ylim=c(y.min2, y.max2), xlim=c(0,sqrt(18)), xaxt='n', lwd=2)
   rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
   matlines(sqrt(time.p), g.exp, type='l', col=rainbow(17, end=0.65), lty=rep(1:2, 9), lwd=2)
   axis(1, at=sqrt(time.p), labels=as.character(time.p))
   text(sqrt(6),y.min2+0.1*(y.max2-y.min2), name.rgn[gene], pos=2, cex=1.5)
   legend('topleft', genotypes, col=rainbow(17, end=0.65), lty=rep(1:2, 9),
          lwd=2, ncol=5, cex=0.9, bg='gray80')
}
par(opar)
dev.off()


######################
##### mock value inference using fls2 dynamics

rm(list=ls())
load('./outputs/flg22.glm.nb.220314.RData')
load('./outputs/flg22.selected.genes.dynamics.and.geno.variation.RData')

### time points
time.p = c(0,1,2,3,5,9,18)
time.pr = rep(time.p, 3)
time.pr.sq = sqrt(time.pr)
time.pr.sq1 = time.pr.sq / max(time.pr.sq) *2 -1
t1.sq = seq(0,sqrt(18),0.02)
t1 = t1.sq^2
t1.sq1 = t1.sq/max(t1.sq) * 2 -1

### fls2 data for genes.dyn
### make it with dev resid

date()
G.corrected.tcv = lapply(glm.coefs[genes.dyn], function(x) {
   est=x$est
   sem = x$sem
   est1 = est[grep('fls2', names(est))]
   sem1 = sem[grep('fls2', names(sem))]
   max.v = max(est1)
   min.v = min(est1)
   v0 = est1[1]
   est1r = rep(est1, 3)
   sem1r = rep(sem1, 3)
   weight1r = 1/sem1r^2
   d.resid = x$d.resid
   d.resid1 = d.resid[grep('fls2', fix.ef)]
   d.vals = est1r + d.resid1
   
   ### fit polynomial up to 4th order
   model.p = lm(d.vals ~ poly(time.pr.sq1, 4),
                weights=1/sem1r^2)
   model.p = step(model.p, trace = 0)

   ### plots for checking
   #ymin = min(d.vals); ymax= max(d.vals) * 1.3
   #plot(time.pr, d.vals, ylim=c(ymin-(ymax-ymin)*0.2, ymax+(ymax-ymin)*0.2))
   #lines(t1, predict(model.p, newdata = data.frame(time.pr.sq1 = t1.sq1)), col='red')
   #text(5.8, ymin-(ymax-ymin)*0.1, paste0(gene,', ', length(coef(model.p))-1), 
   #     pos=2, col='blue' )

   ### the values for time 0:6 for GUS, Ed
   fls2.v = fitted(model.p)[1:7]
   names(fls2.v) = time.p
   
   ### the associated sem
   fls2.ci = predict(model.p, interval = 'c')[1:7,]
   fls2.ci = fls2.ci[,'upr'] - fls2.ci[,'fit']
   fls2.sem = fls2.ci/ qnorm(0.975)
   names(fls2.sem) = time.p
   
   ######## 
   #### adjust for 0hr of each genotype
   est.mat = matrix(est, nrow=17)
   est.0h.mf0 = est.mat[,1] - fls2.v[1]  # this is a constant for each genotype
   est.mat.comp = sapply(1:7, function(x) {
     est.mat[,x] - fls2.v[x] - est.0h.mf0
   })
   
   sem.mat = matrix(sem, nrow=17)
   sem.mat.comp = sapply(1:7, function(x) {
     sqrt(sem.mat[1:16, x]^2 + fls2.sem[x]^2)
   })
   
   colnames(est.mat.comp) = colnames(sem.mat.comp) = names(fls2.v) 
   rownames(est.mat.comp) = genotypes
   rownames(sem.mat.comp) = genotypes[1:16]
   
   return(list(flg22.vsf = est.mat.comp[1:16,], flg22.vsf.se = sem.mat.comp, 
               fls2.vsf = est.mat.comp[17,],
               mock.7tp = fls2.v, mock.7tp.se = fls2.sem))
   
}) # 30 sec

date()
save(G.corrected.tcv, file='./outputs/flg22.corrected.fls2.220314.RData')

###############################
###### visualization by heatmap

rm(list=ls())
#### load data
load('./outputs/flg22.corrected.fls2.220314.RData')
load('./outputs/flg22.hclust.genes.dyn.all.samp.pearson.complete.RData')

#### compare flg22.vsf with all.est, subtracted 0 time point
#### heatmap for 12537 genes.dyn x 119 genotype/time

#### (1) scaled for each gene
#### all est minus 0h values
all.0h.cols = grep('\\.0$', colnames(all.est))
all.0h.cols = rep(all.0h.cols, ea=7)
all.est.0h = all.est[,all.0h.cols]
all.est.vs0h = all.est - all.est.0h

#### all est minus mock simulated
all.est.vsf = t(sapply(G.corrected.tcv, 
                       function(x) c(as.numeric(t(x$flg22.vsf)), x$fls2.vsf)))
time.p = names(G.corrected.tcv[[1]]$fls2.vsf)
time.p.17 = rep(time.p, 17)
ams.1 = all.est.vsf

## remove 0 and 1 hr: 0h all 0; 1 hr not reliable
ams.1 = ams.1[,-grep('0|1$', time.p.17)]

row_distances_subset = 1 - cosine(t(ams.1)) # > 10 min
###
distances = as.dist(row_distances_subset, diag = FALSE, upper = FALSE) 
hclust.all.est.vsf = hclust( distances, method = "complete" ) 
ord = hclust.all.est.vsf$order

data_reordered_full = ams.1[ord,]
all.est.reord = all.est.vs0h[ord,]

## parameters for the heat map
repeat_no1 = 4  # width of each column for mock

geno_maxes = apply(abs(cbind(all.est.reord, data_reordered_full) ), 1, max, na.rm=TRUE) 
geno_maxes_matrix = matrix( rep(geno_maxes, times = ncol(data_reordered_full)), ncol = ncol(data_reordered_full) )
data_reordered_rescaled = data_reordered_full / geno_maxes_matrix
upr.limit = 1
lwr.limit = -1

fit_coef_rep = t(matrix( rep( c(t(data_reordered_rescaled)), each = repeat_no1), ncol = nrow(data_reordered_rescaled)))

gap.pos = 5*1:16*repeat_no1
data_reordered_rescaled_with_gaps = fit_coef_rep
separator=NA
for (coln in rev(gap.pos)) {
  data_reordered_rescaled_with_gaps = cbind(data_reordered_rescaled_with_gaps[,1:coln], 
                                            rep(separator, nrow(data_reordered_rescaled_with_gaps)),
                                            data_reordered_rescaled_with_gaps[,(coln+1):ncol(data_reordered_rescaled_with_gaps)])
} 

colnames(data_reordered_rescaled_with_gaps) = NULL
## add on an extra row: the bottom row is getting chopped off for mysterious reasons. This fixes the problem. (Note: there was no problem with the pdf, just with the console visualization.)
data_reordered_rescaled_with_gaps = rbind(data_reordered_rescaled_with_gaps, data_reordered_rescaled_with_gaps[nrow(data_reordered_rescaled_with_gaps),])

longData <- melt(data_reordered_rescaled_with_gaps)
time_labels <- factor(longData$Var1, levels = rev(unique(longData$Var1)))
sigalloc_labels <- factor(longData$Var2, levels = unique(longData$Var2))
value_labels = longData$value
value_labels = round(value_labels, digits = 2)
value_labels[value_labels==0] = ''
longData = cbind(longData, value_labels)
alloc_vis <- ggplot(longData,
                    aes(y=time_labels, x = sigalloc_labels, fill = value, label=value_labels))
alloc_vis <- alloc_vis + geom_tile() + ggtitle('') + xlab('') + ylab('') + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
bl <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")(200)
re <- colorRampPalette(brewer.pal(9, "Reds"), space="Lab")(200)
alloc_vis <- alloc_vis + scale_fill_gradientn(colours=c(bl,"white", re), limits = c(-upr.limit, upr.limit), na.value="black")
#print(alloc_vis)

png(file = './outputs/flg22.vsf.x.genes.dyn.scaled2.v220315.png', width = 15*300, height = 6.5*300, res=300)
print(alloc_vis) 
dev.off() #15 sec


###################
#### (2) not scaled

## parameters for the heat map
repeat_no1 = 4  # width of each column for mock

data_reordered_rescaled = data_reordered_full
upr.limit = max(abs(data_reordered_rescaled))

fit_coef_rep = t(matrix( rep( c(t(data_reordered_rescaled)), each = repeat_no1), ncol = nrow(data_reordered_rescaled)))

gap.pos = 5*1:16*repeat_no1
data_reordered_rescaled_with_gaps = fit_coef_rep
separator=NA
for (coln in rev(gap.pos)) {
  data_reordered_rescaled_with_gaps = cbind(data_reordered_rescaled_with_gaps[,1:coln], 
                                            rep(separator, nrow(data_reordered_rescaled_with_gaps)),
                                            data_reordered_rescaled_with_gaps[,(coln+1):ncol(data_reordered_rescaled_with_gaps)])
} 

colnames(data_reordered_rescaled_with_gaps) = NULL
## add on an extra row: the bottom row is getting chopped off for mysterious reasons. This fixes the problem. (Note: there was no problem with the pdf, just with the console visualization.)
data_reordered_rescaled_with_gaps = rbind(data_reordered_rescaled_with_gaps, data_reordered_rescaled_with_gaps[nrow(data_reordered_rescaled_with_gaps),])
longData <- melt(data_reordered_rescaled_with_gaps)
time_labels <- factor(longData$Var1, levels = rev(unique(longData$Var1)))
sigalloc_labels <- factor(longData$Var2, levels = unique(longData$Var2))
value_labels = longData$value
value_labels = round(value_labels, digits = 2)
value_labels[value_labels==0] = ''
longData = cbind(longData, value_labels)
#dev.new(width = 15, height = 7)
alloc_vis <- ggplot(longData,
                    aes(y=time_labels, x = sigalloc_labels, fill = value, label=value_labels))
alloc_vis <- alloc_vis + geom_tile() + ggtitle('') + xlab('') + ylab('') + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
bl <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")(100)
re <- colorRampPalette(brewer.pal(9, "Reds"), space="Lab")(200)
alloc_vis <- alloc_vis + scale_fill_gradientn(colours=c(bl,"white", re), limits = c(-upr.limit/2 * 0.8, upr.limit * 0.8), na.value="black")
#print(alloc_vis)

png(file = './outputs/flg22.vsf.x.genes.dyn.nonscaled.v220315.png', width = 15*300, height = 6.5*300, res=300)
print(alloc_vis) 
dev.off()


############################
#########################
##### fitting gamma pdf to the corrected flg22 data

rm(list=ls())
#### load data
load('./outputs/flg22.corrected.fls2.220314.RData')

genotypes = rownames(G.corrected.tcv[[1]]$flg22.vsf)
time.p = colnames(G.corrected.tcv[[1]]$flg22.vsf)

#### gene selection v210906
#### by q < 0.05, > 3-fold (same direction) at least two consecutive time points, after 1hour, in at least one genotype
pvals.flg22.dyn = lapply(G.corrected.tcv, function(x) {
   pvals = t(sapply(genotypes, function(geno) {
      est = x$flg22.vsf[geno,3:7]
      sem = x$flg22.vsf.se[geno,3:7]
      2* pnorm(abs(est/sem), lower.tail = F)
   }))
   rownames(pvals) = genotypes
   colnames(pvals) = time.p[3:7]
   return(pvals)
})
all.pvals = unlist(pvals.flg22.dyn)
hist(all.pvals, breaks=50) # very good
all.qvals = qvalue(all.pvals)$qvalues
th.pval = min(all.pvals[all.qvals > 0.05])
th.pval
#[1] 0.1178119  a lot of positives

zero.mat = matrix(0, nrow=16, ncol=7)
genes.flg22.dyn = sapply(names(G.corrected.tcv), function(x) {
   est = G.corrected.tcv[[x]]$flg22.vsf
   pvals = pvals.flg22.dyn[[x]]
   est.m = zero.mat
   est.m[est >= log2(3)] = 1
   est.m[est <= -log2(3)] = -1
   pvals = pvals < th.pval
   e.p = est.m[,3:7] * pvals
   e.pt = e.p[,1:4] * e.p[,2:5] # two consecutive time points from 2 to 18 hours
   e.pt.s = apply(e.pt, 1, function(x) sum(x > 0))
   
   # positive or negative - consistent
   est.pn = apply(est, 1, function(x) {
     max.v = max(x)
     min.v = min(x)
     pn.ok = T
     if (max.v * min.v < 0) {
       pn.ok = abs(log2(abs(max.v) / abs(min.v))) > 2 # at least 4 times when posi nega
     }
     return(pn.ok)
   })

   sum(e.pt.s * est.pn) > 7  # at least in 8 genotypes
})
genes.flg22.dyn = names(G.corrected.tcv)[genes.flg22.dyn]
length(genes.flg22.dyn)
#[1] 2946 genes

flg22.dyn.dat = t(sapply(G.corrected.tcv[genes.flg22.dyn], function(x) {
   Ed.vals = x$flg22.vsf
   Ed.vals = as.numeric(t(Ed.vals))
   return(Ed.vals)
}))

colnames(flg22.dyn.dat) = paste0(rep(genotypes, ea=7), '.', rep(time.p, 16), 'h')
save(flg22.dyn.dat, file='./outputs/flg22.dyn.dat.v220315.RData')

##########
##### heatmap
ams.1 = flg22.dyn.dat

### cosine correlation among genes
row_distances_subset = 1 - cosine(t(ams.1)) 
###
distances = as.dist(row_distances_subset, diag = FALSE, upper = FALSE) 
hclust.flg22.dyn = hclust( distances, method = "complete" ) 
save(flg22.dyn.dat, hclust.flg22.dyn, file='./outputs/hclust.flg22.dyn.pearson.complete.RData')
ord = hclust.flg22.dyn$order
data_reordered_full = ams.1[ord,]

## parameters for the heat map
repeat_no1 = 4  # width of each column 

data_reordered_rescaled = data_reordered_full
upr.limit = max(abs(data_reordered_rescaled))

fit_coef_rep = t(matrix( rep( c(t(data_reordered_rescaled)), each = repeat_no1), ncol = nrow(data_reordered_rescaled)))

gap.pos = 7*1:15*repeat_no1
data_reordered_rescaled_with_gaps = fit_coef_rep
separator=NA
for (coln in rev(gap.pos)) {
   data_reordered_rescaled_with_gaps = cbind(data_reordered_rescaled_with_gaps[,1:coln], 
                                             rep(separator, nrow(data_reordered_rescaled_with_gaps)),
                                             data_reordered_rescaled_with_gaps[,(coln+1):ncol(data_reordered_rescaled_with_gaps)])
}
colnames(data_reordered_rescaled_with_gaps) = NULL
## add on an extra row: the bottom row is getting chopped off for mysterious reasons. This fixes the problem. (Note: there was no problem with the pdf, just with the console visualization.)
data_reordered_rescaled_with_gaps = rbind(data_reordered_rescaled_with_gaps, data_reordered_rescaled_with_gaps[nrow(data_reordered_rescaled_with_gaps),])
longData <- melt(data_reordered_rescaled_with_gaps)
time_labels <- factor(longData$Var1, levels = rev(unique(longData$Var1)))
sigalloc_labels <- factor(longData$Var2, levels = unique(longData$Var2))
value_labels = longData$value
value_labels = round(value_labels, digits = 2)
value_labels[value_labels==0] = ''
longData = cbind(longData, value_labels)
#dev.new(width = 15, height = 7)
alloc_vis <- ggplot(longData,
                    aes(y=time_labels, x = sigalloc_labels, fill = value, label=value_labels))
alloc_vis <- alloc_vis + geom_tile() + ggtitle('') + xlab('') + ylab('') + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
bl <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")(200)
re <- colorRampPalette(brewer.pal(9, "Reds"), space="Lab")(200)
alloc_vis <- alloc_vis + scale_fill_gradientn(colours=c(bl,"white", re), limits = c(-upr.limit * 0.7, upr.limit * 0.7), na.value="black")
#print(alloc_vis)

png(file = './outputs/flg22.vsf.dyn.nonscaled.v220315.png', width = 15*300, height = 6.5*300, res=300)
print(alloc_vis) 
dev.off()

########################################
##### gamma pdf fit

rm(list=ls())

#### function definition
### gamma pdf defined with shape in log2 (l2a) and peak time (pt), and peak height (ph)
g.distr.pt0.l2a = function(x, ph = 1, l2a = 1.5, pt = 2, t0 = 0) {
   ptx = pt-t0; xn0 = x[x > t0] - t0; am1 = 2^l2a -1
   out = rep(0, length(x))
   out[x > t0] = ph * (exp(1)/ptx)^am1 * xn0^am1 * exp(-am1 * xn0/ptx) 
   return(out)
} 

### AICc
AICc = function(aic, n, k) {
   aic + 2*k*(k+1) / (n-k-1)
}

### expression data visualization genotypes x time
g.exp.plot = function(g.exp, geno, g.name) {
   g.exp = matrix(g.exp, ncol=length(geno))
   colnames(g.exp) = geno
   y.min = min(g.exp)
   y.max = max(g.exp)
   y.max2 = y.min + 1.7 *(y.max-y.min)
   y.min2 = y.min - 0.3 *(y.max-y.min)
   matplot(0:6, g.exp,
           xlab=NA, ylab=NA,
           ylim=c(y.min2, y.max2), xlim=c(0,6), xaxt='n', lwd=2)
   rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
   matlines(0:6, g.exp, type='l', col=rainbow(length(geno), end=0.65), 
            lty=rep(1:2, 20)[1:length(geno)], lwd=2)
   axis(1, at=0:6)
   text(6,y.min2+0.1*(y.max2-y.min2), g.name, pos=2, cex=1.5)
   legend('topleft', geno, col=rainbow(length(geno), end=0.65), 
          lty=rep(1:2, 20)[1:length(geno)],
          lwd=2, ncol=5, cex=0.9, bg='gray80')
}


#### load workspaces
load('./outputs/flg22.glm.nb.220314.RData')
load('./outputs/flg22.dyn.dat.v220315.RData')
load('./outputs/flg22.corrected.fls2.220314.RData')

genes.flg22.dyn = rownames(flg22.dyn.dat)
glm.coefs.gfdyn = glm.coefs[genes.flg22.dyn]
G.c.tcv.dfdyn = G.corrected.tcv[genes.flg22.dyn]
rm(list=c('glm.coefs'))  # remove large unnecessary lists

#### focus on only 16 genotyps
gfdyn.samp = !grepl('fls2|\\.m\\.', as.character(fix.ef))
fix.ef.gfdyn = factor(as.character(fix.ef)[gfdyn.samp])

#### organize the observation and weights values
obs.wei = lapply(genes.flg22.dyn, function(x) {
  d.res = glm.coefs.gfdyn[[x]]$d.resid[gfdyn.samp]
  d.res = matrix(d.res, nrow=16)
  d.res = rbind(d.res[,1:7], d.res[,8:14], d.res[,15:21])
  d.res = as.numeric(t(d.res))
  est = rep(flg22.dyn.dat[x,], 3)
  obs = est + d.res
  sem = rep(as.numeric(t(G.c.tcv.dfdyn[[x]]$flg22.vsf.se)), 3)
  wei = 1/sem^2
  return(list(obs=obs, weights=wei))
})
names(obs.wei) = genes.flg22.dyn


################ use actual observations (d.resid)
##### model related functions

#### gamma time course model fitting
gamma.m.fit = function(gene, u.pt.lim=25) {
  ed.t = as.numeric(flg22.dyn.dat[gene,])
  ed.t = matrix(ed.t, nrow=16, byrow = T)
  obs.g = obs.wei[[gene]]$obs
  obs.g = matrix(obs.g, ncol=7, byrow = T)
  weights = obs.wei[[gene]]$weights
  weights = matrix(weights,ncol=7, byrow = T)
  rownames(obs.g) = rownames(weights) = rep(genotypes[1:16],3)
  rownames(ed.t) = genotypes[1:16]
  
  ### return NA for a badly behaving gene
  ## if max of abs value is at 1h
  geno.1h = apply(ed.t, 1, function(x) {
    which.max(abs(x)) == 2 
  })
  ## if max of abs value sign is different in any of the genot
  geno.sign = apply(ed.t, 1, function(x) sign(x[which.max(abs(x))]))
  ## if max of abs value is at 1h in > 2 genotypes OR 
  ## the sign for max(abs) value in any genotypes different 
  ## then return NA
  if (sum(geno.1h) > 2 | (sum(geno.sign > 0) > 0 & sum(geno.sign < 0) > 0)) {
    return(NA)
  } else {
    
    ### estimate starting parameter values
    ## whether 18h point should be dropped for each genotype
    ed.t = t(apply(ed.t, 1, function(x) {
      x1 = x
      if (abs(x[5]) - abs(x[6]) > 0.1 * abs(x[5]) & abs(x[7]) - abs(x[6]) > 0 &
          abs(x[7] - x[6]) > abs(x[7]) * 0.1 & 
          sign(x[5]) == sign(x[6]) & sign(x[6]) == sign(x[7])) {
        x1[7] = NA
      }
      return(x1)
    }))
    
    if (sum(is.na(ed.t[,7])) > 1) ed.t[,7] = obs.g[,7] = NA
    
    ## remove 1h data
    ed.t[,2] = obs.g[,2] = NA
    
    ## up or down regulation? ud
    ud = sign(as.numeric(ed.t)[which.max(abs(as.numeric(ed.t)))])
    ## pt
    max.t = apply(ed.t, 1, function(x) which.max(ud * x))
    pt.geno = time.p[max.t]
    names(pt.geno) = names(max.t)
    ## ph
    ph.geno = sapply(names(max.t), function(x) {
      ed.t[x, max.t[x]]
    })
    names(ph.geno) = names(max.t)
    
    ### preliminary fit for each genotype
    t = time.p
    l.ph.r = 0.7; u.ph.r = 3
    l.l2a = log2(3); u.l2a = 3.5
    l.pt = 2.4; u.pt = u.pt.lim - 2
    pre.paras = c()
    pre.paras.se = c()
    pre.paras.df = c()
    for (genot1 in names(max.t)) {
      ed.tg = ed.t[genot1,]
      ph1 = as.numeric(ph.geno[genot1])
      pt1 = as.numeric(pt.geno[genot1])
      l2a1 = 3
      if (ud < 0) {
        l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
      } else {
        l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
      }
      if (pt1 < 18) u.pt = 10
      pre.gm = nlsLM(ed.tg ~ g.distr.pt0.l2a(t, ph=ph, l2a=l2a, pt=pt, t0=0),
                     start=c(ph=ph1, l2a=l2a1, pt=pt1),
                     lower=c(ph=l.ph1, l2a=l.l2a, pt=l.pt),
                     upper=c(ph=u.ph1, l2a=u.l2a, pt=u.pt))
      pre.paras = rbind(pre.paras, coef(pre.gm))
      pre.paras.se = rbind(pre.paras.se, summary(pre.gm)$coef[,2])
      pre.paras.df = c(pre.paras.df, summary(pre.gm)$df[2])
    }
    rownames(pre.paras) = rownames(pre.paras.se) = names(pre.paras.df) = names(max.t)
    
    ## fit to "data"
    l.ph.r = l.l2a.r = l.pt.r = 0.7 
    u.ph.r = u.l2a.r = u.pt.r = 1.4
    tx = rep(t, each=3)
    
    pre.paras1 = c()
    pre.paras1.se = c()
    pre.paras1.df = c()
    skipped.genot = c()
    for (genot1 in names(max.t)) {
      genot.rows = which(rownames(obs.g) == genot1)
      obs.tg = as.numeric(obs.g[genot.rows,])
      weights.tg = as.numeric(weights[genot.rows,])
      ph1 = as.numeric(pre.paras[genot1, 'ph'])
      pt1 = as.numeric(pre.paras[genot1, 'pt'])
      l2a1 = as.numeric(pre.paras[genot1, 'l2a'])
      if (ud < 0) {
        l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
      } else {
        l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
      }
      l.l2a = l2a1 * l.l2a.r; u.l2a = l2a1 * u.l2a.r
      if (l.l2a < log2(3)) l.l2a = log2(3)
      if (u.l2a > 3.5) u.l2a = 3.5
      l.pt = pt1 * l.pt.r; u.pt = pt1 * u.pt.r
      if (u.pt > u.pt.lim - 1) u.pt=u.pt.lim - 1
      pre.gm = tryCatch( nlsLM(obs.tg ~ g.distr.pt0.l2a(tx, ph=ph, l2a=l2a, pt=pt, t0=0),
                               start=c(ph=ph1, l2a=l2a1, pt=pt1),
                               lower=c(ph=l.ph1, l2a=l.l2a, pt=l.pt),
                               upper=c(ph=u.ph1, l2a=u.l2a, pt=u.pt),
                               weights = weights.tg), error = function(e){} )
      if (is.null(pre.gm)) {
        coef.pre.gm = pre.paras[genot1, ]
        se.pre.gm = pre.paras.se[genot1, ]
        df.pre.gm = pre.paras.df[genot1]
        skipped.genot = c(skipped.genot, genot1)
      } else {
        coef.pre.gm = coef(pre.gm)
        se.pre.gm = summary(pre.gm)$coef[,2]
        df.pre.gm = summary(pre.gm)$df[2]
      }
      pre.paras1 = rbind(pre.paras1, coef.pre.gm)
      pre.paras1.se = rbind(pre.paras1.se, se.pre.gm)
      pre.paras1.df = c(pre.paras1.df, df.pre.gm)
    }
    rownames(pre.paras1) = rownames(pre.paras1.se) = names(pre.paras1.df) = names(max.t)
    
    if (length(skipped.genot) <= 3) { 
     return(list(para.vals = pre.paras1,
                  se.vals = pre.paras1.se,
                  df.vals = pre.paras1.df))
    } else {
      return(NA)
    }
  }
}

###### time course plotting for data and model for all genotypes per gene

### data and model time course plots for all genotypes per gene
geno.g.tc.plot = function(gene, tx1=0:6, tx1.sq=0:6, ud='up') {
  ### model values
  para.v = gm.fit[[gene]]$para.vals
  m.v = c()
  for (genot1 in rownames(para.v)) {
    para.geno = para.v[genot1,]
    mg.v = g.distr.pt0.l2a(tx1, ph=para.geno['ph'], l2a=para.geno['l2a'], pt=para.geno['pt'],
                           t0 = 0)
    m.v = cbind(m.v, mg.v)
  }
  
  ### plot y-range
  y.min1 = min(m.v)
  y.max1 = max(m.v)
  
  g.exp.Ed = t(G.corrected.tcv[[gene]]$flg22.vsf)[-2,]
  y.min2 = min(min(g.exp.Ed), y.min1)
  y.max2 = max(max(g.exp.Ed), y.max1)
  
  ### data plot
  matplot(time.p.sq[-2], g.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
          xlab=NA, ylab=NA,
          ylim=c(y.min2, y.max2), xlim=range(time.p.sq), xaxt='n')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  matlines(time.p.sq[-2], g.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=1)
  axis(1, at=time.p.sq[-2], labels=time.p[-2])
  if (ud == 'up') {
    text(0,y.max2-0.1*(y.max2-y.min2), paste(gene,'data'), pos=4, cex=0.8)
  } else {
    text(0,y.min2+0.1*(y.max2-y.min2), paste(gene,'data'), pos=4, cex=0.8)
  }
  
  ### model plot
  matplot(tx1.sq, m.v, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
          xlab=NA, ylab=NA,
          ylim=c(y.min2, y.max2), xlim=range(time.p.sq), xaxt='n')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  matlines(tx1.sq, m.v, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=1)
  axis(1, at=time.p.sq[-2], labels = time.p[-2])
  if (ud == 'up') {
    text(0,y.max2-0.1*(y.max2-y.min2), paste(gene,'model'), pos=4, cex=0.8)
  } else {
    text(0,y.min2+0.1*(y.max2-y.min2), paste(gene,'model'), pos=4, cex=0.8)
  }
}

#############
######## MAIN

###########
#### fit the model 
genotypes = rownames(G.corrected.tcv[[1]]$flg22.vsf)
time.p = colnames(G.corrected.tcv[[1]]$flg22.vsf)
time.p = as.numeric(time.p)

gm.fit = list()
u.pt.lim = 25

date()
for (x in genes.flg22.dyn) {
  if (x == "AT3G02480") {  # manually remove difficult to fit genes
    gm.fit[[x]] = NA
  } else {
    gm.fit[[x]] = gamma.m.fit(gene=x, u.pt.lim = u.pt.lim)
    if (!is.na(gm.fit[[x]][1])) {
      pt.vals2 = sort(gm.fit[[x]]$para.vals[,'pt'])[c(2, 15)]
      if (pt.vals2[2]/pt.vals2[1] > 3) gm.fit[[x]] = NA
    }
  }
} # 100 sec
date()

genes.flg22.dyn[sapply(gm.fit, function(x) is.na(x[1]))]
#462 not-modeled genes out of 
length(gm.fit)
# [1] 2946 genes

#### variables corresponding to the observations and weights.
time.p = as.numeric(colnames(G.c.tcv.dfdyn[[1]]$flg22.vsf))
time.p.sq = sqrt(time.p)
ti = rep(time.p, 16)
ti.r = rep(time.p, each=16, times=3)
geno1 = rep(genotypes[-17], each=7)
geno.r = rep(genotypes[-17], 7*3)
a.indx = rep(1:16, each=7)
ar.indx = rep(1:16, 7*3)


## up or down regulation
ud.v = sapply(gm.fit, function(x){
  if (is.na(x[1])) ph.s = NA else {
    ph.v = x$para.vals[,'ph']
    ph.m = ph.v[which.max(abs(ph.v))]
    if (ph.m < 0) {ph.s = 'down'} else {ph.s = 'up'}
  }
  return(ph.s)
})

## sort the genes according to ph (up down) and pt (early to late) of JEPS
ud.rgu = names(ud.v)[ud.v == 'up' & !is.na(ud.v)]
ud.rgd = names(ud.v)[ud.v == 'down'& !is.na(ud.v)]

pt.v = sapply(gm.fit[c(ud.rgu, ud.rgd)], function(x) {
  x$para.vals['JEPS', 'pt']
})

ud.rgu = ud.rgu[order(pt.v[ud.rgu])]
ud.rgd = ud.rgd[order(pt.v[ud.rgd])]
length(ud.rgu); length(ud.rgd)
#[1] 959
#[1] 1525

save(gm.fit, ud.rgu, ud.rgd, file='./outputs/gm.fit.flg22.220320.RData')


##### plotting
tx1.sq = seq(0, sqrt(18), 0.05)
tx1 = tx1.sq^2

pdf('./outputs/FigS5.data.gammamodel.pdf', height = 10, width = 7.5)
opar=par(mfrow=c(6,4), mar=c(2,2,0.5,0.5))
for (gene in c(ud.rgu, ud.rgd)) {
  geno.g.tc.plot(gene, tx1, tx1.sq, ud = ud.v[gene])
} 
par(opar)
dev.off()

#### by heatmap
### model values at 0:6 hpi
tx1 = c(0, 2, 3, 5, 9, 18)
model.v.0.18 = t(sapply(gm.fit[c(ud.rgu, ud.rgd)], function(x) {
  para.v = x$para.vals
  m.v = c()
  for (genot1 in rownames(para.v)) {
    para.geno = para.v[genot1,]
    mg.v = g.distr.pt0.l2a(tx1, ph=para.geno['ph'], l2a=para.geno['l2a'], pt=para.geno['pt'])
    m.v = cbind(m.v, mg.v)
  }
  as.numeric(m.v)
}))
### data
data.v.0.18 = t(sapply(G.corrected.tcv[c(ud.rgu, ud.rgd)], function(x) t(x$flg22.vsf[,-2])))
### discrepancy
discr.v.0.18 = data.v.0.18 - model.v.0.18

### heatmap
geno.n = rownames(gm.fit[[1]]$para.vals)
ams.1 = cbind(data.v.0.18, model.v.0.18, discr.v.0.18)

col.anno1 = rbind(matrix('',ncol=48,nrow=3), geno.n, matrix('',ncol=48,nrow=2))
col.anno1 = as.character(col.anno1)
colnames(ams.1) = col.anno1

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))

bot_anno4x = HeatmapAnnotation(foo = anno_text(colnames(ams.1), just = 'right', rot=60,
                                               gp = gpar(fontsize = 7), location = 1),
                               show_annotation_name = F)
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

ht_4x = Heatmap(ams.1, col = col, cluster_rows = F, cluster_columns = FALSE,
                name = 'model.fit', show_row_dend = F,
                column_gap = unit(4,"mm"), 
                column_split = c(rep('_Data',96), rep('_Model', 96), rep('Data - Model', 96)),
                row_title_gp = gpar(fontsize = 20),
                show_row_names = F, show_column_names = F, 
                column_title_gp = gpar(fontsize=10), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(240,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno4x,
                height = unit(120, 'mm'))

jpeg("./outputs/FigS4.model.fit2.v230327.jpeg",height = 160, width = 270, units = "mm",res = 300)
draw(ht_4x, annotation_legend_list = lgd)
decorate_heatmap_body('model.fit', column_slice = 1, {
  grid.lines(c(1/16,1/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(2/16,2/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(3/16,3/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(4/16,4/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(5/16,5/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(6/16,6/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(7/16,7/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(8/16,8/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(9/16,9/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(10/16,10/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(11/16,11/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(12/16,12/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(13/16,13/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(14/16,14/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(15/16,15/16),c(0,1), gp=gpar(col='gray40'))
})
decorate_heatmap_body('model.fit', column_slice = 2, {
  grid.lines(c(1/16,1/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(2/16,2/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(3/16,3/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(4/16,4/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(5/16,5/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(6/16,6/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(7/16,7/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(8/16,8/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(9/16,9/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(10/16,10/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(11/16,11/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(12/16,12/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(13/16,13/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(14/16,14/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(15/16,15/16),c(0,1), gp=gpar(col='gray40'))
})
decorate_heatmap_body('model.fit', column_slice = 3, {
  grid.lines(c(1/16,1/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(2/16,2/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(3/16,3/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(4/16,4/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(5/16,5/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(6/16,6/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(7/16,7/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(8/16,8/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(9/16,9/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(10/16,10/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(11/16,11/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(12/16,12/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(13/16,13/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(14/16,14/16),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(15/16,15/16),c(0,1), gp=gpar(col='gray40'))
})
grid.text('A', x = unit(5,'mm'),y = unit(138,'mm'),gp = gpar(fontface = 'bold'))
grid.text('B', x = unit(90,'mm'),y = unit(138,'mm'),gp = gpar(fontface = 'bold'))
grid.text('C', x = unit(170,'mm'),y = unit(138,'mm'),gp = gpar(fontface = 'bold'))

dev.off()


#########
##### TableS2 all parameter values and fitted values
### fitted value matrix model.v.0.6
genotypes = rownames(gm.fit[[1]]$para.vals)
gt = rep(genotypes, each=6)
ti = rep(c(0,2,3,5,9,18), 16)
ti = paste0(ti, 'h')
cna = paste(gt,ti, sep="_")
colnames(model.v.0.18) = cna

### model parameter values
par.v = t(sapply(gm.fit[rownames(model.v.0.18)], function(x) {
  p.mat = x$para.vals
  p.line = as.numeric(t(p.mat))
}))
gt = rep(genotypes, each=3)
ti = rep(c('A','log2s','pt'), 16)
cna = paste(gt,ti, sep="_")
colnames(par.v) = cna

ts.tab = cbind(par.v, model.v.0.18)
write.csv(ts.tab, file='./outputs/TableS2.para.vals.fitted.vals.csv')


########################################################
########################################################
############ averaging model for time point data
############ compare with flg22-PTI

##### load workspace
rm(list = ls())
load('./outputs/flg22.corrected.fls2.220314.RData')

##### function definitions
load('./outputs/NRAM.algorithm.w.mean.se.RData')


###### MAIN ##########

library(R.utils)
ave.m.out.all = lapply(G.corrected.tcv, function(v.out) {
  m.est.geno.h = v.out$flg22.vsf
  m.se.geno.h = v.out$flg22.vsf.se
  # select 16 genotypes * 2-6 hours
  m.est.geno.h = m.est.geno.h[1:16, 3:7]
  m.se.geno.h = m.se.geno.h[1:16, 3:7]
  ave.m.out = c()
  for (hpi in colnames(m.est.geno.h)) {
    m.est.geno = m.est.geno.h[,hpi]
    m.se.geno = m.se.geno.h[,hpi]
    names(m.est.geno) = names(m.se.geno) =  geno.names[names(m.est.geno)]
    ave.m.out = cbind(ave.m.out, run.ave.m(m.est.geno, m.se.geno, 
                                           th.1.method = '', th.2.method = ''))
  }
  colnames(ave.m.out) = colnames(m.est.geno.h)
  return(ave.m.out)
})  # 1min 45 sec
detach('package:R.utils', unload=T)

save(ave.m.out.all, file='./outputs/flg22.ave.m.out.all.220321.RData')

sig.nr = sapply(ave.m.out.all, function(x) {
  sum(x[-16,] !=0) 
})
sum(sig.nr > 0)
#[1] 11299 genes out of 12537 genes
hist(sig.nr[sig.nr > 0], breaks=0:50)

########################
######## visualizatoin by heatmap, first in the order of the ETI gene sets

#### well behaving genes
load('./outputs/flg22.dyn.dat.v220315.RData')
dim(flg22.dyn.dat)
#[1] 2946  112

#### PTI value matrices
PTI.wt.mat = t(sapply(G.corrected.tcv, function(x) {
  x$flg22.vsf['JEPS', 3:7]
}))
dim(PTI.wt.mat)
#[1] 12537     5

PTI.rem.mat = t(sapply(ave.m.out.all, function(x){
  x['remainder',]
}))
dim(PTI.rem.mat)
#[1] 12537     5

PTI.NRAM.mat = t(sapply(ave.m.out.all, function(x){
  as.numeric(t( x[1:15,]))
}))
dim(PTI.NRAM.mat)
#[1] 12537    75
sum(apply(PTI.NRAM.mat,1, function(x) sum(x != 0) ) > 0)
#[1] 11299, not bad

###### hierarchical clustering by cosine correlation and complete linkage
#### max(abs()) value distribution for each matrix
quantile(apply(PTI.wt.mat, 1, function(x) max(abs(x))))
#       0%       25%       50%       75%      100% 
#0.1154262 0.8965036 1.3147473 1.9438809 9.7562373   generally lower than ETI

quantile(apply(PTI.rem.mat, 1, function(x) max(abs(x))))
#      0%      25%      50%      75%     100% 
#0.0000000 0.7887365 1.1742251 1.7641864 9.6036540   generally lower than ETI

quantile(apply(PTI.NRAM.mat, 1, function(x) max(abs(x))))
#       0%       25%       50%       75%      100% 
#0.0000000 0.4046819 0.9013417 1.3040485 7.6178899   a bit higher than ETI 

# since ETI.NRAM.mat has multiple times for each allocation, probably no weighting
# is needed among the matrices.

load('./outputs/ETI.gfa.sorted.genes.v230327.RData')

## compare with the ETI gene set
PE.common = intersect(c(gfa.names1.nega, gfa.names1.posi), names(G.corrected.tcv))
length(PE.common)
#[1] 2949 genes out of 
length(c(gfa.names1.nega, gfa.names1.posi))
#[1] 3262 ETI-responsive modeled genes, not bad

time.p=c(0,1,2,3,5,9,18)
PTI.wt.mat.E = matrix(NA, nrow = length(c(gfa.names1.posi, gfa.names1.nega)),
                      ncol=5)
PTI.rem.mat.E = matrix(NA, nrow = length(c(gfa.names1.posi, gfa.names1.nega)),
                      ncol=5)
PTI.NRAM.mat.E = matrix(NA, nrow = length(c(gfa.names1.posi, gfa.names1.nega)),
                        ncol=75)
rownames(PTI.wt.mat.E) = rownames(PTI.rem.mat.E) = rownames(PTI.NRAM.mat.E) =
  c(gfa.names1.posi, gfa.names1.nega)
colnames(PTI.wt.mat.E) = colnames(PTI.rem.mat.E) = time.p[3:7]

PTI.wt.mat.E[PE.common,] = PTI.wt.mat[PE.common,]
PTI.rem.mat.E[PE.common,] = PTI.rem.mat[PE.common,]
PTI.NRAM.mat.E[PE.common,] = PTI.NRAM.mat[PE.common,]


ams.1 = cbind(PTI.wt.mat.E, PTI.rem.mat.E, PTI.NRAM.mat.E)

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

bot_anno1 = HeatmapAnnotation(foo = anno_text(rep(as.character(c(2,3,5,9,18)), 2), rot = 0, 
                                              just = 'center',gp = gpar(fontsize = 7), location = unit(1,'mm')),
                              annotation_label = 'Time (hpi)',
                              show_annotation_name = F)
cnames = rownames(ave.m.out.all[[1]])[-16]
cnames2 = cbind(matrix(rep('', 2*15), ncol=2), matrix(cnames,ncol=1),matrix(rep('', 2*15), ncol=2))
cnames3 = as.character(t(cnames2))
bot_anno2 = HeatmapAnnotation(foo = anno_text(cnames3, rot = 0, 
                                              just = 'center',gp = gpar(fontsize = 7), location = unit(1,'mm')),
                              show_annotation_name = F)
ht_1 = Heatmap(ams.1[,1:10] , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'PTI.wt.rem', show_row_dend = F,
               column_gap = unit(4,"mm"), row_title_gp = gpar(fontsize = 20),
               column_split = rep(c("PTI, WT","PTI, remainder"),each=5),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=8), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               use_raster = F,width = unit(40,'mm'), show_heatmap_legend = F,
               bottom_annotation = bot_anno1,
               height = unit(150, 'mm'))
ht_2 = Heatmap(ams.1[,11:85] , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'PTI.NRAM', show_row_dend = F,
               column_split = rep('NRAM on the mRNA level response at each time point', 75),
               show_row_names = F, show_column_names = F, 
               column_title_gp = gpar(fontsize=8), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               use_raster = F,width = unit(130,'mm'), show_heatmap_legend = F,
               bottom_annotation = bot_anno2,
               height = unit(150, 'mm'))

hm_PTI.wt.rem = grid.grabExpr(draw(ht_1))
hm_PTI.NRAM = grid.grabExpr(draw(ht_2,
                                 annotation_legend_list = lgd))

jpeg("./outputs/PTI.heatmap.E.v230327.jpeg",height = 185, width = 200, units = "mm",res = 300)
grid.arrange(hm_PTI.wt.rem, hm_PTI.NRAM, layout_matrix = rbind(c(1,1,2,2,2,2,2),c(1,1,2,2,2,2,2)))
grid.text('Time (hpi)', x = unit(29,'mm'),y = unit(10,'mm'),gp=gpar(fontsize=8) )
decorate_heatmap_body('PTI.NRAM', {
  grid.lines(c(1/15,1/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(2/15,2/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(3/15,3/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(4/15,4/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(5/15,5/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(6/15,6/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(7/15,7/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(8/15,8/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(9/15,9/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(10/15,10/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(11/15,11/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(12/15,12/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(13/15,13/15),c(0,1), gp=gpar(col='gray80'))
  grid.lines(c(14/15,14/15),c(0,1), gp=gpar(col='gray80'))
})
dev.off()


##############################################################################
############ analysis of the gamma model parameters

### modeled genes among PE.common
load('./outputs/gm.fit.flg22.220320.RData')
gm.fit.early = gm.fit[names(gm.fit) %in% PE.common]
## remove NA genes
gm.fit.early = gm.fit.early[sapply(gm.fit.early, function(x) !is.na(x[1]))]
length(gm.fit.early)
#[1] 1196 genes with models out of 2949 PE.common genes

library(R.utils)
###### NRAM for pt
pt.nram = t(sapply(gm.fit.early, function(x) {
  pt.est = x$para.vals[,'pt']
  pt.se = x$se.vals[,'pt']
  names(pt.est) = names(pt.se) = geno.names[names(pt.est)]
  ave.m.out = run.ave.m(pt.est, pt.se)
}))
sum(apply(pt.nram, 1, function(x) sum(x[-16] !=0)) > 0)
#[1] 805 genes out of 1196 genes

###### NRAM for ph
ph.nram = t(sapply(gm.fit.early, function(x) {
  ph.est = x$para.vals[,'ph']
  ph.se = x$se.vals[,'ph']
  names(ph.est) = names(ph.se) = geno.names[names(ph.est)]
  ave.m.out = run.ave.m(ph.est, ph.se)
}))
sum(apply(ph.nram, 1, function(x) sum(x[-16] !=0)) > 0)
#[1] 883 genes out of 1196 genes

###### NRAM for log shape (l2a)
l2a.nram = t(sapply(gm.fit.early, function(x) {
  l2a.est = x$para.vals[,'l2a']
  l2a.se = x$se.vals[,'l2a']
  names(l2a.est) = names(l2a.se) = geno.names[names(l2a.est)]
  ave.m.out = run.ave.m(l2a.est, l2a.se)
}))
sum(apply(l2a.nram, 1, function(x) sum(x[-16] !=0)) > 0)
#[1] 344 genes out of 1196 genes
detach('package:R.utils', unload=T)

save(ph.nram, pt.nram, l2a.nram, file='./outputs/gammamodel.para.NRAM.flg22.v230327.RData')

#####################
##### Fig 5 heatmaps of ph, pt NRAM for PTI

### load data
rm(list=ls())
load('./outputs/flg22.dyn.dat.v220315.RData') # corrected fls2 subtracted
load('./outputs/gm.fit.flg22.220320.RData')
load('./outputs/gammamodel.para.NRAM.flg22.v230327.RData')
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered, this is recursive

##### JEPS and xxxx time courses for references
common.g = rownames(flg22.dyn.dat) %in% c(gfa.names1.posi, gfa.names1.nega)
sum(common.g)
#[1] 1470  #genes out of
length(c(gfa.names1.posi, gfa.names1.nega))
#[1] 3262
common.ug = rownames(flg22.dyn.dat) %in% c(gfa.names1.posi)
sum(common.ug)
#[1] 820
common.dg = rownames(flg22.dyn.dat) %in% c(gfa.names1.nega)
sum(common.dg)
#[1] 650

quad = flg22.dyn.dat[common.g, 106:112]
wt  = flg22.dyn.dat[common.g, 1:7]
wt.quad = wt - quad

#### compile and order the models
all.gm.fit = gm.fit[common.g]
sum(is.na(all.gm.fit))
#[1] 274  NA
### remove NA
all.gm.fit = all.gm.fit[!is.na(all.gm.fit)]
length(all.gm.fit)
#[1] 1196

#### ph, log-ratio, JEPS/xxxx
ph.wt.quad = sapply(all.gm.fit, function(x) {
  x$para.vals['JEPS','ph'] - x$para.vals['xxxx','ph']
})


#### pt, log-ratio, JEPS/xxxx
pt.wt.quad = sapply(all.gm.fit, function(x) {
  x$para.vals['JEPS','pt'] - x$para.vals['xxxx','pt']
})

##### NRAM ph, pt, repeat each column to increase the width, rownames correspond to the references
r.n = 7 # repeat numb
ph.nram1 = as.numeric(apply(cbind(ph.wt.quad[rownames(ph.nram)], ph.nram[,-16]), 2, function(x) {
  rep(x, r.n)
}))
ph.nram1 = matrix(ph.nram1, nrow=nrow(ph.nram[,-16]))
rownames(ph.nram1) = rownames(ph.nram)
dim(ph.nram1)
#[1] 1196  112

###
pt.nram1 = as.numeric(apply(cbind(pt.wt.quad[rownames(pt.nram)], pt.nram[,-16]), 2, function(x) {
  rep(x, r.n)
}))
pt.nram1 = matrix(pt.nram1, nrow=nrow(pt.nram[,-16]))
rownames(pt.nram1) = rownames(pt.nram)
dim(pt.nram1)
#[1] 1196  112

###
emp.mat = matrix(NA, nrow = length(c(gfa.names1.posi, gfa.names1.nega)), ncol = ncol(ph.nram1))
rownames(emp.mat) = c(gfa.names1.posi, gfa.names1.nega)

ph.nram2 = pt.nram2 = emp.mat
ph.nram2[rownames(ph.nram1),] = ph.nram1
pt.nram2[rownames(pt.nram1),] = pt.nram1

sum(!is.na(ph.nram2[1:1971,1]))
#[1] 602 upregulated genes were subjected to NRAM
sum(!is.na(ph.nram2[1972:3261,1]))
#[1] 594 downregulated genes were subjected to NRAM

emp.mat = matrix(NA, nrow = length(c(gfa.names1.posi, gfa.names1.nega)), ncol = ncol(quad))
rownames(emp.mat) = c(gfa.names1.posi, gfa.names1.nega)
wt1 = quad1 = wt.quad1 = emp.mat
wt1[rownames(wt),] = wt
quad1[rownames(quad),] = quad
wt.quad1[rownames(wt.quad), ] = wt.quad

## remove 1 and 18 hpi
wt1 = wt1[,-c(2,7)]
quad1 = quad1[,-c(2,7)]
wt.quad1 = wt.quad1[,-c(2,7)]

ams.1 = cbind(quad1, wt.quad1, ph.nram2, pt.nram2)

col.anno1 = rbind(matrix('',ncol=2,nrow=2), 
                  c('xxxx', 'JEPS/xxxx'),
                  matrix('',ncol=2,nrow=2))
col.anno1 = as.character(col.anno1)


col.anno2 = rbind(matrix('',ncol=32,nrow=3), 
                  c('PA,JEPS/xxxx', colnames(ph.nram)[-16], 'PT,JEPS-xxxx', colnames(pt.nram)[-16]),
                  matrix('',ncol=32,nrow=3))
col.anno2 = as.character(col.anno2)
colnames(ams.1) = c(col.anno1, col.anno2)

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))
col.pt = colorRamp2(c(-2, -1, -0.4, 0, 0.4, 1, 2), c("forestgreen","green2", "olivedrab3","white","purple3", "magenta","violetred4"))

bot_anno4a = HeatmapAnnotation(foo = anno_text(colnames(ams.1)[1:122], just = 'right', rot=60,
                                               gp = gpar(fontsize = 7), location = 1),
                               show_annotation_name = F)
bot_anno4b = HeatmapAnnotation(foo = anno_text(colnames(ams.1)[123:234], just = 'right', rot=60,
                                               gp = gpar(fontsize = 7), location = 1),
                               show_annotation_name = F)
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))
lgd.pt = Legend( at=seq(-2,2,0.5),col_fun = col.pt, title = "hour",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
                 title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

ht_4a = Heatmap(ams.1[,1:122] , col = col, cluster_rows = F, cluster_columns = FALSE,
                name = 'ph.nram', show_row_dend = F,
                column_gap = unit(4,"mm"), 
                column_split = c(rep(' PTI-FC',10), rep('NRAM on PTI Peak Amplitude', 112)),
                row_title_gp = gpar(fontsize = 20),
                show_row_names = F, show_column_names = F,
                column_title_gp = gpar(fontsize=10), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(85,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno4a,
                height = unit(120, 'mm'))

ht_4b = Heatmap(ams.1[,123:234] , col = col.pt, cluster_rows = F, cluster_columns = FALSE,
                name = 'pt.nram', show_row_dend = F,
                row_title_gp = gpar(fontsize = 20),
                show_row_names = F, show_column_names = F, 
                column_title = 'NRAM on PTI Peak Time',
                column_title_gp = gpar(fontsize=10), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(75,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno4b,
                height = unit(120, 'mm'))

hm_PTI.ph.NRAM = grid.grabExpr(draw(ht_4a,
                                    annotation_legend_list = lgd))
hm_PTI.pt.NRAM = grid.grabExpr(draw(ht_4b,
                                    annotation_legend_list = lgd.pt))

jpeg("./outputs/Fig5.PTI.ph.pt.heatmap.v230203.jpeg",height = 160, width = 200, units = "mm",res = 300)
grid.arrange(hm_PTI.ph.NRAM, hm_PTI.pt.NRAM,
             layout_matrix = rbind(c(1,1,1,1,1,1,2,2,2,2,2),c(1,1,1,1,1,1,2,2,2,2,2)))
decorate_heatmap_body('ph.nram', column_slice=2, {
  grid.lines(c(1/16,1/16),c(0,1), gp=gpar(col='gray10'))
  grid.lines(c(2/16,2/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(3/16,3/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(4/16,4/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(5/16,5/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(6/16,6/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(7/16,7/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(8/16,8/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(9/16,9/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(10/16,10/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(11/16,11/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(12/16,12/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(13/16,13/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(14/16,14/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(15/16,15/16),c(0,1), gp=gpar(col='gray70'))
})
decorate_heatmap_body('pt.nram', {
  grid.lines(c(1/16,1/16),c(0,1), gp=gpar(col='gray10'))
  grid.lines(c(2/16,2/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(3/16,3/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(4/16,4/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(5/16,5/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(6/16,6/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(7/16,7/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(8/16,8/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(9/16,9/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(10/16,10/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(11/16,11/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(12/16,12/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(13/16,13/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(14/16,14/16),c(0,1), gp=gpar(col='gray70'))
  grid.lines(c(15/16,15/16),c(0,1), gp=gpar(col='gray70'))
})
decorate_heatmap_body('ph.nram', column_slice=1,  {
  grid.lines(c(1/2,1/2),c(0,1), gp=gpar(col='gray10'))
})
dev.off()

#############################################

#### stats for upregulated genes
up.mod.genes = gfa.names1.posi[gfa.names1.posi %in% rownames(ph.nram)]
length(up.mod.genes)
#[1] 602

### ph.nram mean
ph.nram.av = apply(ph.nram[up.mod.genes,], 2, mean)[1:15]
col.b = rep('gold3',15)
col.b[ph.nram.av < 0] = 'aquamarine3'

pdf('./outputs/PTI.ph.nram.av.v230327.pdf')
barplot(rev(ph.nram.av), horiz = T, las=1,
        col= rev(col.b), xlim=c(-0.23,0.15),
        xlab=expression('Mean of the term in NRAM on PTI Peak Height (log'[2]*')'))
box()
abline(v=0, col='gray', lty=3, lwd=2)
dev.off()

### Fig 5C pt.nram mean
pt.nram.av = apply(pt.nram[up.mod.genes,], 2, mean)[1:15]
col.b = rep('purple3',15)
col.b[pt.nram.av < 0] = 'olivedrab3'

pdf('./outputs/Fig5C.PTI.pt.nram.av.v230327.pdf')
# pt.nram.av sign is flipped to show acceleration of peak time
barplot(rev(-pt.nram.av), horiz = T, las=1,
        col= rev(col.b), xlim=c(-0.15,0.2),
        xlab=expression('Mean Contribution to Peak Time Acceleration (hour)'),
        main='NRAM on PTI Peak Time')
box()
abline(v=0, col='gray', lty=3, lwd=2)
dev.off()

#### Fig 5D make an equivalent plot with the PC1, 5hpi
### genotype names
c.names = colnames(flg22.dyn.dat)
c.names = strsplit(c.names, '\\.')
c.names = sapply(c.names, '[', 1)
genotypes = unique(c.names)

### gene set
pti.eti.c.genes = rownames(quad)
pti.exp.dat = flg22.dyn.dat[pti.eti.c.genes,]

## remove 1 and 18 hpi
pti.exp.dat1 = pti.exp.dat[,-c(2+0:15*7, 7+0:15*7)]

#### PC1 values, 5hpi
pca.pti.dyn = prcomp(t(pti.exp.dat1), center = T, scale. = F)
pca.results = pca.pti.dyn$x
pca.var = (pca.pti.dyn$sdev)^2
pca.var.p = pca.var/sum(pca.var)

#### PC1-PC2 plot
pc1.geno.time = pca.results[,1]
pc2.geno.time = pca.results[,2]

col.all = rep(rainbow(18, end=0.65)[1:16], each = 5)
text.all = as.character(rep(c(0,2,3,5,9), 16))
text.all1 = paste0('bold(', text.all, ')')

pdf('./outputs/Fig5D.PCA.16geno.res.tc.numb.v230327.pdf', width=10, height = 7.1)
plot(pc1.geno.time, pc2.geno.time, 
     xlab=paste0('PC1 (', round(pca.var.p[1]*100, digits = 1),' %)'),
     ylab=paste0('PC2 (', round(pca.var.p[2]*100, digits = 1),' %)'),
     xaxt='n', yaxt='n', type='n',
     ylim=c(-29,35))
axis(side=1, at=seq(-80, 40,10), las=1)
axis(side=2, at=seq(-30, 20,10), las=1)
text(pc1.geno.time, pc2.geno.time, parse(text=text.all1), cex=1.2, col=col.all)
legend(-54,35,genotypes, ncol = 6, col=rainbow(18, end=0.65)[1:16],
       pch=15)
dev.off()


pc1.geno.time = matrix(pc1.geno.time, nrow=5)
pc1.gt.3h = pc1.geno.time[3,]
names(pc1.gt.3h) = genotypes


