##### glm.nb analysis of Ed-AvrRpt2 RNA-seq data
##### starting with read counts per gene data

#### load packages
library(MASS)
library(qvalue)
library(lsa)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap) 
library(circlize)
library(grid)
library(gridExtra)
library(ggvenn)
library(stringr)
library(minpack.lm)
library(VennDiagram)
library(plotrix)

###########
###########
#### NRAM algorithm
#### inputs = estimates and standard errors

##### function definitions
#### function to make the averaging matrix, called in run.ave.m function
make.rec.mx = function(all.g = LETTERS[1:4], sub.n = LETTERS[1:4], treat.c='y') {
  dec.n = 2:(2^length(sub.n)) -1
  bin.n = intToBin(dec.n)
  bin.n = strsplit(bin.n,'')
  coef.ei = sapply(bin.n, function(x) {
    paste(sub.n[as.logical(as.integer(x))], collapse = ';')
  })
  coef.ei = sort(coef.ei)
  coef.ei = coef.ei[order(nchar(coef.ei))]  # variables in the model
  
  add.m = sapply(coef.ei, function(y) {
    y1 = unlist(strsplit(y,';'))
    if (length(y1) == 1) {
      y1 = as.integer(grepl(y1, all.g))
    } else {
      y1 = sapply(y1, function(z) {
        as.integer(grepl(z, all.g))
      })
      y1 = apply(y1, 1, prod)
    }
    return(y1)
  })
  rownames(add.m) = all.g
  # add.m is the matrix for additive model
  
  rec.mx = add.m
  co.count = sapply(colnames(rec.mx), str_count, pattern=';')
  for (c.c in 0:max(co.count)) {
    ap.cols = co.count == c.c
    if (sum(ap.cols) == 1) next
    rec.mx1 = rec.mx[,ap.cols]
    rec.mx[,ap.cols] = t(apply(rec.mx1, 1, function(x) {
      one.c = sum(x)
      row.v = x
      if (one.c > 1) {
        non.0 = 1/sum(x)
        row.v[x==1] = non.0
      } 
      return(row.v)
    }))
  }
  if (treat.c == 'y'){
    rec.mx = cbind(rec.mx, remainder=1) 
  }
  return(rec.mx)  ## rec.mx is the matrix for averaging model
}

#### function to find significant genes, called in run.ave.m function
sig.genes = function(gene.n, mean.est, se.v, th.1 = 0.05, th.1.method='') {
  geno.names = names(mean.est)
  p.val.mat = c()
  for (g.n in gene.n) {
    geno.nA = geno.names[grep(g.n, geno.names)]
    str1 = paste0('^(.*)',g.n,'(.*)$', collapse = '')
    str2 = paste0('\\1',tolower(g.n),'\\2', collapse = '')
    geno.na = sub(str1, str2, geno.nA, perl=T)
    geno.Aa = rbind(geno.nA, geno.na)
    p.vals = apply(geno.Aa, 2, function(x) {
      mean.diff = mean.est[x[1]] - mean.est[x[2]]
      semd = sqrt(se.v[x[1]]^2 + se.v[x[2]]^2)
      z.val = abs(mean.diff)/semd
      p.val = 2*pnorm(z.val, lower.tail = F)
      return(p.val)
    })
    p.val.mat = rbind(p.val.mat, p.vals)
  }
  if (th.1.method !='') {
    p.val.mat = matrix(p.adjust(p.val.mat, method = th.1.method), 
                       nrow=nrow(p.val.mat), )
  }
  rownames(p.val.mat) = gene.n
  keep.genes = apply(p.val.mat, 1, min) < th.1  # rejection rate
  return(keep.genes)
}

#### function to run averaging model, full model
### mean.est - a vector of mean estimate for each genotype
### se.v - a vector of standard error for the mean estimate (each genotype)
run.ave.m = function(mean.est, se.v, 
                     th.1=0.05, th.1.method='', th.2=0.05, th.2.method='') {
  all.g = names(mean.est)
  gene.n = unlist(strsplit(all.g[1],''))
  keep.genes = sig.genes(gene.n, mean.est, se.v, th.1, th.1.method)
  act.gene.numb = sum(keep.genes)
  ce.vec = emp.vec
  if (act.gene.numb == 0) {
    rd = lm(mean.est ~ 1, weights = 1/(se.v^2) )
    ce.vec['remainder'] = coef(rd)
  } else {
    sub.n = names(keep.genes)[keep.genes]
    treat.c='y'
    rec.mx = make.rec.mx(all.g=all.g, sub.n=sub.n, treat.c=treat.c)
    # rec.mx, matrix for averaging model
    # fit glm.nb
    m. = rec.mx[all.g,]
    lm.ave = lm(mean.est ~ m. -1, weights = 1/(se.v^2) )
    w.mat = diag(1/(se.v^2))
    resid.sd = sd(resid(lm.ave))
    sol.mat = solve(t(m.) %*% w.mat %*% m.) %*% (t(m.) %*% w.mat)  # matrix to solve for fitting
    se.coef = sqrt(diag(sol.mat %*% diag(se.v^2 + resid.sd^2) %*% t(sol.mat)))
    coef.lm.ave = coef(lm.ave)
    names(coef.lm.ave) = substr(names(coef.lm.ave), 3, nchar(names(coef.lm.ave)))
    p.vals = 2* pnorm(abs(coef.lm.ave)/se.coef, lower.tail = F)
    rem.pos = which(names(p.vals) == 'remainder')
    if (th.2.method != '') {
      p.vals[-rem.pos] = p.adjust(p.vals[-rem.pos], method=th.2.method)
    }
    coef.lm.ave = coef.lm.ave[p.vals < th.2]
    ce.vec[names(coef.lm.ave)] = coef.lm.ave
  }
  return(ce.vec)
}

#####

##### needed objects to run NRAM function
emp.vec = rep(0, 16)
av.m.ce.n = c('J','E','P','S','J;E','J;P','J;S','E;P','E;S','P;S',
              'J;E;P','J;E;S','J;P;S','E;P;S','J;E;P;S','remainder')
names(emp.vec) = av.m.ce.n

####
all.geno = c("JEPS", "xEPS", "JxPS", "JExS", "JEPx", 
              "xxPS", "xExS", "xEPx", "JxxS", "JxPx", "JExx", 
              "xxxS", "xxPx", "xExx", "Jxxx", "xxxx", "r2r1", "GUS")
geno.names = c()
for (geno.n in all.geno[1:16]) {
  if (!grepl('J',geno.n)) gn = 'j' else gn = 'J'
  if (!grepl('E',geno.n)) gn = c(gn, 'e') else gn = c(gn, 'E')
  if (!grepl('P',geno.n)) gn = c(gn, 'p') else gn = c(gn, 'P')
  if (!grepl('S',geno.n)) gn = c(gn, 's') else gn = c(gn, 'S')
  geno.names = c(geno.names, paste0(gn, collapse = ''))
}
names(geno.names) = all.geno[1:16]

save(make.rec.mx, sig.genes, run.ave.m, emp.vec, geno.names, file='./outputs/NRAM.algorithm.w.mean.se.RData')

#######################
###### gamma-pdf model
#### gamma-pdf model with peak hight (ph), time (pt) and log2shape(l2a)
g.distr.pt0.l2a = function(x, ph = 1, l2a = 1.5, pt = 2, t0 = 0) {
  ptx = pt-t0; xn0 = x[x > t0] - t0; am1 = 2^l2a -1
  out = rep(0, length(x))
  out[x > t0] = ph * (exp(1)/ptx)^am1 * xn0^am1 * exp(-am1 * xn0/ptx) 
  return(out)
} 
save(g.distr.pt0.l2a, file='./outputs/gamma.pdf.model.RData')

####################


####################################
####################################

#### import read-counts-per-gene dataset
exp.dat = read.delim("./data/final_selected_ETI_libraries_gene_counts.R_input_matrix.txt", header = T, row.names = 1)
dim(exp.dat)
#[1] 33603   540

### 90 percentile normalization offset for each library,
cval.90 = apply(exp.dat, 2, function(x) {
  quantile(x, probs = 0.9)
})

### filter out the genes with many unreliable measurements
## remove genes that are 
## more than half of the data points are 0
## AND the 15th largest value is less than 25.
numb.0.v = apply(exp.dat, 1, function(x) {
  return( ( sum(x==0) > 0.5*ncol(exp.dat) ) &
    ( sort(x, decreasing = T)[15] < 25 ) )
} )
exp.dat = exp.dat[!numb.0.v,]
nrow(exp.dat)
#[1] 18978

### Pseudocounts addition
### to avoid glm fitting problem, for each genotype:time with all rep 0s, 
### Add pseudo count of 1 scaled to the 90%-tile value
min(cval.90)
#[1] 55  # 1 pseudocount for this minimum library
ps.c = round(cval.90/min(cval.90))

table(ps.c)
#ps.c
# 1   2   3   4   5   6 
#28 203 208  83  15   3 

ps.c.mat = matrix(rep(ps.c, nrow(exp.dat)), nrow=nrow(exp.dat), byrow = T)
exp.one.d = exp.dat + ps.c.mat

sum(exp.dat==0)
#[1] 415804
sum(exp.one.d==0)
#[1] 0  # good
hist(as.matrix(log(exp.one.d)))  # looking reasonable

### adjust cval.90 with pseudo counts
cval.90 = cval.90 + ps.c
offset.cval.90 = log(cval.90/500)  # make 90%-tile read count = 500

#### fit glm.nb model
### set up fixed effects
genotypes = c("JEPS", "xEPS", "JxPS", "JExS", "JEPx", 
              "xxPS", "xExS", "xEPx", "JxxS", "JxPx", "JExx", 
              "xxxS", "xxPx", "xExx", "Jxxx", "xxxx", "r2r1", "GUS")
#"GUS" was Ed-GUS transgene in wildtype background

Ed.time = as.vector(paste(c(0:6),"h",sep=""))
Ed.time = as.vector(paste("Ed", Ed.time, sep="."))
#[1] "Ed.0h" "Ed.1h" "Ed.2h" "Ed.3h" "Ed.4h" "Ed.5h"
mock.time = as.vector(paste(c(0,2,5),"h",sep=""))
mock.time = as.vector(paste("m", mock.time, sep="."))
#[1] "m.0h" "m.2h" "m.5h"

Ed.fix.ef = as.vector(outer(genotypes,Ed.time,paste,sep="."))
mock.fix.ef = as.vector(outer(genotypes, mock.time, paste, sep="."))
fix.ef = c(mock.fix.ef, Ed.fix.ef, mock.fix.ef, Ed.fix.ef, mock.fix.ef, Ed.fix.ef)

rep.numb = 3
fix.ef = factor(fix.ef, levels=fix.ef[1:(length(fix.ef)/rep.numb)])
save(exp.one.d, offset.cval.90, fix.ef, cval.90, ps.c, genotypes, file='./outputs/exp.one.d.RData')

#### FigS1. visualize all data, genotype/time/rep for all genes

### load data
rm(list=ls())
load('./outputs/exp.one.d.RData')

### log2 transformed read count data, between-libraries normalized
exp.od.log2 = sapply(1:ncol(exp.one.d), function(x) {
  log2(exp.one.d[,x]) - log2(exp(offset.cval.90[x]))
})
ff.n1 = paste(as.character(fix.ef), rep(c('r1','r2','r3'), each=180), sep = '.')
dimnames(exp.od.log2) = list(rownames(exp.one.d), ff.n1)

### reorder columns
ff.n1 = as.character(matrix(ff.n1, nrow=3, byrow = T))
ff.n = strsplit(ff.n1, '\\.')
genot.n = sapply(ff.n, '[', 1)
genot.nu = unique(genot.n)
genot.nf = factor(genot.n, levels = genot.nu)
treat.n = sapply(ff.n, '[', 2)
time.n = sapply(ff.n, '[', 3)
samp.or = order(time.n)
genot.nf = genot.nf[samp.or]
samp.or = samp.or[order(as.numeric(genot.nf))]
treat.n = treat.n[samp.or]
samp.or = samp.or[order(treat.n, decreasing=T)]

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
repeat_no1 = 7  # width of each column for mock
repeat_no2 = 3  # width of each column for Ed

data_reordered_rescaled = data_reordered_full # no scaling
upr.limit = quantile(data_reordered_full, probs = 0.95) # 95 percentile as max color
lwr.limit = quantile(data_reordered_full, probs = 0.05) # 5 percentile as max color

fit_coef_rep = cbind( t(matrix( rep( c(t(data_reordered_rescaled[,1:162])), each = repeat_no1), ncol = nrow(data_reordered_rescaled))),
                      t(matrix( rep( c(t(data_reordered_rescaled[,163:540])), each = repeat_no2), ncol = nrow(data_reordered_rescaled))))

gap.pos = c(3*3*1:18*repeat_no1, 3*3*18*repeat_no1+1:17*7*3*repeat_no2)
data_reordered_rescaled_with_gaps = fit_coef_rep
separator=NA
for (coln in rev(gap.pos)) {
  data_reordered_rescaled_with_gaps = cbind(data_reordered_rescaled_with_gaps[,1:coln], 
                                            rep(separator, nrow(data_reordered_rescaled_with_gaps)),
                                            data_reordered_rescaled_with_gaps[,(coln+1):ncol(data_reordered_rescaled_with_gaps)])
} # 6sec
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

png(file = './outputs/FigS1.all.540lib.x.18978genes.all.not.scaled2.png', width = 15*900, height = 6.5*600, res=300)
print(alloc_vis) # 2 min
dev.off()

#####################################
#### fit glm.nb model for mean estimates
### load data
rm(list=ls())
load('./outputs/exp.one.d.RData')

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
# ~40 min

save(glm.coefs, genotypes, fix.ef, gene.names, cval.90, ps.c,
     file='./outputs/Ed-AvrRpt2_DI_210223.RData')


### non-endogenous gene
gnames = names(glm.coefs)
ne.genes = gnames[!grepl('^AT.G\\d{5}', gnames)]
ne.genes
#[1] "ATtransgene_Inducible"     "ATtransgene_XVE_fusion_UTR_and_terminator"
# These are AvrRpt2/GUS and the XVE fusion
### four sector hub genes - to be removed as they could have direct mutation effects
four.genes = c('AT5G42650', 'AT5G03280', 'AT3G52430', 'AT1G74710')
names(four.genes) = c('DDE2','EIN2','PAD4','SID2')
### in case of comparison with r2r1, RPS2 gene should be removed
r.genes = c('AT4G26090', 'AT3G07040')
names(r.genes) = c('RPS2','RPM1')

## checking
c(ne.genes, four.genes, r.genes) %in% gnames
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
save(ne.genes, four.genes, r.genes, file='./outputs/genes.to.be.removed.RData')

#######################################
###
#### select dynamical genes compared to t=0, only with Ed treatment
#### q < 0.05, mean.diff > 1, at least two consecutive time points (> 1 hr)
#### including r2r1 and GUS
#### return significant genotypes

fef.Ed = levels(fix.ef)[grep('Ed.', levels(fix.ef))] 
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
}) # 7 sec

### qval cutoff
all.pvals.vs0h = t(sapply(geno.p.md.by.dyn.genes, function(x) {
   as.numeric(sapply(x, function(y) y$p.vals))
}))
hist(all.pvals.vs0h, freq=F) # looking good
all.qvals.vs0h = qvalue(all.pvals.vs0h)$qvalues
upr.pval = min(sort(all.pvals.vs0h)[sort(all.qvals.vs0h) >= 0.05]) # find boundary pval for qval < 0.05
upr.pval
#[1] 0.04813947 # pval should be smaller than this

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
})) # 3 sec

### any single genotype
sum(apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 0))
#[1] 10833
### any two genotypes
sum(apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 1))
#[1] 9272
### any three genotypes
sum(apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 2))
#[1] 8252
# it doesn't decrease quickly. 

### Take any single genotype
genes.dyn = gene.names[apply(sig.geno.by.dyn.gene, 1, function(x) sum(x) > 0)]

#############
##### compare each genotype to r2r1 with same treatment at same time point
geno.p.md.by.genes.var.geno = lapply(glm.coefs, function(x) {
   est = x$est
   sem = x$sem
   est = matrix(est, nrow=length(genotypes))
   sem = matrix(sem, nrow=length(genotypes))
   est.dr2r1 = est[c(1:16,18),] - matrix(rep(est[17,],each=17), nrow=17)
   sem.dr2r1 = sqrt(sem[c(1:16,18),]^2 + matrix(rep(sem[17,],each=17), nrow=17) ^2)
   pval.dr2r1 = 2*pnorm(abs(est.dr2r1)/sem.dr2r1, lower.tail = F)
   return(list(mean.diff=est.dr2r1, p.vals=pval.dr2r1))
}) # 2 sec

### qval cutoff
all.pvals.vsr = t(sapply(geno.p.md.by.genes.var.geno, function(x) {
   as.numeric(x$p.vals)
}))
hist(all.pvals.vsr, freq=F) # looking good
all.qvals.vsr = qvalue(all.pvals.vsr)$qvalues
upr.pval = min(sort(all.pvals.vsr)[sort(all.qvals.vsr) >= 0.05]) # find boundary pval for qval < 0.05
upr.pval
#[1] 0.01185475 # pval should be smaller than this

sig.geno.by.gene.var.geno = t(sapply(geno.p.md.by.genes.var.geno, function(x) {
   p.val.pass = x$p.vals < upr.pval
   md.pass = abs(x$mean.diff) > 1
   geno.p = apply(p.val.pass * md.pass, 1, sum)
})) # 2 sec

### any single genotype
sum(apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 0))
#[1] 12786
### any two genotypes
sum(apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 1))
#[1] 11139
### any three genotypes
sum(apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 2))
#[1] 10104
# it doesn't decrease quickly. 

### Take any single genotype  - this criterion was not used later
genes.var.geno = gene.names[apply(sig.geno.by.gene.var.geno, 1, function(x) sum(x) > 0)]

save(genes.dyn,genes.var.geno, file='./outputs/selected.genes.dynamics.and.geno.variation.RData')
#######

################# 
##### JEPS,xxxx,r2r1,GUS x Ed or mock X time
#### (1) by PCA
### genotype colors: JEPS,xxxx,r2r1,GUS - red, orange, green, blue
### mock pch, 0,2,5hpi: 0, 1, 5
### Ed pch, 0,1,2,3,4,5,6 hpi: 7,12,10,13,8,9,11

#### load data
rm(list=ls())
load('./outputs/Ed-AvrRpt2_DI_210223.RData')
load('./outputs/genes.to.be.removed.RData')

glm.c.slim = glm.coefs[!names(glm.coefs) %in% c(ne.genes, four.genes, r.genes)]
glm.c.slim = sapply(glm.c.slim, function(x) {
  est = x$est
  est.JEPS = est[grep('JEPS',names(est))]
  est.xxxx = est[grep('xxxx',names(est))]
  est.r2r1 = est[grep('r2r1',names(est))]
  est.GUS = est[grep('GUS',names(est))]
  c(est.JEPS, est.xxxx, est.r2r1, est.GUS)
})
dim(glm.c.slim)
#[1]    40 18970
#pch.all = rep(c(0,1,5,7,12,10,13,8,9,11), 4)
text.all = rep(c(paste0('m',c(0,2,5)), paste0('E', 0:6)), 4)
text.all1 = paste0('bold(', text.all, ')')
#col.set = c('red','royalblue','darkolivegreen3','seagreen2')
#col.all = rep(col.set, each = 10)
#col.set1 = brewer.pal(4, 'Set1')
col.set1 = c('red','blue','green','purple')
col.set2 = brewer.pal(10,'Set3')[c(4,5,7,10)]
col.all = as.character(t(matrix(c(rep(col.set2, 3), rep(col.set1, 7)), nrow=4)))

##### Fig 1A PCA
pca.slim = prcomp(glm.c.slim, center=T, scale. =F)
pca.results = pca.slim$x
var.val = (pca.slim$sdev)^2
var.val.p = var.val/sum(var.val)
pdf('./outputs/Fig1A.PCA.JEPS.xxxx.r2r1.GUS.r1.pdf', width=11, height = 6.6)
plot(pca.results[,1], pca.results[,2], col=col.all, cex=1.3,
     xlab=paste0('PC1 (', round(var.val.p[1]*100, digits = 1),' %)'),
     ylab=paste0('PC2 (', round(var.val.p[2]*100, digits = 1),' %)'),
     xaxt='n', yaxt='n', type='n')
text(pca.results[,1], pca.results[,2], parse(text=text.all1), cex=1.2, col=col.all)
axis(side=1, at=seq(-50, 175,25), las=1)
axis(side=2, at=seq(-75, 50,25), las=1)
legend('topright',c('JEPS.Ed','JEPS.mock','xxxx.Ed','xxxx.mock',
                       'r2r1.Ed','r2r1.mock','GUS.Ed','GUS.mock'), 
       inset=0.02,
       ncol = 4, col=as.character(rbind(col.set1, col.set2)),
       pch=rep(c('E','m'),4), pt.cex=1.4)
dev.off()

pdf('./outputs/Supple.Fig1A.PCA23.JEPS.xxxx.r2r1.GUS.r1.pdf', width=11, height = 9)
plot(pca.results[,2], pca.results[,3], col=col.all, cex=1.3,
     xlab=paste0('PC2 (', round(var.val.p[2]*100, digits = 1),' %)'),
     ylab=paste0('PC3 (', round(var.val.p[3]*100, digits = 1),' %)'),
     xaxt='n', yaxt='n', type='n')
text(pca.results[,2], pca.results[,3], parse(text=text.all1), cex=1.2, col=col.all)
axis(side=1, at=seq(-50, 175,25), las=1)
axis(side=2, at=seq(-75, 50,25), las=1)
legend('topleft',c('JEPS.Ed','JEPS.mock','xxxx.Ed','xxxx.mock',
                    'r2r1.Ed','r2r1.mock','GUS.Ed','GUS.mock'), 
       inset=0.02,
       ncol = 4, col=as.character(rbind(col.set1, col.set2)),
       pch=rep(c('E','m'),4), pt.cex=1.4)
dev.off()

pdf('./outputs/Supple.Fig1A.PCA13.JEPS.xxxx.r2r1.GUS.r1.pdf', width=11, height = 5)
plot(pca.results[,1], pca.results[,3], col=col.all, cex=1.3,
     xlab=paste0('PC1 (', round(var.val.p[1]*100, digits = 1),' %)'),
     ylab=paste0('PC3 (', round(var.val.p[3]*100, digits = 1),' %)'),
     xaxt='n', yaxt='n', type='n')
text(pca.results[,1], pca.results[,3], parse(text=text.all1), cex=1.2, col=col.all)
axis(side=1, at=seq(-50, 175,25), las=1)
axis(side=2, at=seq(-75, 50,25), las=1)
legend('topright',c('JEPS.Ed','JEPS.mock','xxxx.Ed','xxxx.mock',
                 'r2r1.Ed','r2r1.mock','GUS.Ed','GUS.mock'), 
       inset=0.02,
       ncol = 4, col=as.character(rbind(col.set1, col.set2)),
       pch=rep(c('E','m'),4), pt.cex=1.4)
dev.off()

save(pca.slim, glm.c.slim, file='./outputs/four.geno.all.pca.RData')


##### leaky AvrRpt2 expression effect? Actually, just the basal level difference across genotypes
#### compare GUS (0h) vs. Col (0h) and r2r1 (0h) vs. Col (0h)
#### collect the data and fix.ef only for 0h data
rm(list=ls())
load('./outputs/exp.one.d.RData')
load('./outputs/genes.to.be.removed.RData')
col.0h = grepl('\\.0h$', fix.ef)
exp.od.0h = exp.one.d[,col.0h]
fix.ef.0h = fix.ef[col.0h]
oc90.0h = offset.cval.90[col.0h]
genot.n = strsplit(as.character(fix.ef.0h), '\\.')
genot.n = sapply(genot.n, '[', 1)
genot.nu = unique(genot.n)
genot.nf = factor(genot.n, levels = genot.nu)

glm.0h.genot = apply(as.matrix(exp.od.0h), 1, function(x) {
  glm.nb.0h = glm.nb(x ~ genot.nf -1 + offset(oc90.0h), link=log)
  return(as.data.frame(summary(glm.nb.0h)$coef))
}) # 2.6 min

#### est and pval calculation for the comparisons
### r2r1 vs. JEPS
glm.0h.r2r1.vs.JEPS = t(sapply(glm.0h.genot, function(x) {
  mean.d = x['genot.nfr2r1','Estimate'] - x['genot.nfJEPS','Estimate']
  semd = sqrt(x['genot.nfr2r1','Std. Error']^2 + x['genot.nfJEPS','Std. Error']^2)
  p.val = 2* pnorm(abs(mean.d)/semd, lower.tail = F)
  return(c(est=mean.d, pval = p.val))
}))

### GUS vs. JEPS
glm.0h.GUS.vs.JEPS = t(sapply(glm.0h.genot, function(x) {
  mean.d = x['genot.nfGUS','Estimate'] - x['genot.nfJEPS','Estimate']
  semd = sqrt(x['genot.nfGUS','Std. Error']^2 + x['genot.nfJEPS','Std. Error']^2)
  p.val = 2* pnorm(abs(mean.d)/semd, lower.tail = F)
  return(c(est=mean.d, pval = p.val))
}))

### r2r1 vs. GUS
glm.0h.r2r1.vs.GUS = t(sapply(glm.0h.genot, function(x) {
  mean.d = x['genot.nfr2r1','Estimate'] - x['genot.nfGUS','Estimate']
  semd = sqrt(x['genot.nfr2r1','Std. Error']^2 + x['genot.nfGUS','Std. Error']^2)
  p.val = 2* pnorm(abs(mean.d)/semd, lower.tail = F)
  return(c(est=mean.d, pval = p.val))
}))

### corrected p-val threshold
all.pvals = c(glm.0h.r2r1.vs.JEPS[,'pval'], glm.0h.GUS.vs.JEPS[,'pval'], glm.0h.r2r1.vs.GUS[,'pval'])
pval.th = min(all.pvals[qvalue(all.pvals)$qvalues > 0.05])
pval.th
#[1] 0.001386738

#### gene selection, q < 0.05, log2FC > 1
### r2r1 vs. JEPS
glm.0h.r2r1.vs.JEPS.sig = 
  glm.0h.r2r1.vs.JEPS[glm.0h.r2r1.vs.JEPS[,'pval'] < pval.th & abs(glm.0h.r2r1.vs.JEPS[,'est']) > 1,]
nrow(glm.0h.r2r1.vs.JEPS.sig)
#[1] 75

### GUS vs. JEPS
glm.0h.GUS.vs.JEPS.sig = 
  glm.0h.GUS.vs.JEPS[glm.0h.GUS.vs.JEPS[,'pval'] < pval.th & abs(glm.0h.GUS.vs.JEPS[,'est']) > 1,]
nrow(glm.0h.GUS.vs.JEPS.sig)
#[1] 11

### r2r1 vs. GUS
glm.0h.r2r1.vs.GUS.sig = 
  glm.0h.r2r1.vs.GUS[glm.0h.r2r1.vs.GUS[,'pval'] < pval.th & abs(glm.0h.r2r1.vs.GUS[,'est']) > 1,]
nrow(glm.0h.r2r1.vs.GUS.sig)
#[1] 151

### union of the genes
r2r1.GUS.JEPS.genes = sort(unique(c(rownames(glm.0h.r2r1.vs.JEPS.sig),
                                    rownames(glm.0h.GUS.vs.JEPS.sig),
                                    rownames(glm.0h.r2r1.vs.GUS.sig))))
## remove ne.genes, r.genes
r2r1.GUS.JEPS.genes = r2r1.GUS.JEPS.genes[!(r2r1.GUS.JEPS.genes %in% c(ne.genes, r.genes))]

length(r2r1.GUS.JEPS.genes)
#[1] 158 genes for Fig 1

r2r1.GUS.JEPS.est = cbind(glm.0h.r2r1.vs.JEPS[,'est'],
                          glm.0h.GUS.vs.JEPS[,'est'],
                          glm.0h.r2r1.vs.GUS[,'est'])[r2r1.GUS.JEPS.genes,]
colnames(r2r1.GUS.JEPS.est) = c('r2r1/JEPS','GUS/JEPS', 'r2r1/GUS')

ams.1 = r2r1.GUS.JEPS.est
### cosine correlation among genes
row_distances_subset = 1 - cosine(t(ams.1))
distances = as.dist(row_distances_subset, diag = FALSE, upper = FALSE)
hclust.genes.0h = hclust( distances, method = "complete" ) 

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = expression(log[2]),legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

bot_anno1a = HeatmapAnnotation(foo = anno_text(colnames(ams.1), rot = 60, 
                                              just = 'right',gp = gpar(fontsize = 7), location = 1),
                              show_annotation_name = F)

ht_1a = Heatmap(ams.1 , col = col, cluster_rows = hclust.genes.0h, cluster_columns = FALSE,
               name = 'r2r1.GUS.JEPS.0h', show_row_dend = F,
               show_row_names = F, show_column_names = F, 
               column_title = NULL,
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               use_raster = F,width = unit(15,'mm'), show_heatmap_legend = F,
               bottom_annotation = bot_anno1a,
               height = unit(70, 'mm'))

jpeg("./outputs/Fig1B.r2r1.GUS.JEPS.comp.jpeg",height = 90, width = 35, units = "mm",res = 300)
draw(ht_1a, annotation_legend_list = lgd)
dev.off()

#### stats of the figure
sum(r2r1.GUS.JEPS.est[,'r2r1/JEPS'] > 0 & r2r1.GUS.JEPS.est[,'r2r1/GUS'] > 0)
#[1] 142  r2r1 > JEPS, GUS
r2r1.up.genes = rownames(r2r1.GUS.JEPS.est)[r2r1.GUS.JEPS.est[,'r2r1/JEPS'] > 0 & 
                                              r2r1.GUS.JEPS.est[,'r2r1/GUS'] > 0]

sum(r2r1.GUS.JEPS.est[,'r2r1/JEPS'] < 0 & r2r1.GUS.JEPS.est[,'r2r1/GUS'] < 0)
#[1] 12   r2r1 < JEPS, GUS
r2r1.down.genes = rownames(r2r1.GUS.JEPS.est)[r2r1.GUS.JEPS.est[,'r2r1/JEPS'] < 0 &
                                                r2r1.GUS.JEPS.est[,'r2r1/GUS'] < 0]

### r2r1/JEPS ~ GUS/JEPS and different, r2r1/GUS similar for weak ETI by leaky expression
## upregulated
sum(r2r1.GUS.JEPS.est[,'r2r1/JEPS'] < -0.1 & r2r1.GUS.JEPS.est[,'GUS/JEPS'] < -0.1)
#[1] 0, note that 2^0.1 ~ 1.07 so it is sufficiently small

## downregulated
sum(r2r1.GUS.JEPS.est[,'r2r1/JEPS'] > 0.1 & r2r1.GUS.JEPS.est[,'GUS/JEPS'] > 0.1)
#[1] 9
r2r1.GUS.JEPS.est[r2r1.GUS.JEPS.est[,'r2r1/JEPS'] > 0.1 & r2r1.GUS.JEPS.est[,'GUS/JEPS'] > 0.1,]
#          r2r1/JEPS   GUS/JEPS   r2r1/GUS
#AT1G54575 1.0073836 0.1008339  0.9065497
#AT1G76780 0.2602259 1.0029939 -0.7427680
#AT3G25770 1.0350507 0.1333504  0.9017003
#AT3G28510 1.1670723 0.4153754  0.7516969
#AT4G22505 0.1053290 1.3367163 -1.2313873
#AT5G10140 0.6457313 1.0332643 -0.3875330  Only this might be but still not that big difference r2r1/JEPS and r2r1/GUS
#AT5G11330 0.1720743 1.1457799 -0.9737056
#AT5G44430 2.3412002 0.1356147  2.2055855
#AT5G65080 1.1579171 0.1423231  1.0155940

# so, no such genes that suggests weak ETI

save(exp.od.0h, genot.nf, r2r1.GUS.JEPS.est, r2r1.up.genes, r2r1.down.genes, glm.0h.genot,
     file='./outputs/r2r1.up.down.genes.RData')

#################
#################
#### Ed/mock transcriptome ratio for gene x genotype/time
#### polynomially modeled mock time course using Ed-GUS 
#### (mock data, 0, 2, 5 hpi to modeled for the entire time)

rm(list=ls())
load('./outputs/Ed-AvrRpt2_DI_210223.RData')
load('./outputs/selected.genes.dynamics.and.geno.variation.RData')

### time points
time.p = 0:6
time.pr = rep(time.p, 3)
time.pr1 = (time.pr-3)/3
time.pr.sq = sqrt(time.pr)
time.pr.sq1 = time.pr.sq / max(time.pr.sq) *2 -1
t1 = seq(0,6,0.02)
t1.1 = (t1-3)/3
t1.sq = sqrt(t1)
t1.sq1 = t1.sq/max(t1.sq) * 2 -1

mock.time.p = c(0,2,5)
mock.time.pr = factor(rep(mock.time.p, each=3*length(genotypes)))

### GUS data for genes.dyn
### make it with dev resid
fef.GUS = as.character(fix.ef[grep('GUS.Ed', fix.ef)])

date()
G.corrected.tcv = lapply(glm.coefs[genes.dyn], function(x) {
  est=x$est
  sem = x$sem
  est1 = est[grep('GUS.Ed', names(est))]
  sem1 = sem[grep('GUS.Ed', names(sem))]
  max.v = max(est1)
  min.v = min(est1)
  v0 = est1[1]
  est1r = rep(est1, 3)
  sem1r = rep(sem1, 3)
  weight1r = 1/sem1r^2
  d.resid = x$d.resid
  d.resid1 = d.resid[grep('GUS.Ed', fix.ef)]
  d.vals = est1r + d.resid1
  
  ### fit polynomial up to 4th order
  model.p = lm(d.vals ~ poly(time.pr.sq1, 4),
               weights=1/sem1r^2)
  model.p = step(model.p, trace = 0)
  
  ### the values for time 0:6 for GUS, Ed
  GUS.Ed.v = fitted(model.p)[1:7]
  names(GUS.Ed.v) = 0:6
  
  ### the associated sem
  GUS.Ed.ci = predict(model.p, interval = 'c')[1:7,]
  GUS.Ed.ci = GUS.Ed.ci[,'upr'] - GUS.Ed.ci[,'fit']
  GUS.Ed.sem = GUS.Ed.ci/ qnorm(0.975)
  names(GUS.Ed.sem) = 0:6
  
  ######## 
  #### adjust for each mock time points 0, 2, 5 hrs
  GUS.Ed.025mean = mean(GUS.Ed.v[c('0','2','5')])
  
  mock.geno.d.resid = d.resid[grep('\\.m\\.', fix.ef)]
  mock.geno.est = rep(est[grep('\\.m\\.', names(est))], 3)
  mock.geno.sem = rep(est[grep('\\.m\\.', names(sem))], 3)
  mock.geno.d.vals = mock.geno.est + mock.geno.d.resid
  geno = factor(rep(genotypes, 9), levels=genotypes)
  
  ### weighted average for each genotype
  mg.lm1 = lm(mock.geno.d.vals ~ geno/mock.time.pr -1, 
              contrasts = list(mock.time.pr = contr.sum),
              weights = 1/mock.geno.sem^2)
  mg.lm = step(mg.lm1, trace=0)
  
  if (length(coef(mg.lm)) == 3*18 ) {
    coef.t.mg.lm = summary(mg.lm)$coef[1:18, ]
  } else {
    if (length(coef(mg.lm)) == 18) {
      coef.t.mg.lm = summary(mg.lm)$coef
    } else {
      coef.t.mg.lm = matrix(rep(summary(mg.lm)$coef, 18), nrow=18, byrow=T)
    }
  } 
  rownames(coef.t.mg.lm) = genotypes
  
  ### simulated mock values for 0:6 time points for each genotype
  mock.geno.sim.v = t(sapply(coef.t.mg.lm[,1], function(x) {
    x - GUS.Ed.025mean + GUS.Ed.v
  }))
  
  ### sem associated with simulated mock values for 0:6 time points for each genotype
  mock.geno.sim.sem = t(sapply(coef.t.mg.lm[,2], function(x) {
    sqrt(x^2 + GUS.Ed.sem^2)
  }))
  
  ### difference from actual 0, 2, 5 hours
  mock.geno.est = est[grep('\\.m\\.', names(est))]
  mock.geno.est = matrix(mock.geno.est, nrow=length(genotypes))
  
  mg.sv025 = mock.geno.sim.v[,c('0','2','5')]
  mg.vsm.geno.est = mock.geno.est - mg.sv025
  
  #### calculate Ed-mock for each genotype at each time point
  Ed.geno.est = est[grep('\\.Ed\\.', names(est))]
  Ed.geno.est = matrix(Ed.geno.est, nrow=length(genotypes))
  rownames(Ed.geno.est) = genotypes
  colnames(Ed.geno.est) = 0:6
  
  Ed.geno.sem = sem[grep('\\.Ed\\.', names(sem))]
  Ed.geno.sem = matrix(Ed.geno.sem, nrow=length(genotypes))
  rownames(Ed.geno.sem) = genotypes
  colnames(Ed.geno.sem) = 0:6
  
  Ed.vsm.geno.est = Ed.geno.est - mock.geno.sim.v
  Ed.vsm.geno.sem = sqrt(Ed.geno.sem^2 + mock.geno.sim.sem^2)
  
  return(list(Ed.vsm = Ed.vsm.geno.est, Ed.vsm.se = Ed.vsm.geno.sem, 
              mock.vsm = mg.vsm.geno.est,
              mock.7tp = mock.geno.sim.v, mock.7tp.se = mock.geno.sim.sem))
  
}) # 3 min
date()

save(G.corrected.tcv, file='./outputs/EdGUS.corrected.mock.v210906.RData')

##################
#### load data
rm(list=ls())
load('./outputs/EdGUS.corrected.mock.v210906.RData')

genotypes = rownames(G.corrected.tcv[[1]]$Ed.vsm)

#### gene selection v210906
#### by q < 0.05, > 3-fold (same direction) at least two consecutive time points, after 2hours, in at least one genotype
pvals.Ed.dyn = lapply(G.corrected.tcv, function(x) {
  pvals = t(sapply(genotypes, function(geno) {
    est = x$Ed.vsm[geno,]
    sem = x$Ed.vsm.se[geno,]
    2* pnorm(abs(est/sem), lower.tail = F)
  }))
  rownames(pvals) = genotypes
  colnames(pvals) = 0:6
  return(pvals)
}) # 3 sec
all.pvals = unlist(pvals.Ed.dyn)
hist(all.pvals, breaks=50) # very good
all.qvals = qvalue(all.pvals)$qvalues
th.pval = min(all.pvals[all.qvals > 0.05])
th.pval
#[1] 0.02116468

zero.mat = matrix(0, nrow=18, ncol=7)
genes.Ed.dyn = sapply(names(G.corrected.tcv), function(x) {
  est = G.corrected.tcv[[x]]$Ed.vsm
  pvals = pvals.Ed.dyn[[x]]
  est.m = zero.mat
  est.m[est >= log2(3)] = 1
  est.m[est <= -log2(3)] = -1
  pvals = pvals < th.pval
  e.p = est.m * pvals
  e.pt = e.p[,1:6] * e.p[,2:7]
  e.pt = e.pt[-c(17,18),4:6] # 3 to 6 hrs
  e.pt.s = apply(e.pt, 1, function(x) sum(x > 0))
  sum(e.pt.s > 0) > 0    ## changed on 230327
})
genes.Ed.dyn = names(G.corrected.tcv)[genes.Ed.dyn]
length(genes.Ed.dyn)
#[1] 3499 genes    #from 3188 w/ change on 230327, i.e., +311

Ed.dyn.dat = t(sapply(G.corrected.tcv[genes.Ed.dyn], function(x) {
  Ed.vals = x$Ed.vsm[genotypes[-(17:18)],]
  Ed.vals = as.numeric(t(Ed.vals))
  return(Ed.vals)
}))

colnames(Ed.dyn.dat) = paste0(rep(genotypes[-(17:18)], ea=7), '.', rep(0:6, 16), 'h')
save(Ed.dyn.dat, file='./outputs/Ed.dyn.dat.v230327.RData')


##################
###################### gamma pdf fit

rm(list=ls())
#### function definition
### gamma pdf defined with shape in log2 (l2a) and peak time (pt), and peak height (ph)
load('./outputs/gamma.pdf.model.RData')

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
load('./outputs/Ed-AvrRpt2_DI_210223.RData')
load('./outputs/Ed.dyn.dat.v230327.RData')
load('./outputs/EdGUS.corrected.mock.v210906.RData')

genes.Ed.dyn = rownames(Ed.dyn.dat)
glm.coefs.gEdyn = glm.coefs[genes.Ed.dyn]
G.c.tcv.dEdyn = G.corrected.tcv[genes.Ed.dyn]
rm(list=c('glm.coefs','G.corrected.tcv'))  # remove large unnecessary lists

#### focus on only 16 genotyps
gEdyn.samp = !grepl('r2r1|GUS|\\.m\\.', as.character(fix.ef))

#### organize the observation and weights values
obs.wei = lapply(genes.Ed.dyn, function(x) {
  d.res = glm.coefs.gEdyn[[x]]$d.resid[gEdyn.samp]
  d.res = matrix(d.res, nrow=16)
  d.res = rbind(d.res[,1:7], d.res[,8:14], d.res[,15:21])
  d.res = as.numeric(t(d.res))
  est = rep(Ed.dyn.dat[x,], 3)
  obs = est + d.res
  sem = rep(as.numeric(t(G.c.tcv.dEdyn[[x]]$Ed.vsm.se[-c(17,18),])), 3)
  wei = 1/sem^2
  return(list(obs=obs, weights=wei))
})
names(obs.wei) = genes.Ed.dyn

#### variables corresponding to the observations and weights.
ti = rep(0:6, 16)
ti.r = rep(0:6, each=16, times=3)
geno1 = rep(genotypes[-(17:18)], each=7)
geno.r = rep(genotypes[-(17:18)], 7*3)
a.indx = rep(1:16, each=7)
ar.indx = rep(1:16, 7*3)

######## visualize early peak genes with 16 genotypes
### ear.p.genes (early.peak genes in the manuscript)
### select the genes that have the maximum response (including repression) earlier than 5 hrs in any genotype
### for the max genotype max-min * f.diff > GUS max-min 
f.diff = 0.1 # the ratio difference in variation
ear.p.genes = sapply(G.c.tcv.dEdyn, function(x) {
  geno.p = apply(x$Ed.vsm[1:16,], 1, function(y) {
    max.t = which.max(abs(y))
    e.peak = max.t < 6 & max.t > 3
    if (sign(max(y) * min(y)) < 0) {
      vari = max(abs(max(y)), abs(min(y)))
    } else {
      vari = abs(max(y)-min(y))
    }
    return(c(e.peak, vari))
  })
  geno.p = apply(geno.p, 2, prod)
  gus.sd = sd(x$Ed.vsm['GUS',])  ## changed in v230327
  sum(geno.p * f.diff > gus.sd )  ## changed in v230327
})
ear.p.genes = genes.Ed.dyn[ear.p.genes > 1]
ear.p.genes = ear.p.genes[-which(ear.p.genes == 'AT2G43820')]
length(ear.p.genes)
#[1] 404 genes   #from 248, w/ change on 230327

##################################
### ear2.p.genes (early.peak2 genes in the manuscript)
### select the genes that have the maximum response (including repression) earlier than 6 hrs in more than 1 genotype
### for the max genotype max-min * f.diff > GUS max-min 
f.diff = 0.1 # the ratio difference in variation
ear2.p.genes = sapply(G.c.tcv.dEdyn, function(x) {
  geno.p = apply(x$Ed.vsm[1:16,], 1, function(y) {
    max.t = which.max(abs(y))
    e.peak = max.t < 7 & max.t > 3
    if (sign(max(y) * min(y)) < 0) {
      vari = max(abs(max(y)), abs(min(y)))
    } else {
      vari = abs(max(y)-min(y))
    }
    return(c(e.peak, vari))
  })
  geno.p = apply(geno.p, 2, prod)
  gus.sd = sd(x$Ed.vsm['GUS',])  ## changed in v230327
  sum(geno.p * f.diff > gus.sd )  ## changed in v230327
})
ear2.p.genes = genes.Ed.dyn[ear2.p.genes > 1]
ear2.p.genes = ear2.p.genes[!(ear2.p.genes %in% ear.p.genes)]
ear2.p.genes = ear2.p.genes[-which(ear2.p.genes == 'AT2G43820')]
length(ear2.p.genes)
#[1] 972 genes    #from 607 genes w/change on 230327

##################### the rest of the genes
### ear3.p.genes (late.peak genes in the manuscript) as the rest of the genes
### select the genes that have the maximum response (including repression) later than 6 hrs in more than 1 genotype
### for the max genotype max-min * f.diff > GUS max-min 
f.diff = 0.15 # the ratio difference in variation
ear3.p.genes = sapply(G.c.tcv.dEdyn, function(x) {
  geno.p = apply(x$Ed.vsm[1:16,], 1, function(y) {
    max.t = which.max(abs(y))
    e.peak = max.t > 3
    if (sign(max(y) * min(y)) < 0) {
      vari = max(abs(max(y)), abs(min(y)))
    } else {
      vari = abs(max(y)-min(y))
    }
    return(c(e.peak, vari))
  })
  geno.p = apply(geno.p, 2, prod)
  gus.sd = sd(x$Ed.vsm['GUS',])  ## changed in v230327
  sum(geno.p * f.diff > gus.sd )  ## changed in v230327
})
ear3.p.genes = genes.Ed.dyn[ear3.p.genes > 1]
ear3.p.genes = ear3.p.genes[!(ear3.p.genes %in% ear.p.genes) & !(ear3.p.genes %in% ear2.p.genes)]
length(ear3.p.genes)
#[1] 1912 genes   #from 1235 genes w/change on 230327

save(ear.p.genes, ear2.p.genes, ear3.p.genes, file='./outputs/genes.for.gam.model.v230327.RData')

################ use actual observations (d.resid)
##### model related function

#### gamma time course model fitting
gamma.m.fit = function(gene, u.pt.lim=8) {
  skipped = c()
  
  ed.t = as.numeric(Ed.dyn.dat[gene,])
  ed.t = matrix(ed.t, nrow=16, byrow = T)
  obs.g = obs.wei[[gene]]$obs
  obs.g = matrix(obs.g, ncol=7, byrow = T)
  weights = obs.wei[[gene]]$weights
  weights = matrix(weights,ncol=7, byrow = T)
  o01= as.numeric(obs.g[,1:2]) 
  w01 = as.numeric(weights[,1:2])
  obs.01m = lm(o01 ~ 1, weights=w01)
  obs.01m = as.numeric(coef(obs.01m))
  ed.t = ed.t - obs.01m
  obs.g = obs.g - obs.01m
  rownames(obs.g) = rownames(weights) = rep(genotypes[1:16],3)
  rownames(ed.t) = genotypes[1:16]
  
  ### estimate starting parameter values
  ## up or down regulation? ud
  ud = sign(as.numeric(ed.t)[which.max(abs(as.numeric(ed.t)))])
  ## pt
  max.t = apply(ed.t, 1, function(x) which.max(ud * x))
  pt.geno = as.integer(max.t) - 1
  names(pt.geno) = names(max.t)
  ## ph
  ph.geno = sapply(names(max.t), function(x) {
    ed.t[x, max.t[x]]
  })
  names(ph.geno) = names(max.t)
  
  ### preliminary fit for each genotype
  t = 0:6
  l.ph.r = 0.7; u.ph.r = 3
  l.l2a = log2(3); u.l2a = 4.5
  l.pt = 3.1; u.pt = u.pt.lim - 1
  l.t0 = 0.4; u.t0 = 1.8
  pre.paras = c()
  pre.paras.se = c()
  pre.paras.df = c()
  for (genot1 in names(max.t)) {
    ed.tg = ed.t[genot1,]
    ph1 = as.numeric(ph.geno[genot1])
    pt1 = as.numeric(pt.geno[genot1])
    l2a1 = 3
    t01 = 1
    if (ud < 0) {
      l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
    } else {
      l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
    }
    pre.gm = nlsLM(ed.tg ~ g.distr.pt0.l2a(t, ph=ph, l2a=l2a, pt=pt, t0=t0),
                   start=c(ph=ph1, l2a=l2a1, pt=pt1, t0=t01),
                   lower=c(ph=l.ph1, l2a=l.l2a, pt=l.pt, t0=l.t0),
                   upper=c(ph=u.ph1, l2a=u.l2a, pt=u.pt, t0=u.t0))
    pre.paras = rbind(pre.paras, coef(pre.gm))
    pre.paras.se = rbind(pre.paras.se, summary(pre.gm)$coef[,2])
    pre.paras.df = c(pre.paras.df, summary(pre.gm)$df[2])
  }
  rownames(pre.paras) = rownames(pre.paras.se) = names(pre.paras.df) = names(max.t)
  
  ## fit to "data"
  l.ph.r = l.l2a.r = l.pt.r = l.t0.r = 0.7 
  u.ph.r = u.l2a.r = u.pt.r = u.t0.r = 1.4
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
    t01 = as.numeric(pre.paras[genot1, 't0'])
    if (ud < 0) {
      l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
    } else {
      l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
    }
    l.l2a = l2a1 * l.l2a.r; u.l2a = l2a1 * u.l2a.r
    if (l.l2a < log2(3)) l.l2a = log2(3)
    if (u.l2a > 4.5) u.l2a = 4.5
    l.pt = pt1 * l.pt.r; u.pt = pt1 * u.pt.r
    if (u.pt > u.pt.lim - 0.2) u.pt=u.pt.lim - 0.2
    l.t0 = t01 * l.t0.r; u.t0 = t01 * u.t0.r
    if (l.t0 < 0.4) l.t0 = 0.4
    if (u.t0 > 1.8) u.t0 = 1.8
    
    pre.gm = tryCatch( nlsLM(obs.tg ~ g.distr.pt0.l2a(tx, ph=ph, l2a=l2a, pt=pt, t0=t0),
                             start=c(ph=ph1, l2a=l2a1, pt=pt1, t0=t01),
                             lower=c(ph=l.ph1, l2a=l.l2a, pt=l.pt, t0=l.t0),
                             upper=c(ph=u.ph1, l2a=u.l2a, pt=u.pt, t0=u.t0),
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
    
    ### scan t0 for the gene
    if (min(pre.paras1[,'t0']) == max(pre.paras1[,'t0'])) {
      t0x = min(pre.paras1[,'t0'])
      t.len = 1
    } else {
      t0x = seq(min(pre.paras1[,'t0']), max(pre.paras1[,'t0']), length.out = 10)
      t.len = 0
    }
    t0x = round(t0x, digits=2)
    
    pre.paras.t0scan = list()
    for (t01 in t0x) {
      l.l2a = log2(3); u.l2a = 4.5
      l.pt.r = 0.8; u.pt.r = 1.2
      
      pre.paras1.t0 = c()
      pre.paras1.t0.se = c()
      pre.paras1.t0.df = c()
      for (genot1 in names(max.t)) {
        genot.rows = which(rownames(obs.g) == genot1)
        obs.tg = as.numeric(obs.g[genot.rows,])
        weights.tg = as.numeric(weights[genot.rows,])
        tx = rep(t, each=3)
        ph1 = as.numeric(pre.paras[genot1, 'ph'])
        pt1 = as.numeric(pre.paras[genot1, 'pt'])
        l2a1 = as.numeric(pre.paras[genot1, 'l2a'])
        if (ud < 0) {
          l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
        } else {
          l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
        }
        l.pt = pt1 * l.pt.r; u.pt = pt1 * u.pt.r
        if (u.pt > u.pt.lim) u.pt = u.pt.lim
        
        obs.tg = obs.tg[tx > 1]; weights.tg = weights.tg[tx > 1] 
        tx = tx[tx > 1]  # only time>=2 for fitting
        
        pre.gm = tryCatch(nlsLM(obs.tg ~ g.distr.pt0.l2a(tx, ph=ph, l2a=l2a, pt=pt, t0=t01),
                                start=c(ph=ph1, l2a=l2a1, pt=pt1),
                                lower=c(ph=l.ph1, l2a=l.l2a, pt=l.pt),
                                upper=c(ph=u.ph1, l2a=u.l2a, pt=u.pt),
                                weights = weights.tg), error = function(e){} )
        if (is.null(pre.gm)) {
          coef.pre.gm = pre.paras1[genot1, ]
          se.pre.gm = pre.paras1.se[genot1, ]
          df.pre.gm = pre.paras1.df[genot1]
          sigma = NA
        } else {
          coef.pre.gm = c(coef(pre.gm), t0=t01)
          se.pre.gm = summary(pre.gm)$coef[,2]
          df.pre.gm = summary(pre.gm)$df[2]
          sigma = sd(resid(pre.gm))
        }
        pre.paras1.t0 = rbind(pre.paras1.t0, c(coef.pre.gm, sigma=sigma))
        pre.paras1.t0.se = rbind(pre.paras1.t0.se, se.pre.gm)
        pre.paras1.t0.df = c(pre.paras1.t0.df, df.pre.gm)
      }
      rownames(pre.paras1.t0) = rownames(pre.paras1.t0.se) = names(pre.paras1.t0.df) = names(max.t)
      pre.paras.t0scan[[as.character(t01)]] = list(coef = pre.paras1.t0, se = pre.paras1.t0.se,
                                                   df = pre.paras1.t0.df)
      
    }
    
    sig.vals = sapply(pre.paras.t0scan, function(x) {
      sig = x$coef[,'sigma']
      sig.all = sqrt(mean(sig^2, na.rm = T))
      return(c(sig.all, max(sig, na.rm = T), sum(is.na(sig))))
    })
    
    if (t.len == 1) {
      pre.paras1.sel = pre.paras.t0scan[[1]]$coef
      pre.paras1.se.sel = pre.paras.t0scan[[1]]$se
      pre.paras1.df.sel = pre.paras.t0scan[[1]]$df
    } else {
      ## first, no NA
      sig.vals.nna = sig.vals[,sig.vals[3,] ==0]
      ## select the one with minimum max(sig)
      sel.n = colnames(sig.vals.nna)[which.min(rank(sig.vals.nna[1,])+rank(sig.vals.nna[2,]))]
      pre.paras1.sel = pre.paras.t0scan[[sel.n]]$coef
      pre.paras1.se.sel = pre.paras.t0scan[[sel.n]]$se
      pre.paras1.df.sel = pre.paras.t0scan[[sel.n]]$df
    }
    
    return(list(para.vals = pre.paras1.sel,
                se.vals = pre.paras1.se.sel,
                df.vals = pre.paras1.df.sel,
                baseline=obs.01m))
  } else {
    return(NA)
  }
}



#############
######## MAIN

###########
#### fit the model for genes with peak time < 5, ear.p.genes
gm.fit = list()
u.pt.lim = 8

for (x in ear.p.genes) {
  if (x == "AT3G03190"|x == "AT3G14870"|x == "AT3G62950"|x == "ATCG00190") next  # these genes are not good, manually removed
        ## see "manually.removed.genes.v230327.r" and its output "manually.removed.genes.v230327.pdf"
  gm.fit[[x]] = gamma.m.fit(gene=x, u.pt.lim = u.pt.lim)
} # 80 sec

sum(sapply(gm.fit, is.na))
#[1] 0  # models were fit for all genes

##########
#### fit the model for genes with peak time 5 to 6, ear2.p.genes
gm.fit2 = list()
u.pt.lim = 8

date()
for (x in ear2.p.genes) {
  if (x == "AT1G21160"|x == "AT2G14610"|x == "AT4G23130"|x == "AT5G02320"|x == "AT5G10140"|x == "AT5G41761"|x == "ATCG01120") next  # these genes are not good, manually removed
  ## see "manually.removed.genes.v230327.r" and its output "manually.removed.genes.v230327.pdf"
  gm.fit2[[x]] = gamma.m.fit(gene=x, u.pt.lim = u.pt.lim)
} # 3 min
date()

sum(sapply(gm.fit2, is.na))
#[1] 0  # models were fit for all genes

##############
#### fit the model for genes with peak time > 6, ear3.p.genes
#### since fitting is difficult for the late time genes, more constrained
#### special model fitting algorithm was used
### since peak time, level inference is not reliable, this set is not used for parameter analysis
### but actual expression value analysis
t0x = seq(1.2, 2.5, length.out = 10)

gm.fit3 = list()
skipped = c()
u.pt.lim = 8

date()
for (gene in ear3.p.genes) {
  ed.t = as.numeric(Ed.dyn.dat[gene,])
  ed.t = matrix(ed.t, nrow=16, byrow = T)
  obs.g = obs.wei[[gene]]$obs
  obs.g = matrix(obs.g, ncol=7, byrow = T)
  weights = obs.wei[[gene]]$weights
  weights = matrix(weights,ncol=7, byrow = T)
  o01= as.numeric(obs.g[,1:2]) 
  w01 = as.numeric(weights[,1:2])
  obs.01m = lm(o01 ~ 1, weights=w01)
  obs.01m = as.numeric(coef(obs.01m))
  ed.t = ed.t - obs.01m
  obs.g = obs.g - obs.01m
  rownames(obs.g) = rownames(weights) = rep(genotypes[1:16],3)
  rownames(ed.t) = genotypes[1:16]
  
  ### estimate starting parameter values
  ## up or down regulation? ud
  ud = sign(as.numeric(ed.t)[which.max(abs(as.numeric(ed.t)))])
  ## pt
  max.t = apply(ed.t, 1, function(x) which.max(ud * x))
  pt.geno = as.integer(max.t) - 1
  names(pt.geno) = names(max.t)
  ## ph
  ph.geno = sapply(names(max.t), function(x) {
    ed.t[x, max.t[x]]
  })
  names(ph.geno) = names(max.t)
  
  ### preliminary fit for each genotype with fixed pt, t0 (pt1 = 7.5, t01=2)
  t = 0:6
  l.ph.r = 0.7; u.ph.r = 3
  l.l2a = log2(3); u.l2a = 4.5
  pt1 = 7.5
  t01 = 2
  pre.paras = c()
  for (genot1 in names(max.t)) {
    ed.tg = ed.t[genot1,]
    ph1 = as.numeric(ph.geno[genot1])
    l2a1 = 3
    if (ud < 0) {
      l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
    } else {
      l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
    }
    pre.gm = nlsLM(ed.tg ~ g.distr.pt0.l2a(t, ph=ph, l2a=l2a, pt=pt1, t0=t01),
                   start=c(ph=ph1, l2a=l2a1),
                   lower=c(ph=l.ph1, l2a=l.l2a),
                   upper=c(ph=u.ph1, l2a=u.l2a))
    pre.paras = rbind(pre.paras, coef(pre.gm))
  }
  rownames(pre.paras) = names(max.t)
  
  pre.paras.t0scan = list()
  for (t01 in t0x) {
    
    ## fit to "data", 1. fix ph
    l.ph.r = l.l2a.r = l.pt.r = 0.7 
    u.ph.r = u.l2a.r = u.pt.r = 1.4
    tx = rep(t, each=3)
    pt1 = 7.5
    
    pre.paras1 = c()
    skipped.genot = c()
    for (genot1 in names(max.t)) {
      genot.rows = which(rownames(obs.g) == genot1)
      obs.tg = as.numeric(obs.g[genot.rows,])
      weights.tg = as.numeric(weights[genot.rows,])
      ph1 = as.numeric(pre.paras[genot1, 'ph'])
      l2a1 = as.numeric(pre.paras[genot1, 'l2a'])
      if (ud < 0) {
        l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
      } else {
        l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
      }
      l.l2a = l2a1 * l.l2a.r; u.l2a = l2a1 * u.l2a.r
      if (l.l2a < log2(3)) l.l2a = log2(3)
      if (u.l2a > 4.5) u.l2a = 4.5
      l.pt = pt1 * l.pt.r; u.pt = u.pt.lim - 0.2
      
      pre.gm = tryCatch( nlsLM(obs.tg ~ g.distr.pt0.l2a(tx, ph=ph1, l2a=l2a, pt=pt, t0=t01),
                               start=c(l2a=l2a1, pt=pt1),
                               lower=c(l2a=l.l2a, pt=l.pt),
                               upper=c(l2a=u.l2a, pt=u.pt),
                               weights = weights.tg), error = function(e){} )
      if (is.null(pre.gm)) {
        coef.pre.gm = c(l2a1,pt1)
        skipped.genot = c(skipped.genot, genot1)
        sigma = NA
      } else {
        coef.pre.gm = coef(pre.gm)
        sigma = sd(resid(pre.gm))
      }
      pre.paras1 = rbind(pre.paras1, c(coef.pre.gm, sigma=sigma))
    }
    rownames(pre.paras1) = names(max.t)
    pre.paras1 = cbind(pre.paras1, ph=pre.paras[,'ph'])
    
    ## fit to "data", 2. fix l2a
    pre.paras2 = c()
    skipped.genot = c()
    for (genot1 in names(max.t)) {
      genot.rows = which(rownames(obs.g) == genot1)
      obs.tg = as.numeric(obs.g[genot.rows,])
      weights.tg = as.numeric(weights[genot.rows,])
      ph1 = as.numeric(pre.paras1[genot1, 'ph'])
      pt1 = as.numeric(pre.paras1[genot1, 'pt'])
      l2a1 = as.numeric(pre.paras1[genot1, 'l2a'])
      if (ud < 0) {
        l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
      } else {
        l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
      }
      l.l2a = l2a1 * l.l2a.r; u.l2a = l2a1 * u.l2a.r
      if (l.l2a < log2(3)) l.l2a = log2(3)
      if (u.l2a > 4.5) u.l2a = 4.5
      l.pt = pt1 * l.pt.r; u.pt = u.pt.lim - 0.1
      
      pre.gm = tryCatch( nlsLM(obs.tg ~ g.distr.pt0.l2a(tx, ph=ph, l2a=l2a1, pt=pt, t0=t01),
                               start=c(ph=ph1, pt=pt1),
                               lower=c(ph=l.ph1, pt=l.pt),
                               upper=c(ph=u.ph1, pt=u.pt),
                               weights = weights.tg), error = function(e){} )
      if (is.null(pre.gm)) {
        coef.pre.gm = c(ph1,pt1)
        skipped.genot = c(skipped.genot, genot1)
        sigma = NA
      } else {
        coef.pre.gm = coef(pre.gm)
        sigma = sd(resid(pre.gm))
      }
      pre.paras2 = rbind(pre.paras2, c(coef.pre.gm, sigma=sigma))
    }
    rownames(pre.paras2) = names(max.t)
    pre.paras2 = cbind(pre.paras2, l2a=pre.paras1[,'l2a'])
    
    if (length(skipped.genot) > 3) {
      skipped = c(skipped, gene)
      next
    }
    
    ## fit to "data", 3. fix pt
    pre.paras3 = c()
    skipped.genot = c()
    for (genot1 in names(max.t)) {
      genot.rows = which(rownames(obs.g) == genot1)
      obs.tg = as.numeric(obs.g[genot.rows,])
      weights.tg = as.numeric(weights[genot.rows,])
      ph1 = as.numeric(pre.paras2[genot1, 'ph'])
      pt1 = as.numeric(pre.paras2[genot1, 'pt'])
      l2a1 = as.numeric(pre.paras2[genot1, 'l2a'])
      if (ud < 0) {
        l.ph1 = ph1 * u.ph.r; u.ph1 = ph1 * l.ph.r
      } else {
        l.ph1 = ph1 * l.ph.r; u.ph1 = ph1 * u.ph.r
      }
      l.l2a = l2a1 * l.l2a.r; u.l2a = l2a1 * u.l2a.r
      if (l.l2a < log2(3)) l.l2a = log2(3)
      if (u.l2a > 4.5) u.l2a = 4.5
      l.pt = pt1 * l.pt.r; u.pt = u.pt.lim
      
      pre.gm = tryCatch( nlsLM(obs.tg ~ g.distr.pt0.l2a(tx, ph=ph, l2a=l2a, pt=pt1, t0=t01),
                               start=c(ph=ph1, l2a=l2a1),
                               lower=c(ph=l.ph1, l2a=l.l2a),
                               upper=c(ph=u.ph1, l2a=u.l2a),
                               weights = weights.tg), error = function(e){} )
      if (is.null(pre.gm)) {
        coef.pre.gm = c(ph1,l2a1)
        skipped.genot = c(skipped.genot, genot1)
        sigma = NA
      } else {
        coef.pre.gm = coef(pre.gm)
        sigma = sd(resid(pre.gm))
      }
      pre.paras3 = rbind(pre.paras3, c(coef.pre.gm, sigma=sigma))
    }
    rownames(pre.paras3) = names(max.t)
    pre.paras3 = cbind(pre.paras3, pt=pre.paras2[,'pt'], t0=t01)
    
    if (length(skipped.genot) > 3) {
      skipped = c(skipped, gene)
      next
    }
    
    pre.paras.t0scan[[as.character(t01)]] = pre.paras3
    
  }
  
  sig.vals = sapply(pre.paras.t0scan, function(x) {
    sig = x[,'sigma']
    sig.all = sqrt(mean(sig^2, na.rm = T))
    return(c(sig.all, max(sig, na.rm = T), sum(is.na(sig))))
  })
  
  ## first, no NA
  sig.vals.nna = sig.vals[,sig.vals[3,] ==0]
  ## select the one with minimum max(sig)
  sel.n = colnames(sig.vals.nna)[which.min(rank(sig.vals.nna[1,])+rank(sig.vals.nna[2,]))]
  pre.paras1.sel = pre.paras.t0scan[[sel.n]]
  
  gm.fit3[[gene]] = list(para.vals = pre.paras1.sel,
                         baseline=obs.01m)
}  #)  ~ 11 min
date()

length(ear3.p.genes); length(gm.fit3)
#[1] 1912
#[1] 1912
# models were fit for all genes

save(gm.fit, gm.fit2, gm.fit3, file='./outputs/gm.fit.v230327.RData')

#################
#### order the genes according to pt, and split by ph sign
rm(list=ls())
load('./outputs/gm.fit.v230327.RData')
load('./outputs/genes.to.be.removed.RData')
gm.fit.all = c(gm.fit,gm.fit2,gm.fit3)
length(gm.fit.all)
#[1] 3278 genes (v230327), I should use this set
gfa.names = names(gm.fit.all)
## remove mutation target genes, transgenes
gfa.names = gfa.names[!gfa.names %in% c(four.genes, ne.genes)]
gfa.names = gfa.names[!grepl('^ATCG', gfa.names)]
length(gfa.names)
#[1] 3262  ## the four genes, transgenes, and ATCG genes are already removed!!!

# order them according to the pt and ph(+/-) in wt
pt.vals = sapply(gm.fit.all[gfa.names], function(x) x$para.vals['JEPS', 'pt'])
ph.vals = sapply(gm.fit.all[gfa.names], function(x) x$para.vals['JEPS', 'ph'])
gfa.names1 = gfa.names[order(pt.vals)]
ph.vals = ph.vals[gfa.names1]
gfa.names1.posi = gfa.names1[ph.vals > 0]
gfa.names1.nega = gfa.names1[ph.vals < 0]
length(gfa.names1.posi); length(gfa.names1.nega)
#[1] 1972
#[1] 1290

save(gfa.names1.posi, gfa.names1.nega, file='./outputs/ETI.gfa.sorted.genes.v230327.RData')
# this workspace was used in the earlier part of the script.

########################
#### Fig 3, fold change gene x genotype/time

### load data
rm(list=ls())
load('./outputs/Ed.dyn.dat.v230327.RData') # corrected mock subtracted
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered, this is recursive

ams.1 = Ed.dyn.dat[c(gfa.names1.posi, gfa.names1.nega),]

diff.wq  = ams.1[,1:7] - ams.1[,106:112]
colnames(diff.wq) = paste('WT-quad_', 0:6, 'h')
ams.1 = cbind(ams.1, diff.wq)

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))

cnames = colnames(Ed.dyn.dat)
cnames = cnames[1:16*7]
cnames = strsplit(cnames, '\\.')
genotypes = unlist(sapply(cnames, '[', 1))
cnames = c(genotypes, 'JEPS/xxxx')
col.anno = rbind(matrix('',ncol=17,nrow=3), cnames, matrix('',ncol=17,nrow=3))
col.anno = as.character(col.anno)

bot_anno2 = HeatmapAnnotation(foo = anno_text(col.anno, just = 'right', rot=60,
                                              gp = gpar(fontsize = 7), location = 1),
                              show_annotation_name = F)
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

ht_2 = Heatmap(ams.1 , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'Ed.vsm', show_row_dend = F,
               column_gap = unit(4,"mm"), 
               column_split = c(rep('A',112), rep('B',7)),
               row_title_gp = gpar(fontsize = 20),
               show_row_names = F, show_column_names = F, 
               column_title = NULL,
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               use_raster = F,width = unit(120,'mm'), show_heatmap_legend = F,
               bottom_annotation = bot_anno2,
               height = unit(120, 'mm'))

jpeg("./outputs/Fig3.corrected.Ed.vsm.v230327.jpeg",height = 160, width = 150, units = "mm",res = 300)
draw(ht_2, annotation_legend_list = lgd)
decorate_heatmap_body('Ed.vsm', column_slice = 1, {
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
})
grid.text('A', x = unit(5,'mm'),y = unit(138,'mm'),gp = gpar(fontface = 'bold'))
grid.text('B', x = unit(117,'mm'),y = unit(138,'mm'),gp = gpar(fontface = 'bold'))
dev.off()

##########
##### Figs 3C and 3D. PC1 coordinates for combinatorial genotype time points

#### load data
rm(list=ls())
load('./outputs/Ed.dyn.dat.v230327.RData') # corrected mock subtracted
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered, this is recursive

#### the data used in Fig 2A
sel.dyn.dat = Ed.dyn.dat[c(gfa.names1.posi, gfa.names1.nega),]

c.names = colnames(sel.dyn.dat)
c.names = strsplit(c.names, '\\.')
c.names = sapply(c.names, '[', 1)
genotypes = unique(c.names)

pch.all = rep(c(7,12,10,13,8,9,11), 16)
col.all = rep(rainbow(18, end=0.65)[1:16], each = 7)
text.all = as.character(rep(0:6, 16))
text.all1 = paste0('bold(', text.all, ')')

##### PCA
pca.dyn = prcomp(t(sel.dyn.dat), center=T, scale. =F)
pca.results = pca.dyn$x
var.val = (pca.dyn$sdev)^2
var.val.p = var.val/sum(var.val)
pdf('./outputs/Fig3C.PCA.16geno.res.tc.v230327.pdf', width=10, height = 5.5)
plot(pca.results[,1], pca.results[,2], col=col.all, pch=pch.all, cex=1.3,
     xlab=paste0('PC1 (', round(var.val.p[1]*100, digits = 1),' %)'),
     ylab=paste0('PC2 (', round(var.val.p[2]*100, digits = 1),' %)'),
     xaxt='n', yaxt='n', type='n')
axis(side=1, at=seq(-75, 100,25), las=1)
axis(side=2, at=seq(-25, 25,25), las=1)
text(pca.results[,1], pca.results[,2], parse(text=text.all1), cex=1.2, col=col.all)
legend(-25,-12,genotypes, ncol = 4, col=rainbow(18, end=0.65)[1:16],
       pch=15)
dev.off()

#### PC1 values
pc1.geno.time = pca.results[,1]
pc1.geno.time = matrix(pc1.geno.time, nrow=7)
pdf('./outputs/Fig3D.PC1.16geno.v230327.pdf', height=5.5, width=6)
matplot(0:6, pc1.geno.time, type='n', 
        xlab='Time (hpi)', ylab='PC1 value', las=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
matlines(0:6, pc1.geno.time, type='l', col=rainbow(18, end=0.65)[1:16], lty=c(1,2,6), lwd=2)
legend('topleft', genotypes, col=rainbow(18, end=0.65)[1:16], lty=c(1,2,6), ncol=2,
       bg='gray', lwd=2)
dev.off()

save(pca.dyn, pc1.geno.time, file='./outputs/pca.dyn.comb.geno.RData')

##################
#### Fig 2. AvrRpt2 mRNA value
#### load workspace
rm(list=ls())
load('./outputs/Ed-AvrRpt2_DI_210223.RData')
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # 3261 modeled gene sets
load('./outputs/pca.dyn.comb.geno.RData')

all.g.est = t(sapply(glm.coefs, function(x) x$est))
dim(all.g.est)
#[1] 18978   180
boxplot(all.g.est, cex=0)
# saved as "all.log2.exp.val.distr.pdf"
# no outliers shown, very consistent across samples
boxplot(as.numeric(all.g.est), cex=0)
# saved as "all.log2.exp.val.together.distr.pdf"
# no outliers shown

### AvrRpt2 expression levels, plot time course for all genotypes
rownames(all.g.est)[18977]
#[1] "ATtransgene_Inducible"
AvrRpt2.exp = all.g.est["ATtransgene_Inducible",]
AvrRpt2.exp = t(matrix(AvrRpt2.exp, nrow=18))
colnames(AvrRpt2.exp) = genotypes
AvrRpt2.exp.Ed = AvrRpt2.exp[4:10,]
pdf('./outputs/Fig2A.AvrRpt2.mRNAlevel.pdf', height = 7, width = 8)
matplot(0:6, AvrRpt2.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
        xlab='Time (hpi)', ylab=expression('AvrRpt2 mRNA level (log'[2]*')'),
        ylim=c(3, 14), xlim=c(0,7), xaxt='n', lwd=2, las=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
matlines(0:6, AvrRpt2.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=2)
axis(1, at=0:6)
boxplot(all.g.est[,1], add=T, at=6.7, col='gray95',lwd=2, cex=0, yaxt='n')
legend('topleft', genotypes, col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
       lwd=2, ncol=6, cex=0.9, bg='gray80')
text(0, 11.7, 'AvrRpt2', cex=1.6, pos=4)
dev.off()

### among the ETI-upregulated genes, take the max time with JEPS
### correlation with AvrRpt2 across the 16 genotypes
geno.v.at.wt.max = t(sapply(glm.coefs[c("ATtransgene_Inducible",gfa.names1.posi)], function(x) {
  est.tab = matrix(x$est, nrow=18)
  est.tab = est.tab[,-(1:3)]
  rownames(est.tab) = genotypes
  val.ac.geno = est.tab[1:16, which.max(est.tab['JEPS',])]
}))
### center by the mean for each gene
geno.rel.v.at.wt.max = t(apply(geno.v.at.wt.max, 1, function(x) {
  x - mean(x)
}))

#### Fig 2B. correlation with AvrRpt2
cor.w.AvrRpt2 = apply(geno.v.at.wt.max[-1,], 1, function(x) {
  cor(x, geno.v.at.wt.max[1,])
})
sum(cor.w.AvrRpt2 > 0)/length(cor.w.AvrRpt2)
#[1] 0.1272819    13% of the upregulated genes with positive correlation
sum(cor.w.AvrRpt2 > 0.5)/length(cor.w.AvrRpt2)
#[1] 0.01470588    1% of the upregulated genes with correlation > 0.05

#### PC1 values, 4 and 6 hpi
pca.results = pca.dyn$x
pc1.geno.time = pca.results[,1]
pc1.geno.time = matrix(pc1.geno.time, nrow=7)
pc1.gt.4h = pc1.geno.time[5,]
pc1.gt.6h = pc1.geno.time[7,]
cor.pc1.gt.4h = cor(pc1.gt.4h, geno.v.at.wt.max[1,])
cor.pc1.gt.6h = cor(pc1.gt.6h, geno.v.at.wt.max[1,])
cor.pc1.gt.4h; cor.pc1.gt.6h
#[1] -0.5384219
#[1] -0.5639828

pdf('./outputs/Fig2B.correlation.wAvrRpt2.cc.v233037.pdf', height = 5, width = 5)
plot(density(cor.w.AvrRpt2),
     xlab='Correlation with AvrRpt2 mRNA level across the genotypes',
     main='', las=1, xlim=c(-1,1))
abline(v=0, col='gray15', lty=2)
#abline(v=cor.pc1.gt.4h, col='gray75', lty=5)
#abline(v=cor.pc1.gt.6h, col='gray65', lty=5)
#text(-0.52, 0.2, 'ETI-PC1, 4hpi', col='gray75', pos=4, srt=90)
#text(-0.71, 0.2, 'ETI-PC1, 6hpi', col='gray65', pos=4, srt=90)

dev.off()

#### Fig 2C xve mRNA level time course
xve.exp = all.g.est["ATtransgene_XVE_fusion_UTR_and_terminator",]
xve.exp = t(matrix(xve.exp, nrow=18))
colnames(xve.exp) = genotypes
xve.exp.Ed = xve.exp[4:10,]
pdf('./outputs/Fig2C.xve.mRNAlevel.pdf', height = 7, width = 7)
matplot(0:6, xve.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
        xlab='Time (hpi)', ylab=expression('XVE mRNA level (log'[2]*')'),
        ylim=c(10, 15), xlim=c(0,6), lwd=2, las=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
matlines(0:6, xve.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=2)
legend('topleft', genotypes, col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
       lwd=2, ncol=6, cex=0.9, bg='gray80')
text(0, 14, 'XVE', cex=1.6, pos=4)
dev.off()

#### correlation between XVE level at 0hpi and the AvrRpt2 induction rate from 1 to 4 hpi
xve.mean = apply(xve.exp.Ed, 2, mean)
avrrpt2.4h.1h.rat = AvrRpt2.exp.Ed[5,]-AvrRpt2.exp.Ed[2,]
cor.t = cor.test(xve.mean, avrrpt2.4h.1h.rat)
cor.t$estimate
#      cor 
#0.7319902
cor.t$p.value
#[1] 0.0005535221

pdf('./outputs/Fig2D.corr.xve.avrrpt2rate.v230327.pdf', height = 5, width = 5)
plot(0,0,
     xlim=c(11, 12.1), ylim = c(2, 5.5),
     xlab=expression('Mean XVE mRNA level (log'[2]*')'),
     ylab=expression('AvrRpt2 mRNA accumulation rate (log'[2]*')'),
     type='n'
)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray75")
abline(lm(avrrpt2.4h.1h.rat ~ xve.mean), col='gray30')
text(avrrpt2.4h.1h.rat ~ xve.mean, labels=names(xve.mean),
     col = rainbow(18, end=0.65), cex=0.6)
text(12.1,2.45, paste0('Pearson Correlation = ', round(cor.t$estimate, digits = 2)), 
     pos=2, cex=0.8)
text(12.1,2.2, paste0('p-value = ', round(cor.t$p.value, digits = 5)), 
     pos=2, cex=0.8)
dev.off()

#### Figs S2, S3
#### visualization of gamma-pdf model fitting data and model 

#### load data
rm(list=ls())
load('./outputs/gm.fit.v230327.RData')
load('./outputs/EdGUS.corrected.mock.v210906.RData')
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered

##### function definitions
#### gamma-pdf model
load('./outputs/gamma.pdf.model.RData')

#### data and model time course plots for all genotypes per gene
geno.g.tc.plot = function(gene, tx1=0:6, ud='up') {
  ### model values
  para.v = all.gm.fit[[gene]]$para.vals
  t01 = all.gm.fit[[gene]]$para.vals[1,'t0']
  m.v = c()
  for (genot1 in rownames(para.v)) {
    para.geno = para.v[genot1,]
    mg.v = g.distr.pt0.l2a(tx1, ph=para.geno['ph'], l2a=para.geno['l2a'], pt=para.geno['pt'],
                           t0 = t01)
    m.v = cbind(m.v, mg.v)
  }
  
  ### plot y-range
  y.min1 = min(m.v)
  y.max1 = max(m.v)
  
  g.exp.Ed = t(G.corrected.tcv[[gene]]$Ed.vsm)
  y.min2 = min(min(g.exp.Ed), y.min1)
  y.max2 = max(max(g.exp.Ed), y.max1)
  
  ### data plot
  matplot(0:6, g.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
          xlab=NA, ylab=NA,
          ylim=c(y.min2, y.max2), xlim=c(0,6), xaxt='n')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  matlines(0:6, g.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=1)
  axis(1, at=0:6)
  if (ud == 'up') {
    text(0,y.max2-0.1*(y.max2-y.min2), paste(gene,'data'), pos=4, cex=0.8)
  } else {
    text(0,y.min2+0.1*(y.max2-y.min2), paste(gene,'data'), pos=4, cex=0.8)
  }
  
  ### model plot
  matplot(tx1, m.v, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
          xlab=NA, ylab=NA,
          ylim=c(y.min2, y.max2), xlim=c(0,6), xaxt='n')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  matlines(tx1, m.v, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=1)
  axis(1, at=0:6)
  if (ud == 'up') {
    text(0,y.max2-0.1*(y.max2-y.min2), paste(gene,'model'), pos=4, cex=0.8)
  } else {
    text(0,y.min2+0.1*(y.max2-y.min2), paste(gene,'model'), pos=4, cex=0.8)
  }
}

###### MAIN
#### compile and order the models
all.gm.fit = c(gm.fit, gm.fit2, gm.fit3)
all.gm.fit = all.gm.fit[c(gfa.names1.posi, gfa.names1.nega)]
length(all.gm.fit)
#[1] 3262

#### plots
tx1 = seq(0,6,0.05)
pdf('./outputs/FigS3.data.gammamodel.v230327.pdf', height = 10, width = 7.5)
opar=par(mfrow=c(6,4), mar=c(2,2,0.5,0.5))
for (gene in gfa.names1.posi) {
  geno.g.tc.plot(gene, tx1, ud = 'up')
}  
for (gene in gfa.names1.nega) {
  geno.g.tc.plot(gene, tx1, ud = 'down')
}  # ~15 sec total for both loops
par(opar)
dev.off()

#### by heatmap
### model values at 0:6 hpi
tx1 = 0:6
model.v.0.6 = t(sapply(all.gm.fit[c(gfa.names1.posi, gfa.names1.nega)], function(x) {
  para.v = x$para.vals
  t01 = x$para.vals[1,'t0']
  m.v = c()
  for (genot1 in rownames(para.v)) {
    para.geno = para.v[genot1,]
    mg.v = g.distr.pt0.l2a(tx1, ph=para.geno['ph'], l2a=para.geno['l2a'], pt=para.geno['pt'],
                           t0 = t01)
    m.v = cbind(m.v, mg.v)
  }
  as.numeric(m.v)
}))
### data
data.v.0.6 = t(sapply(G.corrected.tcv[c(gfa.names1.posi, gfa.names1.nega)], function(x) t(x$Ed.vsm[1:16,])))
### discrepancy
discr.v.0.6 = data.v.0.6 - model.v.0.6

### heatmap
geno.n = rownames(all.gm.fit[[1]]$para.vals)
ams.1 = cbind(data.v.0.6, model.v.0.6, discr.v.0.6)

col.anno1 = rbind(matrix('',ncol=48,nrow=3), geno.n, matrix('',ncol=48,nrow=3))
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
                column_split = c(rep('_Data',112), rep('_Model', 112), rep('Data - Model', 112)),
                row_title_gp = gpar(fontsize = 20),
                show_row_names = F, show_column_names = F, 
                column_title_gp = gpar(fontsize=10), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(240,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno4x,
                height = unit(120, 'mm'))

jpeg("./outputs/FigS2.model.fit2.v230327.jpeg",height = 160, width = 270, units = "mm",res = 300)
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
##### TableS1 all parameter values and fitted values
### fitted value matrix model.v.0.6
genotypes = rownames(all.gm.fit[[1]]$para.vals)
gt = rep(genotypes, each=7)
ti = rep(0:6, 16)
ti = paste0(ti, 'h')
cna = paste(gt,ti, sep="_")
colnames(model.v.0.6) = cna

### model parameter values
par.v = t(sapply(all.gm.fit[rownames(model.v.0.6)], function(x) {
  p.mat = x$para.vals
  p.line = as.numeric(t(p.mat[,1:3]))
  p.line = c(p.line, p.mat[1,4], x$baseline)
}))
gt = rep(genotypes, each=3)
ti = rep(c('A','log2s','pt'), 16)
cna = paste(gt,ti, sep="_")
cna = c(cna, 't0','B')
colnames(par.v) = cna

ts.tab = cbind(par.v, model.v.0.6)
write.csv(ts.tab, file='./outputs/TableS1.para.vals.fitted.vals.csv')

###################
####### NRAM of the gamma-pdf model parameters

##### load data
rm(list=ls())
load('./outputs/gm.fit.v230327.RData')
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered

##### function definitions
load('./outputs/NRAM.algorithm.w.mean.se.RData')

##### MAIN
#### compile and order the models
all.gm.fit = c(gm.fit, gm.fit2, gm.fit3)
all.gm.fit = all.gm.fit[c(gfa.names1.posi, gfa.names1.nega)]

## select genes that have pt >= 7.5 at most 3 genotypes within gm.fit, gm.fit2
gm.fit.early = c(gm.fit, gm.fit2)
well.modeled = sapply(gm.fit.early, function(x) {
  pts = x$para.vals[,'pt']
  (sort(pts, decreasing = T))[4] < 7.5
})
sum(well.modeled)
#[1] 1102 well-modeled genes #from 726 w/change on 230327
gm.fit.early = gm.fit.early[well.modeled]

## sort according to gfa.names1
gm.fit.early = gm.fit.early[c(gfa.names1.posi, gfa.names1.nega)[
  c(gfa.names1.posi, gfa.names1.nega) %in% names(gm.fit.early)]]
length(gm.fit.early)
#[1] 1100 well-modeled genes, meaning well-defined early peaks

###### NRAM
library(R.utils)
###### NRAM for ph
ph.nram = t(sapply(gm.fit.early, function(x) {
  ph.est = x$para.vals[,'ph']
  ph.se = x$se.vals[,'ph']
  names(ph.est) = names(ph.se) = geno.names[names(ph.est)]
  ave.m.out = run.ave.m(ph.est, ph.se)
}))
sum(apply(ph.nram, 1, function(x) sum(x[-16] !=0)) > 0)
#[1] 928 genes out of 1100 genes have non-0

##### NRAM for pt
pt.nram = t(sapply(gm.fit.early, function(x) {
  pt.est = x$para.vals[,'pt']
  pt.se = x$se.vals[,'pt']
  names(pt.est) = names(pt.se) = geno.names[names(pt.est)]
  ave.m.out = run.ave.m(pt.est, pt.se)
}))
sum(apply(pt.nram, 1, function(x) sum(x[-16] !=0)) > 0)
#[1] 398 genes out of 1100 genes have non-0

detach('package:R.utils', unload=T)

save(pt.nram, ph.nram, file='./outputs/gammamodel.para.NRAM.v230327.RData')

#####################
##### Fig 4 heatmaps of ph, pt NRAM

### load data
rm(list=ls())
load('./outputs/Ed.dyn.dat.v230327.RData') # corrected mock subtracted
load('./outputs/gm.fit.v230327.RData')
load('./outputs/gammamodel.para.NRAM.v230327.RData')
load('./outputs/pca.dyn.comb.geno.RData')
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered, this is recursive

##### JEPS and xxxx time courses for references
quad = Ed.dyn.dat[c(gfa.names1.posi, gfa.names1.nega), 106:112]
wt  = Ed.dyn.dat[c(gfa.names1.posi, gfa.names1.nega), 1:7]
wt.quad = wt - quad

#### compile and order the models
all.gm.fit = c(gm.fit, gm.fit2, gm.fit3)
all.gm.fit = all.gm.fit[c(gfa.names1.posi, gfa.names1.nega)]

#### ph, log-ratio, JEPS/xxxx
ph.wt.quad = sapply(all.gm.fit, function(x) {
  x$para.vals['JEPS','ph'] - x$para.vals['xxxx','ph']
})

#### pt, JEPS-xxxx
pt.wt.quad = sapply(all.gm.fit, function(x) {
  x$para.vals['JEPS','pt'] - x$para.vals['xxxx','pt']
})

#### analysis of ph, log2, JEPS/xxx for 866up, 234down genes (pt, selected)
ph.diff = ph.wt.quad[rownames(ph.nram)]
quantile(ph.diff)
#        0%        25%        50%        75%       100% 
#-2.6075436 -0.3812585  0.3387189  0.8320481  5.7429718
ph.for.ph.idff = sapply(all.gm.fit[names(ph.diff)], function(x) x$para.vals['JEPS','ph'])
plot(ph.for.ph.idff, ph.diff)
# no strong pattern
## what about up and down genes separately?
quantile(ph.diff[names(ph.diff) %in% gfa.names1.posi])
#         0%         25%         50%         75%        100%
#-2.60754360  0.07470066  0.50524756  0.94369882  5.74297179
quantile(ph.diff[names(ph.diff) %in% gfa.names1.nega])
#        0%        25%        50%        75%       100% 
#-2.5084648 -1.0619790 -0.7044681 -0.3129273  1.6734395
# both amplitude gets larger in JEPS

#### analysis of pt, JEPS-xxx for 866up, 234down genes (pt, selected)
pt.diff = pt.wt.quad[rownames(pt.nram)]
quantile(pt.diff)
#        0%        25%        50%        75%       100% 
#-4.4326722 -1.5679143 -0.9117306 -0.4732928  3.2609463
pt.for.pt.idff = sapply(all.gm.fit[names(pt.diff)], function(x) x$para.vals['JEPS','pt'])
plot(pt.for.pt.idff, pt.diff)
# no strong correlation other than near the boundary values
## what about up and down genes separately?
quantile(pt.diff[names(pt.diff) %in% gfa.names1.posi])
#        0%        25%        50%        75%       100% 
#-4.4326722 -1.4688488 -0.9211064 -0.5426458  2.2492703
quantile(pt.diff[names(pt.diff) %in% gfa.names1.nega])
#        0%        25%        50%        75%       100% 
#-3.7132351 -2.1318210 -0.7729948 -0.1054881  3.2609463
# as expected, the delay is smaller with downregulated

##### NRAM ph, pt, repeat each column to increase the width, rownames correspond to the references
r.n = 7 # repeat numb
ph.nram1 = as.numeric(apply(cbind(ph.wt.quad[rownames(ph.nram)], ph.nram[,-16]), 2, function(x) {
  rep(x, r.n)
}))
ph.nram1 = matrix(ph.nram1, nrow=nrow(ph.nram[,-16]))
rownames(ph.nram1) = rownames(ph.nram)
dim(ph.nram1)
#[1] 1100  112

###
pt.nram1 = as.numeric(apply(cbind(pt.wt.quad[rownames(pt.nram)], pt.nram[,-16]), 2, function(x) {
  rep(x, r.n)
}))
pt.nram1 = matrix(pt.nram1, nrow=nrow(pt.nram[,-16]))
rownames(pt.nram1) = rownames(pt.nram)
dim(pt.nram1)
dim(pt.nram1)
#[1] 1100  112

###
emp.mat = matrix(NA, nrow = nrow(quad), ncol = ncol(ph.nram1))
rownames(emp.mat) = rownames(quad)

ph.nram2 = pt.nram2 = emp.mat
ph.nram2[rownames(ph.nram1),] = ph.nram1
pt.nram2[rownames(pt.nram1),] = pt.nram1

sum(!is.na(ph.nram2[1:1972,1]))
#[1] 866 upregulated genes were subjected to NRAM
sum(!is.na(ph.nram2[1973:3262,1]))
#[1] 234 downregulated genes were subjected to NRAM

ams.1 = cbind(quad, wt.quad, ph.nram2, pt.nram2)

col.anno2 = rbind(matrix('',ncol=34,nrow=3), 
                  c('xxxx', 'JEPS/xxxx', 'PA,JEPS/xxxx', colnames(ph.nram)[-16], 'PT,JEPS-xxxx', colnames(pt.nram)[-16]),
                  matrix('',ncol=34,nrow=3))
col.anno2 = as.character(col.anno2)
colnames(ams.1) = col.anno2

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))
col.pt = colorRamp2(c(-2, -1, -0.4, 0, 0.4, 1, 2), c("forestgreen","green2", "olivedrab3","white","purple3", "magenta","violetred4"))

bot_anno4a = HeatmapAnnotation(foo = anno_text(colnames(ams.1)[1:126], just = 'right', rot=60,
                                              gp = gpar(fontsize = 7), location = 1),
                              show_annotation_name = F)
bot_anno4b = HeatmapAnnotation(foo = anno_text(colnames(ams.1)[127:238], just = 'right', rot=60,
                                               gp = gpar(fontsize = 7), location = 1),
                               show_annotation_name = F)
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))
lgd.pt = Legend( at=seq(-2,2,0.5),col_fun = col.pt, title = "hour",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
                 title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

ht_4a = Heatmap(ams.1[,1:126] , col = col, cluster_rows = F, cluster_columns = FALSE,
               name = 'ph.nram', show_row_dend = F,
               column_gap = unit(4,"mm"), 
               column_split = c(rep('ETI-FC',14), rep('NRAM on ETI Peak Amplitude', 112)),
               row_title_gp = gpar(fontsize = 20),
               show_row_names = F, show_column_names = F,
               column_title_gp = gpar(fontsize=10), 
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               use_raster = F,width = unit(85,'mm'), show_heatmap_legend = F,
               bottom_annotation = bot_anno4a,
               height = unit(120, 'mm'))

ht_4b = Heatmap(ams.1[,127:238] , col = col.pt, cluster_rows = F, cluster_columns = FALSE,
                name = 'pt.nram', show_row_dend = F,
                row_title_gp = gpar(fontsize = 20),
                show_row_names = F, show_column_names = F, 
                column_title = 'NRAM on ETI Peak Time',
                column_title_gp = gpar(fontsize=10), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(75,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno4b,
                height = unit(120, 'mm'))

hm_ETI.ph.NRAM = grid.grabExpr(draw(ht_4a,
                                    annotation_legend_list = lgd))
hm_ETI.pt.NRAM = grid.grabExpr(draw(ht_4b,
                                    annotation_legend_list = lgd.pt))

### ph proportion explained


jpeg("./outputs/Fig4.pa.pt.heatmap.v230327.jpeg",height = 160, width = 200, units = "mm",res = 300)
grid.arrange(hm_ETI.ph.NRAM, hm_ETI.pt.NRAM,
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
#grid.text('A', x = unit(-107,'mm'),y = unit(144,'mm'),gp = gpar(fontface = 'bold'))
#grid.text('B', x = unit(0,'mm'),y = unit(144,'mm'),gp = gpar(fontface = 'bold'))
dev.off()

### Fig 4D ph.nram mean
up.mod.genes = gfa.names1.posi[gfa.names1.posi %in% rownames(ph.nram)]
length(up.mod.genes)
#[1] 866

ph.nram.av = apply(ph.nram[up.mod.genes,], 2, mean)[1:15]
col.b = rep('gold3',15)
col.b[ph.nram.av < 0] = 'aquamarine3'

pdf('./outputs/Fig4D.pa.nram.av.v230327.pdf')
barplot(rev(ph.nram.av), horiz = T, las=1,
        col= rev(col.b), xlim=c(-0.225,0.28),
        xlab=expression('Mean Contribution to Peak Amplitude (log'[2]*')'),
        main='NRAM on ETI Peak Amplitude')
box()
abline(v=0, col='gray', lty=3, lwd=2)
dev.off()

### Fig 4E pt.nram mean
pt.nram.av = apply(pt.nram[up.mod.genes,], 2, mean)[1:15]
col.b = rep('purple3',15)
col.b[pt.nram.av < 0] = 'olivedrab3'

pdf('./outputs/Fig4E.pt.nram.av1.v230327.pdf')
# pt.nram.av sign is flipped to show acceleration of peak time
barplot(rev(-pt.nram.av), horiz = T, las=1,
        col= rev(col.b), xlim=c(-0.1,0.16),
        xlab=expression('Mean Contribution to Peak Time Acceleration (hour)'),
        main='NRAM on ETI Peak Time')
box()
abline(v=0, col='gray', lty=3, lwd=2)
dev.off()

#### make an equivalent plot with the PC1, 4hpi
### genotype names
c.names = colnames(Ed.dyn.dat)
c.names = strsplit(c.names, '\\.')
c.names = sapply(c.names, '[', 1)
genotypes = unique(c.names)

#### PC1 values, 4hpi
pca.results = pca.dyn$x
pc1.geno.time = pca.results[,1]
pc1.geno.time = matrix(pc1.geno.time, nrow=7)
pc1.gt.4h = pc1.geno.time[5,]
names(pc1.gt.4h) = genotypes

##### function definitions
load('./outputs/NRAM.algorithm.w.mean.se.RData')

############
##### Fig 4F sampling from upregulated genes for PC1
length(intersect(rownames(pt.nram), gfa.names1.posi))
#[1] 866 genes to sample from 
length(gfa.names1.posi)
#[1] 1972 genes

rep.numb = 500; set.seed(7)
comp.nram.pc1 = c()
library(R.utils)
for (i in 1:rep.numb) {
  sel.ugenes = sample(gfa.names1.posi, replace = T) #bootstrap
  sel.udata = Ed.dyn.dat[sel.ugenes,]
  sel.udata.n = apply(sel.udata, 1, function(x) x - mean(x))  # centering
  rot.udata = pca.dyn$rotation[sel.ugenes, ]
  pca.results = sel.udata.n %*% rot.udata 
  pc1.geno.time = pca.results[,1]
  pc1.geno.time = matrix(pc1.geno.time, nrow=7)
  pc1.gt.4h = pc1.geno.time[5,]
  names(pc1.gt.4h) = genotypes
  
  ### NRAM
  pc1.est = pc1.gt.4h
  pc1.se = rep(0.0001, 16)
  names(pc1.est) = names(pc1.se) = geno.names[names(pc1.est)]
  ave.m.out = run.ave.m(pc1.est, pc1.se , th.1=1, th.2=1)
  comp.nram.pc1 = rbind(comp.nram.pc1, ave.m.out[1:15])
}  # 14 sec
detach('package:R.utils', unload=T)

quant.comp.nram.pc1 = apply(comp.nram.pc1, 2, function(x) {
  q.v = quantile(x, probs=c(0.025, 0.975))
  m.v = mean(x)
  return(c(q.v, m.v))
})

col.b = rep('steelblue3', 15)
col.b[quant.comp.nram.pc1[3,] > 0] = 'palevioletred2'

pdf('./outputs/Fig4F.eti.pc1.nram.samp.v230327.pdf')
mp = barplot(rev(quant.comp.nram.pc1[3,]), col=rev(col.b), horiz = T,
        xlim=range(quant.comp.nram.pc1)*1.04, las=1,
        xlab=expression('Contribution to PC1 (log'[2]*')'),
        main='NRAM on ETI-PC1 of ETI-upregulated genes at 4 hpi')
box()
arrows(quant.comp.nram.pc1[1,],rev(mp), quant.comp.nram.pc1[2,], rev(mp), 
       angle=90, code=3, length=0.075)
dev.off()

##### Fig 4G bootstrapping from downregulated genes for PC1
length(intersect(rownames(pt.nram), gfa.names1.nega))
#[1] 234 this is too few to see the trend
length(gfa.names1.nega)
#[1] 1290 genes

rep.numb = 500; set.seed(7)
comp.nram.pc1 = c()
pc1.var.per = c()
library(R.utils)
for (i in 1:rep.numb) {
  sel.ugenes = sample(gfa.names1.nega, replace = T)  # bootstrap
  sel.udata = Ed.dyn.dat[sel.ugenes,]
  sel.udata.n = apply(sel.udata, 1, function(x) x - mean(x))  # centering
  rot.udata = pca.dyn$rotation[sel.ugenes, ]
  pca.results = sel.udata.n %*% rot.udata 
  pc1.geno.time = pca.results[,1]
  pc1.geno.time = matrix(pc1.geno.time, nrow=7)
  pc1.gt.4h = pc1.geno.time[5,]
  names(pc1.gt.4h) = genotypes
  
  ### NRAM
  pc1.est = pc1.gt.4h
  pc1.se = rep(0.0001, 16)
  names(pc1.est) = names(pc1.se) = geno.names[names(pc1.est)]
  ave.m.out = run.ave.m(pc1.est, pc1.se , th.1=1, th.2=1)
  comp.nram.pc1 = rbind(comp.nram.pc1, ave.m.out[1:15])
}  # 10 sec
detach('package:R.utils', unload=T)

quant.comp.nram.pc1 = apply(comp.nram.pc1, 2, function(x) {
  q.v = quantile(x, probs=c(0.025, 0.975))
  m.v = mean(x)
  return(c(q.v, m.v))
})

col.b = rep('steelblue3', 15)
col.b[quant.comp.nram.pc1[3,] > 0] = 'palevioletred2'

pdf('./outputs/Fig4G.eti.pc1.nram.nega.boot.v230327.pdf')
mp = barplot(rev(quant.comp.nram.pc1[3,]), col=rev(col.b), horiz = T,
             xlim=range(quant.comp.nram.pc1)*1.04, las=1,
             xlab=expression('Contribution to PC1 (log'[2]*')'),
             main='NRAM on ETI-PC1 of ETI-downregulated genes at 4 hpi')
box()
arrows(quant.comp.nram.pc1[1,],rev(mp), quant.comp.nram.pc1[2,], rev(mp), 
       angle=90, code=3, length=0.075)
dev.off()


#######################################
##############
###### basal level NRAM

##### glm.nb for 0h data
#### load data
rm(list=ls())
load('./outputs/exp.one.d.RData')

#### reorganize fixed effect factor
fix.ef2= as.character(fix.ef)
genotypes = unique(substr(fix.ef2, 1, 4))
genotypes = genotypes[1:16]
fix.ef2[grep('^JEPS\\..+\\.0h$|^GUS\\..+\\.0h', fix.ef2)] = 'JEPS.m.0h'
for (i in c(genotypes[-1], 'r2r1') ) {
  patt = paste0("^", i,"\\..+\\.0h$")
  rep = paste0(i,".m.0h")
  fix.ef2[grep(patt, fix.ef2)] = rep
}
fix.ef2 = factor(fix.ef2, levels=unique(fix.ef2))

#### fit glm.nb model for mean estimates
gene.names = rownames(exp.one.d)
date()
glm.coefs2 = apply(exp.one.d, 1, function(exp.d){
  glmnb1 = glm.nb(exp.d ~ -1 + fix.ef2 + offset(offset.cval.90),
                  link=log)
  glm.c = summary(glmnb1)$coef[,1:2]
  glm.c = glm.c / log(2)  # to convert the base of log to 2 from e.
  est = glm.c[,1]
  sem = glm.c[,2]
  fitted = glmnb1$fitted.values
  d.resid = residuals(glmnb1, type = 'deviance') / sqrt(fitted) /
    log(2) # to convert the base of log to 2 from e
  return(list(est=est, sem=sem, d.resid=d.resid))
})  # 33 min
date()

save(fix.ef, fix.ef2, offset.cval.90, glm.coefs2, file='./outputs/Ed.AvrRpt2.focus.0h.RData')

################
###### basal level NRAM

#### load data
rm(list=ls())
load('./outputs/Ed.AvrRpt2.focus.0h.RData')
genotypes = unique(substr(fix.ef2, 1, 4))
genotypes = genotypes[1:16]

#### function definitions
load('./outputs/NRAM.algorithm.w.mean.se.RData')

#### Ed-AvrRpt2 glm results, organize for 0h
ETI.glm = glm.coefs2

est.names = names(ETI.glm[[1]]$est)
est.names = substr(est.names, 8, nchar(est.names))
names.0h = grepl('\\.0h$', est.names)

ETI.est.0h = lapply(ETI.glm, function(x) {
  ETI.0h = (x$est[names.0h])
  ETI.0h.se = (x$sem[names.0h])
  
  names(ETI.0h) = names(ETI.0h.se) = genotypes
  return(list(est = ETI.0h[-c(17,18)],
              se = ETI.0h.se[-c(17,18)]))
  
})

#### run NRAM for each gene at 0h
library(R.utils)
ave.m.0h.ETI = t(sapply(ETI.est.0h, function(v.out) {
  m.est.geno = v.out$est
  m.se.geno = v.out$se
  names(m.est.geno) = names(m.se.geno) = geno.names
  run.ave.m(m.est.geno, m.se.geno, 
            th.1.method = '', th.2.method = 'BH')
}))  # 30 sec
detach('package:R.utils', unload=T)

#### select significant genes
sel.genes = apply(ave.m.0h.ETI,1, function(x) {
  sum(x[1:15] != 0) > 0 
})
sum(sel.genes)
#[1] 4562
sel.ave.m.0h.ETI = ave.m.0h.ETI[sel.genes,]

save(sel.ave.m.0h.ETI, ETI.est.0h, file='./outputs/sel.ave.m.0h.ETI.RData')

#################### 
####### Fig 6AB and SuppFig
####### heatmap with 0h.NRAM for 10833 genes.dyn x 180 treatment/genotype/time

#### load data
rm(list=ls())
load('./outputs/Ed.dyn.dat.v230327.RData') # corrected mock subtracted
load('./outputs/sel.ave.m.0h.ETI.RData')
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered

dim(sel.ave.m.0h.ETI)
#[1] 4562   16

#### heatmap
ams.1 = sel.ave.m.0h.ETI[,1:15]
### cosine correlation among genes
row_distances_subset = 1 - cosine(t(ams.1)) # 40 sec
###
distances = as.dist(row_distances_subset, diag = FALSE, upper = FALSE)
hclust.genes.dyn = hclust( distances, method = "complete" ) 
 
#### reorganize according to the gamma-modeled genes
length(gfa.names1.posi)
#[1] 1972
length(intersect(rownames(sel.ave.m.0h.ETI), gfa.names1.posi))
#[1] 765
length(gfa.names1.nega)
#[1] 1290 
length(intersect(rownames(sel.ave.m.0h.ETI), gfa.names1.nega))
#[1] 383

posi.genes = rownames(sel.ave.m.0h.ETI) %in% gfa.names1.posi
nega.genes = rownames(sel.ave.m.0h.ETI) %in% gfa.names1.nega
resp.genes = rep(0,nrow(sel.ave.m.0h.ETI))
resp.genes[posi.genes] = 4
resp.genes[nega.genes] = -4

ams.1 = cbind('ETI-responsive genes'=resp.genes, sel.ave.m.0h.ETI[,1:15])
col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))
bot_anno5x = HeatmapAnnotation(foo = anno_text(colnames(ams.1), just = 'right', rot=60,
                                               gp = gpar(fontsize = 7), location = 1),
                               show_annotation_name = F)

ht_5x = Heatmap(ams.1 , col = col, cluster_rows = hclust.genes.dyn, cluster_columns = FALSE,
               name = 'ETI.0h.NRAM', show_row_dend = F,
               column_gap = unit(2,"mm"), 
               column_split = c('A', rep('B',15)),
               row_title_gp = gpar(fontsize = 20),
               show_row_names = F, show_column_names = F, 
               column_title = NULL,
               border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
               use_raster = F,width = unit(80,'mm'), show_heatmap_legend = F,
               bottom_annotation = bot_anno5x,
               height = unit(150, 'mm'))

jpeg("./outputs/Fig6B.ETI.0h.NRAM.all.v230327.jpeg",height = 220, width = 120, units = "mm",res = 300)
draw(ht_5x, annotation_legend_list = lgd)
decorate_heatmap_body('ETI.0h.NRAM', column_slice = 2, {
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


sum(ams.1[,1] > 0 & ams.1[,'P;S'] > 0)
#[1] 533 out of 765 (sum(ams.1[,1] > 0) common upregulated genes

########## Fig 6A 0hpi NRAM heatmap
##### JEPS and xxxx time courses for references
quad = Ed.dyn.dat[c(gfa.names1.posi, gfa.names1.nega), 106:112]
wt  = Ed.dyn.dat[c(gfa.names1.posi, gfa.names1.nega), 1:7]
wt.quad = wt - quad

### JEPS/xxxx 0h
basal.diff = sapply(ETI.est.0h[c(gfa.names1.posi, gfa.names1.nega)], function(x) {
  x$est['JEPS'] - x$est['xxxx']
})

#### NRAM values in the standard gene order
na.mat = matrix(0, ncol=15, nrow=length(c(gfa.names1.posi, gfa.names1.nega)))
dimnames(na.mat) = list(c(gfa.names1.posi, gfa.names1.nega), colnames(sel.ave.m.0h.ETI)[1:15])
comm.names = c(gfa.names1.posi, gfa.names1.nega)[c(gfa.names1.posi, gfa.names1.nega) %in% rownames(sel.ave.m.0h.ETI)]
na.mat[comm.names,] = sel.ave.m.0h.ETI[comm.names, 1:15]
ams.0 = cbind('basal, JEPS/xxxx'=basal.diff, na.mat)

r.n = 7
ams.1 = matrix(as.numeric(apply(ams.0, 2, function(x) {
  rep(x, r.n)
})), nrow=nrow(ams.0))

ams.1 = cbind(quad, wt.quad, ams.1)

col.anno2 = rbind(matrix('',ncol=18,nrow=(r.n-1)/2), 
                  c('xxxx', 'JEPS/xxxx', colnames(ams.0)),
                  matrix('',ncol=18,nrow=(r.n-1)/2))
col.anno2 = as.character(col.anno2)
colnames(ams.1) = col.anno2

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))

bot_anno6a = HeatmapAnnotation(foo = anno_text(colnames(ams.1), just = 'right', rot=60,
                                               gp = gpar(fontsize = 7), location = 1),
                               show_annotation_name = F)
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

ht_6a = Heatmap(ams.1, col = col, cluster_rows = F, cluster_columns = FALSE,
                name = 'nram.0h', show_row_dend = F,
                column_gap = unit(4,"mm"), 
                column_split = c(rep('ETI-FC',14), rep('NRAM on Basal mRNA level', 112)),
                row_title_gp = gpar(fontsize = 20),
                show_row_names = F, show_column_names = F,
                column_title_gp = gpar(fontsize=10), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(85,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno6a,
                height = unit(120, 'mm'))

jpeg("./outputs/Fig6A.0h.nram.heatmap.v230327.jpeg",height = 160, width = 200, units = "mm",res = 300)
draw(ht_6a, annotation_legend_list = lgd)
decorate_heatmap_body('nram.0h', column_slice=2, {
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
dev.off()


############
###### Fig 7 in Discussion
#### overlap with ETI-upregulated by venn diagram
### load data
rm(list=ls())
load('./outputs/ETI.gfa.sorted.genes.v230327.RData') # modeled gene sets, up and down, pt ordered
load('./outputs/r2r1.up.down.genes.RData')
load('./outputs/sel.ave.m.0h.ETI.RData')

#### flip r2r1_up, r2r1_down to R2R1_down, R2R1_up for easier interpretations
ETI.posi.r2r1.genes = unique(c(gfa.names1.posi,r2r1.up.genes,r2r1.down.genes))
df.posi = data.frame('ETI_up' = ETI.posi.r2r1.genes %in% gfa.names1.posi,
                     'R2R1_up' = ETI.posi.r2r1.genes %in% r2r1.down.genes,
                     'R2R1_down' = ETI.posi.r2r1.genes %in% r2r1.up.genes)

ETI.nega.r2r1.genes = unique(c(gfa.names1.nega,r2r1.up.genes,r2r1.down.genes))
df.nega = data.frame('ETI_down' = ETI.nega.r2r1.genes %in% gfa.names1.nega,
                     'R2R1_up' = ETI.nega.r2r1.genes %in% r2r1.down.genes,
                     'R2R1_down' = ETI.nega.r2r1.genes %in% r2r1.up.genes)

p1 = ggvenn(df.posi, show_percentage=F, fill_color = c('red','orange','green'),
            set_name_size = 4)
p2 = ggvenn(df.nega, show_percentage=F, fill_color = c('blue','orange','green'),
            set_name_size = 4)
pdf('./outputs/Fig7AB.vennD.ETI.r2r1.v230327.pdf', height = 3, width = 5)
grid.arrange(p1,p2, nrow=1)
dev.off()

## R2R1 down strongly associated with P;S at 0hpi > 0
### venn diagram for this
basal.p.s.up.genes = rownames(sel.ave.m.0h.ETI)[sel.ave.m.0h.ETI[,'P;S'] > 0]
length(basal.p.s.up.genes)
#[1] 738

venn.diagram(list('R2R1_down' = r2r1.up.genes, 
                  'P;S,0hpi_positive' = basal.p.s.up.genes),
             filename='./outputs/Fig7C.r2r1.p.s.venn.r1.tiff', 
             width = 1800, height = 1800,
             cex = 1.4,
             cat.dist = c(0.02,0.02),
             cat.pos = c(2,10),
             cat.cex = 1.4,
             fill = c('green','orange'),
             alpha= c(0.5,0.5))



######################
######### For Text S1, details of the PCA in Fig 1A
######### PC1, ETI
######### PC3, infiltration shock
#### GO analysis of the genes that have abs(correlation) > 0.9 for each PC

#### load data
rm(list=ls())
load('./outputs/four.geno.all.pca.RData')

### center glm.c.slim
glm.c.slim.cen = apply(glm.c.slim, 2, function(x) x-mean(x))

#### PC1-correlated genes
pc1.coord = pca.slim$x[,1]
pc1.coord = pc1.coord / sqrt(sum(pc1.coord^2))

gene.pc1.coord = apply(glm.c.slim.cen, 2, function(x) sum(x * pc1.coord))

cor.w.pc1 = apply(glm.c.slim.cen, 2, function(x) cor(x, pc1.coord))
cor.w.pc1.sorted = sort(cor.w.pc1, decreasing = T)
hist(cor.w.pc1) # highly skewed toward negative, but small peak near 1.
PC1.pcor.09.genes = names(cor.w.pc1.sorted[cor.w.pc1.sorted > 0.9])
length(PC1.pcor.09.genes)
#[1] 640
PC1.ncor.09.genes = rev(names(cor.w.pc1.sorted[cor.w.pc1.sorted < -0.9]))
length(PC1.ncor.09.genes)
#[1] 101

write.csv(PC1.pcor.09.genes, './outputs/PC1.09.pcor.genes.csv')
write.csv(PC1.ncor.09.genes, './outputs/PC1.09.ncor.genes.csv')

#### PC2-correlated genes
pc2.coord = pca.slim$x[,2] 
pc2.coord = pc2.coord / sqrt(sum(pc2.coord^2))

gene.pc2.coord = apply(glm.c.slim.cen, 2, function(x) sum(x * pc2.coord))

cor.w.pc2 = apply(glm.c.slim.cen, 2, function(x) cor(x, pc2.coord))
cor.w.pc2.sorted = sort(cor.w.pc2, decreasing = T)
hist(cor.w.pc2) # relatively symmetric
PC2.pcor.09.genes = names(cor.w.pc2.sorted[cor.w.pc2.sorted > 0.9])
length(PC2.pcor.09.genes)
#[1] 15
PC2.ncor.09.genes = rev(names(cor.w.pc2.sorted[cor.w.pc2.sorted < -0.9]))
length(PC2.ncor.09.genes)
#[1] 117

write.csv(PC2.pcor.09.genes, './outputs/PC2.09.pcor.genes.csv')
write.csv(PC2.ncor.09.genes, './outputs/PC2.09.ncor.genes.csv')

#### PC3-correlated genes
pc3.coord = pca.slim$x[,3] 
pc3.coord = pc3.coord / sqrt(sum(pc3.coord^2))

gene.pc3.coord = apply(glm.c.slim.cen, 2, function(x) sum(x * pc3.coord))

cor.w.pc3 = apply(glm.c.slim.cen, 2, function(x) cor(x, pc3.coord))
cor.w.pc3.sorted = sort(cor.w.pc3, decreasing = T)
hist(cor.w.pc3) # pretty symmetric
PC3.pcor.09.genes = names(cor.w.pc3.sorted[cor.w.pc3.sorted > 0.9])
length(PC3.pcor.09.genes)
#[1] 11
PC3.ncor.09.genes = rev(names(cor.w.pc3.sorted[cor.w.pc3.sorted < -0.9]))
length(PC3.ncor.09.genes)
#[1] 2

rad.cor.w.pc3 = acos(cor.w.pc3)
hist(rad.cor.w.pc3 / pi, breaks=50); abline(v=0.5, col='red')
# not that symmetric - not that random

write.csv(PC3.pcor.09.genes, './outputs/PC3.09.pcor.genes.csv')
write.csv(PC3.ncor.09.genes, './outputs/PC3.09.ncor.genes.csv')

### gene by gene visualization for four genotypes
genox = c('JEPS', 'xxxx', 'r1r2', 'GUS')

gene.plot = function(dat, title=rep('',3)) {
  ymax = max(dat)
  ymin = min(dat)
  opar = par(mfrow=c(3, 4), mar=c(4,4,3,0.5) )
  for (j in 1:3) {
    for (i in 0:3) {
      mock.dat = dat[j, 1:3 + 10*i]
      ed.dat = dat[j, 4:10 + 10*i]
      plot(0,0, type='n', 
           xlim = c(0,6), ylim=c(ymin, ymax),
           xlab = 'Time (hpi)', ylab='Relative mRNA level (log2)')
      lines(c(0,2,5), mock.dat, col='blue', type='o')
      lines(0:6, ed.dat, col='red', type='o')
      text(0, 0.95*(ymax-ymin) + ymin, genox[i+1], pos=4)
      if (i == 0) {
        mtext(title[j], line = 0.8)
      }
    }
  }
  par(opar)
}

pdf('./outputs/Fig.for.TextS1.PC1.PC2.genes.pdf', width = 8, height=7)
dat = rbind(pc1.coord, pc2.coord, (pc2.coord - pc1.coord)/sqrt(2))
gene.plot(dat, c('PC1 gene','PC2 gene','PC2 gene - PC1 gene'))
dev.off()

## ETI-up and downregulated
load('./outputs/ETI.gfa.sorted.genes.v230327.RData')

pdf('./outputs/Fig.for.TextS1.PC12gene.plot.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
dev.off()

col.eti = rep(NA, length(gene.pc1.coord))
col.eti[names(gene.pc1.coord) %in% gfa.names1.posi] = 'red'
col.eti[names(gene.pc1.coord) %in% gfa.names1.nega] = 'blue'
pdf('./outputs/Fig.for.TextS1.PC12gene.plot2.v230327.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
points(gene.pc1.coord, gene.pc2.coord, cex=0.4, pch=16, col=col.eti)
dev.off()

col.1plus2 = rep(NA, length(gene.pc1.coord))
col.1plus2[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) > 2] = 'orange'
col.1plus2[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) < -2.5] = 'turquoise'
pdf('./outputs/Fig.for.TextS1.PC12gene.plot3.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
points(gene.pc1.coord, gene.pc2.coord, cex=0.4, pch=16, col=col.1plus2)
dev.off()
pc1p2.up = names(gene.pc1.coord[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) > 2])
pc1p2.down = names(gene.pc1.coord[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) < -2.5])
length(pc1p2.up); length(pc1p2.down)
#[1] 1877
#[1] 967

write.csv(pc1p2.up, './outputs/pc1p2.up.genes.csv')
write.csv(pc1p2.down, './outputs/pc1p2.down.genes.csv')

## for comparison, PC1 only - not used
pc1.up = names(gene.pc1.coord[gene.pc1.coord > 2])
length(pc1.up)
#[1] 2204
write.csv(pc1.up, './outputs/pc1.up.genes.csv')
length(intersect(pc1p2.up, pc1.up))
#[1] 1649 highly overlapping, of course
# look into specific ones
pc1p2.up.m.pc1.up = setdiff(pc1p2.up, pc1.up)
pc1.up.m.pc1p2.up = setdiff(pc1.up, pc1p2.up)
write.csv(pc1.up.m.pc1p2.up, './outputs/pc1.up.minus.pc1p2.up.genes.csv')
write.csv(pc1p2.up.m.pc1.up, './outputs/pc1p2.up.minus.pc1.up.genes.csv')

#### PC2-PC1
col.2minus1 = rep(NA, length(gene.pc1.coord))
col.2minus1[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) <= 2 &
              (gene.pc1.coord + gene.pc2.coord)/sqrt(2) >= -2.5 &
              (gene.pc2.coord - gene.pc1.coord)/sqrt(2) > 2.5] = 'magenta'
col.2minus1[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) <= 2 &
              (gene.pc1.coord + gene.pc2.coord)/sqrt(2) >= -2.5 &
              (gene.pc2.coord - gene.pc1.coord)/sqrt(2) <= -2.5] = 'gold'
pdf('./outputs/Fig.for.TextS1.PC12gene.plot4.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
points(gene.pc1.coord, gene.pc2.coord, cex=0.4, pch=16, col=col.2minus1)
dev.off()
pc2m1.up = names(gene.pc1.coord[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) <= 2 &
                                  (gene.pc1.coord + gene.pc2.coord)/sqrt(2) >= -2.5 &
                                  (gene.pc2.coord - gene.pc1.coord)/sqrt(2) > 2.5])
pc2m1.down = names(gene.pc1.coord[(gene.pc1.coord + gene.pc2.coord)/sqrt(2) <= 2 &
                                    (gene.pc1.coord + gene.pc2.coord)/sqrt(2) >= -2.5 &
                                    (gene.pc2.coord - gene.pc1.coord)/sqrt(2) <= -2.5])
length(pc2m1.up); length(pc2m1.down)
#[1] 1071
#[1] 458
write.csv(pc2m1.up, './outputs/pc2.minus.pc1.up.genes.csv')
write.csv(pc2m1.down, './outputs/pc2.minus.pc1.down.genes.csv')

#### Yang et al. 2020 Circadian and Diurnal data
#### https://doi.org/10.1074/jbc.RA120.013513
circ.g = read.csv('./data/Yang.2020.circadian.159584_1_supp_509367_q8lbct.csv', header=T, row.names=1)
diur.g = read.csv('./data/Yang.2020.diurnal.159584_1_supp_509370_q8lbct.csv', header=T, row.names=1)
### remove transcript entries, leave only gene entries
circ.g.g = unique(circ.g[,'Gene_ID'])
circ.g.g1 = rownames(circ.g)[rownames(circ.g) %in% circ.g.g]
diur.g.g = unique(diur.g[,'Gene_ID'])
diur.g.g1 = rownames(diur.g)[rownames(diur.g) %in% diur.g.g]
length(circ.g.g1); length(diur.g.g1)
#[1] 6020
#[1] 6679
# these numbers are the same as in the text of Yang et al. Good.
circ.gs = circ.g[circ.g.g1,]
diur.gs = diur.g[diur.g.g1,]
# since the number is pretty high, apply more stringent "integrated p" cutoff
circ.gs1 = circ.gs[circ.gs$meta2d_BH.Q < 0.0001,]
diur.gs1 = diur.gs[diur.gs$meta2d_BH.Q < 0.0001,]

circ.gs1.g = rownames(circ.gs1)
diur.gs1.g = rownames(diur.gs1)
length(circ.gs1.g); length(diur.gs1.g)
#[1] 1517
#[1] 967

circ.gs1.g.c = circ.gs1.g[circ.gs1.g %in% names(gene.pc1.coord)]
diur.gs1.g.c = diur.gs1.g[diur.gs1.g %in% names(gene.pc1.coord)]
cd.gs1.gcc = circ.gs1.g.c[circ.gs1.g.c %in% diur.gs1.g.c]  # intersect
length(circ.gs1.g.c); length(diur.gs1.g.c); length(cd.gs1.gcc)
#[1] 1465
#[1] 917
#[1] 370
# pretty good for visualization

### circadian overlay
col.cd = rep(NA, length(gene.pc1.coord))
names(col.cd) = names(gene.pc1.coord)
col.cd[circ.gs1.g.c] = 'magenta'

pdf('./outputs/Fig.for.TextS1.PC12gene.circadian.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
points(gene.pc1.coord, gene.pc2.coord, cex=0.4, pch=16, col=col.cd)
text(20,18,'Circadian', pos=2, cex=2)
dev.off()

### diurnal overlay
col.cd = rep(NA, length(gene.pc1.coord))
names(col.cd) = names(gene.pc1.coord)
col.cd[diur.gs1.g.c] = 'green'

pdf('./outputs/Fig.for.TextS1.PC12gene.diurnal.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
points(gene.pc1.coord, gene.pc2.coord, cex=0.4, pch=16, col=col.cd)
text(20,18,'Diurnal', pos=2, cex=2)
dev.off()

### circadian, diurnal phase
hour.col = rainbow(25)
col.crain = rep(NA, length(gene.pc1.coord))
names(col.crain) = names(gene.pc1.coord)
col.crain[circ.gs1.g.c] = hour.col[1+round(circ.gs1[circ.gs1.g.c,'ARS_adjphase'])]

pdf('./outputs/Fig.for.TextS1.PC12gene.circ.phase.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
points(gene.pc1.coord, gene.pc2.coord, cex=0.4, pch=16, col=col.crain)
text(20,18,'Circadian', pos=2, cex=2)
## color legend
hour.tick = paste0(0:24, 'h')
hour.tick[0:24 %% 6 != 0] = ''
color.legend(-12,-13, -10.5,-3,hour.tick, hour.col, gradient='y')
dev.off()


col.crain = rep(NA, length(gene.pc1.coord))
names(col.crain) = names(gene.pc1.coord)
col.crain[diur.gs1.g.c] = hour.col[1+round(diur.gs1[diur.gs1.g.c,'ARS_adjphase'])]

pdf('./outputs/Fig.for.TextS1.PC12gene.diur.phase.pdf', width=8, heigh=8)
plot(gene.pc1.coord, gene.pc2.coord, cex=0.4, col='gray25',
     xlab='PC1.gene', ylab='PC2.gene',
     xlim=c(-14,19), ylim=c(-14,19))
abline(v=0, col='gray30')
abline(h=0, col='gray30')
abline(0,-1, col='green')
points(gene.pc1.coord, gene.pc2.coord, cex=0.4, pch=16, col=col.crain)
text(20,18,'Diurnal', pos=2, cex=2)
## color legend
hour.tick = paste0(0:24, 'h')
hour.tick[0:24 %% 6 != 0] = ''
color.legend(-12,-13, -10.5,-3,hour.tick, hour.col, gradient='y')
dev.off()


##############
###### Fig S6. mock model fit heatmap

#### load data
rm(list=ls())
load('./outputs/EdGUS.corrected.mock.v210906.RData')

#### data and model time course plots for all genotypes per gene
geno.g.tc.plot = function(gene, tx1=0:6, ud='up') {
  ### model values
  para.v = all.gm.fit[[gene]]$para.vals
  t01 = all.gm.fit[[gene]]$para.vals[1,'t0']
  m.v = c()
  for (genot1 in rownames(para.v)) {
    para.geno = para.v[genot1,]
    mg.v = g.distr.pt0.l2a(tx1, ph=para.geno['ph'], l2a=para.geno['l2a'], pt=para.geno['pt'],
                           t0 = t01)
    m.v = cbind(m.v, mg.v)
  }
  
  ### plot y-range
  y.min1 = min(m.v)
  y.max1 = max(m.v)
  
  g.exp.Ed = t(G.corrected.tcv[[gene]]$Ed.vsm)
  y.min2 = min(min(g.exp.Ed), y.min1)
  y.max2 = max(max(g.exp.Ed), y.max1)
  
  ### data plot
  matplot(0:6, g.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
          xlab=NA, ylab=NA,
          ylim=c(y.min2, y.max2), xlim=c(0,6), xaxt='n')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  matlines(0:6, g.exp.Ed, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=1)
  axis(1, at=0:6)
  if (ud == 'up') {
    text(0,y.max2-0.1*(y.max2-y.min2), paste(gene,'data'), pos=4, cex=0.8)
  } else {
    text(0,y.min2+0.1*(y.max2-y.min2), paste(gene,'data'), pos=4, cex=0.8)
  }
  
  ### model plot
  matplot(tx1, m.v, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6),
          xlab=NA, ylab=NA,
          ylim=c(y.min2, y.max2), xlim=c(0,6), xaxt='n')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  matlines(tx1, m.v, type='l', col=rainbow(18, end=0.65), lty=rep(c(1,2,6), 6), lwd=1)
  axis(1, at=0:6)
  if (ud == 'up') {
    text(0,y.max2-0.1*(y.max2-y.min2), paste(gene,'model'), pos=4, cex=0.8)
  } else {
    text(0,y.min2+0.1*(y.max2-y.min2), paste(gene,'model'), pos=4, cex=0.8)
  }
}

###### MAIN
genotypes = rownames(G.corrected.tcv[[1]]$mock.7tp)

#### compile mock model and discrepancy, -> mock mean estimates
mock.model = t(sapply(G.corrected.tcv, function(x) {
  mock.m = x$mock.7tp[,c('0','2','5')]
  as.numeric(t(mock.m))
}))

mock.disc = t(sapply(G.corrected.tcv, function(x) {
  mock.d = x$mock.vsm
  as.numeric(t(mock.d))
}))

mock.est = mock.disc + mock.model

mean(mock.model);mean(mock.est)  # grand means
#[1] 6.875299
#[1] 6.953372

mock.model = mock.model - mean(mock.model)
mock.est = mock.est - mean(mock.est)

#### by heatmap
### heatmap
ams.1 = cbind(mock.est, mock.model, mock.disc)

col.anno1 = rbind(matrix('',ncol=18,nrow=1), genotypes, matrix('',ncol=18,nrow=1))
col.anno1 = as.character(col.anno1)
colnames(ams.1) = rep(col.anno1, 3)

col = colorRamp2(c(-5, -2.6, -1.3, 0, 1.3, 2.6, 5), c("blue4","dodgerblue2", "aquamarine3","white","gold3", "tan1","red4"))

bot_anno4x = HeatmapAnnotation(foo = anno_text(colnames(ams.1), just = 'right', rot=60,
                                               gp = gpar(fontsize = 7), location = 1),
                               show_annotation_name = F)
lgd = Legend( at=seq(-6,6,2),col_fun = col, title = "log2",legend_height = unit(25,"mm"),grid_width = unit(3,"mm"),
              title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6),title_gap = unit(3,"mm"))

ht_4x = Heatmap(ams.1, col = col, cluster_rows = F, cluster_columns = FALSE,
                name = 'mock.model.fit', show_row_dend = F,
                column_gap = unit(4,"mm"), 
                column_split = c(rep('_Data',54), rep('_Model', 54), rep('Data - Model', 54)),
                row_title_gp = gpar(fontsize = 20),
                show_row_names = F, show_column_names = F, 
                column_title_gp = gpar(fontsize=10), 
                border = "gray20", row_title_rot = 0, na_col = "#E0E0E0", 
                use_raster = F,width = unit(240,'mm'), show_heatmap_legend = F,
                bottom_annotation = bot_anno4x,
                height = unit(120, 'mm'))

jpeg("./outputs/FigS6.mock.model.fit.jpeg",height = 160, width = 270, units = "mm",res = 300)
draw(ht_4x, annotation_legend_list = lgd)
decorate_heatmap_body('mock.model.fit', column_slice = 1, {
  grid.lines(c(1/18,1/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(2/18,2/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(3/18,3/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(4/18,4/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(5/18,5/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(6/18,6/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(7/18,7/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(8/18,8/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(9/18,9/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(10/18,10/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(11/18,11/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(12/18,12/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(13/18,13/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(14/18,14/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(15/18,15/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(16/18,16/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(17/18,17/18),c(0,1), gp=gpar(col='gray40'))
})
decorate_heatmap_body('mock.model.fit', column_slice = 2, {
  grid.lines(c(1/18,1/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(2/18,2/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(3/18,3/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(4/18,4/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(5/18,5/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(6/18,6/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(7/18,7/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(8/18,8/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(9/18,9/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(10/18,10/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(11/18,11/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(12/18,12/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(13/18,13/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(14/18,14/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(15/18,15/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(16/18,16/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(17/18,17/18),c(0,1), gp=gpar(col='gray40'))
})
decorate_heatmap_body('mock.model.fit', column_slice = 3, {
  grid.lines(c(1/18,1/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(2/18,2/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(3/18,3/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(4/18,4/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(5/18,5/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(6/18,6/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(7/18,7/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(8/18,8/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(9/18,9/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(10/18,10/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(11/18,11/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(12/18,12/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(13/18,13/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(14/18,14/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(15/18,15/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(16/18,16/18),c(0,1), gp=gpar(col='gray40'))
  grid.lines(c(17/18,17/18),c(0,1), gp=gpar(col='gray40'))
})
grid.text('A', x = unit(7,'mm'),y = unit(136,'mm'),gp = gpar(fontface = 'bold'))
grid.text('B', x = unit(89,'mm'),y = unit(136,'mm'),gp = gpar(fontface = 'bold'))
grid.text('C', x = unit(170,'mm'),y = unit(136,'mm'),gp = gpar(fontface = 'bold'))

dev.off()



