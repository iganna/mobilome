library('foreach')
library(doParallel)
library("optparse")
source("synteny_infer.R")

# Rscript cmp_genomes.R --path.work ../../ --ref.acc 8236


args = commandArgs(trailingOnly=TRUE)


option_list = list(
  make_option(c("--path.work"), type="character", default=NULL, 
              help="path to working directory", metavar="character"),
  # make_option(c("-q", "--path.query"), type="character", default=NULL, 
  #             help="path to query chromosomes", metavar="character"),
  make_option(c("--path.base"), type="character", default=NULL, 
              help="path to base chromosomes", metavar="character"),
  make_option(c("--path.aln.pref"), type="character", default=NULL, 
              help="path with alignments", metavar="character"),
  make_option(c("--alternative.ref"), type="character", default=NULL, 
              help="alternative references", metavar="character"),
  
  # make_option(c("-c", "--file.chr.len.acc"), type="character", default=NULL, 
  #             help="file with lengths of chromosomes of query accessions", metavar="character"),
  make_option(c("--file.chr.len.ref"), type="character", default=NULL, 
              help="file with lengths of chromosomes of reference accessions", metavar="character"),
  make_option(c("--n.cores"), type="character", default=NULL, 
              help="numer of cores: 10 max", metavar="character"),
  
  
  make_option(c("--n.chr.ref"), type="character", default=NULL, 
              help="number of chromosomes in the reference genome", metavar="character"),
  make_option(c("--n.chr.acc"), type="character", default=NULL, 
              help="number of chromosomes in the accessions", metavar="character"),
  
  make_option(c("--ref.acc"), type="character", default=NULL, 
              help="reference accessions", metavar="character"),
  
  make_option(c("--all.vs.all"), type="character", default=NULL, 
              help="alignment of all chromosomes vs all or not: T/F", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

print(opt)
# return()

if (!is.null(opt$path.work)) path.work <- opt$path.work
# if (!is.null(opt$path.query)) path.query <- opt$path.query
if (!is.null(opt$path.base)) path.base <- opt$path.base
if (!is.null(opt$ref.acc)) base.acc.ref <- opt$ref.acc

if (!is.null(opt$alternative.ref)) alternative.ref <- opt$alternative.ref
if (!is.null(opt$path.aln.pref)) path.pref <- opt$path.aln.pref

# if (!is.null(opt$file.chr.len.acc)) file.chr.len.acc <- opt$file.chr.len.acc
if (!is.null(opt$file.chr.len.ref)) file.chr.len.ref <- opt$file.chr.len.ref

if (!is.null(opt$n.chr.ref)) n.chr.ref <- opt$n.chr.ref
if (!is.null(opt$n.chr.acc)) n.chr.acc <- opt$n.chr.acc
if (!is.null(opt$all.vs.all)) all.vs.all <- as.logical(opt$all.vs.all)

if (!is.null(opt$n.cores)) {
  n.cores <- min(10, as.numeric(opt$n.cores))
} else {
  n.cores = 10
}

# ===========================================================================
# ===========================================================================

# 10 CORES MAXIMUM!!!!!!!!!!!
myCluster <- makeCluster(n.cores, # number of cores to use           # 10 CORES MAXIMUM!!!!!!!!!!!
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

# ------------------------------------------------------------------

# acc.ref = c('0', '10002', '10015','10024', '1741', '22005', '6046', '6124', '6244', '8236', '9537') 
# 0,10002,10015,10024,1741,22005,6046,6124,6244,8236,9537
print(path.pref)
acc.ref = strsplit(alternative.ref, ',')[[1]]
print(acc.ref)

path.aln = paste(path.pref, acc.ref, '/', sep = '')
path.aln.ref =  paste(path.pref, base.acc.ref, '/', sep = '')

# ------------------------------------------------------------------
# Get accession names
aln.pattern = '_postgap3.rds'
accessions = c()
for(p in path.aln){
  accessions = c(accessions, unique(sapply(list.files(p, pattern=paste('*', aln.pattern, sep = '')), function(s) strsplit(s, '_')[[1]][1])))
}
acc.tbl = table(accessions)
accessions = names(acc.tbl)

print(accessions)

# -------------------
# Combinations of chromosomes query-base to chreate the alignments
chromosome.pairs = c()
for(i.query in 1:length(accessions)){
  for(query.chr in 1:n.chr.acc){
    if(!all.vs.all){
      if(query.chr > n.chr.ref) next
      chromosome.pairs = rbind(chromosome.pairs, c(i.query, query.chr, query.chr))
      next
    }
    for(base.chr in 1:n.chr.ref){
      chromosome.pairs = rbind(chromosome.pairs, c(i.query, query.chr, base.chr))
    }
  }
}

# -------------------
# Length of reference chromosomes
if(!file.exists(file.chr.len.ref)){
  chr.len = c()
} else {
  chr.len = readRDS(file.chr.len.ref)
}

add.flag = F
for(i.chr in 1:n.chr.acc){
  for(acc in unique(c(acc.ref, base.acc.ref))){
    # print(c(i.chr, acc))
    if(!is.null(chr.len)){
      idx = (chr.len$acc == acc) & (chr.len$chr == i.chr)
      if(sum(idx) > 0) next
    }
    add.flag = T
    g = seqinr::read.fasta(paste(path.base, acc, '_chr', i.chr, '.fasta', sep = ''))[[1]]
    chr.len = rbind(chr.len, data.frame(acc = acc, chr = i.chr, len = length(g)))
  }
}
if(add.flag) saveRDS(chr.len, file.chr.len.ref)



# ------------------------------------------------------------------


print(base.acc.ref)

if(!dir.exists(path.work)) system(paste('mkdir ', path.work, sep = ''))

# path.work = paste(path.work, 'cmp_', base.acc.ref, '/', sep = '')
# print(path.work)
# if(!dir.exists(path.work)) system(paste('mkdir ', path.work, sep = ''))

path.cmp = paste(path.work, 'cmp_', base.acc.ref, '/', sep = '')
print(path.cmp)
if(!dir.exists(path.cmp)) system(paste('mkdir ', path.cmp, sep = ''))

flag.for = F
tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs))  %dopar% {  # which accession to use
  
# flag.for = T
# for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  acc = accessions[chromosome.pairs[i.chr.pair, 1]]
  query.chr = chromosome.pairs[i.chr.pair, 2]
  base.chr = chromosome.pairs[i.chr.pair, 3]
  
  
  print(acc)
  file.out = paste(path.cmp, acc, '_', query.chr, '_', base.chr, '_cmp_to_ref_consensus.rds', sep = '')
  if (file.exists(file.out)) {
    if(flag.for){
      next
    } else {
      return(NULL)  
    }
  }

  cmp.all = c()
  for (i.alt in 1:length(path.aln)){
    
    p = path.aln[i.alt]
    base.acc.alt = acc.ref[i.alt]
    if(base.acc.alt == acc) next
    
    f.acc2ref = paste(path.aln.ref, acc, '_', query.chr, '_', base.chr, aln.pattern, sep = '')
    f.acc2alt = paste(p, acc, '_', query.chr, '_', base.chr, aln.pattern, sep = '')
    f.alt2ref = paste(path.aln.ref, base.acc.alt, '_', query.chr, '_', base.chr, aln.pattern, sep = '')
    if(!file.exists(f.acc2ref) | !file.exists(f.acc2alt) | !file.exists(f.alt2ref)) next
    
    print(p)      
    t.acc2ref = readRDS(f.acc2ref)
    t.acc2alt = readRDS(f.acc2alt)
    t.alt2ref = readRDS(f.alt2ref)
    
    base.len.ref = chr.len[(chr.len[,1] == base.acc.ref) & (chr.len[,2] == base.chr),3]
    base.len.alt = chr.len[(chr.len[,1] == base.acc.alt) & (chr.len[,2] == base.chr),3]
    
    
    corr.alt2ref = getCorresponding2BasePositions(t.alt2ref, base.len.ref)
    corr.acc2ref = getCorresponding2BasePositions(t.acc2ref, base.len.ref)
    corr.acc2alt = getCorresponding2BasePositions(t.acc2alt, base.len.alt)
    
    
    corr.alt2ref[corr.alt2ref == 0] = NA
    corr.acc2ref2 = corr.acc2alt[corr.alt2ref]
    
    x = cbind(corr.acc2ref, corr.acc2ref2)
    x[is.na(x)] = 0
    
    if(i.alt == 1){
      cmp.all = cbind(cmp.all, x)  
    } else {
      cmp.all = cbind(cmp.all, x[,2])
    }
    
    rm(t.alt2ref)
    rm(t.acc2ref)
    rm(t.acc2alt)
    rm(corr.alt2ref)
    rm(corr.acc2ref)
    rm(corr.acc2alt)
    gc()
    
    
  }
  
  
  if(is.null(cmp.all)){
    if(flag.for){
      next
    } else {
      return(NULL)
    }
  }
  
  # Find all confusing rows
  good.rows = rep(T, nrow(cmp.all))
  for(i in 1:ncol(cmp.all)){
    for(j in 1:ncol(cmp.all)){
      if(i >= j) next
      tmp = (rowSums(cmp.all[,c(i,j)] == 0) > 0) | (cmp.all[,i] == cmp.all[,j])
      good.rows = good.rows & tmp
    }
  }
  print(sum(!good.rows))
  
  confusing.values = unique(c(cmp.all[!good.rows,]))
  confusing.values = sort(setdiff(confusing.values, 0))    
  confusing.rows = F
  for(i in 1:ncol(cmp.all)){
    confusing.rows = confusing.rows | ((cmp.all[,i] %in% confusing.values))
  }
  sum(confusing.rows)
  
  # Mean value
  val = rowSums(cmp.all) / rowSums(cmp.all != 0)
  val[confusing.rows] = 0
  val[rowSums(cmp.all) == 0] = 0
  
  # Remove repeats
  val.tbl = table(val)
  vals.remove = as.numeric(names(val.tbl))[val.tbl != 1]
  val[val %in% vals.remove] = 0
  
  saveRDS(val, file.out, compress = F)
  
  rm(cmp.all)
  rm(val)
  gc()

}

chromosome.pairs = unique(chromosome.pairs[,c(2,3)])

flag.for = F
tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs))  %dopar% {  # which accession to use
# flag.for = T
# for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  query.chr = chromosome.pairs[i.chr.pair, 1]
  base.chr = chromosome.pairs[i.chr.pair, 2]
  
  val.all = c()
  val.names = c()
  for(acc in accessions){
    print(acc)
    
    file.out = paste(path.cmp, acc, '_', query.chr, '_', base.chr, '_cmp_to_ref_consensus.rds', sep = '')
    if (!file.exists(file.out)) {
      next
    } 
    
    val = readRDS(file.out)
    val.all = cbind(val.all, val)
    val.names = c(val.names, acc)
    rm(val)
  }
  
  if(length(val.names) == 0){
    if(flag.for){
      next
    } else {
      return(NULL)  
    }
  }
  colnames(val.all) = val.names
  val.all = cbind(1:nrow(val.all),val.all)
  colnames(val.all)[1] = base.acc.ref
  
  file.out.consensus = paste(path.work, 'consensus_', query.chr, '_', base.chr, '_', base.acc.ref, '.rds', sep = '')
  saveRDS(val.all, file.out.consensus, compress = F)
  rm(val.all)
  gc()
  
}





