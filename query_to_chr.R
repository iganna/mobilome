suppressMessages(library("optparse"))
suppressMessages(library(Biostrings))
suppressMessages(library('seqinr'))  # read.fasta
suppressMessages(library('foreach'))
suppressMessages(library(doParallel))


message('=====================================================')
message('           Genomes into chromosomes  ')
message('-----------------------------------------------------')

# Rscript query_to_chr.R -n 5 -t fasta --path.in ../pb_genomes/ --path.out ../pb_chromosomes/

# Rscript query_to_chr.R -n 8 -t fasta --path.in ../lyrata/ --path.out ../ly_chromosomes/    
# Rscript query_to_chr.R -n 1 -t fasta --path.in ../rhizobia/ --path.out ../rhiz_chromosomes/ -s T

myCluster <- makeCluster(30, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

# len.parts = 5000  # lengths of parts
# n.chr = 5  # number of chromosomes
# query.type = 'fasta'


args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-n", "--n.chr"), type="character", default=NULL, 
              help="number of chromosomes", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of fasta files", metavar="character"),
  make_option(c("-i", "--path.in"), type="character", default=NULL, 
              help="pathway to the input directory", metavar="character"),
  make_option(c("-o", "--path.out"), type="character", default=NULL, 
              help="pathway to the output directory", metavar="character"),
  make_option(c("-s", "--sort"), type="character", default=NULL, 
              help="sort chromosomes by lengths or not", metavar="character"),
  make_option(c("-b", "--bp"), type="character", default=NULL, 
              help="number of base pairs in the part file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)
# print(opt$s)
#stop()

if (!is.null(opt$bp)) len.parts <- as.numeric(opt$bp)
if (!is.null(opt$path.in)) path.query <- opt$path.in
if (!is.null(opt$path.out)) path.chr <- opt$path.out
if (!is.null(opt$type)) query.type <- opt$type
if (!is.null(opt$n.chr)) n.chr <- as.numeric(opt$n.chr)
if(!dir.exists(path.chr)) dir.create(path.chr)
sort.by.lengths = F
if (!is.null(opt$sort)) sort.by.lengths <- opt$sort
message(paste('sort_chr_by_length', sort.by.lengths, sep = ''))

# ------------------------------------------------------------------------------------
msg = '!!!  Please be sure that all chromosomes in files are sorted in the same order or use \"-s T\" flag  !!!'
# message(paste0(rep('-', nchar(msg)), collapse = ''))
message(msg)
# message(paste0(rep('-', nchar(msg)), collapse = ''))
# ------------------------------------------------------------------------------------

message(paste('Directory with genomes:', path.query))
files.query = list.files(path = path.query, pattern = paste0('\\.', query.type, '$', collapse = '') )
query.name = gsub(paste0('*.', query.type, collapse = ''), "" ,files.query)
message(paste0(c('Names of genomes:', query.name), collapse = ' '))

for.flag = F
tmp = foreach(acc = query.name, .packages=c('stringr','Biostrings', 'seqinr')) %dopar% {
#for.flag = T
#for(acc in query.name[1:length(query.name)]){
  
  print(acc)
  q.fasta = read.fasta(paste0(path.query, acc, '.', query.type, collapse = ''))
  
  if(length(q.fasta) < n.chr){
    return(NULL)
  }
  
  if(sort.by.lengths){
    q.len = sapply(q.fasta, length)
    q.fasta = q.fasta[order(-q.len)]
  }
  
  message(paste0(c('Chromosomes', names(q.fasta)[1:(n.chr)], 'will be processed'),
                 sep = ' ', collapse = ''))
  
  if(length(q.fasta) > n.chr){
    message(paste0(c('Chromosomes', names(q.fasta)[(n.chr+1):length(q.fasta)], 'will NOT be processed'), 
                   sep = ' ', collapse = ''))
  }
  
  for(i.chr in 1:n.chr){
    acc.s = gsub('_', '-', acc)
    file.out = paste0(path.chr, acc.s, '_chr', i.chr, '.fasta', collapse = '')
    if(file.exists(file.out)){
      if(for.flag){
        next
      } else {
        return(NULL)
      }
    }
    print(file.out)
    q.tmp = paste0(q.fasta[[i.chr]], collapse = '')
    
    s = toupper(q.tmp)
    
    
    write(paste('>', names(q.fasta)[i.chr] , sep=''), file=file.out, append=F)
    write(s, file=file.out, append=T)
    write('\n', file=file.out, append=T)

  }
}






