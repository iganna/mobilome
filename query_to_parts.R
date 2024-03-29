suppressMessages(library("optparse"))
suppressMessages(library(Biostrings))
suppressMessages(library('seqinr'))  # read.fasta
suppressMessages(library('foreach'))
suppressMessages(library(doParallel))

message('=====================================================')
message('             Chromosomes into parts  ')
message('-----------------------------------------------------') 

#Rscript query_to_parts.R -n 5 -t fasta --path.chr ../pb_chromosomes/ -b 5000 --path.parts ../pb_parts/
# Rscript query_to_parts.R -n 8 -t fasta --path.chr ../ly_chromosomes/ -b 5000 --path.parts ../ly_parts/
# Rscript query_to_parts.R -n 1 -t fasta --path.chr ../rhiz_chromosomes/ -b 5000 --path.parts ../rhiz_parts/ 

myCluster <- makeCluster(30, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

# len.parts = 5000  # lengths of parts
# n.chr = 5  # number of chromosomes
# query.type = 'fasta'
# path.chr = '../pb_databases/'
# path.parts = '../pb_parts/'

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-b", "--bp"), type="character", default=NULL, 
              help="number of base pairs in the part file", metavar="character"),
  make_option(c("-n", "--n.chr"), type="character", default=NULL, 
              help="number of chromosomes", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of fasta files", metavar="character"),
  make_option(c("-i", "--path.chr"), type="character", default=NULL, 
              help="pathway to the input directory", metavar="character"),
  make_option(c("-o", "--path.parts"), type="character", default=NULL, 
              help="pathway to the output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

if (!is.null(opt$bp)) len.parts <- as.numeric(opt$bp)
if (!is.null(opt$path.chr)) path.chr <- opt$path.chr
if (!is.null(opt$path.parts)) path.parts <- opt$path.parts
if (!is.null(opt$type)) query.type <- opt$type
if (!is.null(opt$n.chr)) n.chr <- as.numeric(opt$n.chr)

if(!dir.exists(path.parts)) dir.create(path.parts)

message(paste('Directory with chromosomes:', path.chr))
files.query = list.files(path = path.chr, pattern = paste0('\\.', query.type, '$', collapse = '') )
query.name = unique(sapply(files.query, function(s) strsplit(s, '_chr')[[1]][1]))
message(paste0(c('Names of genomes:', query.name), collapse = ' '))

if(length(query.name) == 0){
  stop('Wrong names of chromosomal files or files are not provided')
}

for.flag = F
tmp = foreach(acc = query.name, .packages=c('stringr','Biostrings', 'seqinr')) %dopar% {

# for.flag = T
# for(acc in query.name){
  
  for(i.chr in 1:n.chr){
    file.in = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')
    
    file.out = paste0(path.parts, acc, '_', i.chr, '.fasta', collapse = '')
    if( file.exists(file.out)) next
    
    q.fasta = read.fasta(file.in)[[1]]
    len.chr = length(q.fasta)
    q.tmp = paste0(q.fasta, collapse = '')
    pos.beg = seq(1, len.chr, len.parts)
    pos.end = c(seq(len.parts, len.chr, len.parts), len.chr)
    s = substring(q.tmp, pos.beg, pos.end)
    s = toupper(s)
    if(sum(sapply(s, nchar)) != len.chr) stop('aa')
    file.out = paste0(path.parts, acc, '_', i.chr, '.fasta', collapse = '')
    write('', file=file.out, append=F)
    for(i in 1:length(s)){
      write(paste0('>acc_', acc, '|chr_', i.chr, '|part_', i, '|', pos.beg[i], collapse=''), file=file.out, append=T)
      write(s[i], file=file.out, append=T)
      write('\n', file=file.out, append=T)
    }
  }
}






