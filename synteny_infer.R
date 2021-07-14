removeLastN <- function(x, n){
  substr(x, 1, nchar(x)-n)
}

removeFirstN <- function(x, n){
  substr(x, n+1, nchar(x))
}

remainLastN <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

remainFirstN <- function(x, n){
  substr(x, 1, n)
}

sumGapLastN <- function(x, n){
  s <- remainLastN(x, n)
  return(sum(strsplit(s, '')[[1]] == '-'))
}

sumGapFirstN <- function(x, n){
  s <- substr(x, 1, n)
  return(sum(strsplit(s, '')[[1]] == '-'))
}

checkCorrespToGenome <- function(t, base.fas.fw, base.fas.bw, query.fas, k = 10){
  for(irow in 1:nrow(t)) {
    s1 = toupper(remainLastN(t[irow, 'V8'], k))
    s1 = gsub("\\-","",s1)
    
    if(nchar(s1) == 0){
      s2 = s1
    } else {
      s2 = toupper(paste0(query.fas[(-(nchar(s1)-1):0) + t[irow, 'V3']], collapse = ''))
    }
    
    
    if(s1 != s2 ){
      print(irow)
      print(s1)
      print(s2)
      stop('Checking the query not passed')
    }
    
    s1 = toupper(remainLastN(t[irow, 'V9'], k))
    s1 = gsub("\\-","",s1)
    if(t$dir[irow] == 0) {
      base.fas = base.fas.fw
    } else {
      base.fas = base.fas.bw
    }
    
    if(nchar(s1) == 0){
      s2 = s1
    } else {
      s2 = toupper(paste0(base.fas[(-(nchar(s1)-1):0) + t[irow, 'V5']], collapse = ''))
    }
    
    
    if(s1 != s2 ){
      print(irow)
      print(s1)
      print(s2)
      stop('Checking the base was not passed')
    }
  }
}

fixDirection <- function(t, base.fas.bw){
  t$dir = c()
  for(irow in 1:nrow(t)) {
    if(t[irow,'V5'] < t[irow,'V4']) {
      
      pos1 = t[irow,'V5']
      pos2 = t[irow,'V4']
      
      t[irow,'V5'] <- length(base.fas.bw) - pos1 + 1
      t[irow,'V4'] <- length(base.fas.bw) - pos2 + 1
      t[irow, 'dir'] = 1
    } else {
      t[irow, 'dir'] = 0
    }
  }
  return(t)
}

showt <- function(t, irow=NULL, chr=F, add = c()){
  idx = c(1:5,7,11, add)
  idx = intersect(idx, 1:ncol(t))
  irow = irow[irow <= nrow(t)]
  
  if(chr) idx = c(idx, 10)
  if(is.null(irow)){
    print(t[1:min(100, nrow(t)),idx])
  } else {
    print(t[irow,idx])
  }
}

glueByThreshold <- function(t, thresholds, base.fas.fw, base.fas.bw, query.fas,
                            file.log = NULL, gap.open = 30, gap.ext = 0.5, maxval = 10^10,
                            base.overlap = T, query.overlap = T, log.append = F) {
  
  if(!is.null(file.log)) {
    write('', file=file.log, append=log.append)
  }
  t = orderT(t)
  for(threshold in thresholds) {
    
    irow = 1
    while(irow < nrow(t)) {
      # print(threshold)
      # print(irow)
      
      # if a base fragment of irow is already within some base fragment
      if(base.overlap){
        idx = (t[,'V4'] <= t[irow,'V4']) & (t[,'V5'] >= t[irow,'V5'])
        idx[t$dir != t$dir[irow]] = F
        idx[irow] = F
        if (sum(idx) > 0) {
          irow <- irow + 1
          next
        }        
      }

      
      # if a query fragment of irow is already within some query fragment
      if(query.overlap){
        idx = (t[,'V2'] <= t[irow,'V2']) & (t[,'V3'] >= t[irow,'V3'])
        idx[t$dir != t$dir[irow]] = F
        idx[irow] = F
        if (sum(idx) > 0) {
          irow <- irow + 1
          next
        }  
      }
      
      
      idx = abs(t[,'V2'] - (t[irow,'V3']))
      # idx = apply(cbind(abs(t[,'V2'] - (t[irow,'V3'])), abs(t[,'V4'] - (t[irow,'V5']))), 1, min)
      idx[irow] = maxval
      idx[t[,'V2'] <= t[irow,'V2']] = maxval
      idx[t[,'V4'] <= t[irow,'V4']] = maxval
      
      # Chromosomes should match
      idx[t[,'V10'] != t[irow,'V10']] = maxval
      
      # If new fragment is in irow fragment by query
      idx.inside = (t[,'V2'] >= t[irow,'V2']) & (t[,'V3'] <= t[irow,'V3'])
      idx[idx.inside] = maxval
      
      # If new fragment is in irow fragment by base
      idx.inside = (t[,'V4'] >= t[irow,'V4']) & (t[,'V5'] <= t[irow,'V5'])
      idx[idx.inside] = maxval
      
      idx[t$dir != t$dir[irow]] <- maxval
      
      idx[abs(t[, 'V2'] - t[irow, 'V3']) > threshold] = maxval
      idx[abs(t[, 'V4'] - t[irow, 'V5']) > threshold] = maxval
      
      if (min(idx) == maxval) {
        irow <- irow + 1
        next
      }
      
      idx = which(idx == min(idx))
      
      if (length(idx) == 0) {
        irow <- irow + 1
        next
      }
      idx = idx[1]
      
      if(threshold > 10^4){
        n.char = 50
        n.char.thresh = 100  
      } else {
        n.char = 9
        n.char.thresh = 50
      }
      
      n.char.thresh = min(min(n.char.thresh, (t[irow, 'V7'] - 2)), (t[idx, 'V7'] - 2) )
      
      d.tmp = 0
      while (d.tmp < n.char.thresh) {
        n.char = n.char + 1
        pos.query.start = max(1, t[irow, 'V3'] - n.char + sumGapLastN(t[irow, 'V8'], n.char))
        pos.query.end = t[idx, 'V2'] + n.char - sumGapFirstN(t[idx, 'V8'], n.char)
        
        
        pos.base.start = max(1, t[irow, 'V5'] - n.char + sumGapLastN(t[irow, 'V9'], n.char))
        pos.base.end = t[idx, 'V4'] + n.char - sumGapFirstN(t[idx, 'V9'], n.char)
        
        if(pos.base.end >= t[idx, 'V5']) break
        if(pos.query.end >= t[idx, 'V3']) break
        if(pos.base.start <= t[irow, 'V4']) break
        if(pos.query.end <= t[irow, 'V2']) break
        
        d.tmp = min(pos.query.end - pos.query.start, pos.base.end - pos.base.start)
      }
      
      if(t$dir[irow] == 0){
        base.fas = base.fas.fw
      } else {
        base.fas = base.fas.bw
      }
      
      
      seq.query <- toupper(query.fas.chr[(pos.query.start+1):(pos.query.end-1)])
      seq.base <- toupper(base.fas[(pos.base.start+1):(pos.base.end-1)])  
      
      globalAlign<- pairwiseAlignment(paste0(seq.query, collapse=''), 
                                      paste0(seq.base, collapse=''), gapOpening = gap.open, 
                                      gapExtension=gap.ext,
                                      type="global")
      globalAlign 
      
      if(!is.null(file.log)) {
        write(paste((pos.query.start+1), (pos.query.end-1), as.character(alignedPattern(globalAlign))),
              file=file.log, append=TRUE)
        write(paste((pos.base.start+1), (pos.base.end-1), as.character(alignedSubject(globalAlign))), 
              file=file.log, append=TRUE)
        write('\n', file=file.log, append=TRUE)
      }
      
      
      seq.query.new = paste(removeLastN(t[irow, 'V8'], n.char),  
                            alignedPattern(globalAlign), 
                            removeFirstN(t[idx, 'V8'], n.char), sep = '')
      
      seq.base.new = paste(removeLastN(t[irow, 'V9'], n.char),  
                           alignedSubject(globalAlign), 
                           removeFirstN(t[idx, 'V9'], n.char), sep = '')
      
      if(nchar(seq.query.new) != nchar(seq.base.new)) stop('aa')
      
      t1 <- t
      
      t[irow, 'V8'] <- seq.query.new
      t[irow, 'V9'] <- seq.base.new
      
      t[irow, 'V3'] <- t[idx, 'V3']
      t[irow, 'V5'] <- t[idx, 'V5']
      t[irow, 'V7'] <- nchar(seq.base.new)
      
      t <- t[-idx,]
      
      # s1 = seq.base.new
      # s2 = paste0(base.fas.fw[t[irow,'V4']:t[irow,'V5']], collapse = '')
      
      s.base = strsplit(seq.base.new,'')[[1]]
      idx.nogap.base = which(s.base != '-')
      if(length(idx.nogap.base) != abs(t[irow,'V5'] - t[irow,'V4']) + 1) stop('wrong in base')
      
      s.query = strsplit(seq.query.new,'')[[1]]
      idx.nogap.query = which(s.query != '-')
      if(length(idx.nogap.query) != abs(t[irow,'V3'] - t[irow,'V2']) + 1) stop('wrong in query')
      
      
      if(idx < irow) {
        irow = irow - 1
      }    
      
      rownames(t) <- NULL
    }
  }
  return(t)
}


glueZero <- function(t){
  t = t[order(t[,'V4']),]
  
  ipos = 1
  while(ipos < nrow(t)) {
    # print(ipos)
    idx = which((t[,'V4'] == (t[ipos,'V5'] + 1)) & (t[,'V2'] == (t[ipos,'V3'] + 1)))
    if (length(idx) == 0){
      ipos = ipos+1
      next
    }
    # if(idx != (ipos+1)) print(idx)  # glue with not the next record
    
    t[ipos, 'V3'] <- t[idx, 'V3']
    t[ipos, 'V5'] <- t[idx, 'V5']
    t[ipos, 'V7'] <- t[ipos, 'V7'] + t[idx, 'V7']
    t[ipos, 'V8'] <- paste(t[ipos, 'V8'], t[idx, 'V8'], sep = '')
    t[ipos, 'V9'] <- paste(t[ipos, 'V9'], t[idx, 'V9'], sep = '')
    t <- t[-idx,]
  }
  
  rownames(t) <- NULL
  return(t)
}

getCoverageBase <- function(t, len){
  idx = rep(F, len)
  for(irow in 1:nrow(t)){
    idx[t[irow,'V4']:t[irow,'V5']] = T
  }
  coverage = sum(idx) /  length(idx)
  print(coverage)
}

getCoverageQuery <- function(t, len){
  idx = rep(F, len)
  for(irow in 1:nrow(t)){
    idx[t[irow,'V2']:t[irow,'V3']] = T
  }
  coverage = sum(idx) /  length(idx)
  print(coverage)
}


removeShortOverlaps <- function(t, echo=F, len.min = 200, gap.thresh = 10000,
                                filter.q = T, filter.b = T){
  t <- t[order(t[,'V2']),]

  t.base = t
  tmp  = base.len - t.base[t.base$dir == 1,'V5'] + 1
  t.base[t.base$dir == 1,'V5'] = base.len - t.base[t.base$dir == 1,'V4'] + 1
  t.base[t.base$dir == 1,'V4'] = tmp

  irow = 1
  while(irow < nrow(t)) {
    
    # if(t$dir[irow] != t$dir[irow+1]) {
    #   irow = irow + 1
    #   next
    # }
    
    gap.q = t[irow + 1,'V2'] - t[irow,'V3'] - 1
    gap.b = t[irow + 1,'V4'] - t[irow,'V5'] - 1
    if(echo) print(c(irow, 
                     t.base[irow + 1,'V2'] - t.base[irow,'V3'] - 1, 
                     t.base[irow + 1,'V4'] - t.base[irow,'V5'] - 1, 
                     t$dir[irow], t$dir[irow+1]))
    
    
    # if((gap.q < 0) && (gap.b < 0) && (gap.b < gap.q)) gap.q = 10
    

    if(filter.q & (gap.q < 0) & (gap.q > -gap.thresh) & 
       (t[irow, 'V2'] < t[irow+1, 'V2']) & (t[irow, 'V3'] < t[irow+1, 'V3'])){
      # overlap from one side
      s.q1 = strsplit(t[irow, 'V8'], '')[[1]]
      s.b1 = strsplit(t[irow, 'V9'], '')[[1]]
      s.i1 = rep(0, t[irow, 'V7'])
      s.i1[s.q1 != '-'] <- t[irow,'V2']:t[irow,'V3']
      s.q1 = s.q1[min(which(s.i1 >= (t[irow,'V3'] + gap.q + 1))) :t[irow, 'V7']]
      s.b1 = s.b1[min(which(s.i1 >= (t[irow,'V3'] + gap.q + 1))) :t[irow, 'V7']]
      
      
      s.q2 = strsplit(t[irow+1, 'V8'], '')[[1]]
      s.b2 = strsplit(t[irow+1, 'V9'], '')[[1]]
      s.i2 = rep(Inf, t[irow+1, 'V7'])
      s.i2[s.q2 != '-'] <- t[irow+1,'V2']:t[irow+1,'V3']
      
      s.q2 = s.q2[1:max(which(s.i2 <= (t[irow+1,'V2'] - gap.q - 1))) ]
      s.b2 = s.b2[1:max(which(s.i2 <= (t[irow+1,'V2'] - gap.q - 1)))]
      
      
      if(sum(c(s.q2, s.b2) == '-') > sum(c(s.q1, s.b1) == '-') ) {
        tmp = 1
      } else if(sum(c(s.q2, s.b2) == '-') < sum(c(s.q1, s.b1) == '-') ) {
        tmp = 2
      } else if (sum(s.q1 != s.b1) > sum(s.q2 != s.b2)){
        tmp = 2 
      } else {
        tmp = 1
      }
    } else if (filter.b & (gap.b < 0) & (gap.b > -gap.thresh) & 
               (t[irow, 'V4'] < t[irow+1, 'V4']) & (t[irow, 'V5'] < t[irow+1, 'V5'])) {
      
      # stop('a')
      # overlap from one side
      s.q1 = strsplit(t[irow, 'V8'], '')[[1]]
      s.b1 = strsplit(t[irow, 'V9'], '')[[1]]
      s.i1 = rep(0, t[irow, 'V7'])
      if(length(s.i1[s.b1 != '-']) != length(t[irow,'V4']:t[irow,'V5'])) stop('aa')
      s.i1[s.b1 != '-'] <- t[irow,'V4']:t[irow,'V5']
      s.q1 = s.q1[min(which(s.i1 >= (t[irow,'V5'] + gap.b + 1))) :t[irow, 'V7']]
      s.b1 = s.b1[min(which(s.i1 >= (t[irow,'V5'] + gap.b + 1))) :t[irow, 'V7']]
      
      
      s.q2 = strsplit(t[irow+1, 'V8'], '')[[1]]
      s.b2 = strsplit(t[irow+1, 'V9'], '')[[1]]
      s.i2 = rep(Inf, t[irow+1, 'V7'])
      s.i2[s.b2 != '-'] <- t[irow+1,'V4']:t[irow+1,'V5']
      s.q2 = s.q2[1:max(which(s.i2 <= (t[irow+1,'V4'] - gap.b - 1))) ]
      s.b2 = s.b2[1:max(which(s.i2 <= (t[irow+1,'V4'] - gap.b - 1)))]
      
      if(sum(c(s.q2, s.b2) == '-') > sum(c(s.q1, s.b1) == '-') ) {
        tmp = 1
      } else if(sum(c(s.q2, s.b2) == '-') < sum(c(s.q1, s.b1) == '-') ) {
        tmp = 2
      } else if (sum(s.q1 != s.b1) > sum(s.q2 != s.b2)){
        tmp = 2 
      } else {
        tmp = 1
      }
    } else {
      tmp = 0
    }
    if(tmp == 1){  # remain a part at the first piece
      # remove a part from the second piece
      t[irow+1, 'V2'] <- t[irow+1, 'V2'] + sum(s.q2 != '-')
      t[irow+1, 'V4'] <- t[irow+1, 'V4'] + sum(s.b2 != '-')
      
      t[irow+1, 'V8'] <- removeFirstN(t[irow+1, 'V8'], length(s.q2))
      t[irow+1, 'V9'] <- removeFirstN(t[irow+1, 'V9'], length(s.b2))
      
      t[irow+1, 'V7'] <- nchar(t[irow+1, 'V8'])
      
    } else if(tmp == 2){ # remain a part at the second piece
      # remove a part from the first piece
      t[irow, 'V3'] <- t[irow, 'V3'] - sum(s.q1 != '-')
      t[irow, 'V5'] <- t[irow, 'V5'] - sum(s.b1 != '-')
      
      t[irow, 'V8'] <- removeLastN(t[irow, 'V8'], length(s.q1))
      t[irow, 'V9'] <- removeLastN(t[irow, 'V9'], length(s.b1))
      t[irow, 'V7'] <- nchar(t[irow, 'V8'])
    }
    irow = irow + 1
  }
  t = t[t[,'V7'] > len.min,]
  return(t)
}

additionalLocalAlignments <- function(t, query.fas.chr, base.fas.fw, base.fas.bw, echo=F, 
                                      n.short=200, gap.max = 15000, gap.open = 30, gap.ext = 0.5, 
                                      file.log = NULL){
  if(!is.null(file.log)) {
    write('', file=file.log, append=F)
  }
  
  t <- t[order(t[,'V2']),]
  t.additional <- c()
  irow = 0
  while(irow < (nrow(t) - 2)) {
    irow = irow+1
    if (t[irow,'dir'] != t[irow+1,'dir']) next
    
    gap.q = t[irow + 1,'V2'] - t[irow,'V3'] - 1
    gap.b = t[irow + 1,'V4'] - t[irow,'V5'] - 1
    
    
    
    if ((gap.q <= 0) || (gap.b <= 0)) next
    if ((gap.q/gap.b > 10) ||(gap.b/gap.q > 10)) next
    if (max(gap.q, gap.b) > gap.max) next
    
    
    s.q = query.fas.chr[(t[irow,'V3']+1): (t[irow + 1,'V2']-1)]
    if(t$dir[irow] == 0) {
      s.b = base.fas.fw[(t[irow,'V5']+1): (t[irow + 1,'V4']-1)]
    } else {
      s.b = base.fas.bw[(t[irow,'V5']+1): (t[irow + 1,'V4']-1)]
    }
    
    s.q <- paste0(s.q, collapse='')
    s.b <- paste0(s.b, collapse='')
    
    localAlign <- pairwiseAlignment(s.q, s.b, 
                                    gapOpening = gap.open, gapExtension=gap.ext,
                                    type="local")
    
    if(!is.null(file.log)) {
      write(paste((t[irow,'V3']+1), (t[irow,'V5']+1), as.character(alignedPattern(localAlign))),
            file=file.log, append=TRUE)
      write(paste((t[irow,'V3']+1), (t[irow,'V5']+1), as.character(alignedSubject(localAlign))), 
            file=file.log, append=TRUE)
      write('\n', file=file.log, append=TRUE)
    }
    
    # localAlign <- pairwiseAlignment(s.q, s.b, type="global")
    
    s.q.range <- localAlign@pattern@range
    s.b.range <- localAlign@subject@range
    
    if(nchar(localAlign@pattern) < n.short) next
    
    if(echo) print(c(irow, gap.q, gap.b, nchar(localAlign@pattern)))
    
    t.tmp = data.frame(V1=paste('add', t[irow,'V1'], sep = '_'),
                       V2=t[irow,'V3'] + s.q.range@start,  # You do not need +1 or -1 here
                       V3=t[irow,'V3'] + s.q.range@start + s.q.range@width - 1,
                       V4=t[irow,'V5'] + s.b.range@start,
                       V5=t[irow,'V5'] + s.b.range@start + s.b.range@width - 1,
                       V6=100, V7=nchar(localAlign@pattern),
                       V8=as.character(alignedPattern(localAlign)), 
                       V9=as.character(alignedSubject(localAlign)), 
                       V10=t[irow,10], dir=t[irow,'dir'])
    t.additional <- rbind(t.additional, t.tmp)
  }
  
  t.additional[, 'V9'] <- as.character(t.additional[, 'V9'])
  
  t <- rbind(t, t.additional)
  
  checkCorrespToGenome(t, query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  t <- t[order(t[,'V2']),]
  
  return(t)
}

removeCompleteOverlaps <- function(t, n.distant = 1000000, n.short = 200, 
                                   filter.query = T, filter.base = T,
                                   centromere.pos = NULL){
  rownames(t) <- NULL
  
  # Get initial positions in base sequence
  t.base = t
  tmp  = base.len - t.base[t.base$dir == 1,'V5'] + 1
  t.base[t.base$dir == 1,'V5'] = base.len - t.base[t.base$dir == 1,'V4'] + 1
  t.base[t.base$dir == 1,'V4'] = tmp
  
  t.base = t.base[rev(order(t.base['V7'])),]
  
  # ------------------------------

  
  # Overlap in query
  if(filter.query){
    irow = 1
    while(irow < nrow(t.base)) {
      # print(nrow(t.base))
      idx = which((t.base[,'V2'] >= t.base[irow,'V2']) & 
                    (t.base[,'V3'] <= t.base[irow,'V3']) )#& (t.base$dir == t.base$dir[irow]))
      idx = setdiff(idx, irow)
      idx = idx[t.base[irow,'V7'] > t.base[idx,'V7']]
      if(length(idx) >= 1) {
        t.base = t.base[-idx,] 
      } else {
        irow = irow + 1
      }
    }    
  }

  
  # Overlap in base
  if(filter.base) {
    irow = 1
    while(irow < nrow(t.base)) {
      # print(nrow(t.base))
      idx = which((t.base[,'V4'] >= t.base[irow,'V4']) & 
                    (t.base[,'V5'] <= t.base[irow,'V5']) )# & (t.base$dir == t.base$dir[irow]))
      idx = setdiff(idx, irow)
      idx = idx[t.base[irow,'V7'] > t.base[idx,'V7']]
      if(length(idx) >= 1) {
        t.base = t.base[-idx,] 
      } else {
        irow = irow + 1
      }
    }    
  }

  
  # idx = rep(F, base.len)
  # for(irow in 1:nrow(t.base)){
  #   idx[t.base[irow,'V4']:t.base[irow,'V5']] = T
  # }
  # 
  # coverage = sum(idx) / base.len
  # print(coverage)
  
  # idx = rep(F, length(query.fas.chr))
  # for(irow in 1:nrow(t.base)){
  #   idx[t.base[irow,'V2']:t.base[irow,'V3']] = T
  # }
  # 
  # coverage = sum(idx) / length(idx)
  # print(coverage)
  
  # if(!is.null(n.short)){
  #   t.base <- t[t.base[,'V7'] > n.short,]    
  # }
  
  # # remove distant outliers
  # if(!is.null(n.distant)){
  #   if(is.null(centromere.pos)){
  #     idx.distant <- c()
  #     for(irow in 1:nrow(t.base)){
  #       if(abs(t.base[irow,'V2'] - t.base[irow,'V4']) > n.distant){
  #         idx.distant <- c(idx.distant, irow)
  #       }
  #     }
  #     if(length(idx.distant) > 1)  t.base <- t.base[-idx.distant,]  
  #   } else {
  # 
  #     t.base = t.base[order(t.base[,'V2']),]
  #     idx = min(which(!((t.base[,'V2'] <= centromere.pos[1]) & (t.base[,'V3'] <= centromere.pos[2])) )) - 1
  #     while(t.base[idx,'V7'] > 10000) idx = idx + 1
  #     cent.start.idx = idx
  #     
  #     idx = max(which(!((t.base[,'V2'] >= centromere.pos[1]) & 
  #                         (t.base[,'V3'] >= centromere.pos[2])) )) + 1
  #     while(t.base[idx,'V7'] > 10000) idx = idx - 1
  #     cent.end.idx = idx
  #     
  #     t.arms = list( t.base[1:(cent.start.idx-1),], 
  #                    t.base[(cent.end.idx+1):nrow(t.base),])
  #     
  #     for(i.tmp in 1:2){
  #       #fit linear model
  #       # lm.res <- summary(lm('V2 ~ offset(1*V4)', t.arms[[i.tmp]][(t.arms[[i.tmp]][,'V7'] > 10^4)&
  #       #                                                     (t.arms[[i.tmp]]$dir == 0),]))
  #       # intercept <- lm.res$coefficients[1,1] + 1.5 * lm.res$coefficients[1,2]
  #       
  #       y = t.arms[[i.tmp]][(t.arms[[i.tmp]][,'V7'] > 10^4)&
  #                             +                         (t.arms[[i.tmp]]$dir == 0),]
  #       intersept <- max(abs(y[,'V2'] - y[,'V4']))
  #       
  #       idx.distant <- c()
  #       for(irow in 1:nrow(t.arms[[i.tmp]])){
  #         if(abs(t.arms[[i.tmp]][irow,'V2'] - t.arms[[i.tmp]][irow,'V4']) > intercept + n.distant){
  #           idx.distant <- c(idx.distant, irow)
  #         }
  #       }
  #       if(length(idx.distant) > 1)  t.arms[[i.tmp]] <- t.arms[[i.tmp]][-idx.distant,]  
  #     }
  #     t.base = t.base[sort(c(rownames(t.arms[[1]]),
  #                       rownames(t.arms[[2]]),
  #                       rownames(t.base[cent.start.idx:cent.end.idx,]))),]
  #     
  #   }
  #   
  # }
  
  t = t[rownames(t.base),]
  
  # plot(t.base[,'V4'], t.base[,'V2'])


  
  # plot(t[,'V4'], t[,'V2'])
  
  t <- t[order(t[,'V2']),]
  rownames(t) <- NULL
  
  return(t)
}


setDir <- function(t, base.len){
  # direction
  t$dir = c()
  idx.dir = t[,'V5'] < t[,'V4']
  pos1 = base.len - t[idx.dir,'V5'] + 1
  pos2 = base.len - t[idx.dir,'V4'] + 1
  
  t[idx.dir,'V5'] <- base.len - t[idx.dir,'V5'] + 1
  t[idx.dir,'V4'] <- base.len - t[idx.dir,'V4'] + 1
  t[, 'dir'] = 0
  t[idx.dir, 'dir'] = 1
  
  return(t)
}



getT <- function(t.file, query.fas.chr, base.fas.fw, base.fas.bw,
                 thresholds = c(100, 500, 1000, 1500, 2000), echo=T){
  # ------- Read blast results -------
  
  if(echo) message(paste0(c('Reading', t.file, '...'), collapse = ' '))
  base.len = length(base.fas.fw)
  t = read.table(t.file, stringsAsFactors = F, header = F)
  
  t <- setDir(t, base.len)
  
  # Get right positions
  start.pos = as.numeric(sapply(strsplit(t[,1], "\\|"), "[", 4)) - 1
  t[,2:3] = t[,2:3] + start.pos
  
  rownames(t) <- NULL
  
  # check
  for(irow in 1:nrow(t)) {
    if(nchar(t[irow,'V8']) != nchar(t[irow,'V9'])) stop('aaa')
  }
  
  source("/Users/anna/OneDrive/pushkin/cryptic/nanopore/scripts/synteny_infer.R")
  checkCorrespToGenome(t, query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  
  
  # ------- Grow up the alignment -------
  
  ## ---- First glue ----
  
  if(echo) message('First glue (exact)...')
  t <- glueZero(t)
  
  ## ---- Pairwise global alignments of gaps ----
  
  if(echo) message('Glue with thresholds')
  t = glueByThreshold(t, thresholds, query.fas = query.fas.chr, 
                      base.fas.fw = base.fas.fw, 
                      base.fas.bw = base.fas.bw, file.log=file.log)
  
  ## ---- Remove complete overlap from both sides (base an query) ----
  
  if(echo) message('Remove complete overlaps')
  t <- removeCompleteOverlaps(t)
  
  ## ---- Remove short overlaps and Local alignments ----
  
  # Remove short overlaps
  if(echo) message('Remove short overlaps')
  t <- removeShortOverlaps(t, echo = F)
  
  # Additional alignments
  if(echo) message('Additional alignments')
  nrow.tmp = nrow(t) - 1
  while(nrow(t) != nrow.tmp) {
    nrow.tmp = nrow(t)
    print(nrow.tmp)
    t <- additionalLocalAlignments(t, query.fas.chr, base.fas.fw, base.fas.bw, echo = F, n.short=100,
                                   file.log = file.log)
    t <- glueZero(t)
    
    checkCorrespToGenome(t, query.fas = query.fas.chr, 
                         base.fas.fw = base.fas.fw, 
                         base.fas.bw = base.fas.bw)
    
    t = glueByThreshold(t, thresholds, query.fas = query.fas.chr, 
                        base.fas.fw = base.fas.fw, 
                        base.fas.bw = base.fas.bw, file.log=file.log)
  }
  
  return(t)
}


plotSyntenyBlocks <-function(t, base.len, idx=NULL){
  df = c()
  for(i in 1:nrow(t)) {
    if(t$dir[i] == 0) {
      df = rbind(df, c(t[i, 'V2'], t[i, 'V4'], i, 0))
      df = rbind(df, c(t[i, 'V3'], t[i, 'V5'], i, 0))
    } else {
      df = rbind(df, c(t[i, 'V2'], base.len - t[i, 'V4'] + 1, i, 1))
      df = rbind(df, c(t[i, 'V3'], base.len - t[i, 'V5'] + 1, i, 1))
    }
  }
  
  df = as.data.frame(df)
  df$clr <- df$V3 %% 2;
  
  if(!is.null(idx)) df$V4[c(idx*2, (idx*2-1))] = 2
  
  p <- ggplot(df, aes(x = V1, y=V2, color = as.factor(V4), group=as.factor(V3)  )) + 
    geom_line(show.legend = FALSE) + theme_bw()
  return(p)
}


orderT <- function(t, order = 'query'){
  if(order == 'query'){
    t = t[order(t[,'V2']),]
  } else if(order == 'base'){
    t = t[order(t[,'V4']),]
  } else {
    message('Not ordered!')
  }
  rownames(t) <- NULL
  return(t)
}

getBlastStat <- function(t.reb, cover.cutoff = 0.9, 
                         flag.max.cover = T,
                         flag.tot.cover = T,
                         flag.pers.cover = T,
                         flag.max.cover.one=F,
                         suff.res = '') {
  
  reb.unique = unique(t.reb[, 'V1'])
  t.reb = t.reb[!duplicated(t.reb[,1:10]),]

  reb.res = c()
  names.res = c()
  
  for(i.reb in 1:length(reb.unique)){
    reb.name = reb.unique[i.reb]
    
    tmp.res = c()
    names.res = c()
    
    len.tmp = as.numeric(strsplit(reb.name, '\\|')[[1]][4])
    n.names.init = length(names.res)
    
    # if record is here
    
    t.tmp = t.reb[t.reb[,'V1'] == reb.name,, drop=F]
    
    
    if(flag.max.cover){
      idx = which(t.tmp[,'V7'] >= cover.cutoff * len.tmp)
      tmp = length(idx)
      if(tmp < 1) {
        tmp = max(t.tmp[,'V7']) / len.tmp
      }
      tmp.res = c(tmp.res, tmp)
      names.res = c(names.res, 'max_coverage')
    }
    
    if(flag.max.cover.one){
      

      tmp = max(t.tmp[,'V3'] - t.tmp[,'V2'] + 1) / len.tmp
      tmp.res = c(tmp.res, tmp)
      names.res = c(names.res, 'max_coverage_one')
    }
    
    if(flag.tot.cover){
      t.coverage = rep(0, len.tmp)
      for(i.tmp in 1:nrow(t.tmp)) {
        t.coverage[t.tmp[i.tmp, 'V2']: t.tmp[i.tmp, 'V3']] = 1
      }
      t.coverage = sum(t.coverage) / len.tmp 
      tmp.res = c(tmp.res, t.coverage)
      names.res = c(names.res, 'tot_coverage')
    }
    
    if(flag.pers.cover){
      t7 = t.tmp[,'V7'] / len.tmp
      t7[t7 > 1] = 1
      tmp = table(c(floor(t7*10), 0:10)) - 1  
      tmp.res = c(tmp.res, rev(tmp))
      names.res = c(names.res, paste('pers', rev(0:10)/10, sep = '_'))
    }
    
    
    if(nchar(suff.res) > 0){
      names.res[(n.names.init+1):length(names.res)] <- 
        paste(names.res[(n.names.init+1):length(names.res)], suff.res, sep = '_')
    }
    
    names(tmp.res) <- names.res
    reb.res = rbind(reb.res, tmp.res)
  }
  rownames(reb.res) <- reb.unique
  
  return(reb.res)
}


# ignore sequences of this size
# thresh - threshold fo remove reduldant sequences
getMobilome <- function(t, base.len, allowed_gap = 5, m.min.length = 15, 
                        thresh = 0.8, m.max.length = Inf, echo=F){
  

  m1 = c()
  m2 = c()
  pos.m1 = c()
  pos.m2 = c()
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    s1 = strsplit(t[irow,'V8'], '')[[1]]
    s2 = strsplit(t[irow,'V9'], '')[[1]]
    
    pos = rep(0, length(s1))
    pos1 = pos; pos1[s1 != '-'] = t[irow,'V2']: t[irow,'V3']
    pos2 = pos; pos2[s2 != '-'] = t[irow,'V4']: t[irow,'V5']
    
    icol = 0
    while(icol < length(s1)){
      icol = icol + 1
      
      if(s1[icol] == '-') {  # start gaps on query
        i.start = icol
        s = c()
        while((sum(s1[icol:min((icol+allowed_gap), length(s1))] == '-') > 0) && (icol <= length(s1))){
          s = c(s, s2[icol])
          icol = icol + 1
        }
        i.end = icol - 1
        
        if((i.end - i.start + 1) < m.min.length) next
        if((i.end - i.start + 1) > m.max.length) next
        
        # Nucleotide frequences
        f = table(s)
        if(max(f) > thresh * length(s)) next
        
        # Di-Nucleotide frequences
        cf = apply(combn(f, 2), 2, sum)
        if(max(cf) > thresh * length(s)) next
        
        if(t[irow, 'dir'] == 1){
          pos2.start <- base.len - pos2[i.end] + 1
          pos2.end <- base.len - pos2[i.start] + 1
        } else {
          pos2.start <- pos2[i.start]
          pos2.end <- pos2[i.end]
        }
        
        if(i.start == 1){
          # print(dim(m2))
          if(length(c(pos2.start, pos2.end, 0, length(s), paste0(s, collapse = ''))) != 5) stop('ddd')
          m2 = rbind(m2, c(pos2.start, pos2.end, 0, length(s), paste0(s, collapse = ''))) 
        } else {
          # print(dim(m2))
          if(length(c(pos2.start, pos2.end, pos1[i.start-1], length(s), paste0(s, collapse = ''))) != 5) stop('eee')
          m2 = rbind(m2, c(pos2.start, pos2.end, pos1[i.start-1], length(s), paste0(s, collapse = ''))) 
        }
        
      } else if(s2[icol] == '-') {  # start gaps on query
        i.start = icol
        s = c()
        while((sum(s2[icol:min((icol+allowed_gap), length(s2))] == '-') > 0) && (icol <= length(s2))){
          s = c(s, s1[icol])
          icol = icol + 1
        }
        i.end = icol - 1
        if((i.end - i.start + 1) < m.min.length) next
        if((i.end - i.start + 1) > m.max.length) next
        
        # Nucleotide frequences
        f = table(s)
        if(max(f) > thresh * length(s)) next
        
        # Di-Nucleotide frequences
        cf = apply(combn(f, 2), 2, sum)
        if(max(cf) > thresh * length(s)) next
        
        if(i.start == 1){
          # print(dim(m1))
          if(length(c(pos1[i.start], pos1[i.end], 0, length(s), paste0(s, collapse = ''))) != 5) stop('bb')
          m1 = rbind(m1, c(pos1[i.start], pos1[i.end], 0, length(s), paste0(s, collapse = '')))  
        } else {
          
          if(t[irow, 'dir'] == 1){
            pos2.start <- base.len - pos2[i.start-1] + 1
          } else {
            pos2.start = pos2[i.start-1]
          }
          # print(dim(m1))
          if(length(c(pos1[i.start], pos1[i.end], pos2.start, length(s), paste0(s, collapse = ''))) != 5) stop('aa')
          m1 = rbind(m1, c(pos1[i.start], pos1[i.end], pos2.start, length(s), paste0(s, collapse = '')))
        }
        
      }
    
      
      if(sum(is.na(m1)) >0) stop('NA1')
      if(sum(is.na(m2)) >0) stop('NA2')
    }
  }
  return(list(query = m1, base = m2))
}

addRowsByNames <- function(reb, s.fasta){
  
  init.names = rownames(reb)
  lost.names = setdiff(names(s.fasta), init.names)
  for(name in lost.names){
    reb = rbind(reb, rep(0, ncol(reb)))
  }
  rownames(reb) <- c(init.names, lost.names)
  reb <- reb[names(s.fasta),, drop=F]
  return(reb)
  
}

showd <- function(t){
  t.base
  links.q = t[2:(nrow(t)),'V2'] - t[1:(nrow(t)-1),'V3'] - 1
  links.b = t[2:(nrow(t)),'V4'] - t[1:(nrow(t)-1),'V5'] - 1
  print(cbind(links.q, links.b, t$dir[2:(nrow(t))], t$dir[1:(nrow(t)-1)]))
}


njColumns <- function(mx){
  n.tmp = ncol(mx)
  d = matrix(0, nrow = n.tmp, ncol = n.tmp, 
             dimnames = list(colnames(mx), colnames(mx)))
  for(i in 1:n.tmp){
    for(j in 1:n.tmp){
      d[i,j] = sum(mx[,i] != mx[,j])
    }
  }
  d1 <- as.dist(d)
  tree <- nj(d1)
  return(tree)
}

hclustColumns <- function(mx){
  n.tmp = ncol(mx)
  d = matrix(0, nrow = n.tmp, ncol = n.tmp, 
             dimnames = list(colnames(mx), colnames(mx)))
  for(i in 1:n.tmp){
    for(j in 1:n.tmp){
      d[i,j] = sum(mx[,i] != mx[,j])
    }
  }
  d1 <- as.dist(d)
  hc <- hclust(d1)
  return(hc)
}

s2blocks <- function(s, block.len, echo=F){
  blocks = c()
  positions = c()
  block.n = floor(length(s) / block.len)
  if(block.n == 0){
    block.n = 1
  }
  for(i in 1:block.n){
    pos.s = 1 + (i-1)*block.len
    if(i == block.n){
      pos.e = length(s)
    } else {
      pos.e = i*block.len
    }
    if(echo) print(c(pos.s, pos.e, pos.e - pos.s + 1))
    positions = rbind(positions, c(pos.s, pos.e))
    blocks = c(blocks, paste0(s[pos.s:pos.e], collapse = ''))
  }
  return(list(b=blocks, p=positions))
}


alignBlocks <- function(t.arm, irow.q, irow.b, 
                        query.fas.chr, base.fas.fw,
                        base.fas.bw, block.len, gap.open=10, gap.ext=0.5,
                        echo = F, diag.cutoff = 3, rev.compl=F){
  
  t.arm.base = t.arm
  tmp  = base.len - t.arm.base[t.arm.base$dir == 1,'V5'] + 1
  t.arm.base[t.arm.base$dir == 1,'V5'] = base.len - t.arm.base[t.arm.base$dir == 1,'V4'] + 1
  t.arm.base[t.arm.base$dir == 1,'V4'] = tmp
  
  
  if(echo) print(c(irow.q, irow.b))
  s.q = query.fas.chr[(t.arm[irow.q,'V3']+1) : (t.arm[irow.q+1,'V2']-1)]
  # if(t.arm$dir[irow.b] == 0){
    s.b = base.fas.fw[(t.arm.base[irow.b,'V5']+1) : (t.arm.base[irow.b+1,'V4']-1)]
  # } else {
    # s.b = base.fas.bw[(t.arm[irow.b,'V5']+1) : (t.arm[irow.b+1,'V4']-1)]
  # }
  
  if(rev.compl){
    s.b = reverseComplement(s.b)
  }
  
  # Split into blocks
  blocks.q <- s2blocks(s.q, block.len)
  blocks.b <- s2blocks(s.b, block.len)
  
  # Alignments
  t.tmp = c()
  mx = matrix(0, nrow = length(blocks.q$b), ncol = length(blocks.b$b))
  for(i.q in 1:length(blocks.q$b)){
    for(i.b in 1:length(blocks.b$b)){
      
      
      if(str_count(blocks.q$b[i.q], "n") / nchar(str_count(blocks.q$b[i.q], "n")) > 0.1) next
      if(str_count(blocks.b$b[i.b], "n") / nchar(str_count(blocks.b$b[i.b], "n")) > 0.1) next
      
      if(length(blocks.q$b) >= length(blocks.b$b)){
        #if((i.q - i.b) > diag.cutoff) next
        #if(-(i.q - i.b) > length(blocks.q$b) + diag.cutoff) next
      }
      # print(c(i.q, i.b))
      
      if(length(blocks.q$b) <= length(blocks.b$b)){
        #if((i.b - i.q) > diag.cutoff) next
        #if((i.b - i.q) > length(blocks.b$b) + diag.cutoff) next
      }
      
      
      # print(c(i.q, i.b))
      
      localAlign <- pairwiseAlignment(blocks.q$b[i.q], blocks.b$b[i.b], 
                                      gapOpening = gap.open, gapExtension=gap.ext,
                                      type="local")
      pos.s.q = localAlign@pattern@range@start + blocks.q$p[i.q] - 1
      pos.e.q = localAlign@pattern@range@start + blocks.q$p[i.q] - 1 +
        localAlign@pattern@range@width - 1
      
      pos.s.b = localAlign@subject@range@start + blocks.b$p[i.b] - 1
      pos.e.b = localAlign@subject@range@start + blocks.b$p[i.b] - 1 +
        localAlign@subject@range@width - 1
      
      s.q.aln <- as.character(alignedPattern(localAlign))
      s.b.aln <- as.character(alignedSubject(localAlign))
      
      if(!rev.compl){
        dir.tmp = 0
      } else {
        dir.tmp = 1
      }
      
      t.tmp = rbind(t.tmp, 
                    data.frame(V1=paste0(c('betw', i.q, i.b), collapse = '|'),
                               V2=pos.s.q, V3=pos.e.q, V4=pos.s.b, V5=pos.e.b, V6=100, V7=nchar(s.q.aln), 
                               V8=s.q.aln, V9=s.b.aln, V10=t.arm[irow.q,'V10'], 
                               dir=dir.tmp, stringsAsFactors = FALSE))
      
      mx[i.q, i.b] = localAlign@pattern@range@width / nchar(blocks.q$b[i.q])
    }
    
  }
  t.tmp1 = t.tmp
  
  if(is.null(t.tmp)) return(list(t = t.tmp, mx = NULL))
  
  # Visializations
  
  # df <- melt(mx)
  # ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  #   geom_tile() + theme_minimal()
  
  
  # Glue up things
  t.tmp = t.tmp1
  
  t.tmp[,c('V2','V3')] <- t.tmp[,c('V2','V3')] + t.arm[irow.q,'V3']
  if(rev.compl == FALSE){
    t.tmp[,c('V4','V5')] <- t.tmp[,c('V4','V5')] + t.arm.base[irow.b,'V5']
  } else {
    t.tmp[,c('V4','V5')] <- t.tmp[,c('V4','V5')]  + base.len - t.arm.base[irow.b,'V5'] - length(s.b)
  }
  
  
  if(nrow(t.tmp) > 1) t.tmp = t.tmp[order(t.tmp[,'V2']),]
  t.tmp <- t.tmp[t.tmp[,'V7'] > 30,]
  # print(nrow(t.tmp))
  # showt(t.tmp)
  if(nrow(t.tmp) == 0)  return(list(t = NULL, mx = mx))
  if(nrow(t.tmp) > 1) {
    t.tmp <- glueZero(t.tmp)
    t.tmp = t.tmp[order(t.tmp[,'V2']),]
    t.tmp = glueByThreshold(t.tmp, block.len/2, query.fas = query.fas.chr, 
                            base.fas.fw = base.fas.fw, 
                            base.fas.bw = base.fas.bw, gap.open = gap.open)
    t.tmp <- removeCompleteOverlaps(t.tmp)
  }
  # showt(t.tmp)
  
  
  # plotSyntenyBlocks(t.tmp)
  

  
  checkCorrespToGenome(t.tmp, query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  return(list(t = t.tmp, mx = mx))
}


isIntersect <- function(x1, x2, y1, y2){
  return(!(( x2 < y1 ) || (x1 > y2)))
}


nIntersect <- function(x1, x2, y1, y2){
  if(!isIntersect(x1, x2, y1, y2)) return(0)
  return(min(x2 - y1 + 1, y2 - x1 + 1))
}

tIntersect <- function(t, i, j){
  return(nIntersect(t[i,'V2'], t[i,'V3'], t[j,'V2'], t[j,'V3']))
}


refineArms <- function(t.arm, block.len = 1000){
  
  # gap.open = 10
  # gap.ext = 0.5
  # irow = 6
  
  
  irow = 0
  while(irow < (nrow(t.arm)-1)){
    irow = irow + 1
    
    for(rev.compl in c(FALSE, TRUE)){
      
      t.arm.base = t.arm
      tmp  = base.len - t.arm.base[t.arm.base$dir == 1,'V5'] + 1
      t.arm.base[t.arm.base$dir == 1,'V5'] = base.len - t.arm.base[t.arm.base$dir == 1,'V4'] + 1
      t.arm.base[t.arm.base$dir == 1,'V4'] = tmp
      
      
      print(c(irow, rev.compl))
      
      irow.q = irow
      irow.b = irow
      
      if(t.arm[irow.q+1,'V2'] - t.arm[irow.q,'V3'] - 1 < block.len) next
      if(t.arm.base[irow.b+1,'V4'] - t.arm.base[irow.b,'V5'] - 1 < block.len) next
      
      if(t.arm[irow.q+1,'V2'] - t.arm[irow.q,'V3'] - 1 > 100 * block.len) next
      if(t.arm.base[irow.b+1,'V4'] - t.arm.base[irow.b,'V5'] - 1 > 100 * block.len) next
      
      
      # Perform Alignment
      aln.blocks <- alignBlocks(t.arm, irow.q, irow.b, 
                                query.fas.chr, base.fas.fw,
                                base.fas.bw, block.len, rev.compl = rev.compl)
      t.tmp = aln.blocks$t
      
      if(is.null(t.tmp)) next
      
      # # Visualisation
      # df <- melt(aln.blocks$mx)
      # ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) +
      #   geom_tile() + theme_minimal()
      
      
      # Split with previous
      t.tmp = rbind(t.arm[c(irow, irow+1),], t.tmp)
      
      t.tmp = glueByThreshold(t.tmp, thresholds, query.fas = query.fas.chr, 
                              base.fas.fw = base.fas.fw, 
                              base.fas.bw = base.fas.bw, file.log=file.log, log.append = T)
      
      t.tmp = t.tmp[t.tmp[,'V7'] > 500,]
      
      if(sum(duplicated(rbind(t.tmp, t.arm[c(irow, irow+1),]))) == 2) next
      
      t.arm = rbind(t.arm[-c(irow, irow+1),], t.tmp)
      t.arm <- t.arm[order(t.arm[,'V2']),]
      checkCorrespToGenome(t.arm, query.fas = query.fas.chr, 
                           base.fas.fw = base.fas.fw, 
                           base.fas.bw = base.fas.bw)
      
      
      showt(t.arm, irow:(irow+5))
      
      # irow = irow - 1 + nrow(t.tmp) - 1
      irow  = irow - 1
      break 
      

      # t.arm = afterRefinement(t.arm)
    }
  }
  return(t.arm)
}


afterRefinement <- function(t, min.len = 500){
  t = t[t[,'V7'] > min.len,]
  t <- removeCompleteOverlaps(t, n.distant = NULL, filter.base = F)
  t = t[order(t[, 'V2']),]
  t <- removeShortOverlaps(t, echo = F, gap.thresh = Inf)
  t = glueZero(t)
  t = t[order(t[,'V2']),]
  
  return(t)
}

findPositionInBase <- function(t, pos.search){
  
  t.base = t
  tmp  = base.len - t.base[t.base$dir == 1,'V5'] + 1
  t.base[t.base$dir == 1,'V5'] = base.len - t.base[t.base$dir == 1,'V4'] + 1
  t.base[t.base$dir == 1,'V4'] = tmp
  
  irows = which((t.base[,'V4'] <= pos.search) & (t.base[,'V5'] >= pos.search))
  
  pos.q = c()
  for(irow in irows){
    s1 = strsplit(t[irow,'V8'], '')[[1]]
    s2 = strsplit(t[irow,'V9'], '')[[1]]
    pos = rep(0, length(s1))
    pos1 = pos; pos1[s1 != '-'] = t[irow,'V2']: t[irow,'V3']
    if(t[irow, 'dir'] == 0){
      pos2 = pos; pos2[s2 != '-'] = t[irow,'V4']: t[irow,'V5']
    } else {
      pos2 = pos; pos2[s2 != '-'] = rev(t.base[irow,'V4']: t.base[irow,'V5'])
    }
    pos.q = c(pos.q, pos1[pos2 == pos.search] * (-1)^t[irow,'dir'])
    
  }
  return(pos.q)  
  
}



