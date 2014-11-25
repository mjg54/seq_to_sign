Args=commandArgs(TRUE)
arg.quadfreq = Args[1]
arg.mast = Args[2]
arg.genfq = Args[3]
arg.name = Args[4]
arg.win = Args[5]
arg.mememotif = Args[6]
arg.pyfile = Args[7]

library(lattice)

total.nuc <- function(str, col, df) {
  length(grep(str, df[,col]))
}

weighted.cuts.quad <- function(df, fasta.freq) {
  df = data.frame(lapply(df, function(v) {
      if (is.character(v)) return(toupper(v))
      else return(v)
  }))
  vec = vector(mode= 'numeric', length = ncol(df))
  for (i in 1:ncol(df)) {
      fasta.vec = vector(mode= 'numeric', length = nrow(fasta.freq))
      for (j in 1:nrow(fasta.freq)) {
          x = total.nuc(fasta.freq[j,1], i, df) * fasta.freq[j,2]
          fasta.vec[j] = x
      }
      vec[i] = sum(fasta.vec)
  }
  return(vec)
}

process.data <- function(functionstring=filepy, mast=FACTORNAME.txt, nm = 'FACTORNAME', freq, fasta = 'hg19.fa', win =20, motif = 1) {
    command=paste('python',functionstring, "-i", mast, "-f", fasta, "-w", win, "-n", nm, '-m', motif, sep=" ")
    cat(command,"\n")
    try(system(command))
    system('echo calculating the seqToSign trace...')
    factor = read.table(paste(nm,'.txt', sep=''),header=F,colClasses = 'character')
    factor.vec = weighted.cuts.quad(factor, freq)
    if (length(factor.vec)%%2 == 1) {
        x.plot = seq(-length(factor.vec)/2 +1, length(factor.vec)/2, 1)        
    } else {
        x.plot=seq(-length(factor.vec)/2 + 1.5, length(factor.vec)/2 + 0.5, 1)
    }
    return(cbind(x.plot, factor.vec))
}


composites <- function(dat, fact = 'DNase', summit = 'Motif', num=20) {
    system(paste('echo generating the figure: composite_', fact, '_signals_', summit, '.pdf', sep=''))
    col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2))
    pdf(paste('composite_', fact, '_signals_', summit, '.pdf', sep=''), width=4, height=4) 
    print(xyplot(dat[,2] ~ dat[,1],
                 type = 'l',
                 scales=list(x=list(cex=1.0, relation = "free"), y =list(cex=1.0, relation="free")),
                 xlim=c(-(as.numeric(as.character(num))),as.numeric(as.character(num))),
                 col = col.lines,
                 auto.key = list(points=F, lines=T),
                 par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines), superpose.line = list(col = col.lines, lwd=3)),
                 cex.axis=1.0,
                 aspect=0.7,
                 lwd=3,
                 ylab = paste(fact," Cut Frequency", sep=''),
                 xlab = paste("Distance from ", summit, "motif center",sep=''),
               ))
  sup = dev.off()
}

#load quad nucleotide cut frequencies
quad = read.table(arg.quadfreq)

#call python functions to extract fasta sequences and then run model
intensities = process.data(functionstring = arg.pyfile, mast = arg.mast, nm = arg.name, quad, fasta = arg.genfq, win = arg.win, motif = arg.mememotif)

#draws the plot
composites(intensities, summit = arg.name, num = arg.win)
