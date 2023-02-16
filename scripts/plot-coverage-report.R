#!/usr/bin/env Rscript

# 16.05.2013 12:56:20 EDT
# Harm van Bakel <hvbakel@gmail.com>
# Divya Kriti <divya.kriti@mssm.edu>
# Ana Gonzalez-Reiche <anasilvia.gonzalez-reiche@mssm.edu>

#############
# ARGUMENTS #
#############

#  Col 1: the long flag name. A multi-character string.
#  Col 2: short flag alias of Column 1. A single-character string.
#  Col 3: 0=no argument, 1=required argument, 2=optional argument.
#  Col 4: Possible values: logical, integer, double, complex, character.
#  Col 5: A brief description of the purpose of the option.

library(getopt)
args = matrix(c('input', 'i', 1, "character", "Tab-delimited input file with genome coverage (name, position, coverage) and variant frequencies",
                'markers', 'm', 2, "character", "Optional tab-delimited input file with positions to draw vertical lines (name, pos)",
                'output', 'o', 2, "character", "Plot output file prefix",
                'lowcov', 'l', 2, "numeric",   "Highlight areas with coverage below this value. Default: 100X",
                'varthr', 't', 2, "numeric",   "Threshold for variant detection. Default: 0.1",
                'rows', 'r', 2, "numeric",   "Number of plot rows per page. Default: 2",
                'cols', 'c', 2, "numeric",   "Number of plot columns per page. Default: 1",
                'scale', 's', 2, "numeric",   "Plot scaling factor. Default: 0.8",
                'help', 'h', 0, "logical",   "Brief help message"
), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$lowcov)   ) { opt$lowcov   = 100 }
if ( is.null(opt$rows)     ) { opt$rows     = 2   }
if ( is.null(opt$cols)     ) { opt$cols     = 1   }
if ( is.null(opt$varthr)   ) { opt$varthr   = 0.1 }
if ( is.null(opt$scale)    ) { opt$scale    = 0.8 }
if ( is.null(opt$output)   ) { opt$output   = "Sample" }


# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$output)) {
  self=sub("^.*/", "", commandArgs()[4], perl=T)
  cat("\n", getopt(args, command=self, usage=T), "\n")
  q(status=1);
}

#############
# LIBRARIES #
#############

library(ggplot2);
library(tools);
library(dplyr);

########
# MAIN #
########

# Read data files
d.var = read.table(opt$input, sep="\t", comment.char="", header=T);
d.var = d.var[d.var$depth>0,];
opt$lowcov = opt$lowcov/1000;
if (!is.null(opt$markers)){
  d.mrk = read.table(opt$markers, colClasses=c("character","integer"), sep="\t", comment.char="#", header=F);
}

# Check input files
if(!any(names(d.var) %in% c("reference","position","reference_base","reference_base_fraction","depth","forward_depth","reverse_detph","A","a","C","c","G","c","T","t")))
  stop("Error: unexpected variance file format\n");

# Get minor variants ### Anotate deletions and insertions if present####
d.var$depth = rowSums(d.var[,13:23])
colnames(d.var)[colnames(d.var) %in% c("insertion", "deletion")] = c("ins", "del")
alt.ratio  = data.frame(A=(d.var$A+d.var$a)/d.var$depth, T=(d.var$T+d.var$t)/d.var$depth, C=(d.var$C+d.var$c)/d.var$depth, G=(d.var$G+d.var$g)/d.var$depth, ins = d.var$ins/d.var$depth, del = d.var$del/d.var$depth ) ;
set.var    = apply(alt.ratio >= opt$varthr & alt.ratio < 0.5, 1, any)
alt.ratio  = alt.ratio[set.var,];
alt.ratio  = alt.ratio[complete.cases(alt.ratio),];
alt.bases  = apply(alt.ratio, 1, function(x, thr, bases){paste(as.character(bases[x>=thr & x<0.5]), collapse="/")}, opt$varthr, names(alt.ratio))

# Annotate only variants present in both forward and reverse reads
set.var[is.na(set.var)] <- FALSE

if (any(set.var)) {
  keep.vars = c()
  if (length(alt.bases)){
    for (i in 1:length(alt.bases)){
      base = alt.bases[[i]]
      row = as.numeric(names(alt.bases)[i])
      base = unlist(strsplit(base,"/"))
      print(base)
      
    }
      # check that base is not the same as pilon base
      if (any(c("ins", "del") %in% base)) {
        keep.vars = c(keep.vars,i)
      } else {
          for (b in base){
              if(d.var[row,b] > 0 & d.var[row,tolower(b)] > 0){
                  keep.vars = c(keep.vars,i)
              }
              if (b == toupper(d.var$pilon_base[row])) {
                  all.bases = colnames(alt.ratio[which(alt.ratio[i,]>opt$varthr)])
                  if (length(all.bases) == 1){
                      alt.bases[[i]] = ""
                  } else {
                      if (all(all.bases %in% b)){
                          alt.bases[[i]] = i 
                      } else {
                          print(all.bases)
                        print("all.bases[!all.bases %in% b]")
                        print(b)
                        print(all.bases %in% b)
                        print(all.bases[!all.bases %in% b])
                        print(alt.bases)
                          #################THIS THROWS ERROR BC IT RETURNS 2. LOOP? IF ALL?)######
                          alt.bases[[i]] = all.bases[!all.bases %in% b]
            }
          }
          
        } 
      }
    }
    alt.bases = alt.bases[keep.vars]
    vars = as.numeric(names(alt.bases))
    alt.change = paste(d.var$position[vars], paste(toupper(d.var$pilon_base[vars]), alt.bases, sep=">"))
    variants   = data.frame(d.var$reference[vars], d.var$position[vars], alt.change);
    variants   = variants[order(variants[,1],variants[,2]),];
  } else {
    cat("None of the detected variants were present in both read pairs. Coverage might be low")
  }
}

# Open plot file
pdf(file=paste(opt$output));
tryCatch( { 
  
  # Make coverage plots
  par(mfrow=c(opt$rows,opt$cols), mar=c(5,4,2,5)+.1, cex=opt$scale);
  
  # Subset data
  d.var[,7] = d.var[,7]/1000;
  
  x.cov    = d.var[,2];
  y.cov    = d.var[,7];
  
  # Plot coverage data
  if (!is.null(opt$ylim)){
    ylim = as.numeric(unlist(strsplit(opt$ylim, split=",")));
    if (length(ylim) != 2) stop("Error: ylim must have two values 'min,max'");
    ylim = sort(ylim);
  }else{
    ylim = c(0, max(y.cov));
  }
  plot(x.cov,y.cov, type="n",xlab="Position (bp)", ylab="Base coverage (x 1000)", ylim=ylim, main='Intra-host variance plot');
  if(opt$lowcov){
    lowcov = x.cov[y.cov<opt$lowcov]
    lines(lowcov,rep(ylim[2],length(lowcov)),type="h", col=rgb(229/255,229/255,229/255));
  }
  if (is.null(opt$bars)){
    polygon(c(min(x.cov), x.cov, max(x.cov)), c(0, y.cov, 0), lwd=1.5, col=rgb(97/255,156/255,255/255,0.3), border=rgb(97/255,156/255,255/255));
  }else{
    lines(x.cov,y.cov,type="h");
  }
  
  # Plot vertical lines
  if (!is.null(opt$markers)){
    #set.mrk     = d.mrk[d.mrk[,1]==seg,];
    for (pos in 1:nrow(d.mrk)){
      abline(v=d.mrk[pos,2], col="blue");
    }
  }
  
  # Plot median coverage
  median.cov = tapply(d.var[,7],d.var[,1], median);
  names(median.cov) = opt$output;
  abline(h=median.cov, col="green");
  text(5000, 1.2*median.cov, paste("median coverage", median.cov, "x 1000", sep = " "), cex = 0.8, col="blue");
  
  # Get available variant data
  d.var[,8] = as.character(d.var[,8])
  d.var[,8] = as.numeric(d.var[,8])
  x.var    = d.var[,2];
  y.var    = 1-d.var[,8];
  # Correct minor variants with higher frequency > 0.5
  y.var[is.na(y.var)] <- 0
  y.var[y.var > 0.5] = 1 - y.var[y.var > 0.5];
  
  # Plot variant data
  par(new=TRUE, cex=opt$scale);
  plot(x.var, y.var, type="h", lwd=1.5, xaxt="n",yaxt="n",xlab="",ylab="", ylim=c(0,0.5), col=rgb(180/255,6/255,205/255));
  grid();
  axis(4);
  mtext("Variant frequency",side=4,line=3, las=0, cex=opt$scale);
  
# Plot variant labels
  if (any(set.var)){
    mv = 0.02;
    text.x  = variants[,2];
    text.y  = y.var[x.var %in% text.x]+0.015;
    xdiff   = c(2500, abs(diff(text.x)));
    ydiff   = c(1, abs(diff(text.y)));
    ymv.set = ydiff<0.05 & xdiff<200;
    if (any(ymv.set)){
      ymv.val = seq(mv, sum(ymv.set)*mv, mv)
      text.y[ymv.set] = text.y[ymv.set] + ymv.val;
    }
    text.y[text.y>0.49] = 0.49;
    text(text.x, text.y, variants[,3], cex=0.7);
  }
}, error=function(e){ # Just print exception message
    message("\n### An error occurred during plotting: possibly because no minor variants were found\n\n");
    print(e);
  #save.image(file = "error.RData");
}
);


# Close plot file
dev=dev.off()
