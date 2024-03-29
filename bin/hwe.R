library(snpStats)
library(hwde)

# Expected arguments: <plink-basename> <individuals-annotation> <outfile-name>

args = commandArgs(TRUE)
basename <- args[1]
dataset <- read.plink(bed=paste(basename,"bed",sep="."),
                      bim=paste(basename,"bim",sep="."),
                      fam=paste(basename,"fam",sep="."))
annotation  <- read.csv(args[2], sep=" ", header=T)
# sort by sample ID
annotation_s <- annotation[order(annotation[,2]),]
outfile = args[3]
ethnicity = args[4]

runHWE <- function(x, ii){
   obs_hom1 <- sum(x[ii]==0, na.rm=T)
   obs_hets <- sum(x[ii]==1, na.rm=T)
   obs_hom2 <- sum(x[ii]==2, na.rm=T)
   nsample <- obs_hom1 + obs_hets +obs_hom2
   if(nsample==0){
      return(list(pval=1, sample=nsample))
   }
   else{
      return(list(pval=hwexact(obs_hom1, obs_hets, obs_hom2), sample=nsample))
   }
}

hwe <-function(x, annotation){
   pval <- c()
   nsample <- c()
   i<-1
   # change the ethnic group to what you want
   eth <- ethnicity # user defined variable

   # calculate HWE across entire collection
   ii <- annotation$ethnicity_predicted ==eth
   ret <- runHWE(x, ii)
   pval[i] <- ret$pval
   nsample[i] <- ret$sample
   i<-i+1
   for (batch in batch){
      # calculate HWE for particular batch
      ii <- annotation$ethnicity_predicted ==eth & annotation$batch==batch
      ret <- runHWE(x, ii)
      pval[i] <- ret$pval
      nsample[i] <- ret$sample
      i<-i+1

      # calculate HWE for all batches excluding particular batch
      ii <- annotation$ethnicity_predicted ==eth & annotation$batch!=batch
      ret <- runHWE(x, ii)
      pval[i] <- ret$pval
      nsample[i] <- ret$sample
      i<-i+1
   }
   return(c(length(pval)*2, c(pval, nsample)))
}

genotypes = as(dataset$genotypes, "numeric")
# sort by sample ID
genotypes_s <- genotypes[order(rownames(genotypes)),]
batch <- sort(unique(annotation$batch))
#   print(paste("Found batches: ", batch))
# for a correct association of the samples from the dataset to the sample in the annotations file, we apply the test on the sorted matrices:
out <- t(apply(genotypes_s, 2, hwe, annotation_s))
out <- out[,-1]
if(file.exists(outfile)) {
   file.remove(outfile)
}
for(row in 1:nrow(out)) {
   cat(dataset$map$chromosome[row], dataset$map$snp.name[row], dataset$map$position[row], dataset$map$allele.1[row], out[row,], "\n", file=outfile, sep="\t", append=TRUE)
}
