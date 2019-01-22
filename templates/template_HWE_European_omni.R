
Rplink <- function(PHENO,GENO,CLUSTER,COVAR){

	library(hwde)
	
        # annotation file that has the affection status, ethnic group and batch information
	annotation  <- read.csv("INDIVIDUALS_ANNOTATION", sep=" ", header=T)  

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
		batch <- sort(unique(annotation$batch))
		
                # change the ethnic group to what you want
		eth <- "European" # user defined variable
			
                # calculate HWE across entire collection
		ii <- annotation$phenotype=="1" & annotation$ethnicity_predicted ==eth
		ret <- runHWE(x, ii)
		pval[i] <- ret$pval
		nsample[i] <- ret$sample
		i<-i+1
		for (batch in batch){
			# calculate HWE for particular batch
			ii <- annotation$phenotype=="1" & annotation$ethnicity_predicted ==eth & annotation$batch==batch
			ret <- runHWE(x, ii)
			pval[i] <- ret$pval
			nsample[i] <- ret$sample
			i<-i+1
				
			# calculate HWE for all batches excluding particular batch
			ii <- annotation$phenotype=="1" & annotation$ethnicity_predicted ==eth & annotation$batch!=batch
			ret <- runHWE(x, ii)
			pval[i] <- ret$pval
			nsample[i] <- ret$sample
			i<-i+1
		}
	return(c(length(pval)*2, c(pval, nsample)))
	}

	apply( GENO , 2 , hwe, annotation)
}
