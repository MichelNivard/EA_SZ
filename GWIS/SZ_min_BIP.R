# read in Bipolar results:
pgc.cross.BIP11.2013.05 <- read.csv("pgc.cross.BIP11.2013-05.txt", sep="")

# read in SCZ results:
daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513 <- read.delim("daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_nogras", stringsAsFactors=FALSE)


# create local objects with simple names, remove SNPs not in the bip file from the scz file:
scz <- daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513[daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513[,2] %in% pgc.cross.BIP11.2013.05$snpid,]
bip <- pgc.cross.BIP11.2013.05

# remove the innitial files from memory:
rm(daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513)
rm(pgc.cross.BIP11.2013.05)

# clean up memory:
gc()

## allign reference alleles between scz and bip:

bip_scz <- merge(bip,scz,by.x=1,by.y=2)

# check for allel flips:
nrow(bip_scz)
sum(bip_scz$a1 == bip_scz$A1)
sum(bip_scz$a2 == bip_scz$A2)

# create objects of the flips, to allow future inspection
flips <- bip_scz[bip_scz$a1 != bip_scz$A1,]
flips2 <- bip_scz[bip_scz$a2 != bip_scz$A2,]

# all flips seem to be genuiine allel flips, not random differences in EAF, therfore I remove them all.

bip_scz2 <- bip_scz[bip_scz$a1 == bip_scz$A1 & bip_scz$a2 == bip_scz$A2,]

# file checks out, overwrite main file:
bip_scz <- bip_scz2

# some cleanup to spare memory:
rm(bip_scz2)
gc()


# add static colms to the file whhich will hold the final output:
scz_res <- bip_scz[,1:5]



### perfrom the GWIS

# cov(X,Y) / Var(X) based on the genetic  covariance between bip and scz, and the genetic variance in BIP

b <- 0.3367 /0.4939 

# frequency terms to account for allele effects:

eaf <- bip_scz$FRQ_U_45670

# in this case function is linear, we can skip accounting for the frequency, at little cost, so we will do so:

beta <- log(bip_scz$OR) - b* log(bip_scz$or)

# check the effect size distribution is reasonable for frequent SNPs with high info:
hist(beta[eaf> 0.05 & eaf < 0.95 & bip_scz$INFO > .9 & bip_scz$INFO < 1.1 & bip_scz$info > .9 & bip_scz$info < 1.1 ],breaks=100)

# object to hold the se:
se <- rep(NA,length(beta))

# correlation matrixt his is required and based on the gcov int from ldscore regression to guard agains inflation due to sample overlap  (:
cor1 <- matrix(c(1,  0.1018 , 0.1018  , 1),2,2)

for( i in 1:nrow(bip_scz)){
  
se[i] <- deltamethod(~x2 - b * x1, mean=log(as.numeric(bip_scz[i,c(6,19)])),  cov = diag(bip_scz[i,c(7,20)]) %*% cor1 %*% diag(bip_scz[i,c(7,20)]), ses=T )

}

# compute p-value
p_val <- pchisq((beta/se)^2,1,lower=F)



# make output, use control frequency, and mean info from BIP and SCZ:

scz_res2 <- cbind(scz_res, bip_scz$FRQ_U_45670,as.numeric(as.character(bip_scz$CEUaf)), (bip_scz$info + bip_scz$INFO)/2, beta,se,p_val)
names(scz_res2)[c(6,7,8)] <- c("EAF" ,"EAF_bip","INFO")

write.table(scz_res2 ,file="scz_res.txt",row.names=F,quote=F)
