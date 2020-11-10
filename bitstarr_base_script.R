library(data.table)
library(plyr)
library(ggplot2)
library(qqman)
library(QuASAR)
#SNPs object contains positional information matched to a rsID from bed file of snps sent to be synthesized
#comes from SNPs <- read.table("~/piquelab/cindy/bitStarr/old/mpraSiteSelection/mpraSNPs_090215.bed", header=T,sep='\t')
load("/wsu/home/fh/fh85/fh8591/piquelab/cindy/bitStarr/bitstarr_objects_dna1_dna2.RData") 
#read in sequencing data counts
bit2 <- read.table("/wsu/home/groups/piquelab/cindy/bitStarr/all_base/dedup.bit2dna.q20.query.txt", sep='\t', header=F)
#cut out unneeded columns
bit2 <- bit2[,-c(5,10,15,20,25,30,35)]
#fix header. header is chr pos1 ref alt refF,refR,altF,altR 
colnames(bit2) <- c("chr","pos1","ref","alt","dna10_1.R.F","dna10_1.R.R","dna10_1.A.F","dna10_1.A.R",
  "dna10_2.R.F","dna10_2.R.R","dna10_2.A.F","dna10_2.A.R",
  "dna10_3.R.F","dna10_3.R.R","dna10_3.A.F","dna10_3.A.R",
  "dna50_1.R.F","dna50_1.R.R","dna50_1.A.F","dna50_1.A.R",
  "dna50_2.R.F","dna50_2.R.R","dna50_2.A.F","dna50_2.A.R",
  "dna50_3.R.F","dna50_3.R.R","dna50_3.A.F","dna50_3.A.R",
  "dna50_pcrcocktail.R.F","dna50_pcrcocktail.R.R","dna50_pcrcocktail.A.F","dna50_pcrcocktail.A.R")
#grab rsID and make sure the SNP was part of the original designed set
bit2_rs <- merge(SNPs, bit2[,-c(3:4)], by=c("chr","pos1"))
#sum counts for each construct 
bit2_rs$DNA_R_fw=tryCatch(rowSums(bit2_rs[ ,grepl("[.]R[.]F", names(bit2_rs))]), error=function(x) bit2_rs[ ,grepl("[.]R[.]F", names(bit2_rs))])
bit2_rs$DNA_A_fw=tryCatch(rowSums(bit2_rs[ ,grepl("[.]A[.]F", names(bit2_rs))]), error=function(x) bit2_rs[ ,grepl("[.]A[.]F", names(bit2_rs))])
bit2_rs$DNA_R_rv=tryCatch(rowSums(bit2_rs[ ,grepl("[.]R[.]R", names(bit2_rs))]), error=function(x) bit2_rs[ ,grepl("[.]R[.]R", names(bit2_rs))])
bit2_rs$DNA_A_rv=tryCatch(rowSums(bit2_rs[ ,grepl("[.]A[.]R", names(bit2_rs))]), error=function(x) bit2_rs[ ,grepl("[.]A[.]R", names(bit2_rs))])
#sum total counts for each SNP
bit2_rs$DNA=bit2_rs$DNA_R_fw+bit2_rs$DNA_R_rv+bit2_rs$DNA_A_fw+bit2_rs$DNA_A_rv
#mpileup produces some duplicate positions, this sums them
i3 <- aggregate(cbind(DNA_R_fw,DNA_A_fw,DNA_R_rv,DNA_A_rv) ~rsID , data=bit2_rs,sum)
#I probably could have put this in the previous line
i3$DNA=i3$DNA_R_fw+i3$DNA_R_rv+i3$DNA_A_fw+i3$DNA_A_rv
#remove low covered SNPs
bit2_dna <- subset(i3, DNA>7)
bit2_dna_F <- subset(bit2_dna[,-c(4:5)], DNA_R_fw>0 & DNA_A_fw>0)
bit2_dna_R <- subset(bit2_dna[,-c(2:3)], DNA_R_rv>0 & DNA_A_rv>0)
#calculate proportion for each SNP, for each direction separately
bit2_dna_F <- transform(bit2_dna_F, DNA_prop_fw=DNA_R_fw/(DNA_R_fw+DNA_A_fw)) 
bit2_dna_R <- transform(bit2_dna_R, DNA_prop_rv=DNA_R_rv/(DNA_R_rv+DNA_A_rv))
bit2_dna_m <- transform(bit2_dna, DNA_prop_fw=ifelse(DNA_R_fw==0, DNA_R_fw+1/(DNA_R_fw+DNA_A_fw+1),DNA_R_fw/(DNA_R_fw+DNA_A_fw)), DNA_prop_rv=ifelse(DNA_R_rv==0, DNA_R_rv+1/(DNA_R_rv+DNA_A_rv+1),DNA_R_rv/(DNA_R_rv+DNA_A_rv)))
#calculate total SNP counts, keeping direction separate, and creating an identifier that includes direction information
bit2_dna_R <- transform(bit2_dna_R, DNA=DNA_A_rv+DNA_R_rv, identifier=paste0(rsID,"_","rv"))
bit2_dna_F <- transform(bit2_dna_F, DNA=DNA_A_fw+DNA_R_fw, identifier=paste0(rsID,"_","fw"))
#put the directions back together for the final list of identifiers present in the library
bit2_dna_all <- rbind(bit2_dna_F[,c("rsID","identifier","DNA")],bit2_dna_R[,c("rsID","identifier","DNA")])

################ RNA ###################
### IMPORTANT NEED FILES IN SINGLE LOCATION ###
#directory containing all sequencing counts data, one file per replicate
myDir <- "~/piquelab/scratch/cindy/UMI/working"
filenames <- list.files(myDir) 
filenames <- filenames[grep("[.]q20.query.txt$", filenames)] 
data_names <- gsub("[.]q20.query.txt", "", filenames) 
shortnames <- gsub("dedup.", "", data_names)
#read in counts files and assign df name based on the replicate identity
for(i in 1:length(filenames)) assign(data_names[i], setnames(fread(file.path(myDir, filenames[i]),header = FALSE, sep='\t', drop = c(3:5)), 
  c("chr","pos1",paste0(shortnames[i],".R.F"),paste0(shortnames[i],".R.R"),paste0(shortnames[i],".A.F"),paste0(shortnames[i],".A.R"))))

df <- lapply(data_names, function(x) get(x))
#recursively merge each replicate
df_merge <- as.data.frame(Reduce(function(x, y) merge(x, y, all=TRUE,by=c("chr","pos1"),all.x=TRUE, all.y=TRUE),df,accumulate=F))
df_merge[is.na(df_merge)] <- 0
#sum counts for each construct 
df_merge$R_fw=tryCatch(rowSums(df_merge[ ,grepl("[.]R[.]F", names(df_merge))]), error=function(x) df_merge[ ,grepl("[.]R[.]F", names(df_merge))])
df_merge$A_fw=tryCatch(rowSums(df_merge[ ,grepl("[.]A[.]F", names(df_merge))]), error=function(x) df_merge[ ,grepl("[.]A[.]F", names(df_merge))])
df_merge$R_rv=tryCatch(rowSums(df_merge[ ,grepl("[.]R[.]R", names(df_merge))]), error=function(x) df_merge[ ,grepl("[.]R[.]R", names(df_merge))])
df_merge$A_rv=tryCatch(rowSums(df_merge[ ,grepl("[.]A[.]R", names(df_merge))]), error=function(x) df_merge[ ,grepl("[.]A[.]R", names(df_merge))])
#grab rsID and make sure the SNP was part of the original designed set
i1 <- merge(SNPs, df_merge, by=c("chr","pos1"))
#make sure SNP was part of DNA library (possibly could skip this since a directional merge is done later)
i2 <- subset(i1, rsID %in% bit2_dna_m$rsID)
#mpileup produces some duplicate positions, this sums them
i3 <- aggregate(cbind(R_fw,A_fw,R_rv,A_rv) ~rsID , data=i2,sum)
#I probably could have put this in the previous line
i3$RNA=i3$R_fw+i3$R_rv+i3$A_fw+i3$A_rv
#remove low covered SNPs
i3_1 <- subset(i3, RNA>1)
bit_F <- subset(i3, R_fw>0 & A_fw>0)
bit_R <- subset(i3, R_rv>0 & A_rv>0)
#make sure SNP was part of DNA library, each direction separate
bit_dna_F_m <- merge(bit_F, bit2_dna_F, by="rsID")
bit_dna_R_m <- merge(bit_R, bit2_dna_R, by="rsID")
bit_dna_m <- merge(bit_dna_F_m, bit_dna_R_m, by=c("rsID","R_fw","R_rv","A_fw","A_rv","DNA","RNA"), all=T)
#run quasar-mpra on each replicate separately
l.df <- lapply(data_names, function(x) get(x))
names(l.df) <- shortnames
test <- lapply(l.df, function(i){ 
  #i <- l.df[[1]]
  i[is.na(i)] <- 0
  i$R_fw=rowSums(i[ ,3]) #ref forward counts
  i$A_fw=rowSums(i[ ,4]) #alt forward counts
  i$R_rv=rowSums(i[ ,5]) #ref reverse counts
  i$A_rv=rowSums(i[ ,6]) #alt reverse counts
  i1 <- merge(SNPs, i, by=c("chr","pos1")) #grab rsID
  i2 <- subset(i1, rsID %in% bit_dna_m$rsID) #remove low count SNPs
  i3 <- aggregate(cbind(R_fw,A_fw,R_rv,A_rv) ~rsID , data=i2,sum) #sum duplicate positions
  data_all <- merge(i3[,-c(4:5)], bit_dna_F_m[,-c(2:6)], by="rsID") #remove info from other direction and add in DNA prop
  data_all_R <- merge(i3[,-c(2:3)], bit_dna_R_m[,-c(2:6)], by="rsID") 
  data_F <- subset(data_all, R_fw>0 | A_fw>0) #make sure SNP has at least one allele with counts
  data_R <- subset(data_all_R, R_rv>0 | A_rv>0)
  colnames(data_F) <- c("rsID","R","A","DNA_R","DNA_A","DNA","DNA_prop") #fix col names to rbind later
  data_F$identifier <- paste(data_F$rsID,"fw",sep='_') #create identifier
  data_F$Direction <- "fw"
  colnames(data_R) <- c("rsID","R","A","DNA_R","DNA_A","DNA","DNA_prop")
  data_R$identifier <- paste(data_R$rsID,"rv",sep='_')
  data_R$Direction <- "rv"
  both <- rbind(data_F, data_R)
  both.res <- fitQuasarMpra(both$R,both$A,both$DNA_prop) #quasar-mpra
  i_both2 <- cbind(both, both.res) #add original information to the quasar results
  i_both2$betas_w <- 1/i_both2$betas_se^2 #calculate beta weights
  i_both2$betas <- i_both2$betas_w*i_both2$betas.beta.binom #recalculate betas using weights
  return(i_both2)
})
#convert list to a df
dat_both <- ldply(test, data.frame)
#preform fixed-effects meta-analysis to combine replicates
ddat_both <- ddply(dat_both, c("identifier","Direction","rsID"), plyr::summarize, #for each identifier
  betas_T=sum(betas)/sum(betas_w), #new beta
  betas_se_comb = sqrt(1/sum(betas_w)), #new beta se
  betas_T_low = betas_T-(1.96*betas_se_comb), #CI interval 
  betas_T_high = betas_T+(1.96*betas_se_comb),
  betas_z_comb = (betas_T - qlogis(median(DNA_prop)))/betas_se_comb, #z score
  DNA_prop = median(DNA_prop), #not really a median, all should be the same number
  p_comb = 2 * pnorm(-abs(betas_z_comb)) # new pvalue
)
ddat_both$padj_comb <- p.adjust(ddat_both$p_comb,method="BH") #BH adjusted pval
