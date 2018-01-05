
###########  binary coding ########
.transBinary <- function(fasta){
  marker <- str_replace_all(string = fasta, pattern = "[^AaTUtuGgCc]", replacement = " 0 0 0 0")
  marker <- str_replace_all(string = marker, pattern = "[Aa]", replacement = " 1 0 0 0")
  marker <- str_replace_all(string = marker, pattern = "[TUtu]", replacement = " 0 1 0 0")
  marker <- str_replace_all(string = marker, pattern = "[Gg]", replacement = " 0 0 1 0")
  marker <- str_replace_all(string = marker, pattern = "[Cc]", replacement = " 0 0 0 1")
  marker <- na.omit(as.numeric(unlist(strsplit(marker, " "))))
  return(marker)
}
###########  K-mers coding ########   !!!!! just for one seq
.Kmercoding <- function(Seq,kmer=2){
  #### 1/0 to 1-4 sequence ###
  coding14 <- c()
  for (n in 1:(length(Seq)/4)){
    coding14[n] <- sum(which(Seq[1:4+4*(n-1)]==1))
  }
  #coding14[which(is.na(coding14))] <- 0
  #### k-mer types ###
  a=1:4
  matchtypes <- switch (kmer,
                        do.call(paste0,expand.grid(a)),
                        do.call(paste0,expand.grid(a,a)),
                        do.call(paste0,expand.grid(a,a,a)),
                        do.call(paste0,expand.grid(a,a,a,a)),
                        do.call(paste0,expand.grid(a,a,a,a,a))
  )
  #### incise sequence ####
  kmerseq <- c()
  for ( i in 1:(length(coding14)-(kmer-1))){
    kmerseq[i] <- paste(as.character( coding14[1:kmer+i-1]),collapse="")
  }
  #### matching kmers ###
  .Kmercoding <- table(kmerseq)[matchtypes]/length(kmerseq)
  .Kmercoding[which(is.na(.Kmercoding))] <- 0
  #print(paste0(kmer,"-mers calculate over"))
  .Kmercoding
}
###########  pse coding ########   !!!!! just for one seq
.phyrna <- function(Seq, lambda = 6, w = 0.9, phyrna_nucleo_scale){

  # if(is.null(phyrna_nucleo_scale)){
  #   phyrna_nucleo_scale <- load(system.file('data/phyrna_nucleo_scale.rda', package = "RNAMAP"))
  # }


  .phyrna_seq <- c()
  #### 1/0 to 1-4 sequence ###
  coding14 <- c()
  for (n in 1:(length(Seq)/4)){
    coding14[n] <- sum(which(Seq[1:4+4*(n-1)]==1))
  }
  coding14 <- coding14[which(coding14!=0)]
  #coding14[which(is.na(coding14))] <- 0
  #### k-mer types ###
  a=1:4
  matchtypes <- do.call(paste0,expand.grid(a,a))
  #### incise sequence ####
  kmerseq <- c()
  kmer <- 2
  for ( i in 1:(length(coding14)-(kmer-1))){
    kmerseq[i] <- paste(as.character( coding14[1:kmer+i-1]),collapse="")
  }
  #### tier
  k_tier <- function(kmerseq=kmerseq,lambda=lambda) {
    first_tier <- c()
    ### v  2.0
    #for(i in 1:(length(kmerseq)-lambda)){
    #  first_tier[i] <- mean((phyrna_nucleo_scale[kmerseq[i],]-phyrna_nucleo_scale[kmerseq[i+min((length(kmerseq)-1),lambda)],])^2)
    #}
    ####v 3.0
    if (lambda >= length(kmerseq)) {
      lambda <-  length(kmerseq)-1
      cat("Warning","lambda",lambda, ">= length(kmerseq)",length(kmerseq))
    }
    for(i in 1:(length(kmerseq)-lambda)){
      first_tier[i] <- mean((phyrna_nucleo_scale[kmerseq[i],]-phyrna_nucleo_scale[kmerseq[i+lambda],])^2)
    }

    first_tier
  }

  mer2 <- scale(.Kmercoding(Seq,2))
  p <- c()
  for (j in 1:lambda){
    p[j] <-  mean(k_tier(kmerseq,lambda=j))
  }
  for(j in 1:16){
    .phyrna_seq[j] <- mer2[j]/(sum(mer2)+w*sum(p))
  }
  for(j in 17:(16+lambda)){
    .phyrna_seq[j] <- (w*p[j-16])/(sum(mer2)+w*sum(p))
  }
  .phyrna_seq
}

###########  feature coding ########   !!!!! just for one seq
#' @export
FeatureExtract <- function(RNAseq, feature = c("Binary", "PCP", "Kmer"), lambda = 6, w = 0.9){
  require(stringr)
  phyrna_nucleo <- matrix(c(-12.2,-13.3,-14.2,-10.2,-7.6,-6.6,-10.2,-5.7,-8.0,-10.5,-12.2,-7.6,-7.6,-8.1,-10.2,-6.6,
                            -29.7,-35.5,-34.9,-26.2,-19.2,-18.4,-26.2,-15.5,-19.4,-27.8,-29.7,-19.2,-19.2,-22.6,-26.2,-18.4,
                            -3.26,-2.35,-3.42,-2.24,-2.08,-0.93,-2.24,-1.10,-2.36,-2.11,-3.26,-2.08,-2.11,-1.33,-2.35,-0.93),
                          16,3)
  colnames(phyrna_nucleo) <- c("Enthalpy(Ka/mol)","Entropy(eU)","FreeEnergy(Ka/mol)")
  rownames(phyrna_nucleo) <- c("33","31","34","32","13","11","14","12","43","41","44","42","23","21","24","22")
  phyrna_nucleo_scale <- apply(phyrna_nucleo, 2, scale)
  rownames(phyrna_nucleo_scale) <- c("33","31","34","32","13","11","14","12","43","41","44","42","23","21","24","22")



  ## binary feature
  cat("start converting sequences to binary features......", "\n")
  featureBinary <- t(sapply(RNAseq, .transBinary))
  class(featureBinary) <- "numeric"

  cat("start calculating kmer-based features......", "\n")
  featureMatKmer3 <- t(apply(featureBinary, 1, .Kmercoding, kmer = 3))
  featureMatKmer2 <- t(apply(featureBinary, 1, .Kmercoding, kmer = 2))
  featureMatKmer1 <- t(apply(featureBinary, 1, .Kmercoding, kmer = 1))

  cat("start calculating physicochemical features......", "\n")
  featureMatPC <- t(apply(featureBinary, 1, .phyrna, lambda = lambda, w = w, phyrna_nucleo_scale = phyrna_nucleo_scale))

  cat("merge features......", "\n")
  featureMat <- cbind(featureBinary, featureMatKmer1, featureMatKmer2, featureMatKmer3,featureMatPC)
  colnames(featureMat) <- paste0("Feature", 1:ncol(featureMat))

  featureMat
}



###########  export sequences around motif ########
#' @export
extra_motif_seq <- function(input_seq_dir,text='c',up=5,end=5) {
  require(seqinr)
  require(reshape2)
  cat("loading sequences...","\n")
  input_seq <- read.fasta(input_seq_dir,as.string = T)

  input_seq_c <- lapply(input_seq, function(x) words.pos(x,pattern = 'c'))

  input_seq_name <- names(input_seq)

  input_seq_c_length <- lengths(input_seq_c)

  input_seq_length <- lengths(lapply(input_seq,  s2c))

  cat("make beds...","\n")
  input_seq_bed <- cbind(unlist(sapply(input_seq_name, function(x)
    rep(x,input_seq_c_length[x]))),unlist(input_seq_c)-up,unlist(input_seq_c)+end)
  rownames(input_seq_bed) <-  paste0(input_seq_bed[,1],"_",as.numeric(input_seq_bed[,2])+up)

  input_seq_bed <- input_seq_bed[which(as.numeric(input_seq_bed[,2])>0),]
  input_seq_bed <- input_seq_bed[which(apply(input_seq_bed, 1, function(x) input_seq_length[x[1]]>=as.numeric(x[3]))),]
  #input_seq_bed <- as.matrix(input_seq_bed)
  cat("make c sequences...","\n")
  seqqq <- lapply(input_seq,  s2c)
  c_seq <- apply(input_seq_bed, 1, function(x) seqqq[[x[1]]][as.numeric(x[2]):as.numeric(x[3])])
  c_seq <- as.data.frame((c_seq))
  c_seq <- lapply(c_seq, as.character)

  cat("over","\n")
  return(c_seq)
}

###########  predicted m5c position  ########
#' @export
predict_m5c <- function(sample_feature){
  cutoff.k <- function (attrs, k) {
    if (dim(attrs)[1] == 0)
      return(character(0))
    if (k < 1)
      stop("k too small")
    if (k > dim(attrs)[1])
      k = dim(attrs)[1]
    sorted_names = rownames(attrs)[order(attrs, decreasing = TRUE)]
    return(sorted_names[1:k])
  }
  require(randomForest)
  cat("loading models...\n")
  load(paste0(system.file(package = "PEAm5c"),"/data/models_m5c_randomforest.Rds"))
  #load(file = "~/R/PEAm5c/models_m5c_randomforest.RData")
  res <- 0
  for(i in 1:10){
    alldata <- sample_feature
    alldata <- as.data.frame(alldata)
    model <-  models[[i]]
    res <- predict(model,alldata[,cutoff.k(weights,50)],type=c("vote"))[,2]/10+res
    cat("predicting part ",i," ...\t")

  }
  res_table <- cbind(str_split(names(res),"_",simplify = T),res,NA)
  res_table[which(res>0.484),4] <- "L_m5C"
  res_table[which(res>0.622),4] <- "N_m5C"
  res_table[which(res>0.765),4] <- "H_m5C"
  res_table[which(res>0.891),4] <- "VH_m5C"

  colnames(res_table) <- c("transcript","position","score","mode")
  cat("over")
  return(res_table)
}

#
# cdna <- read.fasta("~/data/arabidopsis/Araport11_genes.201606.cdna.fasta")
# sample_cdna_1000 <- cdna[sample(names(cdna),1000)]
# write.fasta(sequences = sample_cdna_1000,names = names(sample_cdna_1000),file.out = "~/data/arabidopsis/sample_cdna_1000.fa")
#
# aaa <- extra_motif_seq(input_seq_dir = "~/work/ArabidopsisMethylation_m5c/newdata/all/cdna1_6.fa",up = 5)
# aaa <- extra_motif_seq(input_seq_dir = "~/data/arabidopsis/sample_cdna_1000.fa",up = 5)
#
# aaa <- lapply(aaa, c2s)
# bbb <- FeatureExtract(aaa)
# ccc <- predict_m5c(bbb)
#
# ddd <- aaa[names(which(ccc[,4]=="VH_m5C"))]
# source("~/scriptlib/seqlogo.R")
# logoplot(ddd)
