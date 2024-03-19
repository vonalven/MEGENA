
count.FD <- function(rho,thresh.vec,direction = "ascending")
{
  # combine threshold and correlation values for efficiency
  n.rho <- length(rho)
  label.vec <- c(rep(1,length(rho)),rep(0,length(thresh.vec)));# label correlation value and threshold value
  rho <- c(rho,thresh.vec)# combine correlation and threhold values for efficiency.
  
  if (direction == "ascending")
  {
    # sort correlation values
    i <- order(rho)
    rho <- rho[i]
    label.vec <- label.vec[i]
    
    # count permuted correlation values above thresholds
    j <- which(label.vec == 0)
    output <- cbind(rho[j],(n.rho - cumsum(label.vec)[j])/n.rho)
    colnames(output) <- c("rho.cutoff","FDR")
    
  }
  if (direction == "descending")
  {
    # sort correlation values
    i <- order(rho,decreasing = TRUE)
    rho <- rho[i]
    label.vec <- label.vec[i]
    
    # count permuted correlation values above thresholds
    j <- which(label.vec == 0)
    output <- cbind(rho[j],(n.rho - cumsum(label.vec)[j])/n.rho)
    colnames(output) <- c("rho.cutoff","FDR")
    
  }
  return(output)
} 

calculate.rho.signed <- function(datExpr,n.perm,FDR.cutoff,estimator = "pearson",
                                 use.obs = "na.or.complete",
                                 direction = "positive",
                                 rho.thresh = NULL,sort.el = TRUE)
{
  if (is.null(rownames(datExpr))) rownames(datExpr) <- paste("g",1:nrow(datExpr),sep = "")
  gid <- rownames(datExpr)
  datExpr <- t(datExpr)
  
  if (direction == "absolute") rho <- abs(cor(datExpr,method = estimator,use = use.obs))
  if (direction != "absolute") rho <- cor(datExpr,method = estimator,use = use.obs)
  
  if (is.null(rho.thresh)) 
  {
    if (direction == "absolute") rho.thresh <- seq(0,1,0.01)
    if (direction != "absolute") rho.thresh <- seq(-1,1,0.01)
  }
  
  #### permute data matrix to calculate FDR
  set.seed(1234)
  nc <- nrow(datExpr)
  perm.ind <- lapply(1:n.perm,function(i,n) sample(1:n,n),n = nc)
  count.out <- vector("list",n.perm)
  for (i in 1:n.perm)
  {
    cat("i = ");cat(i);cat("\n");
    if (direction == "absolute") random.rho <- abs(cor(datExpr,datExpr[perm.ind[[i]],],method = estimator,use = use.obs))
    if (direction != "absolute") random.rho <- cor(datExpr,datExpr[perm.ind[[i]],],method = estimator,use = use.obs)
    
    random.rho <- as.vector(random.rho[upper.tri(random.rho)]);
    
    if (direction == "absolute" | direction == "positive") count.out[[i]] <- count.FD(random.rho,rho.thresh,direction = "ascending")
    if (direction == "negative") count.out[[i]] <- count.FD(random.rho,rho.thresh,direction = "descending")
    
    rm(random.rho)
  }
  
  if (direction == "absolute" | direction == "positive") PR = count.FD(as.vector(rho[upper.tri(rho)]),rho.thresh,direction = "ascending");
  if (direction == "negative") PR = count.FD(as.vector(rho[upper.tri(rho)]),rho.thresh,direction = "descending");
  
  PR = PR[,2]
  FPR = Reduce("+",lapply(count.out,function(x) x[,2]))/n.perm;FPR[1] <- 1;
  FDR = FPR/PR;FDR[which(FPR == 0)] <- 0;FDR[which(FDR > 1)] <- 1;
  
  ### apply constraint that higher threshold yields less FDR than lower thresholds
  mx = FDR[1]
  for (i in 2:length(FDR))
  {
    if (FDR[i] > mx) 
    {
      FDR[i] = mx 
    }else{
      mx = FDR[i]
    }
  }
  
  if (direction == "absolute" | direction == "positive") FDR.table <- data.frame(cut.off = rho.thresh,FPR = FPR,PR = PR,FDR = FDR)
  if (direction == "negative") FDR.table <- data.frame(cut.off = rev(rho.thresh),FPR = FPR,PR = PR,FDR = FDR)
  
  
  # choose threshold respect to FDR.cutoff
  if (direction == "absolute" | direction == "positive") 
  {
    rho.cutoff <- min(FDR.table$cut.off[FDR.table$FDR < FDR.cutoff])
    
    # coerce into edgelist
    ij <- which(rho > rho.cutoff & upper.tri(rho),arr.ind = T)
    w <- apply(ij,1,function(xy,m) m[xy[1],xy[2]],m = rho)
    ijw <- cbind(ij,w);colnames(ijw) <- c("row","col","rho")
    
    ijw <- as.data.frame(ijw)
    ijw[[1]] <- gid[ijw[[1]]];ijw[[2]] <- gid[ijw[[2]]]
    
    if (sort.el) ijw <- ijw[order(ijw[,3],decreasing = T),]
    
    output <- list(signif.ijw = ijw,FDR = FDR.table)
    
  }
  if (direction == "negative")
  {
    rho.cutoff <- max(FDR.table$cut.off[FDR.table$FDR < FDR.cutoff])
    
    # coerce into edgelist
    ij <- which(rho <= rho.cutoff & upper.tri(rho),arr.ind = T)
    w <- apply(ij,1,function(xy,m) m[xy[1],xy[2]],m = rho)
    ijw <- cbind(ij,w);colnames(ijw) <- c("row","col","rho")
    
    ijw <- as.data.frame(ijw)
    ijw[[1]] <- gid[ijw[[1]]];ijw[[2]] <- gid[ijw[[2]]]
    
    if (sort.el) ijw <- ijw[order(ijw[,3],decreasing = F),]
    
    output <- list(signif.ijw = ijw,FDR = FDR.table)
    
  }
  return(output)
}



calculate.rho.twoMat <- function(data.mat1,data.mat2,n.perm,FDR.cutoff,estimator = "pearson",rho.thresh = NULL,sort.el = TRUE)
{
	if (ncol(data.mat1) != ncol(data.mat2)) stop("columns of data.mat1 and data.mat2 do not agreen.")
	if (is.null(rownames(data.mat1))) rownames(data.mat1) <- paste("A",1:nrow(data.mat1),sep = "")
	if (is.null(rownames(data.mat2))) rownames(data.mat2) <- paste("B",1:nrow(data.mat2),sep = "")
	gid1 <- rownames(data.mat1);gid2 <- rownames(data.mat2)
	
	data.mat1 <- t(data.mat1);data.mat2 <- t(data.mat2);
	
	rho <- abs(cor(x = cbind(data.mat1),y = cbind(data.mat2),method = estimator))
	
	if (is.null(rho.thresh)) rho.thresh <- seq(0,1,0.01)
	
	#### permute data matrix to calculate FDR
	nc <- nrow(data.mat1)
	perm.ind <- lapply(1:n.perm,function(i,n) sample(1:n,n),n = nc)
	count.out <- vector("list",n.perm)
	for (i in 1:n.perm)
	{
		cat("i = ");cat(i);cat("\n");
		random.rho <- abs(cor(x = data.mat1,y = data.mat2[perm.ind[[i]],],method = estimator))
		random.rho <- as.vector(random.rho);
		
		count.out[[i]] <- count.FD(random.rho,rho.thresh)
		rm(random.rho)
	}
	PR = count.FD(as.vector(rho),rho.thresh);PR = PR[,2]
	FPR = Reduce("+",lapply(count.out,function(x) x[,2]))/n.perm;FPR[1] <- 1;
	FDR = FPR/PR;FDR[which(FPR == 0)] <- 0;FDR[which(FDR > 1)] <- 1;
	FDR.table <- data.frame(cut.off = rho.thresh,FPR = FPR,PR = PR,FDR = FDR)

	# choose threshold respect to FDR.cutoff
	rho.cutoff <- min(FDR.table$cut.off[FDR.table$FDR < FDR.cutoff])

	# coerce into edgelist
	ij <- which(rho > rho.cutoff & upper.tri(rho),arr.ind = T)
	w <- apply(ij,1,function(xy,m) m[xy[1],xy[2]],m = rho)
	ijw <- cbind(ij,w);colnames(ijw) <- c("row","col","rho")

	ijw <- as.data.frame(ijw)
	ijw[[1]] <- gid1[ijw[[1]]];ijw[[2]] <- gid2[ijw[[2]]]

	if (sort.el) ijw <- ijw[order(ijw[,3],decreasing = T),]

	output <- list(signif.ijw = ijw,FDR = FDR.table)

	return(output)
}

################################################

test.pairwiseCor <- function(data.mat1,data.mat2 = NULL,alternative = "two.sided",method = "pearson",use.obs = "na.or.complete")
{
 cat("##### Pairwise Correlation Analysis ######\n")
 
 if (is.null(data.mat2))
 {
  cor.pairs <- do.call(rbind,lapply(1:(nrow(data.mat1)-1),function(i,n) cbind(rep(i,(n-i)),(i+1):n),n = nrow(data.mat1)))
  data.mat2 <- data.mat1;
 }else{
  if (ncol(data.mat1) != ncol(data.mat2)) stop("data.mat1 and data.mat2 do not have same number of columns.");
  cor.pairs <- do.call(rbind,lapply(1:nrow(data.mat1),function(i,n) cbind(rep(i,n),1:n),n = nrow(data.mat2)))
 }
 
 cat(paste("method=",method,"\nN1=",nrow(data.mat1),"\nN2=",nrow(data.mat2),"\nTotal pairs=",nrow(cor.pairs),"\n",sep = ""))
 
 row.names <- rownames(data.mat1);
 col.names <- rownames(data.mat2);
 cat("Calculating correlation coefficient and respective p-value...\n")
 output <- apply(cor.pairs,1,function(ij,mat1,mat2,alternative,method,ub) {
                                      out <- cor.test(x = mat1[ij[1],],y = mat2[ij[2],],alternative = alternative,method = method,use = ub);
							       	  out <- c(out$estimate,out$p.value);
									  names(out) <- c("rho","p.value")
									  return(out)},mat1 = data.mat1,mat2 = data.mat2,alternative = alternative,method = method,ub = use.obs)
 if (nrow(output) != nrow(cor.pairs)) output <- t(output)
 output <- as.data.frame(output)
 output <- data.frame(row = cor.pairs[,1],col = cor.pairs[,2],output)
 
 output <- list(correlation = output,row.names = row.names,col.names = col.names);
 return(output)
}

test.pairwiseCor.par <- function(data.mat1,data.mat2 = NULL,n.cores,
                                 alternative = "two.sided",method = "pearson",use.obs = "na.or.complete")
{
 cat("##### Pairwise Correlation Analysis ######\n")
 
 if (is.null(data.mat2))
 {
  cor.pairs <- do.call(rbind,lapply(1:(nrow(data.mat1)-1),function(i,n) cbind(rep(i,(n-i)),(i+1):n),n = nrow(data.mat1)))
  data.mat2 <- data.mat1;
 }else{
  if (ncol(data.mat1) != ncol(data.mat2)) stop("data.mat1 and data.mat2 do not have same number of columns.");
  cor.pairs <- do.call(rbind,lapply(1:nrow(data.mat1),function(i,n) cbind(rep(i,n),1:n),n = nrow(data.mat2)))
 }
 
 cat(paste("method=",method,"\nN1=",nrow(data.mat1),"\nN2=",nrow(data.mat2),"\nTotal pairs=",nrow(cor.pairs),"\n",sep = ""))
 
 row.names <- rownames(data.mat1);
 col.names <- rownames(data.mat2);
 cat("Calculating correlation coefficient and respective p-value...\n")
 # split pairs into n.cores
 dn <- ceiling(nrow(cor.pairs)/n.cores)
 fact <- factor(do.call(c,lapply(1:ceiling(nrow(cor.pairs)/dn),function(n,dn) rep(n,dn),dn = dn))[1:nrow(cor.pairs)])
 split.pairs <- lapply(split(1:nrow(cor.pairs),fact),function(ii,m) m[ii,],m = cor.pairs);rm(cor.pairs,fact,dn)
 
 output <- foreach(cpair = split.pairs) %dopar% {
 
                  out <- apply(cpair,1,function(ij,mat1,mat2,alternative,method,use.obs) {
                                                    out <- cor.test(x = mat1[ij[1],],y = mat2[ij[2],],
                                                                    alternative = alternative,
                                                                    method = method,use = use.obs);
							                    	out <- c(out$estimate,out$p.value);
													names(out) <- c("rho","p.value")
													return(out)},mat1 = data.mat1,mat2 = data.mat2,alternative = alternative,method = method,use.obs = use.obs)
				  if (nrow(out) != nrow(cpair)) out <- t(out)
				  out <- as.data.frame(out)
                  out <- data.frame(row = cpair[,1],col = cpair[,2],out)
				  return(out)
		   }

 output <- do.call("rbind.data.frame",output)
 output <- list(correlation = output,row.names = row.names,col.names = col.names);
 return(output)
}

calculate.correlation <- function(datExpr,doPerm = 100,doPar = FALSE,num.cores = 8,method = "pearson",
                                  use.obs = "na.or.complete",
                                  FDR.cutoff = 0.05,n.increment = 100,is.signed = FALSE,
                                  output.permFDR = TRUE,output.corTable = TRUE,saveto = NULL)
{
 # Input
 # datExpr = expression matrix (row = probe,column = sample)
 # doPerm = number of permutations
 # num.cores = number of cores to call when parallel computing 
 # FDR.cutoff = FDR cut-off
 # n.increment = number of increments to segment (rho_min,rho_max).
 # is.signed = FALSE (unsigned correlation), TRUE (signed correlation)
 # output.permFDR = TRUE (output permutation indices into .txt file)
 # output.corTable = TRUE (output final correlation into .txt file)
 # saveto = designated save folder
 if (doPerm == 0)
 {
  if (!doPar)
  {
   cor.output <- test.pairwiseCor(datExpr,method = method,use.obs = use.obs);
  }else{
   cor.output <- test.pairwiseCor.par(datExpr,n.cores = num.cores,method = method,use.obs = use.obs);
  }
  # extract significant correlation
  vertex.names <- cor.output$row.names
  edgelist <- cor.output[[1]];
  edgelist <- cbind(edgelist,fdr.q.value = p.adjust(edgelist$p.value,"fdr"));
  edgelist <- edgelist[which(edgelist$fdr.q.value < FDR.cutoff),]
  
  if (is.signed)
  {
   edgelist <- edgelist[order(edgelist[,3],decreasing = T),]
   edgelist <- data.frame(row = vertex.names[edgelist[[1]]],col = vertex.names[edgelist[[2]]],as.data.frame(edgelist[,3:ncol(edgelist)]));
  }else{
   sign <- rep("",nrow(edgelist));sign[which(edgelist[,3] < 0)] <- "negative";sign[which(edgelist[,3] > 0)] <- "positive";
   edgelist[,3] <- abs(edgelist[,3])
   edgelist <- edgelist[order(edgelist[,3],decreasing = T),]
   edgelist <- data.frame(row = vertex.names[edgelist[[1]]],col = vertex.names[edgelist[[2]]],as.data.frame(edgelist[,3:ncol(edgelist)]),sign = sign);
   
   # output results to files
   
  }
  
  if (output.corTable)
  {
    cat("- outputting correlation results...\n");
    if (!is.null(saveto))
    {
      write.table(edgelist,file = paste(saveto,"Data_Correlation.txt",sep = "/"),sep = "\t",row.names = F,col.names = T,quote = F)
    }else{
      write.table(edgelist,file = "Data_Correlation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
    }
  }
  
 }else{
  
  rho.output <- calculate.rho(datExpr,n.perm = doPerm,FDR.cutoff = FDR.cutoff,estimator = method,use.obs = use.obs,
                              rho.thresh =  seq(0,1,1/n.increment),sort.el = TRUE)
  
  if (output.permFDR)
  {
   if (!is.null(saveto))
   {
    write.table(rho.output$FDR,file = paste(saveto,"correlation_FDR_table.txt",sep = "/"),sep = "\t",row.names = F,col.names = T,quote = F)
   }else{
    write.table(rho.output$FDR,file = "correlation_FDR_table.txt",sep = "\t",row.names = F,col.names = T,quote = F)
   }
  }
  
  if (output.corTable)
  {
    cat("- outputting correlation results...\n");
    if (!is.null(saveto))
    {
      write.table(rho.output$signif.ijw,file = paste(saveto,"Data_Correlation.txt",sep = "/"),sep = "\t",row.names = F,col.names = T,quote = F)
    }else{
      write.table(rho.output$signif.ijw,file = "Data_Correlation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
    }
  }
  edgelist <- rho.output$signif.ijw;
 }
 
 
 
 
 return(edgelist)
}

count.obs_new <- function(rho_vec, 
                          rho_thresh_vec = seq(from = 0, to = 1, length.out = 101), 
                          direction      = "ascending"){
  
  if(direction == "ascending"){
    nb_obs_above <- unlist(lapply(rho_thresh_vec, function(x) sum(rho_vec > x)))
  } else if(direction == "descending"){
    nb_obs_above <- unlist(lapply(rho_thresh_vec, function(x) sum(rho_vec < x)))
  }
  out.df      <- data.frame(rho.cutoff   = rho_thresh_vec,
                            nb.obs.above = nb_obs_above,
                            nb.obs.total = length(rho_vec),
                            stringsAsFactors = F)
  
  out.df
  
}

calculate.correlation_new <- function(datExpr,
                                      doPerm            = 100,
                                      num.cores         = 8,
                                      method            = "pearson",
                                      use.obs           = "pairwise.complete.obs",
                                      FDR.cutoff        = 0.05,
                                      direction         = "absolute",
                                      rho.thresh.step   = 0.0025,
                                      output.permFDR    = TRUE,
                                      output.corTable   = TRUE,
                                      saveto            = NULL){
  
  
  # Notes: using direction = positive or negative corresponds to a one-side test (alternative = greater or less respectively)
  
  if(!(direction %in% c("absolute", "positive", "negative"))){
    stop(">> direction must be one of the following: absolute, positive, or negative.")
  }
  if(!is.null(saveto) & !dir.exists(saveto)){
    dir.create(saveto, recursive = T)
  }
  
  if(doPerm == 0){
    
    cat("Computing correlation and pValue...\n")

    # Use fast correlation and pValue implementation (the difference with the normal correlation is in the 1e-15/1e-16 order of magnitude)
    if(direction == "absolute"){
      cor.output <- WGCNA::corAndPvalue(x = t(datExpr), use = use.obs, alternative = "two.sided", method = method)
    } else if(direction == "positive"){
      cor.output <- WGCNA::corAndPvalue(x = t(datExpr), use = use.obs, alternative = "greater", method = method)
    } else if(direction == "negative"){
      cor.output <- WGCNA::corAndPvalue(x = t(datExpr), use = use.obs, alternative = "less", method = method)
    }
    
    cor.df                                 <- cor.output$cor
    cor.df[upper.tri(cor.df, diag = T)]    <- NA
    cor.df                                 <- reshape2::melt(cor.df)
    pval.df                                <- cor.output$p
    pval.df[upper.tri(pval.df, diag = T)]  <- NA
    pval.df                                <- reshape2::melt(pval.df)
    colnames(cor.df)                       <- c("col", "row", "rho")
    colnames(pval.df)                      <- c("col", "row", "p.value")
    # match rows instead of merging - faster
    stopifnot(identical(cor.df$row, pval.df$row) & identical(cor.df$col, pval.df$col))
    
    cat("Formatting objects...\n") 
    edgelist             <- type.convert(as.data.frame(cbind(cor.df, pval.df[, -c(1:2), drop = F])), as.is = T)
    edgelist             <- edgelist[!is.na(edgelist$rho), ]
    # edgelist             <- edgelist[edgelist$row != edgelist$col, ]
    edgelist$fdr.q.value <- p.adjust(edgelist$p.value, "fdr")
    edgelist             <- edgelist[edgelist$fdr.q.value < FDR.cutoff, ]
    edgelist             <- dplyr::select(edgelist, c("row", "col", "rho", "p.value", "fdr.q.value"))
    rownames(edgelist)   <- NULL
    
    if(direction == "absolute"){
      edgelist$rho <- abs(edgelist$rho)
    }
    
    if(direction %in% c("absolute", "positive")){
      edgelist <- edgelist[order(edgelist$rho, decreasing = T), ]
    } else if(direction == "negative"){
      edgelist <- edgelist[order(edgelist$rho, decreasing = F), ]
    }
    FDR.table <- NULL
    
    if(output.corTable){
      if(!is.null(saveto)){
        write.table(edgelist, file = paste(saveto, "Data_Correlation.txt", sep = "/"), sep = "\t", row.names = F, col.names = T, quote = F)
      } else{
        write.table(edgelist, file = "Data_Correlation.txt", sep = "\t", row.names = F, col.names = T, quote = F)
      }
    }
    
  } else{
    
    # Notes:
    # We operate here under the assumption that 
    # (1) the observed "real" correlation are a combination of False Positive 
    #     (FP) and True Positive (FP), i.e. FP + TP
    #     to optimize computations we split this computations in N permutations
    # (2) the permuted correlations are FP
    # (3) combining this, we compute the global False Discovery Rate for each 
    #     correlation threhsold in rho.thresh as:
    #     FDR = FP / (FP + TP)
    #     Ideally, a similar permutation approach could be used to compute an 
    #     empirical pValue for each edge (gene-gene correlation pair). This was
    #     the original adapted implementation here, available in a previous 
    #     commit of this forked MEGENA branch, but it is too computationally
    #     expensive.
    
    # build correlation cutoffs vector
    # round rho.thresh.step to divide the range [0, 1] into equal bins
    rho.thresh.step <- 1/floor(1/rho.thresh.step)
    if(direction == "absolute"){
      rho.thresh <- seq(from = 0, to = 1, by = rho.thresh.step)
    } else{
      rho.thresh <- seq(from = -1, to = 1, by = rho.thresh.step)
    }
    
    # Set random seed inside lapply to guarantee maximum reproducibility
    nc        <- ncol(datExpr)
    perm.ind  <- lapply(1:doPerm, function(i){
      set.seed(i)
      sample(1:nc, nc)
    })
    
    # Bootstrap
    # Use fast correlation implementation (the difference with the normal correlation is in the 1e-15/1e-16 order of magnitude)
    # Note: count false discovery numbers without computing a rate since the total number of observations is always the same
    cat("Computing bootstrap correlation and counting number of observations above correlation thresholds (assumed to be FP - using shuffled data, i.e. FP)...\n")
    perm.corr.counts <- parallel::mclapply(mc.cores = num.cores, X = perm.ind, FUN = function(x){
      # x <- perm.ind[[1]]
      
      if(direction == "absolute"){
        random.cor <- abs(WGCNA::cor(t(datExpr), t(datExpr[, x]), use = use.obs, method = method))
      } else{
        random.cor <- WGCNA::cor(t(datExpr), t(datExpr[, x]), use = use.obs, method = method)
      }
      
      # Note: the upper.tri section correspond to the correlation with the real gene vector (row name) and the shuffled gene vector (column name)
      random.cor <- as.vector(random.cor[upper.tri(random.cor)])
      
      if(direction %in% c("absolute", "positive")){
        count.obs_new(random.cor, rho.thresh, "ascending")
      } else if(direction == "negative"){
        out <- count.obs_new(random.cor, rho.thresh, "descending")
      }
    })
    
    # Compute observed correlation
    # Use fast correlation implementation (the difference with the normal correlation is in the 1e-15/1e-16 order of magnitude)
    cat("Computing observed correlation and counting number of observations above correlation thresholds (assumed to be total positive observations, i.e. FP + TN)...\n")
    if(direction == "absolute"){
      obs_corr <- abs(WGCNA::cor(t(datExpr), use = use.obs, method = method))
    } else{
      obs_corr <- WGCNA::cor(t(datExpr), use = use.obs, method = method)
    }
    
    obs_corr.vec  <- as.vector(obs_corr[upper.tri(obs_corr)])
    
    if(direction %in% c("absolute", "positive")){
      obs.counts.df <- count.obs_new(obs_corr.vec, rho.thresh, "ascending")
    } else if(direction == "negative"){
      obs.counts.df <- count.obs_new(obs_corr.vec, rho.thresh, "descending")
    }
    
    cat("Computing FDR and subsetting edges...\n")
    
    # Compute false positive (FP) rate
    count.fp.tot  <- apply(do.call(cbind, lapply(perm.corr.counts, function(x) x[, 2])), 1, sum)
    count.obs.tot <- apply(do.call(cbind, lapply(perm.corr.counts, function(x) x[, 3])), 1, sum)
    fpr.vec       <- count.fp.tot / count.obs.tot
    
    if(direction %in% c("absolute", "positive")){
      fpr.vec[1] <- 1
    } else if(direction == "negative"){
      fpr.vec[length(fpr.vec)] <- 1
    }
    
    
    # Compute true positive + false positive (TP + FP) rate
    tp.fp.r.vec   <- obs.counts.df$nb.obs.above / obs.counts.df$nb.obs.total
    
    # Compute empirical FDR pValue
    fdr.vec                   <- fpr.vec / tp.fp.r.vec
    fdr.vec[fpr.vec == 0]     <- 0
    fdr.vec[tp.fp.r.vec == 0] <- 0
    
    # adjust monotonicity and fdr range
    fdr.vec.monotonic                         <- fdr.vec
    fdr.vec.monotonic[fdr.vec.monotonic > 1]  <- 1
    
    # apply constraint that higher correlation threshold yields lower FDR than 
    # lower correlation thresholds (<< if corr increase, FDR decrease >>)
    if(direction %in% c("absolute", "positive")){
      
      mx <- fdr.vec.monotonic[1]
      for(i in 2:length(fdr.vec.monotonic)){
        if(fdr.vec.monotonic[i] > mx){
          fdr.vec.monotonic[i] <- mx 
        } else{
          mx <- fdr.vec.monotonic[i]
        }
      }
      
    } else if(direction == "negative"){
      
      mx <- fdr.vec.monotonic[length(fdr.vec.monotonic)]
      for(i in (length(fdr.vec.monotonic) - 1):1){
        if(fdr.vec.monotonic[i] > mx){
          fdr.vec.monotonic[i] <- mx 
        } else{
          mx <- fdr.vec.monotonic[i]
        }
      }
      
    }
    
    FDR.table <- data.frame(cut.off       = rho.thresh,
                            FPR           = fpr.vec,
                            PR            = tp.fp.r.vec,
                            FDR           = fdr.vec,
                            FDR_monotonic = fdr.vec.monotonic)
    
    
    # choose threshold respect to FDR.cutoff
    if(direction %in% c("absolute", "positive")){
      
      rho.cutoff <- min(FDR.table$cut.off[FDR.table$FDR_monotonic < FDR.cutoff])
      
      edgelist                 <- obs_corr
      edgelist[lower.tri(edgelist)] <- NA
      edgelist                 <- reshape2::melt(edgelist, na.rm = T)
      edgelist$Var1            <- as.character(edgelist$Var1)
      edgelist$Var2            <- as.character(edgelist$Var2)
      edgelist                 <- edgelist[!is.na(edgelist$value), ]
      edgelist                 <- edgelist[edgelist$Var1 != edgelist$Var2, ]
      colnames(edgelist)       <- c("row", "col", "rho")
      edgelist                 <- edgelist[edgelist$rho >= rho.cutoff, ]
      rownames(edgelist)       <- NULL
      edgelist                 <- edgelist[order(edgelist$rho, decreasing = T), ]
      
      
    list(signif.ijw = ijw,FDR = FDR.table)
      
    } else if(direction == "negative"){
      
      rho.cutoff <- max(FDR.table$cut.off[FDR.table$FDR < FDR.cutoff])
      
      edgelist                 <- obs_corr
      edgelist[lower.tri(edgelist)] <- NA
      edgelist                 <- reshape2::melt(edgelist, na.rm = T)
      edgelist$Var1            <- as.character(edgelist$Var1)
      edgelist$Var2            <- as.character(edgelist$Var2)
      edgelist                 <- edgelist[!is.na(edgelist$value), ]
      edgelist                 <- edgelist[edgelist$Var1 != edgelist$Var2, ]
      colnames(edgelist)       <- c("row", "col", "rho")
      edgelist                 <- edgelist[edgelist$rho <= rho.cutoff, ]
      rownames(edgelist)       <- NULL
      edgelist                 <- edgelist[order(edgelist$rho, decreasing = F), ]
      
    }
    
    if (output.permFDR){
      if (!is.null(saveto)){
        write.table(FDR.table, file = paste(saveto, "correlation_FDR_table.txt", sep = "/"), sep = "\t", row.names = F, col.names = T, quote = F)
      } else{
        write.table(FDR.table, file = "correlation_FDR_table.txt", sep = "\t", row.names = F, col.names = T, quote = F)
      }
    }
    
    if (output.corTable){
      cat("- outputting correlation results...\n");
      if (!is.null(saveto)){
        write.table(edgelist, file = paste(saveto,"Data_Correlation.txt", sep = "/"), sep = "\t", row.names = F, col.names = T, quote = F)
      }else{
        write.table(edgelist, file = "Data_Correlation.txt", sep = "\t", row.names = F, col.names = T, quote = F)
      }
    }
  }
  
  return(list(signif.ijw = edgelist,
              FDR        = FDR.table))
}

##########################
