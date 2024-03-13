
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



calculate.correlation_new <- function(datExpr,
                                      doPerm          = 100,
                                      doPar           = FALSE,
                                      num.cores       = 8,
                                      method          = "pearson",
                                      use.obs         = "na.or.complete",
                                      FDR.cutoff      = 0.05,
                                      n.increment     = 100,
                                      is.signed       = FALSE,
                                      output.permFDR  = TRUE,
                                      output.corTable = TRUE,
                                      saveto          = NULL){
  
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
  
  if(doPerm == 0){
    
    cor.output                             <- WGCNA::corAndPvalue(x = t(datExpr), use = use.obs, alternative = "two.sided", method = method)
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
    
    edgelist             <- type.convert(as.data.frame(cbind(cor.df, pval.df[, -c(1:2), drop = F])), as.is = T)
    edgelist             <- edgelist[!is.na(edgelist$rho), ]
    # edgelist             <- edgelist[edgelist$row != edgelist$col, ]
    edgelist$fdr.q.value <- p.adjust(edgelist$p.value, "fdr")
    edgelist             <- edgelist[edgelist$fdr.q.value < FDR.cutoff, ]
    edgelist             <- dplyr::select(edgelist, c("row", "col", "rho", "p.value", "fdr.q.value"))
    rownames(edgelist)   <- NULL
    edgelist             <- edgelist[order(edgelist$rho, decreasing = T), ]
    
    if(!is.signed){
      edgelist$rho <- abs(edgelist$rho)
    }
    
    if(output.corTable){
      if(!is.null(saveto)){
        write.table(edgelist, file = paste(saveto, "Data_Correlation.txt", sep = "/"), sep = "\t", row.names = F, col.names = T, quote = F)
      } else{
        write.table(edgelist, file = "Data_Correlation.txt", sep = "\t", row.names = F, col.names = T, quote = F)
      }
    }
    
  } else{
    
    # cat("Notice: is.signed parameter not applied, returning absolute value correlations...\n")
    
    nc        <- ncol(datExpr)
    perm.ind  <- lapply(1:doPerm, function(i){
      set.seed(i)
      sample(1:nc, nc)
    })
    
    # Compute observed correlation
    obs_corr <- cor(t(datExpr))
    obs_corr[lower.tri(obs_corr)] <- NA
    obs_corr <- type.convert(reshape2::melt(obs_corr), as.is = T)
    obs_corr <- obs_corr[!is.na(obs_corr$value), ]
    
    # Bootstrap
    boot_corrs <- parallel::mclapply(mc.cores = num.cores, X = perm.ind, FUN = function(x){
      out <- cor(t(datExpr), t(datExpr[, x]), method = method, use = use.obs)
      # Note: the upper.tri section correspond to the correlation with the real gene vector (row name) and the shuffled gene vector (column name)
      out[lower.tri(out)] <- NA
      out <- type.convert(reshape2::melt(out), as.is = T)
      out <- out[!is.na(out$value), ]
      out
    })
    # double-check that the dimnames of the data.frames are identical
    tmp.df <- do.call(cbind, lapply(boot_corrs, function(x) x$Var1))
    stopifnot(all(apply(tmp.df, 1, function(x) length(unique(x)) == 1)))
    tmp.df <- do.call(cbind, lapply(boot_corrs, function(x) x$Var2))
    stopifnot(all(apply(tmp.df, 1, function(x) length(unique(x)) == 1)))
    stopifnot(all(boot_corrs[[1]]$Var1 == obs_corr$Var1))
    stopifnot(all(boot_corrs[[1]]$Var2 == obs_corr$Var2))
    
    boot_corrs[2:length(boot_corrs)] <- lapply(boot_corrs[2:length(boot_corrs)], function(x) x[, -c(1:2), drop = F])
    boot_corrs                       <- do.call(cbind, boot_corrs)

    
    # # Compute p-value based on the alternative hypothesis
    # if(alternative == "two.sided") {
    #   p_value <- (sum(abs(boot_corrs) >= abs(obs_corr)) + 1) / (n_boot + 1)
    # } else if (alternative == "less") {
    #   p_value <- (sum(boot_corrs <= obs_corr) + 1) / (n_boot + 1)
    # } else if (alternative == "greater") {
    #   p_value <- (sum(boot_corrs >= obs_corr) + 1) / (n_boot + 1)
    # } else {
    #   stop("Invalid alternative hypothesis. Choose from 'two.sided', 'less', or 'greater'")
    # }
    
    # add index to make sure that mclapply doesn't switch the initial row orders in any operative system
    boot_corrs <- cbind(idx = 1:nrow(boot_corrs), boot_corrs)
    pVal_df    <- parallel::mclapply(mc.cores = num.cores, X = 1:nrow(boot_corrs), FUN = function(x){
      # x <- 1
      # + 1 to the numerator:   ensures that even if the observed correlation is not present in the bootstrapped samples, it still receives a non-zero probability, avoiding a p-value of 0. 
      # + 1 to the denominator: ensures that the p-value remains within the range [0, 1], even when the observed statistic is extreme.
      data.frame(idx = boot_corrs$idx[x], p.value = (sum(abs(boot_corrs[x, -c(1:3), drop = F]) >= abs(obs_corr$value[x])) + 1) / (doPerm + 1))
    })
    pVal_df <- as.data.frame(data.table::rbindlist(pVal_df))
    pVal_df <- pVal_df[match(boot_corrs$idx, pVal_df$idx), ]
    stopifnot(identical(pVal_df$idx, boot_corrs$idx))
    pVal_range <- base::range(pVal_df$p.value)
    stopifnot(pVal_range[1] >= 0 & pVal_range[2] <= 1)
    
    boot_corrs$p.value <- pVal_df$p.value
    boot_corrs$rho     <- obs_corr$value
    
    # format edgelist
    edgelist                <- dplyr::select(boot_corrs, c("Var1", "Var2", "rho", "p.value"))
    edgelist                <- edgelist[edgelist$Var1 != edgelist$Var2, ]
    edgelist$fdr.q.value    <- p.adjust(edgelist$p.value, "fdr")
    colnames(edgelist)[1:2] <- c("row", "col")
    edgelist                <- edgelist[edgelist$fdr.q.value < FDR.cutoff, ]
    rownames(edgelist)      <- NULL
    edgelist                <- edgelist[order(edgelist$rho, decreasing = T), ]
    
    if(!is.signed){
      edgelist$rho <- abs(edgelist$rho)
    }
    
    # format boot_corrs
    boot_corrs <- boot_corrs[, !(colnames(boot_corrs) %in% c("idx", "p.value", "rho"))]
    colnames(boot_corrs)[grepl("value", colnames(boot_corrs))] <- paste0("rho.perm", 1:sum(grepl("value", colnames(boot_corrs))))
    colnames(boot_corrs)[1:2] <- c("row", "col")
    boot_corrs                <- boot_corrs[boot_corrs$row != boot_corrs$col, ]
    
    if (output.permFDR){
      if (!is.null(saveto)){
        write.table(boot_corrs, file = paste(saveto, "correlation_FDR_table.txt", sep = "/"), sep = "\t", row.names = F, col.names = T, quote = F)
      } else{
        write.table(boot_corrs, file = "correlation_FDR_table.txt", sep = "\t", row.names = F, col.names = T, quote = F)
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
  
  return(edgelist)
}


##########################
