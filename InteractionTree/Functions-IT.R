# ####################################################################
# SUBGROUP ANALYSIS VIA INTERACTION TREES 
# ============================================
# WRITTEN BY XIAOGANG SU, PH.D. 
# EMAIL: XIAOGANGSU@GMAIL.COM
# WEBSITE: https://sites.google.com/site/xgsu00/
# ####################################################################

options(warn=-1)

# ==============================================
# FUNCTION rdat() SIMULATES A SIMPLE DATA SET
# ==============================================

rdat <- function(n=100, beta=c(2, 2, 2, 2, 2, 2), p= 5, 
	sigma=1, cut0=.5, digits=2, add.nominal=TRUE, n.LETTERS = 5)
{
    if (p < (length(beta)-2)/2 ) stop("Your p is too small.")  
    trt <- sample(0:1, n, replace=T)
    #### Generate Covariates
    X <- matrix(round(runif(n=n*p, min=0, max=1), digits = digits), 
                nrow = n, ncol=p)
    ### Generate Responses  
    mu <- beta[1] + beta[2]*trt + beta[3]*sign(X[,1]<=cut0) + 
      beta[4]*sign(X[,2]<=cut0) + beta[5]*sign(X[,1]<=cut0)*trt + 
      beta[6]*sign(X[,2]<=cut0)*trt
    y <- mu + rnorm(n, mean=0, sd=sigma) 
    
    ##### Output
    dat <- data.frame(cbind(y=y, trt=trt, X))
    if (add.nominal) {
      x0 <- sample(LETTERS[1:n.LETTERS], size=n, replace=TRUE)
      dat$x0 <- x0
    }
    names(dat) <- c("y", "trt", paste("x", 1:(NCOL(dat)-2), sep=""))
    return(dat)
}


# ------------------------------------------------------------------
# FUNCTION int2binary CONVERTS N (NATURAL NUMBER) INTO BINARY
# -------------------------------------------------------------------
N2binary <- function(x) {
  y <- rev(as.integer(intToBits(x))) 
  id.rm <- seq_len(match(1,y,length(y))-1)
  x <- paste(y[-id.rm], collapse = "")
  return(x)
}
# N2binary(123)



# ===============================================================================
# FUNCTION ttest() COMPUTES THE t TEST FOR INTERACTION WITH CONTINUOUS RESPONSE
# ===============================================================================

ttest <- function(dat, z, n0)
{  
    y <- dat$y; trt <- dat$trt; score <- NA; n <- nrow(dat)
    if (length(y)!=length(z)) stop("the vector z must have the same length as data.")       
    y11 <- y[trt==0 & z==1]; y12 <- y[trt==0 & z==0]
    y21 <- y[trt==1 & z==1]; y22 <- y[trt==1 & z==0]
    n11 <- length(y11); n12 <- length(y12); n21 <- length(y21); n22 <- length(y22); 
    # print(cbind(k, n11, n12, n21, n22))
    t <- NA
    if (min(n11, n12, n21, n22) >= n0) {    
        ybar11 <- mean(y11); ybar12 <- mean(y12); ybar21 <- mean(y21); ybar22 <- mean(y22); 
        dif <- (ybar11 + ybar22) - (ybar12 + ybar21)
        mse <- (sum((y11-ybar11)^2) + sum((y12-ybar12)^2)  + sum((y21-ybar21)^2) + sum((y22-ybar22)^2))/(n-4) 
        t <- dif/sqrt(mse *(1/n11 + 1/n12 + 1/n21 + 1/n22))
        score <- t^2
    }
   score
}




# ===========================================================================
# THE power.set() FUNCTION PROVIDES THE POWER SET FOR A CATEGORICAL VARIABLE
# ===========================================================================
          
power.set <- function(x) {
  if(length(x) == 0) return(vector(mode(x), 0))
  x <- sort(unique(x)); n <- length(x); K <- NULL
  for(m in x) K <- rbind(cbind(K, FALSE), cbind(K, TRUE))
  out <- apply(K, 1, function(x, s) s[x], s = x)
  out <- out[-c(1, length(out))]
  l <- length(out); i <- 1
  out[!sapply(out, length)>=ceiling(n/2+.5)]
}


# THIS FUNCTION as.numeric.factor() CONVERTS FACTOR INTO NUMERIC
as.numeric.factor <- function(x){as.numeric(levels(x))[x]}


# ================================================================
# ONE SINGLE SPLIT USING THE t TEST STATISTICS FOR INTERACTION
# ================================================================

# WHEN USING FOR RANDOM FORESTS, SET test=NULL. 
# min.ndsz= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR CLAIMING A TERMINAL NODE
# no= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR (n11, n10, n01, n00)
# split.var= ASKS FOR THE COLUMNS OF SPLITTING VARIABELS, INCLUDING BOTH CONTINUOUS AND CATEGORICAL ONES
# ctg= SPECIFIES THE COLUMNS OF CATEGORICAL VARIABLES
# max.depth= SPECIFIED THE MAXIMUM HEIGHT OF A TREE (ANOTHER WAY OF STOPPING THE GROWTH).
# mtry= SPECIFIES THE NUMBER OF COVARIATES IN THE RANDOMLY SELECTED SUBSETS FOR SPLITTING

partition.INT <- function(dat, test=NULL, name="1", min.ndsz=20, n0=5, 
	split.var, ctg=NULL, 
      max.depth=15, mtry=length(split.var))
{   
    # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES. 
    call <- match.call(); out <- match.call(expand = F)
    out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
    name.l <- paste(name, 0, sep=""); name.r <- paste(name, 1, sep="")
    n <- nrow(dat); 
    if (!is.null(test)) {n.test <- NA; score.test <- NA;}  ########## TEST SAMPLE ########  
    var <- vname <- NA; cut <- NA; max.score <- -1e20;   
    trt <- dat$trt; y <- dat$y; vnames <- colnames(dat)
    # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
    trt.effect <- NA; n.1 <- sum(trt==1); n.0 <- n-n.1
    if (min(n.1, n.0) >0) {trt.effect <- mean(y[trt==1]) - mean(y[trt==0])}
    # CONTROL THE MAX TREE DEPTH
    depth <- nchar(name) 
    if (depth <= max.depth && n >= min.ndsz) {
        m.try <- ifelse(is.null(mtry), length(split.var), mtry)  
        for(i in sample(split.var, size=m.try, replace=F)) {
            x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));  
            if(length(temp) > 1) { 
                if (is.element(i,ctg)) zcut <- power.set(temp)                                                  ############################ CLASS VARIABLE
                else zcut <- temp[-length(temp)]
                # print(i); print(temp); print(zcut)
                for(j in zcut) {
                    score <- NA
                    if (is.element(i,ctg)) {grp <- sign(is.element(x, j)); cut1 <- paste(j, collapse=" ")}      ############################ CLASS VARIABLE
                    else  {grp <- sign(x <= j); cut1 <- as.character(j)}                          
                    score <- ttest(dat, z=grp, n0)
                    # print(cbind(var=i, cut=j, score=score))
                    if (!is.na(score) && score >= max.score) {max.score <- score; var <- i; vname <- v.name; cut <- cut1; best.cut<-j} 
    }}}}
    if (!is.null(test)) { 
	n.test <- nrow(test); score.test <- NA;
    	if (!is.na(var)) {
        if (is.element(var,ctg)) grp.test <- sign(is.element(test[,var], best.cut))                              ############################
        else  grp.test <- sign(test[,var] <= best.cut)   
        score.test <- ttest(test, z=grp.test, n0=(n0/2))
    	  if (!is.na(score.test)){
        	out$name.l <- name.l; out$name.r <- name.r
        	out$left.test <- test[grp.test==1,  ]
        	out$right.test <- test[grp.test==0,  ]
        	if (is.element(var,ctg)) {                                                                               ############################                      
            	out$left  <- dat[is.element(dat[,var], best.cut),]      
            	out$right <- dat[!is.element(dat[,var], best.cut), ]}
        	else {
            	out$left  <- dat[dat[,var]<= best.cut,]
            	out$right <- dat[dat[,var]> best.cut, ]
		}
	  }
    	  else {var <- NA; vname <- NA; cut <- NA;  max.score <- NA}
    	  out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
		var = var, vname=vname, cut= cut, score=ifelse(max.score==-1e10, NA, max.score), 
		score.test, size.test=n.test)
	}
	else {
		out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
			var = NA, vname=NA, cut= NA, score=NA, 
			score.test=NA, size.test=n.test)
	}
    }	else {
	if (!is.na(var)) {
	   out$name.l <- name.l; out$name.r <- name.r
	   if (is.element(var,ctg)) {                                                                               ############################                      
            out$left  <- dat[is.element(dat[,var], best.cut),]      
            out$right <- dat[!is.element(dat[,var], best.cut), ]}
         else {
            out$left  <- dat[dat[,var]<= best.cut,]
            out$right <- dat[dat[,var]> best.cut, ]
	   }
	   out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
				var = var, vname=vname, cut= cut, score=ifelse(max.score==-1e10, NA, max.score))
	}
   	else {
	   out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
		var=NA, vname=NA, cut=NA, score=NA)
 	}
   }
   out 
}

# TEST CODE 
# split <- partition.INT(dat, test=dat, name="1", min.ndsz=20, n0=5, 
#	split.var, ctg=NULL, 
#      max.depth=15, mtry=length(split.var))
# split$info






# =================================================
# THE grow.INT() FUNCTION CONSTRUCTS A LARGE TREE 
# =================================================

grow.INT <- function(data, test=NULL, min.ndsz=20, n0=5, 
	split.var, ctg=NULL, 
	max.depth=15, mtry=length(split.var))
{
    out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
    list.nd <- list(data); 
    if (!is.null(test)) list.test <- list(test)
    name <- 1
    while (length(list.nd)!=0) {    
      for (i in 1:length(list.nd)){
	  # print(i)
        if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
		test0 <- NULL
		if (!is.null(test)) test0 <- list.test[[i]]
		split <- partition.INT(list.nd[[i]], test0, name[i], min.ndsz=min.ndsz, 
                n0=n0, split.var=split.var, ctg=ctg, max.depth=max.depth, mtry=mtry)
        	out <- rbind(out, split$info)
  	  	if (!is.null(split$left) && is.null(test)) {
			temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
			temp.name <- c(temp.name, split$name.l, split$name.r)
		
	  	} else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
			temp.list <- c(temp.list, list(split$left, split$right))
			temp.name <- c(temp.name, split$name.l, split$name.r)
			temp.test <- c(temp.test, list(split$left.test, split$right.test))
		}
	  }
	}
      list.nd <- temp.list; list.test <- temp.test; name <- temp.name
      temp.list <- temp.test <- temp.name <- NULL
    }   
    out$node <- as.character(out$node)
    out <- out[order(out$node), ]
    out
}

# grow.INT(data=dat, test=dat, min.ndsz=20, n0=5, 
#	split.var, ctg=NULL, 
#	max.depth=5)



# ==========================================================
# FUNCTION de() FINDS ALL THE DESCENDENTS OF NODE x IN tree
# ==========================================================


de <- function(x, tree)
{
    if(length(x) != 1) stop("The length of x in function de must be 1.")    
    y <- tree$node;  de <- NA
    if(sum(match(x, y), na.rm = T) != 0) {
        temp <- 1:length(y)
        start <- match(x, y) + 1    
        end <- length(y)
        if(start <= length(y) & nchar(y[start]) > nchar(x)) {
            temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
            if(!is.na(temp1)) end <- temp1
            de <- y[start:end]
    }}
    de
}



# ############################################################################
# Pruning and Size Selection Based on LeBlanc and Crowley (JASA, 1992)
# ############################################################################


# =================================
# METHOD I: THE TEST SAMPLE METHOD
# =================================

# -----------------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE TREE SIZE SELECTION VIA TEST SAMPLE METHOD
# -----------------------------------------------------------------------------

prune.size.testsample <- function(tre)
{
     out <- as.list(NULL)     
     ntest <- as.numeric(tre[1, ncol(tre)])
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre$var)); subtree <- 1            
     a <- cbind(Ga.2=2, Ga.3=3, Ga.4=4, Ga.BIC=log(ntest))
     max.Ga <- rep(-1e20, 4); size <- rep(0, 4); btree <-as.list(1:4) 
     while (n.tmnl > 1 ) {
            # print(tre)
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
            G <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T); 
            Ga <- G - a*l 
            for (k in 1:4){if (Ga[k] > max.Ga[k]) {max.Ga[k] <- Ga[k]; size[k] <- n.tmnl; btree[[k]] <- tre}}                        
            result <- rbind(result, c(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, Ga))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score", "score.test")] <- NA
            n.tmnl <- sum(is.na(tre$cut))
            if (n.tmnl ==1) {for (k in 1:4){if (0 > max.Ga[k]) {max.Ga[k] <- 0; size[k] <- 1; btree[[k]] <- tre}}}
            subtree <- subtree + 1          
      }
      # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
      result <- rbind(result, c(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                    size.tmnl=1, alpha=9999, G=0, Ga.2=0, Ga.3=0, Ga.4=0, Ga.BIC=0))     
      result <- data.frame(result)
      names(result) <- c("subtree", "node.rm", "size.tree", "size.tmnl", "alpha", 
                          "G", "Ga.2", "Ga.3", "Ga.4", "Ga.BIC") 
      out$result <- result; out$size <- size; out$btree <- btree
      out 
}







# =================================
# METHOD II: THE BOOTSTRAP METHOD
# =================================

# -----------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE BOOTSTRAP TREE SIZE SELECTION METHOD
# -----------------------------------------------------------------------

prune.size <- function(tre)
{
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre$var)); subtree <- 1            
     while (n.tmnl > 1 ) {
            # if (n.tmnl==5) {btre <- tre; print(btre)}
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
            G <- sum(as.numeric(as.vector(tre$score)), na.rm=T);
            G.test <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T)
            result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, G.test=G.test))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score", "score.test")] <- NA
            n.tmnl <- sum(is.na(tre$cut))
            subtree <- subtree + 1          
      }
      # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                    size.tmnl=1, alpha=9999, G=0, G.test=0))    
      result <- as.data.frame(result)
      result
}


# ==========================================================================================  #
#  FUNCTIONS RELATED TO THE PRUNING AND THEN BOOTSTRAP FOR TREE SIZE SELECTION
# ==========================================================================================  #

# OPTION LeBlanc IS TO APPLY THE WHOLE SAMPLE (TRUE) OR THE OUT-OF-BAD SAMPLE (FALSE) IN THE BOOTSTRAP PROCEDURE
# OPTION min.boot.tree.size IS TO MAKE SURE A NON-NULL TREE IS GROWN AT EACH BOOTSTRAP 11/1/2007

boottrap.grow.prune <- function(B=30, data, N0=20, n0=5, 
	split.var, ctg=NULL, max.depth=10, 
	mtry=length(split.var), LeBlanc=TRUE, min.boot.tree.size=1)  
{
    call <- match.call(); out <- match.call(expand = F)
    out$boot.tree <- out$boot.prune <- out$... <- NULL
    time.start <- date()
    tree0 <- grow.INT(data, data, min.ndsz=N0, n0=n0, split.var=split.var, ctg=ctg, 
            max.depth=max.depth, mtry=mtry);  
	if (NROW(tree0)==1) {print(tree0); stop("Your intital tree is a NULL tree! There is no need for tree size selection.")}
    # print(tree0);  
    prune0 <- prune.size(tree0); 
    boot.tree <- list(tree0); boot.prune <- list(prune0) 
    b <- 1
    while (b <= B) {
        print(paste("###################### b = ", b, " ###########################", sep=""))
        # SAMPLING OBSERVATION
        samp <- sample(1:nrow(data), size=nrow(data), replace=T) 
        dat <- data[samp, ];     
        dat.oob <- data[-unique(samp),]
        n.oob <- nrow(dat.oob); # print(n.oob)        
        if (LeBlanc) {tre <- grow.INT(dat, data, min.ndsz=N0, n0=n0, split.var=split.var, 
                            ctg=ctg, max.depth=max.depth, mtry=mtry)}
        else {tre <- grow.INT(dat, dat.oob, min.ndsz=N0, n0=n0, split.var=split.var, ctg=ctg, 
                max.depth=max.depth, mtry=mtry)}
        print(tre)        
        if (nrow(tre)> min.boot.tree.size) {
            boot.tree <- c(boot.tree, list(tre)); 
            prune <- prune.size(tre); # print(prune)
            boot.prune <- c(boot.prune, list(prune));
            b <- b+1
        }
    }
    time.end <- date(); 
    print(paste("The Start and End time for ", B, "bootstrap runs is:"))
    print(rbind(time.start, time.end))
    out$boot.tree <- boot.tree
    out$boot.prune <- boot.prune
    # THE INITIAL LARGE TREE
    out$initial.tree <- tree0;    
    out
}   





# ========================================================================
# FUNCTION obtain.btree() OBTAINS THE BEST SUBTREE WITH KNOW SIZE bsize=
# ========================================================================

obtain.btree <- function(tre, bsize=6)
{
	btre <- NULL
	if (bsize==1) { btre <- tre[1,]; btre[, c("var", "cut", "score", "score.test")] <- NA}  
    	else if (bsize <1) stop("THE BEST TREE SIZE bsize= MUST BE >=1!")
	else {
     		n.tmnl <- sum(is.na(tre$cut)); indicator <- T  
		if (bsize > n.tmnl) stop("THE BEST TREE SIZE bsize PROVIDED IS LARGER THAN THE FULL TREE THAT YOU HAVE GROWN.")          
     		while (n.tmnl >= bsize && indicator ==T) {
			# print(tre); print(cbind(n.tmnl, bsize))
            	internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            	r.value <- 1:l
            	for(i in 1:l) {
                		branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                		score <- as.numeric(as.vector(branch$score))
                		r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            	}
            	alpha <- min(r.value)
            	nod.rm <- internal[r.value == alpha]; 
            	tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            	tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score", "score.test")] <- NA
            	n.tmnl <- sum(is.na(tre$cut))      
			# print(cbind(n.tmnl, bsize))  
			# if (n.tmnl==bsize) {btre <- tre; print(btre); indicator <- F}  
    	if (n.tmnl==bsize) {btre <- tre; indicator <- F}  
      	}
	}
	if (is.null(btre)) print(paste("The optimally-pruned subtree sequence does not have a subtree of bsize = ", bsize, sep="")) 
      return(btre)
}



# ===============================================================
# SELECT THE BEST SUBTREE SIZE AND OBTAIN THE BEST SUBTREE MODEL
# ===============================================================

bootstrap.size <- function(boot.prune, n, plot.it=TRUE, filename=NULL, horizontal=T)
{   
    OUT <- as.list(NULL)
    #  COMPUTE THE ALPHA PRIME'S
    prune0 <- boot.prune[[1]] 
    n.subtree <- nrow(prune0)
    alpha <- as.numeric(as.vector(prune0$alpha));
    # temp <- c(alpha[1], alpha[-length(alpha)])  	##############
    temp <- c(0, alpha[-length(alpha)])  			############## CHANGE FIRST VALUE OF ALPHA TO 0
    alpha.prime <- sqrt(alpha*temp)  
    # cbind(alpha,  alpha.prime=prune0$alpha.prime)
    b <- length(boot.prune)
    G <- as.numeric(as.vector(prune0$G)); 
    size.tmnl <- as.numeric(as.vector(prune0$size.tmnl)); 
    subtree <- as.numeric(as.vector(prune0$subtree)); 
    # tree.penalty <- log(nrow(teeth))
    G.a <- matrix(0, n.subtree, 5)
    penalty <- c(0, 2:4, log(n))
    for (i in 1:n.subtree) {
        a.i <- alpha.prime[i]
        bias <- 0
        for (j in 2:b){
            prune.bs <- boot.prune[[j]]
            alpha.bs <- as.numeric(as.vector(prune.bs$alpha)); 
            g <- as.numeric(as.vector(prune.bs$G)); 
            g.test <- as.numeric(as.vector(prune.bs$G.test)); 
            indx <- 1
            if (sum(alpha.bs <= a.i)>0) {          
                temp1 <- which.max(which(alpha.bs<=a.i))
                indx <- ifelse(is.null(temp1), 1, temp1)
            }
            temp2 <- (g-g.test)[indx]
            bias <- bias + temp2 
            # print(cbind(i, a.i, j, bias, indx, temp2))
        }
        G.honest <- G[i] - bias/(b-1) 
        G.a[i,] <- G.honest - penalty*(size.tmnl[i]-1)
    }
    out <- data.frame(cbind(size.tmnl, G.a))
    colnames(out) <- c("tmnl.size", "G", "G.2", "G.3", "G.4", "G.log(n)")
    G.a <- out

    # PLOT THE G.a WITH DIFFERENT CHOICES OF a
    if (plot.it) {
	if (!is.null(filename)) postscript(file=filename, horizontal=horizontal)
	par(mfrow=c(1, 1), mar=rep(4, 4))   ##################### SET THE PLOTTING PARAMETERS
	n.subtrees <- nrow(G.a)
	subtree.size <- G.a[,1]
	min.x <- min(subtree.size); max.x <- max(subtree.size)
	min.y <- min(G.a$"G.log(n)"); max.y <- max(G.a$G.2)
	plot(x=c(min.x, max.x), y=c(min.y, max.y), type="n", xlab="tree size", ylab="G(a)")
	for (j in 3:6) lines(subtree.size, G.a[,j], lty=j-1, col=j-1, lwd=2)
	legend(x=min.x, y=(max.y+min.y)/2, lty=2:5, col=2:5, legend=c("G(2)", "G(3)", "G(4)", "G(ln(n))")) 
	if (!is.null(filename)) dev.off()
    }   
    # OBTAIN THE BEST TREE SIZE 
    bsize <- btree <- as.list(NULL)
    Ga.cols <- c("G.2", "G.3", "G.4", "G.log(n)")
    for (j in Ga.cols) {
	best.size <- subtree.size[which.max(G.a[,j])]
	bsize[[j]] <- best.size     
      btree[[j]] <- obtain.btree(tree0, bsize=best.size)
    }
    OUT$G.a <- G.a; OUT$bsize <- bsize; OUT$btree <- btree
    return(OUT)
}	


# ==============================================================
# FUNCTION send.down() RUNS A TREE STRUCTURE DOWN A DATA SET
# ==============================================================

send.down <- function(data, tre, char.var=1000)
{
    call <- match.call(); out <- match.call(expand = F)
    out$tree <- out$data <- out$... <- NULL
    dat <- cbind(data, node=1); tre <- cbind(tre, n.test=NA)
    cut.point <- as.vector(tre$cut); 
    split.var <- as.numeric(as.vector(tre$var)); 
    for (i in 1:nrow(tre)){
        in.node <- (dat$node)==(tre$node[i]);
        tre$n.test[i] <- sum(in.node)                       
        if (!is.na(split.var[i])){
            # print(cbind(i, var=tre$var[i], cut=tre$cut[i]))
            var.split <- dat[,split.var[i]]; 
            cut <- cut.point[i]
            if (!is.element(split.var[i], char.var)) { 
                cut1 <- as.numeric(cut)    
                l.nd <- dat$node[in.node & var.split <= cut1] 
                r.nd <- dat$node[in.node & var.split > cut1]
                dat$node[in.node & var.split <= cut1] <- paste(l.nd, 0, sep="")
                dat$node[in.node & var.split >  cut1] <- paste(r.nd, 1, sep="")  
            }
            else {
                var.split <- as.character(var.split)
                cut1 <- unlist(strsplit(as.character(cut), split=" ")) #####################
                l.nd <- dat$node[in.node & is.element(var.split, cut1)] 
                r.nd <- dat$node[in.node & !is.element(var.split, cut1)]                  
                dat$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 0, sep="")  
                dat$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 1, sep="")}                   
    }}
    # print(tre)
    out$data <- dat
    out$tree <- tre
    out 
}


# ================================================================
# PLOTTING IT TREE STRUCTURE, MODIFIED FROM PETER CALHOUN'S CODE
# ================================================================

plot.tree <- function(tree, cols.nominal=NULL, color0=NULL, 
                      leaf.label=FALSE, pch.leaf=21, plot.score=TRUE, 
                      textDepth=10, lines="rectangle", digits=4)
{
  depth <- max(nchar(tree[,1]))
  par(xaxs='i')
  par(mar=c(1,1,1,1))
  par(xpd=TRUE)
  plot(1, type="n", xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), axes=FALSE,xaxs="i",yaxs="i")
  nodes <- tree$node
  # nodesBin <- gsub("1", "0", nodesBin)
  # nodesBin <- gsub("2", "1", nodesBin)
  nodesBin <- nodes
  lastObs <- nchar(nodesBin)
  nodesBin <- substr(nodesBin,2,lastObs)
  var <- tree$vname 
  col.var <- tree$var   
  suppressWarnings(cut <- as.character(tree$cut)) 
  size <- tree$size  
  m <- tree$trt.effect
  if (plot.score) suppressWarnings(score <- round(as.numeric(as.character(tree$score)), digits = digits))
  if (is.null(color0)) {
    if (!is.null(tree$color))  color0 <- tree$color
    else color0 <- rep("chocolate3", length(nodesBin))
  } else if (length(color0) != NROW(tree)) 
    stop("color0 must have the same length as nrow(tree).")
  
  for(i in 1:length(nodesBin)){
    nChar <- nchar(nodesBin[i])
    # print(cbind(i, nChar))
    if(!is.na(var[i])){
      j <- ifelse(i==1, i+1, i)
      if(lines=="rectangle"){
        # print(cbind(i, nodesBin[i], strtoi(nodesBin[i])))
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      } else if(lines=="triangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      }         
      if(nChar <= textDepth){ 
        # print(cbind(i=i, var[i], cut[i]))
        if (plot.score) {
          score.i <- as.character(score[i])
          text(x=(1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2), 
               y=(depth-nChar)/(depth+1)+1/(depth+20)-0.03, 
               labels=substitute(paste(bold(score.i))), cex=1, col="gray50")
        }
        # HANDLE SPLITS WITH NOMINAL OR ORDINAL X
        if (is.na(as.numeric(as.character(cut[i]))) | 
            (!is.null(cols.nominal) && is.element(col.var[i], cols.nominal))) {
          cutpoint <- unlist(strsplit(as.character(cut[i]), split=" "))
          cutpoint <- paste("{", paste(cutpoint, collapse=" "), "}", sep="")
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2), 
               (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(var[i]))%in%.(cutpoint)),cex=1, col="deepskyblue4")
        } else {
          cutpoint <- round(as.numeric(cut[i]), digits=digits)
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2), 
               (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(var[i]))<=.(cutpoint)),cex=1, col="deepskyblue4")
        }
      } 
    } else {
      if(nChar <= textDepth){
        m.i <- round(as.numeric(m[i]), digits=digits)
        points((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+0.01, 
               pch=pch.leaf, col=color0[i], bg=color0[i], cex=3)
        # points((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+0.01, 
        #       pch=21, col="gray65", bg="gray65", cex=3)
        
        if (leaf.label) {
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.02, 
               nodes[i],cex=1, offset=1, col=color0[i])
        } else {
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.02, 
               paste("n=", size[i], sep=""),cex=1, offset=1, col="dimgray")
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.05, 
               paste("d=", m.i, sep=""), cex=1, offset=1, col="dimgray")
        }
      }
    }
  }
}






# =================================================================
# Generating LaTeX CODES FOR PLOTTING TREES WITH TreeTex PACKAGE 
# =================================================================

latex.plot <- function(file="tree-code.tex", tree, group=rep("I", n.node), cat.var=c(3))
{
  n.node <- nrow(tre)
  sink(file=file)
#  cat("\n \\begin{Tree} \n")
  for (i in n.node:1) {
    node <- tre[i, ] 
    if (is.na(node$var))  {
      cat(paste("\\node{\\external\\type{circle}\\cntr{", node$size, "}\\lft{\\textbf{", group[i], "}}} \n", sep="")) 
    } else {
      if (!is.element(node$var, cat.var)) { 
        cat(paste("\\node{\\type{frame}\\cntr{$\\texttt{V", node$vname, "} \\leq ", node$cut, "$}} \n", sep=""))
      } else {
        # cut <- node$cut
        cut1 <- unlist(strsplit(as.character(node$cut), split=" "))
        cat(paste("\\node{\\type{frame}\\cntr{\\texttt{V", node$vname, "} $ \\in \\{", paste(cut1, collapse=","), "\\} $}} \n", sep="")) 
    }}
  }
  sink()  
} 
# latex.plot(file="tree-code.tex", tree=tree0, group=rep("I", n.node), cat.var=c(3))


# =========================================================================================
# FUNCTION plot.tree.latex() GENERATES LATEX CODES FOR PLOTTING TREE IN PACKAGE pstricks
# =========================================================================================
# \usepackage{pstricks,pst-node,pst-tree}  % Plotting Trees

plot.tree.latex <- function(tree, file="tree-code.tex", digits=5)
{
	n.node <- nrow(tree)
	sink(file=file)
	# DEFINE SYMBOLS
	cat("\\begin{figure} \\centering \n")
	cat("\\newcommand{\\Tgcircle}{\\Tcircle[linecolor=black, fillcolor=gray]} \n")
	cat("\\newcommand{\\Tgoval}{\\Toval[linecolor=black, fillcolor=gray]} \n")
	cat("\\newcommand{\\Tgframe}[1]{\\Tr{\\psframebox[linecolor=black, fillstyle=solid, fillcolor=orange]{#1}}} \n")
	# OPTION
	cat("\\psset{nodesep=0.7pt, linecolor=black, treesep=1.2cm, levelsep=1.8cm} \n")
	I0 <- i0 <- NULL
	for (i in 1:n.node) {
		node.i <- as.character(tree[i, 1])
		de.i <- de(node.i, tree)
		blanks <- paste(rep(" ", (nchar(node.i)-1)*8), sep="")  # 8 SPACES IN ONE TAB		
 		n.i <- tree$size[i]
		lable.i <- ifelse(!is.null(tree$trt.effect), tree$trt.effect[i], " "); 	##### THIS MEASURE MAY BE DIFFERENT FOR OTHER TYPES OF TREES
		if (!is.na(de.i[1])) {	# INTERNAL NODE
			if (nchar(node.i)==1 ||  substr(node.i, nchar(node.i), nchar(node.i))=="2") 
				cat(blanks, "\\pstree{\\Tgcircle{~~}} \n", sep = "")
			else cat(blanks, "\\pstree{\\Tgcircle{~~} \\tlput{\\color{blue}", rule.i, "\\hspace{-.6in}}} \n", sep = "")				
			cat(blanks, "{ \n", sep = "") 
			I0 <- c(I0, i)
			i0 <- c(i0, i + length(de.i))
			# UPDATE THE SPLITTING RULE
			vname.i <- tree$vname[i]; cut.i <- strtrim(as.character(tree$cut[i]), width=digits); 
			rule.i <- paste("\\texttt{", vname.i, "}", "$\\leq ", cut.i, "$", sep="")
		} else if (substr(node.i, nchar(node.i), nchar(node.i))=="1") { # TERMINAL NODE
			cat(blanks, "\\Tgframe{",  lable.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{", n.i, "}}",
              		"\\tlput{\\color{blue} ", rule.i, " \\hspace{-.3in}} \n", sep = "")
		} else cat(blanks, "\\Tgframe{",  lable.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{",  n.i, "}} \n", sep = "")
		if (is.element(i, i0)) {
			rep0 <- rep("}", sum(i0==i))
			node.i0 <- as.character(tree[I0[i0==i][1] , 1])
			blanks0 <- paste(rep(" ", (nchar(node.i0)-1)*8), sep="")  		
			cat(blanks0, rep0, "\n", sep = "") 
		}
	}
	cat("\\end{figure} \n")
	sink()  
}
# plot.tree.latex(btre.SSS, file="tree-code.tex", digits=5)








# ===============================================
# AMALAGAMATION FOLLOWING CIAMPI ET AL.'S IDEA 
# ===============================================

#  char.var= asks for columns of categorical or character covariates
#  n.group= is the number of final groups desired

amalagamation <- function(dat, btree, char.var =NULL, n.groups=3, details=T)
{
	OUT <- as.list(NULL)
	dat1 <- send.down(data=dat, tre=as.data.frame(btree), 
		char.var=char.var)$data
	out <- NULL
	nd <- sort(unique(dat1$node)); 
	i0 <- j0 <- 1
	step <- 1; action <- "Merged"
	if (length(nd) <= n.groups) stop("# Final groups > # Terminal Nodes. No amalgamation necessary!")
	while (length(nd) > n.groups) {
    		dat1$node[dat1$node==j0] <- i0
    		nd <- as.numeric(sort(unique(dat1$node))); 
		if (details) print(nd)
		n.nodes <- length(nd)
    		i0 <- j0 <- 1
    		min.score <- 1e20
    		for (i in nd) {
        		for (j in nd){
            		if (i<j) {
                			if (details) print(cbind(i, j))
                			dat.tmp <- dat1[dat1$node==i | dat1$node==j, c("y", "trt", "node")] 
                			z <- ifelse(dat.tmp$node==i, 1, 0); 
                			score <- ttest(dat.tmp, z, n0=2)                 
                			if (details) print(cbind(i, j, score)) 
                			if (score < min.score) {min.score <- score; i0 <- i; j0 <- j} 
            		}
        		}
    		}
		if (n.nodes==n.groups) action <- "Not Merged"
    		out <- rbind(out, cbind(step=step, n.nodes=n.nodes, i0, j0, action, min.score))
    		if (details) print(cbind(i0, j0, min.score))
		step <- step +1
	}
	out <- as.data.frame(out)
	OUT$merge.info <- out
	OUT$data <- dat1
	return(OUT)
}


# =====================================================================
# MATRIX OF PAIRWISE t TEST FOR INTERACTION BETWEEN EVERY TWO GROUPS
# =====================================================================

paired.interaction <- function(dat, log.worth=T, col.group, n0=2){
	colname.group <- colnames(dat)[col.group]
	subgroup <- sort(unique(dat[, col.group])); 
	n.groups <- length(subgroup)
	P0 <- matrix(1, nrow=n.groups, ncol=n.groups) 
	for (i in 1:n.groups) {
		for (j in 1:n.groups) {
			if (i <j) {
				grp.i <- subgroup[i]; grp.j <- subgroup[j]; 
				dat.tmp <- dat[is.element(dat[, col.group], c(grp.i, grp.j)), c("y", "trt", colname.group)] 
				z <- ifelse(dat.tmp[, colname.group]==grp.i, 1, 0); 
                		score <- ttest(dat.tmp, z, n0=n0)         
				df.t <- nrow(dat.tmp) - 4
				p.value <- 2*pt(sqrt(score), df=df.t, lower.tail=F) 
				P0[i,j] <- P0[j, i] <- p.value
	}}}
	row.names(P0) <- colnames(P0) <- subgroup 
	if (log.worth) P0 <- - log10(P0)
	return(P0)
}



# ====================================================
# SUMMARIZE THE FINAL SUBGROUPS (AFTER AMALGAMATION)
# ====================================================

# THE ARGUMENT dat= COMES FROM OUTPUT OF amalgamation(). IT SHOULD HAS ONE COLUMN CALLED "subgroup".

summary.subgroups <- function(dat, var.equal=T, digits=3, col.grp){
	subgroup <- dat[, col.grp]
	levels.subgroup <- sort(unique(subgroup)); 
	n.groups <- length(levels.subgroup)
	out <- NULL
	# SUMMARY STAT AT EACH SUBGROUP
	for (i in levels.subgroup) {
    		n1 <- sum((subgroup==i) & (dat$trt==1)); n0 <- sum((subgroup==i) & (dat$trt==0))
    		y1.bar <- mean(dat$y[subgroup==i & dat$trt==1], na.rm=T)
    		y0.bar <- mean(dat$y[subgroup==i & dat$trt==0], na.rm=T) 
		diff.mean <- y1.bar - y0.bar 		
		fit <- t.test(y~trt, data=dat, subset=(dat[ ,col.grp]==i), var.equal = var.equal)
		two.sample.ttest <- fit$statistic
		DF <- as.numeric(fit$parameter)
		p.value <- fit$p.value
    		out <- rbind(out, cbind(subgroup=i, n1, y1.bar, n0, y0.bar, diff.mean, two.sample.ttest, DF, p.value))
	}
	# COMPUTE THE P-VALUES FROM PAIRWISE t TEST OF INTERACTION
	P0 <- paired.interaction(dat, log.worth=F, col.group=col.grp, n0=2)
	# PREPARE THE OUTPUT
	out <- cbind(out, P0)
	out <- as.data.frame(out); row.names(out) <- NULL 
	for (j in c(3, 5:7, 9:(9+n.groups)))  out[,j] <- round(as.numeric(as.vector(out[,j])), digits=digits)
	out <- out[order(out$diff.mean), ]
	return(out)
}





# =========================================== END ===============================================  #
