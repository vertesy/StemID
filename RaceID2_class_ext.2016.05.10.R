############################################################
# RaceID2 Extended via Mauro 2016.05.10
############################################################

if ( !exists("Extensions")) Extensions = F

#### load required packages ####

require(amap)
require(cluster)
require(colorspace)
require(flexmix)
require(fpc)
require(lattice)
require(locfit)
require(MASS)
require(mclust)
require(pheatmap)
require(proxy)
require(RColorBrewer)
require(tsne)


#### Extensions ####

if (Extensions) {
  require(caTools)
  require(epiR)
  require(igraph)
  require(MASS)
  require(network)
  require(parallel)
  require(som)
  require(TSP)
  require(vegan)
}

#### Class definition  ####

SCseq <- setClass("SCseq", slots = c(expdata = "data.frame", ndata = "data.frame", fdata = "data.frame", distances = "matrix", tsne = "data.frame", cluster = "list", background = "list", out = "list", cpart = "vector", fcol = "vector", filterpar = "list", clusterpar = "list", outlierpar ="list" ))

setValidity("SCseq",
            function(object) {
              msg <- NULL
              if ( ! is.data.frame(object@expdata) ){
                msg <- c(msg, "input data must be data.frame")
              }else if ( nrow(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
                msg <- c(msg, "NAs are not allowed in input data")
              }else if (sum( apply( object@expdata,1,min ) ) < 0 ){
                msg <- c(msg, "negative values are not allowed in input data")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )

setMethod("initialize",
          signature = "SCseq",
          definition = function(.Object, expdata ){
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
          )

#### Filtering and Downsampling methods and functions ####

setGeneric("filterdata", function(object, mintotal=1000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000, dsversion = 'JCB') standardGeneric("filterdata"))

setMethod("filterdata",
          signature = "SCseq",
          definition = function(object,mintotal,minexpr,minnumber,maxexpr,downsample,dsn,rseed,dsversion) {
            if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
            if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
            if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
            if ( ! ( is.numeric(downsample) | is.logical(downsample) ) ) stop( "downsample has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
            if ( ! dsversion %in% c('DG', 'JCB') ) stop ("dsversion must be either DG or JCB")
            object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=downsample, dsn=dsn, dsversion = dsversion)
            object@ndata <- object@expdata[,apply(object@expdata,2,sum,na.rm=TRUE) >= mintotal]
            if ( downsample ) {
              set.seed(rseed)
              object@ndata <- downsample(object@expdata, n=mintotal, dsn=dsn, dsversion = dsversion)
            } else {
              x <- object@ndata
              object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) + .1 )
            }
            x <- object@ndata
            object@fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
            x <- object@fdata
            object@fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]
            return(object)
          }
          )

downsample <- function(x, n, dsn, dsversion){
    if (dsversion == "DG") {  # Dominic's downsampler
        x <- round(x,0)
        x <- x[apply(x,2,sum,na.rm=TRUE) >= n]
        for ( j in 1:dsn ){
            z  <- data.frame(GENEID=rownames(x), row.names = rownames(x))
            for ( i in 1:ncol(x) ){
                y <- aggregate(rep(1,n), by = list(sample(rep(rownames(x),x[,i]),n)), FUN = sum)
                na <- colnames(x[i])
                colnames(y) <- c("GENEID", na)
                rownames(y) <- y$GENEID
                z[,na] <- rep(0,nrow(z))
                k <- intersect(rownames(z),y$GENEID)
                z[k,na] <- y[k,na]
                z[is.na(z[,na]),na] <- 0
            }
            ds <- if ( j == 1 ) z[,-1] else ds + z[,-1]
        }
        ds <- ds/dsn + .1
        
    } else if (dsversion == "JCB") {  # Jean-Charles' downsampler
        rnd <- round(x, 0)
        rnd <- rnd[, which(colSums(rnd) >= n), drop = FALSE]
        ds <- rnd
        ds[,] = 0
        pool <- vector(mode = 'list', length = ncol(rnd))
        for (i in 1:ncol(rnd)){pool[[i]] <- rnd[rnd[,i]!=0,i,drop=FALSE]}
        for (j in 1: dsn) {
            poollists <- lapply(pool, function(y){table(sample(rep(row.names(y),y[,1]), size =
                n,replace = FALSE))})
            for (i in 1:ncol(rnd)){
                ds[match(names(poollists[[i]]),rownames(rnd)),i] <- as.numeric(poollists[[i]]) + ds[match(names(poollists[[i]]),rownames(rnd)),i]
            }
        }
        ds <- ds/dsn + .1
    }
    
    return(ds)
}


#### Clustering methods and functions ####

setGeneric("clustexp", function(object, clustnr=20, bootnr=50, metric="pearson", do.gap=TRUE, sat=FALSE, SE.method="Tibs2001SEmax", SE.factor=.25, B.gap=50, cln=0, rseed=17000, FUNcluster="kmeans") standardGeneric("clustexp"))

setMethod("clustexp",
signature = "SCseq",
definition = function(object,clustnr,bootnr,metric,do.gap,sat,SE.method,SE.factor,B.gap,cln,rseed,FUNcluster) {
    if ( ! is.numeric(clustnr) ) stop("clustnr has to be a positive integer") else if ( round(clustnr) != clustnr | clustnr <= 0 ) stop("clustnr has to be a positive integer")
    if ( ! is.numeric(bootnr) ) stop("bootnr has to be a positive integer") else if ( round(bootnr) != bootnr | bootnr <= 0 ) stop("bootnr has to be a positive integer")
    if ( ! ( metric %in% c( "spearman","pearson","kendall","euclidean","maximum","manhattan","canberra","binary","minkowski") ) ) stop("metric has to be one of the following: spearman, pearson, kendall, euclidean, maximum, manhattan, canberra, binary, minkowski")
    if ( ! ( SE.method %in% c( "firstSEmax","Tibs2001SEmax","globalSEmax","firstmax","globalmax") ) ) stop("SE.method has to be one of the following: firstSEmax, Tibs2001SEmax, globalSEmax, firstmax, globalmax")
    if ( ! is.numeric(SE.factor) ) stop("SE.factor has to be a non-negative integer") else if  ( SE.factor < 0 )  stop("SE.factor has to be a non-negative integer")
    if ( ! ( is.numeric(do.gap) | is.logical(do.gap) ) ) stop( "do.gap has to be logical (TRUE/FALSE)" )
    if ( ! ( is.numeric(sat) | is.logical(sat) ) ) stop( "sat has to be logical (TRUE/FALSE)" )
    if ( ! is.numeric(B.gap) ) stop("B.gap has to be a positive integer") else if ( round(B.gap) != B.gap | B.gap <= 0 ) stop("B.gap has to be a positive integer")
    if ( ! is.numeric(cln) ) stop("cln has to be a non-negative integer") else if ( round(cln) != cln | cln < 0 ) stop("cln has to be a non-negative integer")
    if ( ! is.numeric(rseed) ) stop("rseed has to be numeric")
    if ( !do.gap & !sat & cln == 0 ) stop("cln has to be a positive integer or either do.gap or sat has to be TRUE")
    if ( ! ( FUNcluster %in% c("kmeans","hclust","kmedoids") ) ) stop("FUNcluster has to be one of the following: kmeans, hclust,kmedoids")
    
    object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,sat=sat,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed,FUNcluster=FUNcluster)
    if ( clustnr < 2) stop("Choose clustnr > 1")

    y <- clustfun(object@fdata,clustnr,bootnr,metric,do.gap,sat,SE.method,SE.factor,B.gap,cln,rseed,FUNcluster)
    
    object@cluster   <- list(kpart=y$clb$result$partition, jaccard=y$clb$bootmean, gap=y$gpr, clb=y$clb)
    object@distances <- as.matrix( y$di )
    set.seed(111111)
    object@fcol <- sample(rainbow(max(y$clb$result$partition)))
    return(object)
}
)

clustfun <- function(x,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,sat=FALSE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids",distances=NULL,link="single")
{
    di <- dist.gen(t(x), method=metric)
    if ( do.gap | sat | cln > 0 ){
        gpr <- NULL
        f <- if ( cln == 0 ) TRUE else FALSE
        if ( do.gap ){
            set.seed(rseed)
            if ( FUNcluster == "kmeans" )   gpr <- clusGapExt(as.matrix(di), FUN = kmeans, K.max = clustnr, B = B.gap, iter.max=100)
            if ( FUNcluster == "kmedoids" ) gpr <- clusGapExt(as.matrix(t(x)), FUN = function(x,k) pam(dist.gen(x,method=metric),k), K.max = clustnr, B = B.gap, method=metric)
            if ( FUNcluster == "hclust" )   gpr <- clusGapExt(as.matrix(di), FUN = function(x,k){ y <- hclusterCBI(x,k,link=link,scaling=FALSE); y$cluster <- y$partition; y }, K.max = clustnr, B = B.gap)
            if ( f ) cln <- maxSE(gpr$Tab[,3],gpr$Tab[,4],method=SE.method,SE.factor)
        }
        if ( sat ){
            if ( ! do.gap ){
                if ( FUNcluster == "kmeans" )   gpr <- clusGapExt(as.matrix(di), FUN = kmeans, K.max = clustnr, B = B.gap, iter.max=100, random=FALSE)
                if ( FUNcluster == "kmedoids" ) gpr <- clusGapExt(as.matrix(t(x)), FUN = function(x,k) pam(dist.gen(x,method=metric),k), K.max = clustnr, B = B.gap, random=FALSE, method=metric)
                if ( FUNcluster == "hclust" )   gpr <- clusGapExt(as.matrix(di), FUN = function(x,k){ y <- hclusterCBI(x,k,link=link,scaling=FALSE); y$cluster <- y$partition; y }, K.max = clustnr, B = B.gap, random=FALSE)
            }
            g <- gpr$Tab[,1]
            y <- g[-length(g)] - g[-1]
            mm <- numeric(length(y))
            nn <- numeric(length(y))
            for ( i in 1:length(y)){
                mm[i] <- mean(y[i:length(y)])
                nn[i] <- sqrt(var(y[i:length(y)]))
            }
            if ( f ) cln <- max(min(which( y - (mm + nn) < 0 )),1)
        }
        if ( cln <= 1 ) {
            clb <- list(result=list(partition=rep(1,dim(x)[2])),bootmean=1)
            names(clb$result$partition) <- names(x)
            return(list(x=x, clb=clb, gpr=gpr, di=di))
        }
        if ( FUNcluster == "kmeans" )   clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=kmeansCBI,krange=cln,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
        if ( FUNcluster == "kmedoids" ) clb <- clusterboot(di,B=bootnr,bootmethod="boot",clustermethod=pamkCBI,k=cln,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
        if ( FUNcluster == "hclust" )   clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=hclusterCBI,k=cln,link=link,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
        return(list(x=x, clb=clb, gpr=gpr, di=di))
    }
}


dist.gen <- function(x, method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x), method=method,...) ) else dist(x, method=method,...)

clusGapExt <-function (x, FUNcluster, K.max, B = 100, verbose = interactive(), method="euclidean",random=TRUE,...)
{
    stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >=
    2, (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
    if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)
    stop("'B' has to be a positive integer")
    if (is.data.frame(x))
    x <- as.matrix(x)
    ii <- seq_len(n)
    W.k <- function(X, kk) {
        clus <- if (kk > 1)
        FUNcluster(X, kk, ...)$cluster
        else rep.int(1L, nrow(X))
        0.5 * sum(vapply(split(ii, clus), function(I) {
            xs <- X[I, , drop = FALSE]
            sum(dist.gen(xs,method=method)/nrow(xs))
        }, 0))
    }
    logW <- E.logW <- SE.sim <- numeric(K.max)
    if (verbose)
    cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. ",
    sep = "")
    for (k in 1:K.max) logW[k] <- log(W.k(x, k))
    if (verbose)
    cat("done\n")
    xs <- scale(x, center = TRUE, scale = FALSE)
    m.x <- rep(attr(xs, "scaled:center"), each = n)
    V.sx <- svd(xs)$v
    rng.x1 <- apply(xs %*% V.sx, 2, range)
    logWks <- matrix(0, B, K.max)
    if (random){
        if (verbose)
        cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n",
        sep = "")
        for (b in 1:B) {
            z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1],
            max = M[2]), nn = n)
            z <- tcrossprod(z1, V.sx) + m.x
            ##z <- apply(x,2,function(m) runif(length(m),min=min(m),max=max(m)))
            ##z <- apply(x,2,function(m) sample(m))
            for (k in 1:K.max) {
                logWks[b, k] <- log(W.k(z, k))
            }
            if (verbose)
            cat(".", if (b%%50 == 0)
            paste(b, "\n"))
        }
        if (verbose && (B%%50 != 0))
        cat("", B, "\n")
        E.logW <- colMeans(logWks)
        SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
    }else{
        E.logW <- rep(NA,K.max)
        SE.sim <- rep(NA,K.max)
    }
    structure(class = "clusGap", list(Tab = cbind(logW, E.logW,
    gap = E.logW - logW, SE.sim), n = n, B = B, FUNcluster = FUNcluster))
}


KmeansCBI <- function (data, krange, k = NULL, scaling = FALSE, runs = 1,
criterion = "ch", method="euclidean",...)
{
    if (!is.null(k))
    krange <- k
    if (!identical(scaling, FALSE))
    sdata <- scale(data, center = TRUE, scale = scaling)
    else sdata <- data
    c1 <- Kmeansruns(sdata, krange, runs = runs, criterion = criterion, method = method,
    ...)
    partition <- c1$cluster
    cl <- list()
    nc <- krange
    for (i in 1:nc) cl[[i]] <- partition == i
    out <- list(result = c1, nc = nc, clusterlist = cl, partition = partition,
    clustermethod = "kmeans")
    out
}

Kmeansruns <- function (data, krange = 2:10, criterion = "ch", iter.max = 100,
runs = 100, scaledata = FALSE, alpha = 0.001, critout = FALSE,
plot = FALSE, method="euclidean", ...)
{
    data <- as.matrix(data)
    if (criterion == "asw")
    sdata <- dist(data)
    if (scaledata)
    data <- scale(data)
    cluster1 <- 1 %in% krange
    crit <- numeric(max(krange))
    km <- list()
    for (k in krange) {
        if (k > 1) {
            minSS <- Inf
            kmopt <- NULL
            for (i in 1:runs) {
                options(show.error.messages = FALSE)
                repeat {
                    kmm <- try(Kmeans(data, k, iter.max = iter.max, method=method,
                    ...))
                    if (class(kmm) != "try-error")
                    break
                }
                options(show.error.messages = TRUE)
                swss <- sum(kmm$withinss)
                if (swss < minSS) {
                    kmopt <- kmm
                    minSS <- swss
                }
                if (plot) {
                    par(ask = TRUE)
                    pairs(data, col = kmm$cluster, main = swss)
                }
            }
            km[[k]] <- kmopt
            crit[k] <- switch(criterion, asw = cluster.stats(sdata,
            km[[k]]$cluster)$avg.silwidth, ch = calinhara(data,
            km[[k]]$cluster))
            if (critout)
            cat(k, " clusters ", crit[k], "\n")
        }
    }
    if (cluster1)
    cluster1 <- dudahart2(data, km[[2]]$cluster, alpha = alpha)$cluster1
    k.best <- which.max(crit)
    if (cluster1)
    k.best <- 1
    km[[k.best]]$crit <- crit
    km[[k.best]]$bestk <- k.best
    out <- km[[k.best]]
    out
}

#### Outlier detection ####

setGeneric("findoutliers", function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.95, version = 2) standardGeneric("findoutliers"))

setMethod("findoutliers",
signature = "SCseq",
definition = function(object,outminc,outlg,probthr,thr,outdistquant, version) {
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before findoutliers")
    if ( ! is.numeric(outminc) ) stop("outminc has to be a non-negative integer") else if ( round(outminc) != outminc | outminc < 0 ) stop("outminc has to be a non-negative integer")
    if ( ! is.numeric(outlg) ) stop("outlg has to be a non-negative integer") else if ( round(outlg) != outlg | outlg < 0 ) stop("outlg has to be a non-negative integer")
    if ( ! is.numeric(probthr) ) stop("probthr has to be a number between 0 and 1") else if (  probthr < 0 | probthr > 1 ) stop("probthr has to be a number between 0 and 1")
    if ( ! is.numeric(thr) ) stop("thr hast to be a vector of numbers between 0 and 1") else if ( min(thr) < 0 | max(thr) > 1 ) stop("thr hast to be a vector of numbers between 0 and 1")
    if ( ! is.numeric(outdistquant) ) stop("outdistquant has to be a number between 0 and 1") else if (  outdistquant < 0 | outdistquant > 1 ) stop("outdistquant has to be a number between 0 and 1")
    if (! version %in% c(1,2) ) stop("version has to be 1 or 2")
    
    object@outlierpar <- list( outminc=outminc,outlg=outlg,probthr=probthr,thr=thr,outdistquant=outdistquant )
    ### calibrate background model
    m <- log2(apply(object@fdata,1,mean))
    v <- log2(apply(object@fdata,1,var))
    f <- m > -Inf & v > -Inf
    m <- m[f]
    v <- v[f]
    mm <- -8
    repeat{
        fit <- lm(v ~ m + I(m^2))
        if( coef(fit)[3] >= 0 | mm >= 3){
            break
        }
        mm <- mm + .5
        f <- m > mm
        m <- m[f]
        v <- v[f]
    }
    object@background <- list()
    object@background$vfit <- fit
    object@background$lvar <- function(x,object) 2**(coef(object@background$vfit)[1] + log2(x)*coef(object@background$vfit)[2] + coef(object@background$vfit)[3] * log2(x)**2)
    object@background$lsize <- function(x,object) x**2/(max(x + 1e-6,object@background$lvar(x,object)) - x)
    
    ### identify outliers
    out   <- c()
    stest <- rep(0,length(thr))
    cprobs <- c()
    for ( n in 1:max(object@cluster$kpart) ){
        if ( sum(object@cluster$kpart == n) == 1 ){ cprobs <- append(cprobs,.5); names(cprobs)[length(cprobs)] <- names(object@cluster$kpart)[object@cluster$kpart == n]; next }
        x <- object@fdata[,object@cluster$kpart == n]
        x <- x[apply(x,1,max) > outminc,]
        z <- t( apply(x,1,function(x){ apply( cbind( pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) , 1 - pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) ),1, min) } ) )
        cp <- apply(z,2,function(x){ y <- p.adjust(x,method="BH"); y <- y[order(y,decreasing=FALSE)]; return(y[outlg]);})
        f <- cp < probthr
        cprobs <- append(cprobs,cp)
        if ( sum(f) > 0 ) out <- append(out,names(x)[f])
        for ( j in 1:length(thr) )  stest[j] <-  stest[j] + sum( cp < thr[j] )
    }
    object@out <-list(out=out,stest=stest,thr=thr,cprobs=cprobs)
    
    ### cluster outliers
    clp2p.cl <- c()
    cols <- names(object@fdata)
    cpart <- object@cluster$kpart
    di   <- as.data.frame(object@distances)
    for ( i in 1:max(cpart) ) {
        tcol <- cols[cpart == i]
        if ( sum(!(tcol %in% out)) > 1 ) clp2p.cl <- append(clp2p.cl,as.vector(t(di[tcol[!(tcol %in% out)],tcol[!(tcol %in% out)]])))
    }
    clp2p.cl <- clp2p.cl[clp2p.cl>0]
    
    cadd  <- list()
    if ( length(out) > 0 ){
        if (length(out) == 1){
            cadd <- list(out)
        }else{
            n <- out
            m <- as.data.frame(di[out,out])
            
            for ( i in 1:length(out) ){
                if ( length(n) > 1 ){
                    o   <- order(apply(cbind(m,1:dim(m)[1]),1,function(x)  min(x[1:(length(x)-1)][-x[length(x)]])),decreasing=FALSE)
                    m <- m[o,o]
                    n <- n[o]
                    f <- m[,1] < quantile(clp2p.cl,outdistquant) | m[,1] == min(clp2p.cl)
                    ind <- 1
                    if ( sum(f) > 1 ) for ( j in 2:sum(f) ) if ( apply(m[f,f][j,c(ind,j)] > quantile(clp2p.cl,outdistquant) ,1,sum) == 0 ) ind <- append(ind,j)
                    cadd[[i]] <- n[f][ind]
                    g <- ! n %in% n[f][ind]
                    n <- n[g]
                    m <- m[g,g]
                    if ( sum(g) == 0 ) break
                    
                }else if (length(n) == 1){
                    cadd[[i]] <- n
                    break
                }
            }
        }
        
        for ( i in 1:length(cadd) ){
            cpart[cols %in% cadd[[i]]] <- max(cpart) + 1
        }
    }
    
    ### determine final clusters
    if (version == 1) {
        for ( i in 1:max(cpart) ){
            d <- object@fdata[,cols[cpart == i]]
            if ( sum(cpart == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
            if ( i == 1 ) dcent <- data.frame(cent) else dcent <- cbind(dcent,cent)
            if ( i == 1 ) tmp <- data.frame(apply(object@fdata[,cols],2,dist.gen.pairs,y=cent,method=object@clusterpar$metric)) else tmp <- cbind(tmp,apply(object@fdata[,cols],2,dist.gen.pairs,y=cent,method=object@clusterpar$metric))
        }
    } else if (version == 2) {
        flag <- 1
        for ( i in 1:max(cpart) ){
            f <- cols[cpart == i]
            if ( length(f) > 0 ){
                d <- object@fdata
                if ( length(f) == 1 ){
                    cent <- d[,f]
                }else{
                    if ( object@clusterpar$FUNcluster == "kmedoids" ){
                        x <- apply(as.matrix(dist.gen(t(d[,f]),method=object@clusterpar$metric)),2,mean)
                        cent <- d[,f[which(x == min(x))[1]]]
                    }else{
                        cent <- apply(d[,f],1,mean)
                    }
                }
                if ( flag == 1 ) dcent <- data.frame(cent) else dcent <- cbind(dcent,cent)
                if ( flag == 1 ) tmp <- data.frame(apply(d,2,dist.gen.pairs,y=cent,method=object@clusterpar$metric)) else tmp <- cbind(tmp,apply(d,2,dist.gen.pairs,y=cent,method=object@clusterpar$metric))
                flag <- 0
            }
        }
    }
    
    cpart <- apply(tmp,1,function(x) order(x,decreasing=FALSE)[1])
    
    for  ( i in max(cpart):1){if (sum(cpart==i)==0) cpart[cpart>i] <- cpart[cpart>i] - 1 }
    
    object@cpart <- cpart
    
    set.seed(111111)
    object@fcol <- sample(rainbow(max(cpart)))
    return(object)
}
)

#### tSNEmap methods and functions ####

setGeneric("comptsne", function(object,rseed=15555,sammonmap=FALSE,initial_cmd=TRUE,...) standardGeneric("comptsne"))

setMethod("comptsne",
signature = "SCseq",
definition = function(object,rseed,sammonmap,...){
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before comptsne")
    set.seed(rseed)
    di <- if ( object@clusterpar$FUNcluster == "kmedoids") as.dist(object@distances) else dist.gen(as.matrix(object@distances))
    if ( sammonmap ){
        object@tsne <- as.data.frame(sammon(di,k=2)$points)
    } else {
        #  ts <- if ( initial_cmd ) tsne(di,initial_config=cmdscale(di,k=2),...) else
        ts <- tsne(di,k=2,...)
        object@tsne <- as.data.frame(ts)
    }
    return(object)
}
)

#### Differential gene expression analysis

#### Plot clustheatmap

setGeneric("clustheatmap", function(object, final=FALSE, hmethod="single") standardGeneric("clustheatmap"))

setMethod("clustheatmap",
signature = "SCseq",
definition = function(object,final,hmethod){
    if ( final & length(object@cpart) == 0 ) stop("run findoutliers before clustheatmap")
    if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before clustheatmap")
    x <- object@fdata
    part <- if ( final ) object@cpart else object@cluster$kpart
    na <- c()
    j <- 0
    for ( i in 1:max(part) ){
        if ( sum(part == i) == 0 ) next
        j <- j + 1
        na <- append(na,i)
        d <- x[,part == i]
        if ( sum(part == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
        if ( j == 1 ) tmp <- data.frame(cent) else tmp <- cbind(tmp,cent)
    }
    names(tmp) <- paste("cl",na,sep=".")
    ld <- if ( object@clusterpar$FUNcluster == "kmedoids" ) dist.gen(t(tmp),method=object@clusterpar$metric) else dist.gen(as.matrix(dist.gen(t(tmp),method=object@clusterpar$metric)))
    if ( max(part) > 1 )  cclmo <- hclust(ld,method=hmethod)$order else cclmo <- 1
    q <- part
    for ( i in 1:max(part) ){
        q[part == na[cclmo[i]]] <- i
    }
    part <- q
    di <-  if ( object@clusterpar$FUNcluster == "kmedoids" ) object@distances else as.data.frame( as.matrix( dist.gen(t(object@distances)) ) )
    pto <- part[order(part,decreasing=FALSE)]
    ptn <- c()
    for ( i in 1:max(pto) ){ pt <- names(pto)[pto == i]; z <- if ( length(pt) == 1 ) pt else pt[hclust(as.dist(t(di[pt,pt])),method=hmethod)$order]; ptn <- append(ptn,z) }
    col <- object@fcol
    mi  <- min(di,na.rm=TRUE)
    ma  <- max(di,na.rm=TRUE)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    ColorRamp   <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    if ( mi == ma ){
        ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
    }
    par(mar = c(3,5,2.5,2))
    image(as.matrix(di[ptn,ptn]),col=ColorRamp,axes=FALSE)
    abline(0,1)
    box()
    
    tmp <- c()
    for ( u in 1:max(part) ){
        ol <- (0:(length(part) - 1)/(length(part) - 1))[ptn %in% names(x)[part == u]]
        points(rep(0,length(ol)),ol,col=col[cclmo[u]],pch=15,cex=.75)
        points(ol,rep(0,length(ol)),col=col[cclmo[u]],pch=15,cex=.75)
        tmp <- append(tmp,mean(ol))
    }
    axis(1,at=tmp,lab=cclmo)
    axis(2,at=tmp,lab=cclmo)
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
    matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
    col=ColorRamp,
    xlab="",ylab="",
    xaxt="n")
    layout(1)
    return(cclmo)
}
)




#### Others ####

dist.gen.pairs <- function(x,y,...) dist.gen(t(cbind(x,y)),...)



binompval <- function(p,N,n){
  pval   <- pbinom(n,round(N,0),p,lower.tail=TRUE)
  pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
  return(pval)
}



setGeneric("plotgap", function(object) standardGeneric("plotgap"))

setMethod("plotgap",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotgap")
            plot(object@cluster$gap)
          }
          )

setGeneric("plotsilhouette", function(object) standardGeneric("plotsilhouette"))

setMethod("plotsilhouette",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotsilhouette")
            if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
            kpart <- object@cluster$kpart
            distances  <- dist.gen(object@distances)
            si <- silhouette(kpart,distances)
            plot(si)
          }
          )

setGeneric("plotjaccard", function(object) standardGeneric("plotjaccard"))

setMethod("plotjaccard",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotjaccard")
            if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no Jaccard's similarity plot")
            barplot(object@cluster$jaccard,names.arg=1:length(object@cluster$jaccard),ylab="Jaccard's similarity")
          }
          )

setGeneric("plotoutlierprobs", function(object) standardGeneric("plotoutlierprobs"))

setMethod("plotoutlierprobs",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotoutlierprobs")
            p <- object@cluster$kpart[ order(object@cluster$kpart,decreasing=FALSE)]
            x <- object@out$cprobs[names(p)]
            fcol <- object@fcol
            for ( i in 1:max(p) ){
              y <- -log10(x + 2.2e-16)
              y[p != i] <- 0
              if ( i == 1 ) b <- barplot(y,ylim=c(0,max(-log10(x + 2.2e-16))*1.1),col=fcol[i],border=fcol[i],names.arg=FALSE,ylab="-log10prob") else barplot(y,add=TRUE,col=fcol[i],border=fcol[i],names.arg=FALSE,axes=FALSE)
  }
            abline(-log10(object@outlierpar$probthr),0,col="black",lty=2)
            d <- b[2,1] - b[1,1]
            y <- 0
            for ( i in 1:max(p) ) y <- append(y,b[sum(p <=i),1] + d/2)
            axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,lab=1:max(p))
            box()
          }
          )

setGeneric("plotbackground", function(object) standardGeneric("plotbackground"))

setMethod("plotbackground",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotbackground")
            m <- apply(object@fdata,1,mean)
            v <- apply(object@fdata,1,var)
            fit <- locfit(v~lp(m,nn=.7),family="gamma",maxk=500)
            plot(log2(m),log2(v),pch=20,xlab="log2mean",ylab="log2var",col="grey")
            lines(log2(m[order(m)]),log2(object@background$lvar(m[order(m)],object)),col="red",lwd=2)
            lines(log2(m[order(m)]),log2(fitted(fit)[order(m)]),col="orange",lwd=2,lty=2)
            legend("topleft",legend=substitute(paste("y = ",a,"*x^2 + ",b,"*x + ",d,sep=""),list(a=round(coef(object@background$vfit)[3],2),b=round(coef(object@background$vfit)[2],2),d=round(coef(object@background$vfit)[3],2))),bty="n")
          }
          )

setGeneric("plotsensitivity", function(object) standardGeneric("plotsensitivity"))

setMethod("plotsensitivity",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotsensitivity")
            plot(log10(object@out$thr), object@out$stest, type="l",xlab="log10 Probability cutoff", ylab="Number of outliers")
            abline(v=log10(object@outlierpar$probthr),col="red",lty=2)
          }
          )

setGeneric("clustdiffgenes", function(object,pvalue=.01) standardGeneric("clustdiffgenes"))

setMethod("clustdiffgenes",
          signature = "SCseq",
          definition = function(object,pvalue){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before clustdiffgenes")
            if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@cpart
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              m <- apply(x,1,mean)
              n <- if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
              no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
              m <- m*no
              n <- n*no
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.all=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
              cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
            }
            return(cdiff)
          }
          )

setGeneric("diffgenes", function(object,cl1,cl2,mincount=5) standardGeneric("diffgenes"))

setMethod("diffgenes",
          signature = "SCseq",
          definition = function(object,cl1,cl2,mincount){
            part <- object@cpart
            cl1 <- c(cl1)
            cl2 <- c(cl2)
            if ( length(part) == 0 ) stop("run findoutliers before diffgenes")
            if ( ! is.numeric(mincount) ) stop("mincount has to be a non-negative number") else if (  mincount < 0 ) stop("mincount has to be a non-negative number")
            if ( length(intersect(cl1, part)) < length(unique(cl1)) ) stop( paste("cl1 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            if ( length(intersect(cl2, part)) < length(unique(cl2)) ) stop( paste("cl2 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            f <- apply(object@ndata[,part %in% c(cl1,cl2)],1,max) > mincount
            x <- object@ndata[f,part %in% cl1]
            y <- object@ndata[f,part %in% cl2]
            if ( sum(part %in% cl1) == 1 ) m1 <- x else m1 <- apply(x,1,mean)
            if ( sum(part %in% cl2) == 1 ) m2 <- y else m2 <- apply(y,1,mean)
            if ( sum(part %in% cl1) == 1 ) s1 <- sqrt(x) else s1 <- sqrt(apply(x,1,var))
            if ( sum(part %in% cl2) == 1 ) s2 <- sqrt(y) else s2 <- sqrt(apply(y,1,var))

            d <- ( m1 - m2 )/ apply( cbind( s1, s2 ),1,mean )
            names(d) <- rownames(object@ndata)[f]
            if ( sum(part %in% cl1) == 1 ){
              names(x) <- names(d)
              x <- x[order(d,decreasing=TRUE)]
            }else{
              x <- x[order(d,decreasing=TRUE),]
            }
            if ( sum(part %in% cl2) == 1 ){
              names(y) <- names(d)
              y <- y[order(d,decreasing=TRUE)]
            }else{
              y <- y[order(d,decreasing=TRUE),]
            }
            return(list(z=d[order(d,decreasing=TRUE)],cl1=x,cl2=y,cl1n=cl1,cl2n=cl2))
          }
          )

plotdiffgenes <- function(z,gene=g){
  if ( ! is.list(z) ) stop("first arguments needs to be output of function diffgenes")
  if ( length(z) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( sum(names(z) == c("z","cl1","cl2","cl1n","cl2n")) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( length(gene) > 1 ) stop("only single value allowed for argument gene")
  if ( !is.numeric(gene) & !(gene %in% names(z$z)) ) stop("argument gene needs to be within rownames of first argument or a positive integer number")
  if ( is.numeric(gene) ){ if ( gene < 0 | round(gene) != gene ) stop("argument gene needs to be within rownames of first argument or a positive integer number") }
  x <- if ( is.null(dim(z$cl1)) ) z$cl1[gene] else t(z$cl1[gene,])
  y <- if ( is.null(dim(z$cl2)) ) z$cl2[gene] else t(z$cl2[gene,])
  plot(1:length(c(x,y)),c(x,y),ylim=c(0,max(c(x,y))),xlab="",ylab="Expression",main=gene,cex=0,axes=FALSE)
  axis(2)
  box()
  u <- 1:length(x)
  rect(u - .5,0,u + .5,x,col="red")
  v <- c(min(u) - .5,max(u) + .5)
  axis(1,at=mean(v),lab=paste(z$cl1n,collapse=","))
  lines(v,rep(mean(x),length(v)))
  lines(v,rep(mean(x)-sqrt(var(x)),length(v)),lty=2)
  lines(v,rep(mean(x)+sqrt(var(x)),length(v)),lty=2)

  u <- ( length(x) + 1 ):length(c(x,y))
  v <- c(min(u) - .5,max(u) + .5)
  rect(u - .5,0,u + .5,y,col="blue")
  axis(1,at=mean(v),lab=paste(z$cl2n,collapse=","))
  lines(v,rep(mean(y),length(v)))
  lines(v,rep(mean(y)-sqrt(var(y)),length(v)),lty=2)
  lines(v,rep(mean(y)+sqrt(var(y)),length(v)),lty=2)
  abline(v=length(x) + .5)
}




setGeneric("plottsne", function(object,final=TRUE) standardGeneric("plottsne"))

setMethod("plottsne",
          signature = "SCseq",
          definition = function(object,final){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( final & length(object@cpart) == 0 ) stop("run findoutliers before plottsne")
            if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before plottsne")
            part <- if ( final ) object@cpart else object@cluster$kpart
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=object@fcol[i],cex=.75,font=4)
            }
          }
          )

setGeneric("plotlabelstsne", function(object,labels=NULL) standardGeneric("plotlabelstsne"))

setMethod("plotlabelstsne",
          signature = "SCseq",
          definition = function(object,labels){
            if ( is.null(labels ) ) labels <- names(object@ndata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotlabelstsne")
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            text(object@tsne[,1],object@tsne[,2],labels,cex=.5)
          }
          )


setGeneric("plotsymbolstsne", function(object,types=NULL) standardGeneric("plotsymbolstsne"))

setMethod("plotsymbolstsne",
          signature = "SCseq",
          definition = function(object,types){
            if ( is.null(types) ) types <- names(object@fdata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotsymbolstsne")
            if ( length(types) != ncol(object@fdata) ) stop("types argument has wrong length. Length has to equal to the column number of object@ndata")
            coloc <- rainbow(length(unique(types)))
            syms <- c()
            par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,col="grey")
            for ( i in 1:length(unique(types)) ){
              f <- types == sort(unique(types))[i]
              syms <- append( syms, ( (i-1) %% 25 ) + 1 )
              points(object@tsne[f,1],object@tsne[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
            }
            legend("topright", ,inset=c(-0.35,0),legend=sort(unique(types)), col=coloc, pch=syms)
          }
          )

setGeneric("plotsymbolstsne2", function(object,types=NULL) standardGeneric("plotsymbolstsne2"))

setMethod("plotsymbolstsne2",
          signature = "SCseq",
          definition = function(object,types){
            if ( is.null(types) ) types <- names(object@fdata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotsymbolstsne")
            if ( length(types) != ncol(object@fdata) ) stop("types argument has wrong length. Length has to equal to the column number of object@ndata")
            col <- rainbow_hcl(length(unique(types)))
            xrange <- 0.05*(max(object@tsne[,1])-min(object@tsne[,1]))
            yrange <- 0.05*(max(object@tsne[,2])-min(object@tsne[,2]))
            par(mar=c(5.5, 5.5, 5.5, 5.5), xpd=TRUE)

            plot(sc@tsne,cex=1.5,col='lightgrey',pch=20,ylim=c(min(object@tsne[,2])-yrange,max(object@tsne[,2])+yrange),xlim=c(min(object@tsne[,1])-xrange,max(object@tsne[,1])+xrange))
            syms <- c()
            #            plot(object@tsne,axes=0,pch=20,col="grey")
            for ( i in 1:length(unique(types)) ){
              f <- types == sort(unique(types))[i]
              #              syms <- append( syms, ( (i-1) %% 25 ) + 1 )
              points(object@tsne[f,1],object@tsne[f,2],col=col[i],pch=20 ,cex=0.8)
            }
            legend("topright",inset=c(-0.20,0),legend=sort(unique(types)), col=col, pch=20,cex=0.75,bty='n')
          }
)


#lennarts version of plottsne2
setGeneric("plottsne2", function(object,final=TRUE) standardGeneric("plottsne2"))


setMethod("plottsne2",
          signature = "SCseq",
          definition = function(object,final){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( final & length(object@cpart) == 0 ) stop("run findoutliers before plottsne")
            if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before plottsne")
            part <- if ( final ) object@cpart else object@cluster$kpart
            col <- sample(rainbow_hcl(max(part)))
            col<-sc@fcol

            #plot(object@tsne,axes=0,pch=20,cex=1.5,col="lightgrey")
            xrange <- 0.05*(max(object@tsne[,1])-min(object@tsne[,1]))
            yrange <- 0.05*(max(object@tsne[,2])-min(object@tsne[,2]))
            plot(sc@tsne,cex=2.5,col='lightgrey',pch=20,ylim=c(min(object@tsne[,2])-yrange,max(object@tsne[,2])+yrange),xlim=c(min(object@tsne[,1])-xrange,max(object@tsne[,1])+xrange))
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) points(sc@tsne[part == i,1],sc@tsne[part == i,2],col=col[i],cex=2,pch=20)
            }
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) text(sc@tsne[part == i,1],sc@tsne[part == i,2],i,font=4,cex=0.5)
            }
          }
)


###lennart plotexptsne
setGeneric("plotexptsne", function(object,g,n="",logsc=FALSE,ax=TRUE) standardGeneric("plotexptsne"))

setMethod("plotexptsne",
          signature = "SCseq",
          definition = function(object,g,n="",logsc=FALSE,ax=TRUE){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotexptsne")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if ( !is.numeric(logsc) & !is.logical(logsc) ) stop("argument logsc has to be logical (TRUE/FALSE)")
            if ( n == "" ) n <- g[1]
            l <- apply(object@ndata[g,] - .1,2,sum) + .1
            if (logsc) {
              f <- l == 0
              l <- log2(l)
              l[f] <- NA
            }
            mi <- min(l,na.rm=TRUE)
            ma <- max(l,na.rm=TRUE)
            ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
            ColorRamp <- rev(terrain.colors(100))
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            v <- round((l - mi)/(ma - mi)*99 + 1,0)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            par(mar = c(3,5,2.5,2))
            xrange <- 0.05*(max(object@tsne[,1])-min(object@tsne[,1]))
            yrange <- 0.05*(max(object@tsne[,2])-min(object@tsne[,2]))
            if(ax){
              plot(object@tsne,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=2,col="grey",ylim=c(min(object@tsne[,2])-yrange,max(object@tsne[,2])+yrange),xlim=c(min(object@tsne[,1])-xrange,max(object@tsne[,1])+xrange))
            }else{
              plot(object@tsne,main=n,pch=20,cex=2,col="grey",xaxt='n',yaxt='n',bty='n',xlab="",ylab="")
            }

            for ( k in 1:length(v) ){
              points(object@tsne[k,1],object@tsne[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
              #points(object@tsne[k,1],object@tsne[k,2],col='lightgrey',bg=ColorRamp[v[k]],pch=21,cex=1.5)
            }
            #points(object@tsne[names(preHSC)[which(preHSC == 1)],],cex=1.5,pch=1,col='red',lwd=0.75)
            #points(object@tsne[names(preHSC)[which(preHSC == 2)],],cex=1.5,pch=1,col='blue',lwd=0.75)
            #legend("bottomright",legend=c("preHSC I","preHSC II"),pch=1,col=c('red','blue'),bty='n')
            par(mar = c(3,2.5,2.5,2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",
                  xaxt="n")
            layout(1)
          }
)

###dominics plotexptsne
# setGeneric("plotexptsne", function(object,g,n="",logsc=FALSE) standardGeneric("plotexptsne"))
#
# setMethod("plotexptsne",
#           signature = "SCseq",
#           definition = function(object,g,n="",logsc=FALSE){
#             if ( length(object@tsne) == 0 ) stop("run comptsne before plotexptsne")
#             if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
#             if ( !is.numeric(logsc) & !is.logical(logsc) ) stop("argument logsc has to be logical (TRUE/FALSE)")
#             if ( n == "" ) n <- g[1]
#             l <- apply(object@ndata[g,] - .1,2,sum) + .1
#             if (logsc) {
#               f <- l == 0
#               l <- log2(l)
#               l[f] <- NA
#             }
#             mi <- min(l,na.rm=TRUE)
#             ma <- max(l,na.rm=TRUE)
#             ColorRamp <-colorRampPalette(brewer.pal(n = 7,name = "YlOrRd"))(100)
#             #ColorRamp <-rev(terrain.colors(100))
#             ColorLevels <- seq(mi, ma, length=length(ColorRamp))
#             v <- round((l - mi)/(ma - mi)*99 + 1,0)
#             layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
#             par(mar = c(3,5,2.5,2))
#             plot(object@tsne,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=0,col="grey")
#             for ( k in 1:length(v) ){
#               points(object@tsne[k,1],object@tsne[k,2],bg=ColorRamp[v[k]],pch=21,cex=1.5,col="black",lwd=0.5)
#             }
#             par(mar = c(3,2.5,2.5,2))
#             image(1, ColorLevels,
#                   matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
#                   col=ColorRamp,
#                   xlab="",ylab="",
#                   xaxt="n")
#             layout(1)
#           }
#           )


plot.2.file <- function(name, type="pdf", res.f=1){
  full.name <- paste(name,type,sep=".")
  if (type == "pdf"){
    pdf(full.name)
  }else if(type == "png"){
    png(full.name,width = 480*res.f, height = 480*res.f, pointsize = 12*res.f)
  }else if (type == "tiff"){
    tiff(full.name,width = 480*res.f, height = 480*res.f, pointsize = 12*res.f)
  }else if (type == "jpeg"){
    jpeg(full.name,width = 480*res.f, height = 480*res.f, pointsize = 12*res.f)
  }
}

id2name <- function(x) sub("\\_\\_chr\\w+","",x)

name2id <- function(x,id) id[sub("\\_\\_chr\\w+","",id) %in% x]

plot.err.bars.y <- function(x, y, y.err, col="black", lwd=1, lty=1, h=0.1){
  arrows(x,y-y.err,x,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y-y.err,x+h,y-y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y+y.err,x+h,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
}

plot.set <- function(x,absexp=TRUE,pname=""){
  if ( absexp ) x <- x/apply(x,1,sum)
  x.a <- apply(x,2,mean,na.rm=TRUE)
  x.s <- sqrt(apply(x,2,var,na.rm=TRUE))
  y.lab <- "Expression"
  y.lim <- if ( absexp ) c(0, max(x.a + x.s) *1.1 ) else c(min(x.a - x.s),max(x.a + x.s) + .1*( max(x.a + x.s) - min(x.a - x.s)))
  mi <- if ( absexp ) 0 else min(x.a)
  plot(1:length(x.a),x.a,cex=0,axes=FALSE,ylim=y.lim,ylab=y.lab,xlim=c(0,length(cols) + 1),xlab="",main=pname)
  box()
  axis(2)
  rect(1:length(cols) - .5, mi, 1:length(cols) + .5,x.a[1:length(cols)],col="red")
  plot.err.bars.y(1:length(x.a),x.a,x.s)
  legend("topleft",legend=substitute(paste("N=",a,sep=""),list(a=dim(x)[1])),bty="n")
}


heat.set <- function(x,absexp=TRUE,pname=""){
  if ( absexp ) x <- x/apply(x,1,sum)
  d.l <- hclust(dist(x))
  g.l <- hclust(dist(t(x)))
  ramp  <- colorRamp(c("darkblue", "yellow"))
  c.col <- rgb( ramp(seq(0, 1, length = 256)), max = 255)
  heatmap(as.matrix(x),Colv=as.dendrogram(g.l),Rowv=as.dendrogram(d.l),labCol=names(x),col=c.col,main=pname,keep.dendro=TRUE,labRow=NA)
}

downsample.binom <- function(x,n) as.data.frame( t(apply(t(t(x)/apply(d,2,sum)),1,function(x,n) tapply(x,1:length(x),function(x,n) rbinom(1,size=n,prob=x), n=n ), n=n )) ) + .1

mindist <- function(g,d,cl1,cl2,stringent=FALSE){
  if ( stringent ){
    d1 <- as.data.frame(d[g,cl1])
    d2 <- as.data.frame(d[g,cl2])
    up <- apply(d1,1,min) - apply(d2,1,max)
    up[up<0] <- 0
    do <- apply(d2,1,min) - apply(d1,1,max)
    do[do<0] <- 0
    sum(apply(data.frame(up,do),1,max))
  }else{
    k <- as.data.frame(as.matrix(dist(t(d[g,]))))
    min(apply(data.frame( k[cl2[! ( cl2 %in% cl1 ) ],cl1] ),1,min))
  }
}

setGeneric("diffgenescomb", function(object,cl1,cl2=NULL,g=NULL,n=1,mincount=5,stringent=FALSE,clmedian=TRUE) standardGeneric("diffgenescomb"))

setMethod("diffgenescomb",
          signature = "SCseq",
          definition = function(object,cl1,cl2,g,n,mincount,stringent,clmedian){
            if ( is.null(g) ) g <- rownames(object@ndata)
            part <- object@cpart
            if ( length(part) == 0 ) stop("run findoutliers before diffgenes")
            cl1 <- c(cl1)
            if ( is.null(cl2) ) cl2 <- sort(unique(part[! ( part %in% cl1 )]))
            cl2 <- c(cl2)
            if ( ! is.numeric(mincount) ) stop("mincount has to be a non-negative number") else if (  mincount < 0 ) stop("mincount has to be a non-negative number")
            if ( ! is.numeric(n) ) stop("n has to be a positive integer") else if (  n <= 0 | round(n) != n ) stop("n has to be a positive integer")
            if ( length(intersect(cl1, part)) < length(unique(cl1)) ) stop( paste("cl1 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            if ( length(intersect(cl2, part)) < length(unique(cl2)) ) stop( paste("cl2 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            if ( length(intersect(g, rownames(object@ndata))) < length(unique(g)) ) stop( "g has to be a subset of rownames of transcript count table" )
            part <- object@cpart[object@cpart %in% c(cl1,cl2)]
            x <- object@ndata[apply(object@ndata[,part %in% c(cl1,cl2)],1,max) > mincount,names(part)]
            g <- intersect(g,rownames(x))
            if ( length(g) < n ) return(NULL)
            f <- names(part)[part %in% cl2]
            m <- if (length(f) > 1 ) apply(x[,f],1,mean) else x[,f]
            y <- ( x - m )/sqrt(object@background$lvar(m,object))
            z <- if (clmedian) aggregate(t(y),by=list(part=part),FUN=median) else aggregate(t(y),by=list(part=part),FUN=mean)
            d <- as.data.frame( t( z[,-1] ) )
            names(d) <- paste("cl",as.vector(z[,1]),sep=".")

            zx <- if (clmedian) aggregate(t(x),by=list(part=part),FUN=median) else aggregate(t(x),by=list(part=part),FUN=mean)
            dx <- as.data.frame( t( zx[,-1] ) )
            names(dx) <- paste("cl",as.vector(zx[,1]),sep=".")
            g <- g[g %in% rownames(d)]
            v <- (1:nrow(d))[rownames(d) %in% g]
            #v <- t(tapply(g,1:length(g),grep,x=rownames(d)))
            sets <- as.data.frame(combs(v,n))
            w <- apply(as.data.frame(sets),1,mindist,d=d,cl1=paste("cl",cl1,sep="."),cl2=paste("cl",cl2,sep="."),stringent=stringent)
            f <- order(w,decreasing=TRUE)
            wn <- t(apply(as.data.frame(as.data.frame(sets)[f,]),1,function(x,d) rownames(d)[x],d=d))
            wn <- if ( n == 1 ) data.frame(V1=t(wn),score=w[f]) else data.frame(wn,score=w[f])
            return(list(marker=wn,expression=dx,zscore=d))
          }
          )


sortbyentropy <- function(object,clset,entr){
  for ( i in 1:length(clset) ){
    k <- names(entr[ object@cpart == clset[i] ])[ order(entr[ object@cpart == clset[i] ],decreasing=TRUE) ]
    n <- if ( i == 1 ) k else append(n,k)
  }
  n
}

sortbysom <- function(object,clset,di=NULL){
  x <- object@ndata[,object@cpart %in% clset]
  if ( is.null(di) ) d <- object@distances[object@cpart %in% clset,object@cpart %in% clset] else d <- di
  sm <- som(d,as.numeric(ncol(d)),1)
  names(x)[order(sm$visual$x,decreasing=FALSE)]
}

plotdynamics <- function(object,x,clset){
  m <- min(as.vector(as.matrix(x)),na.rm=TRUE)
  for ( i in 1:ncol(x) ) x[is.na(x[,i]),i] <- m
  for ( i in 1:ncol(x) ){ x[,i] <- x[,i]/x[i,i]}
  d <- rep(1,sum(part == clset[1]))
  co <- rep(object@fcol[clset[1]],sum(part == clset[1]))
  for ( i in 1:(length(clset)-1) ){
    d <- append(d,rep(x[clset[i+1],clset[i]] - 1 + d[length(d)],sum(part == clset[i+1])))
    co <- append(co,rep(object@fcol[clset[i + 1]],sum(part == clset[i+1])))
  }
  plot(1:length(d),d-1,pch=20,xlab="Pseudotime",ylab="Relative Transcriptome Change",cex=0)
  lines(1:length(d),d-1,pch=20,col="grey")
  for ( i in 1:length(d) ) points(i,d[i]-1,col=co[i],pch=20)
}

plotentropy <- function(object,entr,n,cluster=TRUE){
  cl <- unique(object@cpart[n])
  plot(1:length(n),entr[n],type="l",col="grey",axes=FALSE,xlab="",ylab="Entropy")
  y <- 0
  for ( i in 1:length(cl) ){
    f <- object@cpart[n] == cl[i]
    points((1:length(n))[f],entr[n[f]],col=object@fcol[cl[i]],pch=20)
    z <- if ( i == 1 ) sum(f)   else append(z, z[i-1] + sum(f))
    x <- if ( i == 1 ) sum(f)/2 else append(x, z[i-1] + sum(f)/2)
    if ( cluster ) abline(v=z[i],col="grey",lty=2)
  }
  axis(2)
  box()
  if ( cluster ) axis(1,at=x,lab=cl)
}

plotcumpdist <- function(object,n,di=NULL,cluster=TRUE){
  cl <- unique(object@cpart[n])
  if ( is.null(di) ){
    di <- as.data.frame(object@distances)
    names(di) <- names(object@ndata)
    rownames(di) <- names(object@ndata)
  }
  pd <- 0
  for ( i in 2:length(n) ){
    pd <- append(pd,pdist(n[1:i],di))
  }
  plot(1:length(n),pd,type="l",col="grey",axes=FALSE,xlab="",ylab="Cumulative Distance")
  y <- 0
  for ( i in 1:length(cl) ){
    f <- object@cpart[n] == cl[i]
    points((1:length(n))[f],pd[f],col=object@fcol[cl[i]],pch=20)
    z <- if ( i == 1 ) sum(f)   else append(z, z[i-1] + sum(f))
    x <- if ( i == 1 ) sum(f)/2 else append(x, z[i-1] + sum(f)/2)
    if ( cluster ) abline(v=z[i],col="grey",lty=2)
  }
  axis(2)
  box()
  if ( cluster ) axis(1,at=x,lab=cl)
}

filterset <- function(object,n,minexpr=5,minnumber=1){
  object@ndata[apply(object@ndata[,n] >= minexpr,1,sum) >= minnumber,n]
}

getsom <- function(x,k=5,nb=50,d=2,logscale=TRUE,locreg=FALSE,alpha=.25){
  if (! d %in% c(1,2) ) stop("The dimension d of the SOM has to be 1 or 2\n")
  if ( locreg ){
    x <- t(apply(x,1,function(x,alpha){ v <- 1:length(x); predict(loess( x ~ v, span=alpha ))},alpha=alpha))
  }else{
    x <- t(apply(x,1,runmean,k=k))
  }
  x <- x/apply(x,1,sum)
  if ( logscale ) x <- log2(x)
  if ( d == 2 ) return( som(x,round(sqrt(nb),0),round(sqrt(nb),0)) )
  if ( d == 1 ) return( som(x,1,nb) )
}

plotexpression <- function(object,g,n,k=5,name=NULL,cluster=TRUE,locreg=FALSE,alpha=.25,types=NULL){
  cl <- unique(object@cpart[n])
  m <- object@ndata
  xlim <- c(1,length(n))
  if ( !is.null(types) ) xlim[1] <- 1.25 * xlim[1]
  y <- if ( length(g) == 1 ) m[g,n] else t(apply(m[g,n],2,sum))
  if ( is.null(name) ) name <- g[1]
  plot(1:length(n),t(y),cex=0,axes=FALSE,xlab="",ylab="Expression",main=name,xlim=xlim)
  ##y <- 0
  if ( ! is.null(types) ){
    coloc <- rainbow(length(unique(types)))
    syms <- c()
    for ( i in 1:length(unique(types)) ){
      f <- types == sort(unique(types))[i]
      syms <- append( syms, ( (i-1) %% 25 ) + 1 )
      points((1:length(n))[f],t(y)[f],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
    }
  }
  for ( i in 1:length(cl) ){
    f <- object@cpart[n] == cl[i]
    if ( is.null(types) ){
      points((1:length(n))[f],t(y)[f],col=object@fcol[cl[i]],pch=20)
    }
    z <- if ( i == 1 ) sum(f)   else append(z, z[i-1] + sum(f))
    x <- if ( i == 1 ) sum(f)/2 else append(x, z[i-1] + sum(f)/2)
    if ( cluster ) abline(v=z[i],col="grey",lty=2)
    u <- 1:length(n)
    if ( locreg ){
      v <- t(y)
      lines(u,predict(loess( v ~ u, span=alpha )))
    }else{
      lines(u,runmean(t(y),k=k))
    }
  }
  if ( !is.null(types) ) legend("topleft", legend=sort(unique(types)), col=coloc, pch=syms)

  axis(2)
  box()
  if ( cluster ) axis(1,at=x,lab=cl)
}

hclusterCBI <- function (data, k, cut = "number", link="ward", method="euclidean", scaling = TRUE, noisecut = 0,
    ...)
{
    if (!identical(scaling, FALSE))
        sdata <- scale(data, center = TRUE, scale = scaling)
    else sdata <- data
    n <- nrow(data)
    noise <- FALSE
    c1 <- hcluster(sdata, link=link, method = method)
    if (cut == "number")
        partition <- cutree(c1, k = k)
    else partition <- cutree(c1, h = k)
    cl <- list()
    nc <- max(partition)
    clsizes <- numeric(0)
    for (i in 1:nc) clsizes[i] <- sum(partition == i)
    ncn <- sum(clsizes > noisecut)
    if (ncn < nc) {
        noise <- TRUE
        newcln <- (1:nc)[clsizes > noisecut]
        nc <- ncn + 1
        newpart <- rep(nc, n)
        for (i in 1:ncn) newpart[partition == newcln[i]] <- i
        partition <- newpart
    }
    for (i in 1:nc) cl[[i]] <- partition == i
    out <- list(result = c1, noise = noise, nc = nc, clusterlist = cl,
        partition = partition, clustermethod = "hclust/cutree")
    out
}






plotheatmap <- function(x,xpart=NULL,xcol=NULL,xlab=TRUE,ypart=NULL,ycol=NULL,ylab=TRUE,xgrid=FALSE,ygrid=FALSE){

  mi  <- min(x,na.rm=TRUE)
  ma  <- max(x,na.rm=TRUE)
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(5,1), heights=c(5,1))
  ColorRamp   <- rev(colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100))
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  if ( mi == ma ){
    ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
  }
  par(mar = c(3,5,2.5,2))
  image(t(as.matrix(x)),col=ColorRamp,axes=FALSE)
  box()
  set.seed(20)
  if ( !is.null(xpart) ){
    tmp <- c()
    for ( u in unique(xpart) ){
      ol <- (0:(length(xpart) - 1)/(length(xpart) - 1))[xpart == u]
     if ( !is.null(xcol) ) points(ol,rep(0,length(ol)),col=xcol[u],pch=15,cex=.75)
      tmp <- append(tmp,mean(ol))
      delta <- .5/(length(xpart) - 1)
      if ( xgrid & max(ol) < 1) abline(v=max(ol) + delta,col="grey",lty=2)
    }
    if ( xlab ) axis(1,at=tmp,lab=unique(xpart))
  }
  set.seed(20)
  if ( !is.null(ypart) ){
    tmp <- c()
    for ( u in unique(ypart) ){
      ol <- (0:(length(ypart) - 1)/(length(ypart) - 1))[ypart == u]
      if ( !is.null(ycol) ) points(rep(0,length(ol)),ol,col=ycol[u + 1],pch=15,cex=.75)
      tmp <- append(tmp,mean(ol))
      delta <- .5/(length(ypart) - 1)
      if ( ygrid & max(ol) < 1) abline(a=max(ol) + delta,b=0,col="grey",lty=2)
    }
    if ( ylab ) axis(2,at=tmp,lab=unique(ypart))
  }
  par(mar = c(20,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  layout(1)
}

pdist <- function(n,di){ x <- 0; for ( i in 1:( length(n) - 1 ) ){ x <- x + di[n[i],n[i + 1]] }; x}

solve_tsp_cells <- function(object,q,di=NULL,seed=NULL){
  if (!is.null(seed) ) set.seed(seed)
  if ( is.null(di) ){
    di <- as.data.frame(object@distances)
    names(di) <- names(object@ndata)
    rownames(di) <- names(object@ndata)
  }
  p <- object@cpart
  n <- names(object@cpart[object@cpart %in% q])
  CP <- c()
  for ( i in 1:length(q) ){
    CP <- append(CP,sample(n[p[n] == q[i]],1))
  }
  MP <- pdist(CP,di)
  n  <- n[! (n %in% CP ) ]
  while ( length(n) > 0 ){
    st <- sample(n,1)

    mpi <- Inf
    cpi <- Inf
    if ( p[st] == p[CP[1]] ){
      cpi <- append(st,CP)
      mpi <- pdist(cpi,di)
    }
    for ( i in 1:(length(CP)-1) ){
      if ( p[st] %in% p[CP[i:(i+1)]] ){
        cp <-c(CP[1:i],st,CP[(i+1):length(CP)])
        mp <- pdist(cp,di)
        if (mp < mpi ){
          cpi <- cp
          mpi <- mp
        }
      }
    }
    if ( p[st] == p[CP[length(CP)]] ){
      cp <-c(CP,st)
      mp <- pdist(cp,di)

      if (mp < mpi ){
        cpi <- cp
        mpi <- mp
      }
    }
    CP <- cpi
    MP <- mpi
    n  <- n[! (n %in% CP ) ]
  }
  return(cp)
}



plot.err.shade.y <- function(x,y,y.err.up,y.err.down,col="black",...){
  k <- as.vector(col2rgb(col))
  m <- rgb(k[1],k[2],k[3],maxColorValue=255,alpha=50)
  polygon(c(x,rev(x)),c(y - y.err.down, rev(y + y.err.up)), col=m, border=NA)
  lines(x,y,col=col,...)
  box()
}

plotprofile <- function(y,object,n,wl,z=NULL,partition=FALSE,ylab="",xlab="",...){
  ramp  <- colorRamp(c("darkblue", "orange"))
  co <- rgb( ramp(seq(0, 1, length = nrow(y))), max = 255)
  ylim <- if ( is.null(z) ) c(min(y,na.rm=TRUE),max(y,na.rm=TRUE)) else c(min(y - z,na.rm=TRUE),max(y + z,na.rm=TRUE))
  plot(1:ncol(y),rep(0,ncol(y)),cex=0,ylim=ylim,xlim=c(1 - wl/2,ncol(y) + wl/2),axes=FALSE,ylab=ylab,xlab=xlab)
  for ( i in 2:ncol(y)){
    lines(1:ncol(y),runmean(t(y[i,]),k=5),col=co[i])
    if ( !is.null(z) ) plot.err.shade.y(1:ncol(y),runmean(t(y[i,]),k=5),runmean(t(z[i,]),k=5),runmean(t(z[i,]),k=5),col=co[i])
  }
  if ( partition ){
    cl <- unique(object@cpart[n])
    for ( i in 1:length(cl) ){
      f <- object@cpart[n] == cl[i]
      m <- if ( i == 1 ) sum(f)   else append(m, m[i-1] + sum(f))
      x <- if ( i == 1 ) sum(f)/2 else append(x, m[i-1] + sum(f)/2)
      abline(v=m[i] - wl/2,col="grey",lty=2)
    }
    axis(1,at=x - wl/2,lab=cl)
  }
  axis(2)
  box()
}

plot3dtsne <- function(object){
  require(rgl)
  di <- dist.gen(as.matrix(object@distances))
  tt <- tsne(di,k=3)
  plot3d(tt[,1], tt[,2], tt[,3], xlab = "Dim 1", ylab = "Dim 2", zlab = "Dim 3", alpha = 0.75, col = "grey", pch="16", type="p", size = 8, point_antialias = TRUE)
  for ( i in sort(unique(object@cpart)) ){ f <- object@cpart == i; text3d(tt[f,1], tt[f,2], tt[f,3], rep(i,sum(f)), font=10, size=9, depth_test = "always", color=object@fcol[i])}
}

compnoise <- function(d,n,ebin,wl){
  rb <- data.frame(bin=0:length(ebin))
  sd <- data.frame(bin=0:length(ebin))
  av <- data.frame(bin=0:length(ebin))
  nb <- data.frame(bin=0:length(ebin))
  for ( i in 1:( length(n) - wl ) ){
    cat(i,"\n")
    nl <- n[i:(i+wl)]
    x <- apply(d[,nl],1,function(x){ w <- which(log2(mean(x)) >= ebin); if ( length(w) == 0 ) return(0) else rev(w)[1] })
    m <- apply(d[,nl],1,mean)
    v <- sqrt(apply(d[,nl],1,var))
    f <- v > 0 & m > 0
    z <- aggregate(log2(v/m)[f],by=list(bin=x[f]),mean,na.rm=TRUE)
    y <- aggregate(log2(v/m)[f],by=list(bin=x[f]),var,na.rm=TRUE)
    e <- aggregate(log2(m)[f],by=list(bin=x[f]),mean,na.rm=TRUE)
    k <- aggregate(rep(1,sum(f)),by=list(bin=x[f]),sum,na.rm=TRUE)
    names(z)[2] <- paste("w",i,sep="")
    names(e)[2] <- paste("w",i,sep="")
    names(y)[2] <- paste("w",i,sep="")
    names(k)[2] <- paste("cl",i,sep="")
    y[,2] <- sqrt(y[,2])
    rb <- merge(rb,z,by="bin",all.x=TRUE)
    av <- merge(av,e,by="bin",all.x=TRUE)
    sd <- merge(sd,y,by="bin",all.x=TRUE)
    nb <- merge(nb,k,by="bin",all.x=TRUE)
  }
  return(list(rb=rb,av=av,sd=sd,nb=nb))
}

compnoisecluster <- function(d,n,set,object,ebin,wl){
  rbc <- data.frame(bin=0:length(ebin))
  avc <- data.frame(bin=0:length(ebin))
  sdc <- data.frame(bin=0:length(ebin))
  nbc <- data.frame(bin=0:length(ebin))
  for ( i in set){
    nl <- names(object@cpart)[object@cpart == i]
    if ( length(nl) > 1 ){
      x <- apply(as.data.frame(d[,nl]),1,function(x){ w <- which(log2(mean(x)) >= ebin); if ( length(w) == 0 ) return(0) else rev(w)[1] })
      m <- apply(as.data.frame(d[,nl]),1,mean)
      v <- sqrt(apply(as.data.frame(d[,nl]),1,var))
      f <- v > 0 & m > 0
      z <- aggregate(log2(v/m)[f],by=list(bin=x[f]),mean,na.rm=TRUE)
      y <- aggregate(log2(v/m)[f],by=list(bin=x[f]),var,na.rm=TRUE)
      e <- aggregate(log2(m)[f],by=list(bin=x[f]),mean,na.rm=TRUE)
      k <- aggregate(rep(1,sum(f)),by=list(bin=x[f]),sum,na.rm=TRUE)
      names(z)[2] <- paste("cl",i,sep=".")
      names(y)[2] <- paste("cl",i,sep="")
      names(e)[2] <- paste("cl",i,sep="")
      names(k)[2] <- paste("cl",i,sep="")
      rbc <- merge(rbc,z,by="bin",all.x=TRUE)
      avc <- merge(avc,e,by="bin",all.x=TRUE)
      sdc <- merge(sdc,y,by="bin",all.x=TRUE)
      nbc <- merge(nbc,k,by="bin",all.x=TRUE)
    }else{
      rbc[,paste("cl",i,sep=".")] <- rep(NA,nrow(rbc))
      avc[,paste("cl",i,sep=".")] <- rep(NA,nrow(avc))
      sdc[,paste("cl",i,sep=".")] <- rep(NA,nrow(sdc))
      nbc[,paste("cl",i,sep=".")] <- rep(NA,nrow(nbc))
    }
  }
  return(list(rbc=rbc,avc=avc,sdc=sdc,nbc=nbc))
}


sigcor <- function(x,y,cthr=.4){
  fit <- lm(x ~ y)
  pv <- as.data.frame(summary(fit)[4])[2,4]
  y <- as.data.frame(summary(fit)[4])[2,1]
  if ( is.na(pv) | is.na(y) ) return( NA )
  z <- sign(y)*sqrt(summary(fit)$r.square)
  if ( is.na(z) ) return(NA)
  if ( pv < .01 & abs(z) >= cthr ) return(z) else return(NA)
}

thrcor <- function(x,y,cthr=.4){
  z <- cor(x,y)
  if ( abs(z) >= cthr ) return(z) else return(NA)
}


plotgenegroup <- function(gr,n,object,k,partition=FALSE){
  y <- if ( length(gr) == 1) t(runmean(t(object@ndata[gr,n]),k=k)) else t(apply(object@ndata[gr,n],1,runmean,k=k))
  h <- length(gr)
  set.seed(10)
  co <- sample(rainbow(h),h)
  plot(1:ncol(y),rep(0,ncol(y)),cex=0,ylim=c(0,max(apply(y,2,max))),xlim=c(1,ncol(y)*1.5),axes=FALSE,ylab="Expression",xlab="")
  for ( i in 1:h ) lines(1:ncol(y),y[i,],col=co[i])
  if ( partition ){
    cl <- unique(object@cpart[n])
    for ( i in 1:length(cl) ){
      f <- object@cpart[n] == cl[i]
      m <- if ( i == 1 ) sum(f)   else append(m, m[i-1] + sum(f))
      x <- if ( i == 1 ) sum(f)/2 else append(x, m[i-1] + sum(f)/2)
      abline(v=m[i] - wl/2,col="grey",lty=2)
    }
    axis(1,at=x - wl/2,lab=cl)
  }
  axis(2)
  box()
  legend("topright",legend=gr,col=co,lty=1,cex=.75,bty="n")
}

getneighbours <- function(sc,g,d,k=1,cthr,method=sigcor){
  d <- d[apply(abs(d)>cthr,1,sum,na.rm=TRUE) > 0,]
  if ( ! g %in% rownames(d) | nrow(d) <= 1 ) return(NULL)
  for ( i in 1:k ){
    g <- if ( length(g) == 1 ) unique(append(g,rownames(d)[!is.na(d[,g])])) else unique(append(g,rownames(d)[apply(!is.na(d[,g]),1,sum) > 0]))
  }
  return(g)
}

plotnetwork <- function(x,minx=-1,maxx=1){
  mcol <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(201)
  lev  <- seq(minx, maxx, length=length(mcol))
  dcc <- t(apply(round(100*(x + 1) + 1,0),1,function(x){y <- c(); for ( n in x ) y <- append(y,mcol[n]); y }))
  layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
  par(mar = c(3,1,1,1))
  plot(network(as.matrix(x)),edge.col=dcc,edge.lwd=2,vertex.cex=1,vertex.col="grey",displaylabels=TRUE,label=id2name(names(x)),label.cex=.75,usearrows=FALSE)
  par(mar = c(20,2.5,2.5,2))
  image(1, lev, matrix(data=lev, ncol=length(lev),nrow=1), col=mcol, xlab="", ylab="", xaxt="n")
  layout(1)
}

plotcvbins <- function(d,ebin){
  av <- apply(d,1,mean)
  sd <- sqrt(apply(d,1,var))
  plot(log2(av),log2(sd/av),pch=20,col="grey",xlab="log2 mean expression",ylab="log2 CV")
  ramp  <- colorRamp(c("darkblue", "orange"))
  co <- rgb( ramp(seq(0, 1, length = length(ebin))), max = 255)
  for ( i in 1:length(ebin) ) abline(v=ebin[i],col=co[i],lty=2)
}

hyper.pval <- function(N.0,n.0,N,n){
  p.val <- phyper(n,n.0,N.0-n.0,N)
  if (p.val > 0.5) p.val <- 1 - p.val
  return(p.val)
}

plotexpmap <- function(object,g,map,mapnames,n="",logsc=FALSE){
  if ( n == "" ) n <- g[1]
  l <- apply(object@ndata[g,mapnames],2,sum)
  if (logsc) {
    f <- l == 0
    l <- log2(l)
    l[f] <- NA
  }
  mi <- min(l,na.rm=TRUE)
  ma <- max(l,na.rm=TRUE)
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  v <- round((l - mi)/(ma - mi)*99 + 1,0)
  layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
  par(mar = c(3,5,2.5,2))
  plot(map,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=0,col="grey")
  for ( k in 1:length(v) ){
    points(map[k,1],map[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
  }
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  layout(1)
}


diffexpnb <- function(x, cells_of_interest, cells_background ,norm=TRUE,DESeq=FALSE,method="per-condition",vfit=NULL,locreg=FALSE){
  if ( ! method %in% c("per-condition","pooled","pooled-CR") ) stop("invalid method")
  x <- x[,c(cells_background,cells_of_interest)]
  if ( DESeq ){
    des <- data.frame( "row.names" = colnames(x),
                       "condition" = c( rep(1,length(cells_background)), rep(2,length(cells_of_interest)) ),
                       "libType" = rep("single-end", dim(x)[2]))
    cds <- newCountDataSet( round(x,0), des$condition )
    cds <- estimateSizeFactors( cds )
    cds <- estimateDispersions( cds, method=method, fitType="local" )
    res <- nbinomTest( cds, 1, 2 )
    rownames(res) <- res$id
    res <- res[,-1]
    list("des"=des, "cds"=cds, "res"=res)
  }else{
    if (norm) x <- as.data.frame( t(t(x)/apply(x,2,sum)) * median(apply(x,2,sum,na.rm=TRUE)) )
    fit = Meanz = Medianz = v = list()
    for ( i in 1:2 ){
      group <- if ( i == 1 ) cells_background else cells_of_interest
      Meanz[[i]]	<- signif(if ( length(group) > 1 ) apply(x[,group],1,mean) else x[,group], digits = 3)
      # Mean_AllSamples =	signif(apply(x, 1, mean), digits=2) # Not used as other functions rely on the old calculation
      Medianz[[i]]	<- if ( length(group) > 1 ) apply(x[,group],1,median) else x[,group]
      v[[i]] <- if ( length(group) > 1 ) apply(x[,group],1,var)  else apply(x,1,var)
      
      if ( method == "pooled"){
        mg <- apply(x,1,mean)
        vg <- apply(x,1,var)
        f <- vg > 0 & mg > .5
        logv <- log2(vg[f])
        logm <- log2(mg[f])
      }else{
        f <- v[[i]] > 0 & Meanz[[i]] > .5
        logv <- log2(v[[i]][f])
        logm <- log2(Meanz[[i]][f])
      }
      
      if ( locreg ){
        f <- order(logm,decreasing=FALSE)
        u <- 2**logm[f]
        y <- 2**logv[f]
        lf <- locfit(y~lp(u,nn=.7), family="gamma", maxk=500)
        fit[[i]] <- approxfun(u, fitted(lf), method = "const")
      }else{
        fit[[i]] <- if ( is.null(vfit) ) lm(logv ~ logm + I(logm^2)) else vfit
      }
    }
    
    if ( locreg ){
      vf  <- function(x,i) fit[[i]](x)
    }else{
      vf  <- function(x,i) 2**(coef(fit[[i]])[1] + log2(x)*coef(fit[[i]])[2] + coef(fit[[i]])[3] * log2(x)**2)
    }
    sf  <- function(x,i) x**2/(max(x + 1e-6,vf(x,i)) - x)
    
    pv <- apply("X" = data.frame(Meanz[[1]],Meanz[[2]]), "MARGIN" = 1, "FUN" = function(x){
      p12 <- dnbinom(	"x" = 0:round(x[1]*length(cells_background) + x[2]*length(cells_of_interest),0),
                      "mu" = mean(x)*length(cells_background),size=length(cells_background)*sf(mean(x),1))*
        dnbinom(	"x" = round(x[1]*length(cells_background) + x[2]*length(cells_of_interest),0):0,
                 "mu" = mean(x)*length(cells_of_interest),size=length(cells_of_interest)*sf(mean(x),2));
      
      sum(p12[p12 <= p12[round(x[1]*length(cells_background),0) + 1]])/sum(p12)} )
    
    if ( sum(pv=="NaN") ) { # Mauro's fix
      pv[which(pv=="NaN")] <-0.000000e+00
      cat(paste(length(which(pv=="NaN")),"NaNs were replaced by 0.\n"))
    }
    
    foldChange = Meanz[[2]]/Meanz[[1]]
    res <- data.frame(	"baseMean" =		signif((Meanz[[1]] + Meanz[[2]])/2, digits = 2),
                       "baseMeanA" =		Meanz[[1]],
                       "baseMeanB" = 		Meanz[[2]],
                       "foldChange" = 		round(foldChange, digits = 2),
                       "log2FoldChange" = 	round(log2(foldChange), digits = 2),
                       "MedianA"=			Medianz[[1]],
                       "MedianB"=			Medianz[[2]],
                       "foldChange_Median"=signif(Medianz[[2]]/Medianz[[1]], digits = 1),
                       "pval" = 			signif(pv, digits = 2),
                       "padj" = 			signif(p.adjust(pv,method="BH"), digits = 1))
    vf1 <- data.frame("m" = Meanz[[1]], "v" = v[[1]], "vm" = vf(Meanz[[1]],1))
    vf2 <- data.frame("m" = Meanz[[2]],"v" = v[[2]], "vm" = vf(Meanz[[2]],2))
    rownames(res) <- rownames(vf1) <- rownames(vf2) <- rownames(x)
    list("vf1"=vf1, "vf2"=vf2, "res"=res)
  }
}

setGeneric("clustdiffgenesnb", function(object,pvalue=.01) standardGeneric("clustdiffgenesnb"))

setMethod("clustdiffgenesnb",
          signature = "SCseq",
          definition = function(object,pvalue){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before clustdiffgenesnb")
            if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            part  <- object@cpart
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              g1 <- names(part)[part == i]
              g2 <- names(part)[part != i]
              y <- diffexpnb(x,g1,g2,vfit=sc@background$vfit)$res
              d <- data.frame(mean.ncl=y$baseMeanB,mean.cl=y$baseMeanA,fc=1/y$foldChange,pv=y$padj)
              rownames(d) <- rownames(x)
              d <- d[order(d$pv,decreasing=FALSE),]
              cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
            }
            return(cdiff)
          }
          )



scshuffle <- function(sc,lp,metric="pearson",all=FALSE){
  fl <- TRUE
  if ( all ){
    z <-  data.frame( t(apply(as.data.frame(sc@fdata),1,sample)) )
    names(z) <- names(sc@fdata)
  }else{
    for ( i in unique(sc@cpart)){
      if ( sum(sc@cpart == i) > 1 ){
        y <- data.frame( t(apply(as.data.frame(sc@fdata[,sc@cpart == i]),1,sample)) )
      }else{
        y <- as.data.frame(sc@fdata[,sc@cpart == i])
      }
      ##y <- as.data.frame(y[sample(1:nrow(y)),])
      names(y) <- names(sc@cpart)[sc@cpart == i]
      rownames(y) <- rownames(sc@fdata)
      z <- if (fl) y else cbind(z,y)
      fl <- FALSE
    }
  }
  di <- as.matrix(dist.gen(t(z[,names(sc@cpart)]),method=metric))
  return(di[names(lp),names(lp)])
}

pdishuffle <- function(pdi,lp,cn,m,all=FALSE){
  if ( all ){
    d <- as.data.frame(pdi)
    ##y <- t(apply( cn[sapply(lp,function(x,m) which(x==m),m=m),],1,function(x) runif(length(x),min=-1,max=1) + x ))
    y <- t(apply(pdi,1,function(x) runif(length(x),min=min(x),max=max(x))))
    ##y <- t(apply(pdi,1,function(x) sample(x)))
    ##y <- t(apply(t(apply(pdi,1,sample)),2,sample))
    names(y)    <- names(d)
    rownames(y) <- rownames(d)
    ##diag(y) <- 0
    return(y)
  }else{
    fl <- TRUE
    for ( i in unique(lp)){
      if ( sum(lp == i) > 1 ){
        y <-data.frame( t(apply(as.data.frame(pdi[,lp == i]),1,sample)) )
      }else{
        y <- as.data.frame(pdi[,lp == i])
      }
      names(y) <- names(lp)[lp == i]
      rownames(y) <- names(lp)
      z <- if (fl) y else cbind(z,y)
      fl <- FALSE
    }
    z <- t(z[,names(lp)])
    return(z)
  }
}

compproj <- function(pdiloc,lploc,cnloc,mloc,d=NULL){
  pd    <- data.frame(pdiloc)
  k     <- paste("X",sort(rep(1:nrow(pdiloc),length(mloc))),sep="")
  pd$k  <- paste("X",1:nrow(pdiloc),sep="")
  pd    <- merge(data.frame(k=k),pd,by="k")

  if ( is.null(d) ){
    cnv   <- t(matrix(rep(t(cnloc),nrow(pdiloc)),nrow=ncol(pdiloc)))
    pdcl  <- paste("X",lploc[as.numeric(sub("X","",pd$k))],sep="")
    rownames(cnloc) <- paste("X",mloc,sep="")
    pdcn  <- cnloc[pdcl,]
    v     <- cnv - pdcn
  }else{
    v    <- d$v
    pdcn <- d$pdcn
  }
  w <- pd[,names(pd) != "k"] - pdcn

  h <- apply(cbind(v,w),1,function(x){
    x1 <- x[1:(length(x)/2)];
    x2 <- x[(length(x)/2 + 1):length(x)];
    x1s <- sqrt(sum(x1**2)); x2s <- sqrt(sum(x2**2)); y <- sum(x1*x2)/x1s/x2s; return( if (x1s == 0 | x2s == 0 ) NA else y ) })

  rma <- as.data.frame(matrix(h,ncol=nrow(pdiloc)))
  names(rma) <- unique(pd$k)
  pdclu  <- lploc[as.numeric(sub("X","",names(rma)))]
  ##pdclp  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else mloc[if (max(x,na.rm=TRUE)>0) which(x == max(x,na.rm=TRUE)) else which(abs(x) == max(abs(x),na.rm=TRUE))])
  ##pdclh  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else x[if (max(x,na.rm=TRUE)>0) which(x == max(x,na.rm=TRUE)) else which(abs(x) == max(abs(x),na.rm=TRUE))])
  pdclp  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else mloc[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdclh  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else x[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdcln  <-  names(lploc)[as.numeric(sub("X","",names(rma)))]
  names(rma) <- pdcln
  rownames(rma) <- paste("X",m,sep="")
  res    <- data.frame(o=pdclu,l=pdclp,h=pdclh)
  rownames(res) <- pdcln
  return(list(res=res[names(lploc),],rma=as.data.frame(t(rma[,names(lploc)])),d=list(v=v,pdcn=pdcn)))
}



compproj_old <- function(pdiloc,lploc,cnloc,mloc){
  pd    <- data.frame(pdiloc)
  k     <- paste("X",sort(rep(1:nrow(pdiloc),length(mloc))),sep="")
  pd$k  <- paste("X",1:nrow(pdiloc),sep="")
  pd    <- merge(data.frame(k=k),pd,by="k")
  cnv   <- t(matrix(rep(t(cnloc),nrow(pdiloc)),nrow=ncol(pdiloc)))
  pdcl  <- lploc[as.numeric(sub("X","",pd$k))]

  pdcli <- apply(data.frame(pdcl),1,function(x,m) which(x == m),m=mloc)
  pdcn  <- cnloc[pdcli,]
  w <- pd[,names(pd) != "k"] - cnloc[pdcli,]
  v <- cnv - cnloc[pdcli,]

  lv <- sqrt(apply(v*v,1,sum))
  lw <- sqrt(apply(w*w,1,sum))
  lv[lv == 0 | lw == 0] <- NA
  h <- apply(v*w,1,sum)/lv/lw

#  pdn <- rep(1:length(mloc),ncol(pdiloc))
#  y <- cbind(pd[,names(pd) != "k"], cnv, pdcn)
#  h <- apply(y,1,function(x){x <- as.numeric(x); u <- x[(2*ncol(cnloc) + 1):(3*ncol(cnloc))]; w <- x[1:ncol(cnloc)] - u; v <- x[(ncol(cnloc) + 1):(2*ncol(cnloc))] - u; lv <- sqrt(sum(v*v)); lw <- sqrt(sum(w*w)); return( if ( min(lv,lw) == 0 ) NA else sum(v*w)/lv/lw ) })
  #rma <- as.data.frame(matrix(h,ncol=ncol(pdiloc)))
  rma <- as.data.frame(matrix(h,ncol=nrow(pdiloc)))
  names(rma) <- unique(pd$k)
  pdclu  <- lploc[as.numeric(sub("X","",names(rma)))]
  ##pdclp  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else mloc[if (max(x,na.rm=TRUE)>0) which(x == max(x,na.rm=TRUE)) else which(abs(x) == max(abs(x),na.rm=TRUE))])
  ##pdclh  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else x[if (max(x,na.rm=TRUE)>0) which(x == max(x,na.rm=TRUE)) else which(abs(x) == max(abs(x),na.rm=TRUE))])
  pdclp  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else mloc[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdclh  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else x[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdcln  <-  names(lploc)[as.numeric(sub("X","",names(rma)))]
  names(rma) <- pdcln
  rownames(rma) <- paste("X",m,sep="")
  res    <- data.frame(o=pdclu,l=pdclp,h=pdclh)
  rownames(res) <- pdcln
  return(list(res=res[names(lploc),],rma=as.data.frame(t(rma[,names(lploc)]))))
}

cellsfromtree <- function(prtr,z,sc,cthr){
  f <- c()
  g <- c()
  ##z <- sapply(z,function(x,m) which( x == m ),m= getm(sc,cthr))
  for ( i in 1:( length(z) - 1 ) ){
    rf <- if ( z[i+1] > z[i] ) FALSE else TRUE
    k <- if ( rf ) paste(z[i + 1],z[i],sep=".") else paste(z[i],z[i+1],sep=".")
    p <- prtr$l[[k]]
    n <- prtr$n[[k]]
    if ( rf ){
      ##h <- p < Inf & p > -Inf
      if ( i == 1 & i + 1 == length(z) ) h <- p < Inf & p > -Inf
      if ( i == 1 & i + 1 <  length(z) ) h <- p < Inf & p >= 0
      if ( i >  1 & i + 1 == length(z) ) h <- p <= 1  & p > -Inf
      if ( i >  1 & i + 1 <  length(z) ) h <- p <= 1  & p >= 0
    }else{
      ##h <- p > -Inf & p <  Inf
      if ( i == 1 & i + 1 == length(z) ) h <- p > -Inf & p <  Inf
      if ( i == 1 & i + 1 <  length(z) ) h <- p > -Inf & p <= 1
      if ( i >  1 & i + 1 == length(z) ) h <- p >= 0   & p <  Inf
      if ( i >  1 & i + 1 <  length(z) ) h <- p >= 0   & p <= 1
    }
    v <- n[h][order(p[h],decreasing=FALSE)]
    if ( rf ) v <- rev(v)
    f <- append(f,v)
    g <- append(g,rep(i,length(v)))
  }
  return(list(f=f,g=g))
}

cellsfromtree2 <- function(rma,z,sc,pthr,prbackn,sig=FALSE){
  z <- sapply(z,function(x,m) which( x == m ),m=getm(sc,pthr))
  l <- list()
  for ( i in 1:length(z) ){
    if ( i < length(z) ) l[[i]] <- list(n=c(),p=c())
    if ( i == 1 ){
      g <- prback$o == z[i] & prback$l == z[i+1] & !is.na(prback$o)
      if ( sum(g) == 0 | ! sig ){
        ql <- 0
        qh <- 0
      }else{
        ql <- quantile(prback[g,"h"],pthr)
        qh <- quantile(prback[g,"h"],1 - pthr)
      }
      x <- rma[sc@cpart[rownames(rma)] == z[i],]
      f <- x[,z[i+1]] > qh | ( x[,z[i+1]] < 0 & x[,z[i+1]] < ql ) | is.na(x[,z[i+1]])
      if ( sum(f,na.rm=TRUE) > 0 ){
        l[[i]]$n <- append(l[[i]]$n,rownames(x)[f])
        l[[i]]$p <- append(l[[i]]$p,x[f,z[i+1]])
      }
    }
    if ( i > 1 & i < length(z)){
      g <- prback$o == z[i] & prback$l == z[i-1] & !is.na(prback$o)
      if ( sum(g) == 0 | ! sig ){
        ql1 <- 0
        qh1 <- 0
      }else{
        ql1 <- quantile(prback[g,"h"],pthr)
        qh1 <- quantile(prback[g,"h"],1 - pthr)
      }
      g <- prback$o == z[i] & prback$l == z[i+1] & !is.na(prback$o)
      if ( sum(g) == 0 | ! sig ){
        ql2 <- 0
        qh2 <- 0
      }else{
        ql2 <- quantile(prback[g,"h"],pthr)
        qh2 <- quantile(prback[g,"h"],1 - pthr)
      }
      x <- rma[sc@cpart[rownames(rma)] == z[i],z[c(i-1,i+1)]]
      y <- apply(x,1,function(x) if (sum(is.na(x))>0) NA else which( max(abs(x),na.rm=TRUE) == abs(x) ))
      f <- ( y == 1 & ( x[,1] > qh1 | ( x[,1] < 0 & x[,1] < ql1 ) ) ) & !is.na(x[,1])
      if ( sum(f,na.rm=TRUE) > 0 ){
        l[[i-1]]$n <- append(l[[i-1]]$n,rownames(x)[f])
        l[[i-1]]$p <- append(l[[i-1]]$p,1 - x[f,1])
      }
      f <- ( y == 2 & ( x[,2] > qh2 | ( x[,2] < 0 & x[,2] < ql2 ) ) ) | is.na(x[,2])
      if ( sum(f,na.rm=TRUE) > 0 ){
        l[[i]]$n <- append(l[[i]]$n,rownames(x)[f])
        l[[i]]$p <- append(l[[i]]$p,x[f,2])
      }
    }
    if ( i == length(z)){
      g <- prback$o == z[i] & prback$l == z[i-1] & !is.na(prback$o)
      if ( sum(g) == 0 | ! sig ){
        ql <- 0
        qh <- 0
      }else{
        ql <- quantile(prback[g,"h"],pthr)
        qh <- quantile(prback[g,"h"],1 - pthr)
      }

      x <- rma[sc@cpart[rownames(rma)] == z[i],]
      f <- ( x[,z[i-1]] > qh | ( x[,z[i-1]] < 0 & x[,z[i-1]] < ql ) ) & !is.na(x[,z[i-1]])
      if ( sum(f,na.rm=TRUE) > 0 ){
        l[[i-1]]$n <- append(l[[i-1]]$n,rownames(x)[f])
        l[[i-1]]$p <- append(l[[i-1]]$p,1-x[f,z[i-1]])
      }
      f <- is.na(x[,z[i-1]])
      if ( sum(f,na.rm=TRUE) > 0 ){
        l[[i-1]]$n <- append(l[[i-1]]$n,rownames(x)[f])
        l[[i-1]]$p <- append(l[[i-1]]$p,rep(1,sum(f)))
      }

    }
  }
  tr <- c()
  for ( i in 1:length(l) ){
    if ( length(l[[i]]$p) > 0 ){
      f <- is.na(l[[i]]$p)
      l[[i]]$p[f] <- 0
      f <- order(l[[i]]$p,decreasing=FALSE)
      l[[i]]$p <- l[[i]]$p[f]
      l[[i]]$n <- l[[i]]$n[f]
      tr <- append(tr,l[[i]]$n)
    }
  }
  return(list(tr=tr,l=l))
}

getm <- function(sc,cthr){
  lp <- sc@cpart
  n <- aggregate(rep(1,length(lp)),list(lp),sum)
  n <- as.vector(n[order(n[,1],decreasing=FALSE),-1])
  return( (1:length(n))[n>cthr] )
}


kmedians <- function(x, n, iter.max=100, verbose=FALSE){
  if ( is.null(dim(n)) ){
    centroids <- x[,sample(1:ncol(x),n)]
  }else{
    centroids <- n
    n <-ncol(centroids)
  }
  ##set.seed(seed)
  clusters  <- rep(0,ncol(x));
  iter <- 0;
  old_clusters <- clusters
  while (iter < iter.max ){
    clusters <- apply(x,2,function(x,y){ u <- apply(y,2,function(y,z) sqrt(sum((y-z)**2)),z=x); which(u == min(u)) },y=centroids)
    for ( i in 1:n ){
      f <- clusters==i;
      centroids[,i] <- smedian(x[,f])
    }
    if ( ( iter > 0 ) & sum( old_clusters != clusters ) == 0 ){
      if ( verbose ) cat("kmedians: Convergence after ",iter," iterations\n")
      break
    }
    old_clusters = clusters
    iter <- iter + 1
  }
  return(list(clusters=clusters,iter.max=iter.max,iter=iter,n=n))
}

smedian <- function(x,iter.max=100,verbose=FALSE){
  ## calculates the spatial (or multivariate) L_1-median of X (data dimension x samples matrix)
  ## it's based on the modified Weiszfeld algorithm of Vardi and Zhang
  ## for simplicity WLOG y~=X(:,i) otherwise slight modification necessary (see Vardi paper)
  ## if verbose is set (1) convergence details are given (default: 0)
  ##
  ## 7.7.2005 Fabian Theis (fabian@theis.name)
  if ( is.null(dim(x)) ) return(x)
  y  <- apply(x,1,mean)
  de <- sqrt(2.2204e-16)
  for ( iter in 1:iter.max ){
    dx <- sqrt(apply((x - y)**2,2,sum))
    a  <- apply(t(t(x)/dx),1,sum);
    b  <- sum(1/dx)
    T <- a/b
    ynew <- T;
    if ( sqrt(sum((y-ynew)**2)) <= de ){
      if ( verbose ) cat("smedian: convergence after ",iter," iterations\n")
      return(ynew)
    }
    y <- ynew
  }
  if ( verbose ) cat("smedian: warning - no convergence after ",iter.max," iterations\n")
}

## rows are clustered !!


setGeneric("plotsaturation", function(object,disp=FALSE) standardGeneric("plotsaturation"))

setMethod("plotsaturation",
          signature = "SCseq",
          definition = function(object,disp){
            if ( length(object@cluster$gap) == 0 ) stop("run clustexp before plotsaturation")
            g <- object@cluster$gap$Tab[,1]
            y <- g[-length(g)] - g[-1]
            mm <- numeric(length(y))
            nn <- numeric(length(y))
            for ( i in 1:length(y)){
              mm[i] <- mean(y[i:length(y)])
              nn[i] <- sqrt(var(y[i:length(y)]))
            }
            cln <- max(min(which( y - (mm + nn) < 0 )),1)
            x <- 1:length(y)
            if (disp){
              x <- 1:length(g)
              plot(x,g,pch=20,col="grey",xlab="k",ylab="log within cluster dispersion")
              f <- x == cln
              points(x[f],g[f],col="blue")
            }else{
              plot(x,y,pch=20,col="grey",xlab="k",ylab="Change in log within cluster dispersion")
              points(x,mm,col="red",pch=20)
              plot.err.bars.y(x,mm,nn,col="red")
              points(x,y,col="grey",pch=20)
              f <- x == cln
              points(x[f],y[f],col="blue")
            }
          }
          )




setGeneric("clustdiffgenes2", function(object,pvalue=.01) standardGeneric("clustdiffgenes2"))

setMethod("clustdiffgenes2",
          signature = "SCseq",
          definition = function(object,pvalue){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before clustdiffgenes")
            if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@cpart
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              m <-  if ( sum(part != i) > 1 ) apply(x[,part != i],1,mean) else x[,part != i]
              n <-  if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
              no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
              m <- m*no
              n <- n*no
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.ncl=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
              cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
            }
            return(cdiff)
          }
          )

setGeneric("plotsilhouette2", function(object) standardGeneric("plotsilhouette2"))

setMethod("plotsilhouette2",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotsilhouette")
            if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
            kpart <- object@cluster$kpart
            distances  <- if ( object@clusterpar$FUNcluster == "kmedoids" ) as.dist(object@distances) else dist.gen(object@distances)
            si <- silhouette(kpart,distances)
            plot(si)
          }
          )


compmedoids <- function(x,part,metric="pearson"){
  m <- c()
  for ( i in sort(unique(part)) ){
    f <- names(x)[part == i]
    if ( length(f) == 1 ){
      m <- append(m,f)
    }else{
      y <- apply(as.matrix(dist.gen(t(x[,f]),method=metric)),2,mean)
      m <- append(m,f[which(y == min(y))[1]])
    }
  }
  m
}




plotnetwork2 <- function(x){
  x <- (x + t(x))/2
  ramp  <- colorRamp(c( "yellow","blue"))
  mcol <- rgb( ramp(seq(0, 1, length = 101)), max = 255)

  #mcol <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(101)
  lev  <- seq(min(x), max(x), length=length(mcol))
  dcc <- t(apply(round(100*(x - min(x))/(max(x) - min(x)),0) + 1,1,function(x){y <- c(); for ( n in x ) y <- append(y,mcol[n]); y }))
  layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
  par(mar = c(3,1,1,1))
  plot(network(as.matrix(x)),edge.col=dcc,edge.lwd=2,vertex.cex=1,vertex.col="grey",displaylabels=TRUE,label=id2name(names(x)),label.cex=.75,usearrows=FALSE)
  par(mar = c(20,2.5,2.5,2))
  image(1, lev, matrix(data=lev, ncol=length(lev),nrow=1), col=mcol, xlab="", ylab="", xaxt="n")
  layout(1)
}


barplotgene <- function(sc,g){
  f <- grep(g,rownames(sc@fdata))
 if ( length(f) == 1 ){
   x <- t(sc@fdata[f,names(sort(sc@cpart))])
 }else{
   x <- apply(sc@fdata[f,names(sort(sc@cpart))],2,sum)
 }
 names(x) <- sort(sc@cpart)
 barplot(x,beside=TRUE,border=FALSE,col="white",main=g,names.arg=rep("",length(x)))
 for ( i in unique(sc@cpart)){ y <- x; y[sort(sc@cpart) != i] <- 0; barplot(y,col=sc@fcol[i],beside=TRUE,add=TRUE,border=FALSE,names.arg=rep("",length(x)),axes=FALSE)}
}

plotsat <- function(g,plot=TRUE){
  y <- g[-length(g)] - g[-1]
  mm <- numeric(length(y))
  nn <- numeric(length(y))
  for ( i in 1:length(y)){
    mm[i] <- mean(y[i:length(y)])
    nn[i] <- sqrt(var(y[i:length(y)]))
  }
  cln <- max(min(which( y < mm + nn  & y > mm - nn  )),1)
  if ( plot ){
    x <- 1:length(y)
    plot(x,y,pch=20,col="grey",xlab="k",ylim=c(min(c(y,mm - nn),na.rm=TRUE),max(c(y,mm + nn),na.rm=TRUE)),ylab="Change in variance of distance ratio")
    points(x,mm,col="red",pch=20)
    plot.err.bars.y(x,mm,nn,col="red")
    points(x,y,col="grey",pch=20)
    f <- x == cln
    points(x[f],y[f],col="blue")
  }
  return(cln)
}



branchcells <- function(s,br){
  cl <- intersect( as.numeric(strsplit(br[[1]],"\\.")[[1]]), as.numeric(strsplit(br[[2]],"\\.")[[1]]))
  n <- list()
  scl <- sco[[s]]
  k <- c()
  for ( i in 1:length(br) ){
    n[[i]] <- names(sco[[s]]@cpart[ prtree[[s]]$n[[br[[i]]]]])[sco[[s]]@cpart[ prtree[[s]]$n[[br[[i]]]]] == cl]
    k <- append(k, max( scl@cpart ) + 1)
    scl@cpart[n[[i]]] <- max( scl@cpart ) + 1
  }
  set.seed(111111)
  scl@fcol <- sample(rainbow(max(scl@cpart)))

  return( list(n=n,scl=scl,k=k) )
}


plotsom <- function(smp,ylim=NULL){
#  if ( is.null(ylim) ) ylim <- c( min(smp$data,na.rm=TRUE),max(smp$data,na.rm=TRUE) )
   plot(c(0,smp$xdim),c(0,smp$ydim),cex=0,xlab="x",ylab="y",main="")
   for ( i.x in 1:smp$xdim){
     for ( i.y in 1:smp$ydim ){
       rect(i.x - 1, i.y - 1, i.x, i.y, col="grey")
     }
   }
   if ( is.null(ylim) ){
     ylim <- c(Inf,-Inf)
     for ( i.x in 1:smp$xdim){
       for ( i.y in 1:smp$ydim ){
         f <- smp$visual$x == ( i.x - 1 ) & smp$visual$y == ( i.y - 1 )
         if ( sum(f) == 0 ) next
         if ( sum(f) != 1 ){
           m  <- apply(smp$data[f,],2,mean)
           sd <- sqrt(apply(smp$data[f,],2,var))
         }else{
           m  <- t(smp$data[f,])
           sd <- rep(0,ncol(smp$data))
         }
         ylim[1] <- min(ylim[1],m-sd)
         ylim[2] <- max(ylim[2],m+sd)
       }
     }
   }

   for ( i.x in 1:smp$xdim){
     for ( i.y in 1:smp$ydim ){
       l <- list()
       f <- smp$visual$x == ( i.x - 1 ) & smp$visual$y == ( i.y - 1 )
       if ( sum(f) != 1 ){
         l[["av"]] <-      apply(smp$data[f,],2,mean)
         sd        <- sqrt(apply(smp$data[f,],2,var))
       }else{
         l[["av"]] <-      t(smp$data[f,])
         sd        <- rep(0,ncol(smp$data))
       }

       l[["up"]] <- l[["av"]] + sd
       l[["lo"]] <- l[["av"]] - sd
       for ( q in c("lo","up","av") ){
         x.c <- ( 1:length(l[[q]]) )/length(l[[q]]) + ( i.x - 1 )
         l[[q]][l[[q]] < ylim[1]] <- ylim[1]
         l[[q]][l[[q]] > ylim[2]] <- ylim[2]
         l[[q]] <- ( l[[q]] - ylim[1] )/abs(ylim[1] - ylim[2]) + ( i.y - 1 )
         if ( q == "av" ){
           lines(x.c,l[[q]],col="black")
           #lines(vec.2.inv(x.c,x.c[2] - x.c[1]),vec.2(l[[q]]),col="black",lwd=.5)
         }else{
           lines(x.c,l[[q]],col="white")
           #lines(vec.2.inv(x.c,x.c[2] - x.c[1]),vec.2(l[[q]]),col="white",lwd=.5)
         }
       }
     }
   }
   for ( i.x in 1:smp$xdim){
     for ( i.y in 1:smp$ydim ){
       rect(i.x - 1, i.y - 1, i.x, i.y)
     }
   }
 }
