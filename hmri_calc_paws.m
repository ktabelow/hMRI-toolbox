function [] = hmri_calc_paws(ESTATICSmodel, mpmData, mask, kstar, patchsize, ladjust)

mscbw = 5;
alpha = 0.025; % this is not needed, could be additional input to the qf() function if desired
wghts = [];

sdim = size(ESTATICSmodel.R2s);

nv = ESTATICSmodel.nv;
if isempty(mask)
    mask= ones(sdim);
end
nvoxel = prod(sdim);

% begin consistency checks
  if(any(dim(mpmESTATICSModel$invCov)!=c(nv,nv,nvoxel))) stop("inconsistent invCov")
  if(any(dim(mask)!=sdim)) stop("inconsistent mask")
  if(any(dim(mpmESTATICSModel$modelCoeff)!=c(nv,nvoxel))) stop("inconsistent parameter length")
  if(!is.null(mpmData)){#2
    # allow for mpmData to be expanded or to only contain data within mask
    if(any(dim(mpmData)[-1]!=nvoxel)){#3
      if(all(dim(mpmData)[-1]==sdim)){#4
        #  reduce mpmData to voxel within mask
        dim(mpmData) <- c(dim(mpmData)[1],prod(sdim))
        mpmData <- mpmData[,mask]
      } else stop("inconsistent mpmData")
    }#3
  }#2
% end consistency checks


% determine a suitable adaptation bandwidth
switch patchsize
    case 1
        corr_fac_patchsize = 1
    case 2
        corr_fac_patchsize = 2.77
    case 3
        corr_fac_patchsize = 3.46
end
% factor 2 (analog to 2 sigma in KL) to have more common values for alpha
% corr_fac_patchsize adjusted using simulated data
lambda = ladjust * 2 * nv * qf(nv, ESTATICSmodel.necho - nv) * corr_fac_patchsize;

% begin write to console
if(verbose) cat("using lambda=", lambda, " patchsize=", patchsize,"\n")
% end write to console

    invCov <- extract(mpmESTATICSModel,"invCov")
  if(mscbw>0){
    rsdx <- extract(mpmESTATICSModel,"rsigma")^2
    rsdhat <- medianFilter3D(rsdx,mscbw,mask)
    rsdhat[!mask] <- mean(rsdhat[mask])
    invCov <- sweep(invCov,3:5,rsdx/rsdhat,"*")
  }
  dim(invCov) <- c(nv,nv,prod(sdim))
  invCov <- invCov[,,mask]
  zobj <- vpawscov2(mpmESTATICSModel$modelCoeff,
                    kstar,
                    invCov,
                    mask,
                    lambda = lambda,
                    wghts = wghts,
                    patchsize = patchsize,
                    data = mpmData,
                    verbose = verbose)
  ## assign values
  obj <- list(modelCoeff = zobj$theta,
              invCov = mpmESTATICSModel$invCov,
              isConv = mpmESTATICSModel$isConv,
              rsigma = mpmESTATICSModel$rsigma,
              bi = zobj$bi,
              smoothPar = c(zobj$lambda, zobj$hakt, alpha, patchsize, mscbw),
              smoothedData = zobj$data,
              sdim = mpmESTATICSModel$sdim,
              nFiles = mpmESTATICSModel$nFiles,
              t1Files = mpmESTATICSModel$t1Files,
              pdFiles = mpmESTATICSModel$pdFiles,
              mtFiles = mpmESTATICSModel$mtFiles,
              model = mpmESTATICSModel$model,
              maskFile = mpmESTATICSModel$maskFile,
              mask = mask,
              sigma = mpmESTATICSModel$sigma,
              L = mpmESTATICSModel$L,
              TR = mpmESTATICSModel$TR,
              TE = mpmESTATICSModel$TE,
              FA = mpmESTATICSModel$FA,
              TEScale = mpmESTATICSModel$TEScale,
              dataScale = mpmESTATICSModel$dataScale)
  class(obj) <- "sESTATICSModel"
  invisible(obj)
  
  ## END function smoothESTATICS()
}

vpawscov2 <- function(y,
                      kstar = 16,
                      invcov = NULL,
                      mask = NULL,
                      scorr = 0,
                      spmin = 0.25,
                      lambda = NULL,
                      ladjust = 1,
                      wghts = NULL,
                      patchsize = 1,
                      data = NULL,
                      verbose = TRUE) {#1
  ##
  ##  this is the version with full size invcov (triangular storage)
  ##  and optional smoothing of vector-valued images supplied in data
  ##  for internal use in package qMRI
  ##  Uses condensed data (voxel within mask only)
  ##  returns a list
  ##
  dy <- dim(y)
  nvec <- dy[1]
  if(nvec>5) stop("limited to 5 parameters")
  indcov <- switch(nvec,1,
                   c(1,2,4),
                   c(1,2,5,3,6,9),
                   c(1,2,6,3,7,11,4,8,12,16),
                   c(1,2,7,3,8,13,4,9,14,19,5,10,15,20,25))
  if(!is.null(data)) nsample <- dim(data)[1]
  dy <- dim(mask)
  d <- length(dy)
  if (d != 3)
    stop("need 3D mask")
  if(is.null(lambda)){#2
    lambda <- 2 * ladjust * qchisq(pchisq(8.82, 1), nvec)
    lambda <- lambda * switch(patchsize+1,1,1.3,1.6)
  }#2
  if (is.null(wghts)) wghts <- c(1, 1, 1)
  wghts <- wghts[1] / wghts[2:3]
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  h0 <- 0
  if (any(scorr > 0)) {#3
    h0 <- numeric(length(scorr))
    for (i in 1:length(h0))
      h0[i] <- geth.gauss(scorr[i])
    if (length(h0) < d)
      h0 <- rep(h0[1], d)
    if(verbose) cat("Corresponding bandwiths for specified correlation:",
                    h0,
                    "\n")
  }#3
  ## create index information for voxel in mask
  nvoxel <- sum(mask)
  position <- array(0,dy)
  position[mask] <- 1:nvoxel
  dim(mask) <- NULL
  dim(y) <- c(nvec,nvoxel)
  dim(invcov) <- c(nvec * nvec,nvoxel)
  hseq <- 1
  zobj <- list(bi = rep(1, nvoxel), theta = y)
  bi <- zobj$bi
  if(verbose) cat("Progress:")
  total <- cumsum(1.25 ^ (1:kstar)) / sum(1.25 ^ (1:kstar))
  mc.cores <- setCores(, reprt = FALSE)
  np1 <- 2 * patchsize + 1
  np2 <- if (n2 > 1) 2 * patchsize + 1 else 1
  np3 <- if (n3 > 1) 2 * patchsize + 1 else 1
  k <- 1
  hmax <- 1.25 ^ (kstar / d)
  lambda0 <- 1e32
  mae <- NULL
  while (k <= kstar) {#4
    hakt0 <- gethani(1, 1.25 * hmax, 2, 1.25 ^ (k - 1), wghts, 1e-4)
    hakt <- gethani(1, 1.25 * hmax, 2, 1.25 ^ k, wghts, 1e-4)
    if(verbose) cat("step", k, "hakt", hakt, "time", format(Sys.time()), "\n")
    hseq <- c(hseq, hakt)
    dlw <- (2 * trunc(hakt / c(1, wghts)) + 1)
    if(k==kstar & !is.null(data)){#5
      dim(data) <- c(nsample,nvoxel)
      zobj <- .Fortran(C_pvawsme,
                       as.double(y),
                       as.double(data), ## data to smooth additionally
                       as.integer(position),
                       as.integer(nvec),
                       as.integer(nvec * (nvec + 1) / 2),
                       as.integer(nsample), ## leading dimension of data
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       as.double(zobj$bi),
                       bi = double(nvoxel), #binn
                       theta = double(nvec * nvoxel),
                       data = double(nsample*nvoxel),
                       as.double(invcov[indcov,]),#
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       double(nsample * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt","data")]
      dim(zobj$data) <- c(nsample, nvoxel)
    } else {#6
      zobj <- .Fortran(C_pvaws2,
                       as.double(y),
                       as.integer(position),
                       as.integer(nvec),
                       as.integer(nvec * (nvec + 1) / 2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       as.double(zobj$bi),
                       bi = double(nvoxel), #binn
                       theta = double(nvec * nvoxel),
                       as.double(invcov[indcov,]),# compact storage
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt")]
    }#6
    x <- 1.25 ^ k
    scorrfactor <- x / (3 ^ d * prod(scorr) * prod(h0) + x)
    lambda0 <- lambda * scorrfactor
    if (verbose & max(total) > 0) {#7
      cat(signif(total[k], 2) * 100, "%  ", sep = "")
      cat("mean(bi)", signif(mean(zobj$bi),3)," ")
    }#7
    k <- k + 1
    gc()
  }
  dim(zobj$theta) <- c(nvec, nvoxel)
  if(verbose) cat("\n")
  list(
    theta=zobj$theta,
    hakt=hakt,
    lambda=lambda,
    hseq = hseq,
    bi = zobj$bi,
    data= if(!is.null(zobj$data)) zobj$data else NULL
  )
}
    end

function qval = qf(df1, df2)
quantile = ones(4, 40);
quantile(1,:)=
quantile(2,:)=
quantile(3,:)=
quantile(4,:)=
qval = quantile(df1,df2);
end

