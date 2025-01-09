function [] = hmri_calc_paws(ESTATICSmodel, mpmData, mask, kstar, patchsize, ladjust)
  
  % DEFINE ALL CONSTANTS
  spmin = 0.25, % the statistical kernel function is a plateau to spmin with linear decrease till 1
  lambda0 = 1e32; % the first iteration step uses this adaptation parameter lambda in order to create a stable non-adaptive first estimate
  mc.cores = 1; % number of cores for OMP parallel 
  corr_fac_patchsize = [1, 2.77, 3.46]; % adjustment factor for adaptation bandwidth for different patchsizes (1, 2, 3), determined using simulated data

  % EXTRACT ALL REQUIRED VALUES FROM INPUT
  nechos = ESTATICSmodel.nechos; % this should be the total number of echos over all contrasts
  [n1, n2, n3] = size(ESTATICSmodel.R2s); % this is the spatial dimension of the data
  nvoxel = n1 * n2 * n3; % this is the total number of voxel
  nvec = ESTATICSmodel.nv; % number of ESTATICS model parameters
  [np1, np2, np3] = deal(2 * patchsize + 1); % spatial dimension of the patches
  hmax = 1.25 ^ (kstar / 3); % maximum spatial bandwidth corresponding to the number of iteration steps kstar in 3D
  lambda = ladjust * 2 * nvec * qf(nv, nechos - nvec) * corr_fac_patchsize(patchsize); % determine the adaptation bandwidth lambda
  wghts = []; % this adjust for non-cubic voxel: if voxel size is 1.2 x 1.2 x 2.4mm wghts should be [1 1 2]
   
    
   
 
   
  %  this is the version with full size invcov (triangular storage)
  %  Uses condensed data (voxel within mask only)
  % create index information for voxel in mask
  nvoxel; % we need the number of voxel within the mask here!
  position = zeros(n1, n2, n3); % this is an array of size dy (spatial size of data)
  position[mask] <- 1:nvoxel 
  % position has the spatial dimensions of the data
  % it is first filled with zeros
  % all the voxels within the brain mask 
  % are assigned numbers 1, 2, 3, ... (number of voxel within the mask) 
  % in this order

  bi = ones(nvoxel);
  theta = modelCoeff;

  k = 1;
  while k <= kstar

    hakt = gethani(1, 1.25 * hmax, 2, 1.25 ^ k, wghts, 1e-4); % This function requires FORTRAN code, take from qMRI package!
    dlw = 2 * floor(hakt ./ [1, wghts]) + 1;

    if k == kstar % use this for the last iteration step
      dim(data) <- c(nechos,nvoxel)
      zobj <- .Fortran(C_pvawsme,
                       as.double(modelCoeff),
                       as.double(mpmData), ## data to smooth additionally
                       as.integer(position),
                       as.integer(nvec),
                       as.integer(nvec * (nvec + 1) / 2),
                       as.integer(nechos), ## leading dimension of data
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       hakt = as.double(hakt),
                       as.double(lambda0),
                       as.double(zobj$theta),
                       as.double(zobj$bi),
                       bi = double(nvoxel), #binn
                       theta = double(nvec * nvoxel),
                       data = double(nechos*nvoxel),
                       as.double(invcov),#
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       double(nechos * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt","data")]
      dim(zobj$data) <- c(nechos, nvoxel)
    else % use this for all but the last iteration step
      zobj <- .Fortran(C_pvaws2,
                       as.double(modelCoeff),
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
                       as.double(invcov),# compact storage
                       as.integer(mc.cores),
                       as.double(spmin),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nvec * mc.cores),
                       as.integer(np1),
                       as.integer(np2),
                       as.integer(np3))[c("bi", "theta", "hakt")]
    end
  
    lambda0 = lambda; % use the adaptation after the first step
    k <- k + 1
  end

  dim(zobj$theta) <- c(nvec, nvoxel)
  if(verbose) cat("\n")
  list(
    theta=zobj$theta,
    hakt=hakt,
    lambda=lambda,
    bi = zobj$bi,
    data= if(!is.null(zobj$data)) zobj$data else NULL
  )

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
end

function qval = qf(df1, df2)
quantile = ones(4, 40);
quantile(1,:)= [647.789011  38.506329  17.443443  12.217863  10.006982   8.813101   8.072669   7.570882   7.209283
6.936728   6.724130   6.553769   6.414254   6.297939   6.199501   6.115127   6.042013   5.978052
5.921631   5.871494   5.826648   5.786299   5.749805   5.716639   5.686366   5.658624   5.633109
5.609564   5.587768   5.567535   5.548702   5.531129   5.514693   5.499288   5.484820   5.471206
5.458372   5.446254   5.434793   5.423937];
quantile(2,:)=[799.500000  39.000000  16.044106  10.649111   8.433621   7.259856   6.541520   6.059467   5.714705
5.456396   5.255889   5.095867   4.965266   4.856698   4.765048   4.686665   4.618874   4.559672
4.507528   4.461255   4.419918   4.382768   4.349202   4.318726   4.290932   4.265483   4.242094
4.220525   4.200572   4.182061   4.164840   4.148779   4.133765   4.119700   4.106496   4.094076
4.082372   4.071326   4.060881   4.050992];
quantile(3,:)= [864.162972  39.165495  15.439182   9.979199   7.763589   6.598799   5.889819   5.415962   5.078119
 4.825621   4.630025   4.474185   4.347178   4.241728   4.152804   4.076823   4.011163   3.953863
3.903428   3.858699   3.818761   3.782886   3.750486   3.721080   3.694273   3.669736   3.647192
3.626408   3.607187   3.589359   3.572778   3.557318   3.542868   3.529334   3.516631   3.504685
3.493429   3.482807   3.472766   3.463260];
quantile(4,:)= [899.583310  39.248418  15.100979   9.604530   7.387886   6.227161   5.522594   5.052632   4.718078
  4.468342   4.275072   4.121209   3.995898   3.891914   3.804271   3.729417   3.664754   3.608344
3.558706   3.514695   3.475408   3.440126   3.408268   3.379359   3.353009   3.328894   3.306741
3.286321   3.267438   3.249925   3.233640   3.218456   3.204267   3.190977   3.178505   3.166777
3.155728   3.145301   3.135445   3.126114];
qval = quantile(df1,df2);
end

