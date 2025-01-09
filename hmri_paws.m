function [denoised_weighted_data] = hmri_paws(weighted_data, params)
    %main function that does PAWS

    mscbw = 5; % bandwidth to smooth the inverse covariance matrix 

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
  
%%inputs
%weighted_data (cell): cell array of available contrasts
%params (str): 3 adaptive denoising parameters kstar = 16, lambda,
%patchsize = [0, 1, 2] make sure in the GUI (default = 1)

%%outputs:
%denoised_weighted_data: adaptive denoised volumes

%populate parameters from arguments

%call hmri_coreg before denoising

%call hmri_calc_R2s with 4 arguments (including design matrix D and SError)
[R2s, extrapolated, SError, D] = hmri_calc_R2s(dataToFit,fit_method);

%from R2s and extrapolates calculate residuals and variance (using design
%matrix)

res_var = 0;
num_echo=0;
for ccon =1:numconn
    for te = data[ccon].TE
        res_var = res_var + (data[cccon].te -  extrapolated{ccon}.*exp(-R2s*te))^2;
    end
    num_echo = num_echo + length(data[ccon].TE);
end
res_var = res_var / (num_echo - numconn + 1);
variance = ones([4 , 4, size(R2s)]);
model_var = inv(D'D);
for i =1:(numconn + 1)
    for j= 1:(numconn +1)
        variance(i, j, :, :,:) = res_var * model_var(i, j);
    end
end
% here comes the reduction to the voxel within the mask
% here comes the reduction for the symmetric part

%call PAWS with data, residuals, params (include here?)

ESTATICSmodel.extrapolated = extrapolated;
ESTATICSmodel.R2s = R2s;
ESTATICSmodel.variance = variance;
ESTATICSmodel.nv = numconn + 1;
ESTATICSmodel.nechos = num_echo;
mask = [];

% CREATE MASK DEFAULT if not given
if isempty(mask)
    mask = ones(n1, n2, n3);
  end

  % extract array nv x nv x nvoxel
  invCov = ESTATICSmodel.variance

  % start smoothing the covariance martix to stabilise
  % we might skip this part because we use the linear model
  % if we use it this should to hmri_paws
  if mscbw > 0
    rsdx <- extract(mpmESTATICSModel,"rsigma")^2
    rsdhat <- medianFilter3D(rsdx,mscbw,mask)
    rsdhat[!mask] <- mean(rsdhat[mask])
    invCov <- sweep(invCov,3:5,rsdx/rsdhat,"*")
  end
  dim(invCov) <- c(nv, nv, prod(sdim))

  % projecting out the voxel within the mask only to save memory
  % mask is TRUE/FALSE
  invCov_d = invConv()
  % Baris mask the R2s directly in hmri_paws.
  invCov <- invCov[,,mask]

  
  % we expect modelCoeff to be a nv x nvoxel_within_mask 
  % we expect invCov to be nv x nv x nvoxel_within_mask


  % the next switch can only be executed if nvec < 5, pls CHECK
  % this has to be done in hmri_paws alreadz when the variance array is
  % created
  switch nvec
    case 1
      indcov = [1];
    case 2
      indcov = [1 2 4];
    case 3
      indcov = [1 2 6 3 7 11 4 8 12 16];
    case 4
      indcov = [1 2 7 3 8 13 4 9 14 19 5 10 15 20 25];
  end
  dim(invcov) = c(nvec * nvec, nvoxel)
  invcov = invcov(indcov, :);
  % end of " this has to de done in ..."


% add ladjust = 1 to the list of defaults
hmri_calc_paws(ESTATICSmodel, dataToFit,mask,  kstar = 16, patchsize = 1, ladjust)

%take out denoised_weighted_data
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end