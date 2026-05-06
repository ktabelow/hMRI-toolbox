function [denoised_weighted_data] = hmri_paws(weighted_data, params)
  %main function that does PAWS

  %%%%%%%%% define constants for manual testing (move the constants out and delete after testing)
  %paws params
  kstar = params.kstar;
  patchsize = params.patchsize;
  ladjust = params.ladjust;
  mask = params.mask{1}; %take as TRUE for testing

  %hmri params
  fit_method = 'OLS';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % DEFINE ALL CONSTANTS
  mscbw = 5; % bandwidth to smooth the inverse covariance matrix 
      
  %%inputs
  %weighted_data (cell): cell array of available contrasts
  %params (str): 3 adaptive denoising parameters kstar = 16, lambda,
  %patchsize = [0, 1, 2] make sure in the GUI (default = 1)

  %%outputs:
  %denoised_weighted_data: adaptive denoised volumes

  % populate parameters from arguments

  %call hmri_coreg before denoising

  %call hmri_calc_R2s with 4 arguments (including design matrix D and SError)
  [R2s, extrapolated, SError, D] = hmri_calc_R2s(dataToFit,fit_method);

  %from R2s and extrapolates calculate residuals and variance (using design
  %matrix)

  res_var = 0;
  num_echo=0;
  for ccon =1:numconn
      for te = data{ccon}.TE
          res_var = res_var + (data{ccon}.te - extrapolated{ccon} .* exp(-R2s * te)).^2;
      end
      num_echo = num_echo + length(data{ccon}.TE);
  end
  res_var = res_var / (num_echo - numconn + 1);
  variance = ones([4 , 4, size(R2s)]);
  model_var = inv(D' * D);
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
  mask = [];

  % CREATE MASK DEFAULT if not given
  if isempty(mask)
      mask = ones(n1, n2, n3);
  end
  
  wghts = []; % this adjust for non-cubic voxel: if voxel size is 1.2 x 1.2 x 2.4mm wghts should be [1 1 2]
 
  % extract array nv x nv x nvoxel
  invCov = ESTATICSmodel.variance

  % start smoothing the covariance martix to stabilise
  % we might skip this part because we use the linear model
  % if we use it this should to hmri_paws
  if mscbw > 0
    rsdx   = rsigma .^ 2;                          % residual variance per voxel
    rsdhat = medianFilter3D(rsdx, mscbw, mask);     % spatially smoothed version
    rsdhat(~mask) = mean(rsdhat(mask));
    invCov = invCov .* (rsdx ./ rsdhat);            % rescale inv covariance
  end
  invCov = reshape(invCov, nv * nv, prod(size(mask)));

  % keep only voxels within mask
  invCov = invCov(:, mask(:));

  
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
      indcov = [1 2 5 3 6 9];
    case 4
      indcov = [1 2 6 3 7 11 4 8 12 16];
  end
  invCov = reshape(invCov, nvec * nvec, nnz(mask));
  invCov = invCov(indcov, :);
  % end of " this has to de done in ..."

  % CONSISTENCY CHECKS (TBD somewhere) TO CHECK
  % size(mask) = [n1, n2, n3]
  % size(modelCoeff) = [nvec, n1, n2, n3]
  % size(mpmData) = [nechos, n1, n2, n3]
  % mpmData = mpmData[, mask];

  % add ladjust = 1 to the list of defaults
  smoothedmpmData = hmri_calc_paws(modelCoeff, dataToFit, invCov, mask, num_echo, numconn + 1, wghts, kstar, patchsize, ladjust);

  %take out denoised_weighted_data
  outputArg1 = inputArg1;
  outputArg2 = inputArg2;
end