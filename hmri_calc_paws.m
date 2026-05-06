function [smoothedmpmData] = hmri_calc_paws(modelCoeff, mpmData, invCov, mask, nechos, nvec, wghts, kstar, patchsize, ladjust)
  % modelCoeff : nvec x nmask double matrix of model parameters (mask-restricted)
  % mpmData    : nechos x nmask double matrix of raw echo data (mask-restricted)
  % invCov     : nvec*(nvec+1)/2 x nmask double matrix, lower-triangular inverse covariance (mask-restricted)
  % mask       : n1 x n2 x n3 logical array
  % nechos     : total number of echoes across all contrasts
  % nvec       : number of model parameters
  % wghts      : 3-element voxel-size ratio vector [wx, wy, wz], or [] for isotropic
  % kstar      : number of smoothing iterations
  % patchsize  : patch radius (0, 1, or 2)
  % ladjust    : adaptation bandwidth adjustment factor (default 1)

  % CONSTANTS
  spmin  = 0.25;
  lambda0 = 1e32;
  corr_fac_patchsize = [1, 2.77, 3.46];

  % wghts: default to isotropic, then convert 3-element [wx,wy,wz] to
  % 2-element ratio vector [wx/wy, wx/wz] as expected by the C functions
  if isempty(wghts)
    wghts = [1, 1, 1];
  end
  wghts = wghts(1) ./ wghts(2:3);

  % SPATIAL DIMENSIONS from mask
  [n1, n2, n3] = size(mask);

  % patch size in each dimension (1 when that dimension is flat)
  np1 = 2 * patchsize + 1;
  if n2 > 1, np2 = 2 * patchsize + 1; else np2 = 1; end
  if n3 > 1, np3 = 2 * patchsize + 1; else np3 = 1; end

  hmax   = 1.25 ^ (kstar / 3);
  lambda = ladjust * 2 * nvec * qf(nvec, nechos - nvec) * corr_fac_patchsize(patchsize + 1);

  % position array: 0 outside mask, 1..nmask inside mask (int32 for MEX)
  position = zeros(n1, n2, n3, 'int32');
  position(mask > 0) = int32(1:nnz(mask));

  % initialise adaptation weights and current estimate
  bi    = ones(1, nnz(mask));
  theta = modelCoeff;

  % SMOOTHING ITERATION
  k = 1;
  while k <= kstar

    hakt    = gethani(1, 1.25 * hmax, 2, 1.25 ^ k, wghts, 1e-4);
    dlw     = 2 * floor(hakt ./ [1, wghts]) + 1;
    lwghtSz = prod(dlw);  % scalar size of the lwght scratch buffer

    if k == kstar && ~isempty(mpmData)

      [bi, theta, smoothedmpmData] = pvawslast(modelCoeff, ...
                                               mpmData, ...
                                               position, ...
                                               nvec, ...
                                               nvec * (nvec + 1) / 2, ...
                                               nechos, ...
                                               n1, n2, n3, ...
                                               hakt, lambda0, ...
                                               theta, bi, invCov, ...
                                               1, ...
                                               spmin, wghts, lwghtSz, ...
                                               np1, np2, np3);

    else

      [bi, theta] = pvaws(modelCoeff, ...
                          position, ...
                          nvec, ...
                          nvec * (nvec + 1) / 2, ...
                          n1, n2, n3, ...
                          hakt, lambda0, ...
                          theta, bi, invCov, ...
                          1, ...
                          spmin, wghts, lwghtSz, ...
                          np1, np2, np3);

    end

    lambda0 = lambda;
    k = k + 1;
  end

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
