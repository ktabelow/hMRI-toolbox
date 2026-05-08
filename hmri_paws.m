function [denoised_weighted_data] = hmri_paws(weighted_data, params)
  %
  % Adaptive smoothing of MPM weighted data using PAWS.
  %
  % Mirrors the two-step R pipeline:
  %   Step 1 — estimateESTATICS : fit the ESTATICS signal model and derive
  %             per-voxel model coefficients and inverse covariance matrix.
  %   Step 2 — smoothESTATICS   : optional spatial stabilisation of invCov,
  %             then adaptive smoothing via hmri_calc_paws (= vpawscov2).
  %
  % R uses non-linear least squares (NLS); here hmri_calc_R2s provides a
  % log-linear OLS approximation. modelCoeff = [R2s, A_T1, A_MT, A_PD] is
  % in signal space, matching R's NLS convention. invCov is built from the
  % NLS Jacobian J'J/sigma2_signal per voxel, which is consistent with
  % signal-space parameters.

  % ── Parameters ────────────────────────────────────────────────────────────
  kstar      = params.kstar;
  patchsize  = params.patchsize;
  ladjust    = params.ladjust;
  mask       = logical(params.mask{1});
  wghts      = params.voxel_size;
  if isfield(params, 'alpha'), alpha = params.alpha; else, alpha = 0.025; end
  fit_method = 'OLS';
  mscbw      = 5;   % bandwidth for median smoothing of residual variance
                    % (mirrors mscbw in smoothESTATICS; set 0 to disable)

  % ══ Step 1: mimic estimateESTATICS ════════════════════════════════════════
  %
  % estimateESTATICS fits S_k = A_k * exp(-R2* * TE) via NLS per voxel and
  % collects modelCoeff [nv x nvoxel], invCov [nv x nv x nvoxel], rsigma.
  % We use hmri_calc_R2s (log-linear OLS) as an approximation, then build
  % modelCoeff and a signal-space invCov from its output.

  [R2s, extrapolated, SError] = hmri_calc_R2s(weighted_data, fit_method);

  % Derived dimensions
  numcon       = numel(weighted_data);      % number of contrasts
  nvec         = numcon + 1;                % R2s + one intercept per contrast
  [n1, n2, n3] = size(R2s);
  nechos       = sum(arrayfun(@(w) numel(w.TE), weighted_data));

  % modelCoeff [nvec x nmask]: R2s in row 1, signal at TE=0 in rows 2..nvec.
  modelCoeff      = zeros(nvec, nnz(mask));
  modelCoeff(1,:) = R2s(mask);
  for w = 1:numcon
    modelCoeff(w+1,:) = extrapolated{w}(mask);
  end

  % ══ Step 2: mimic smoothESTATICS ══════════════════════════════════════════
  %
  % Build invCov from the NLS signal-model Jacobian, consistent with
  % signal-space modelCoeff = [R2s, A_1, ..., A_ncon].
  %
  % Jacobian of S_k = A_{c(k)} * exp(-R2s * TE_k) w.r.t. [R2s, A_1,...]:
  %   J[k, 1]   = -A_{c(k)} * TE_k * exp(-R2s * TE_k)
  %   J[k, w+1] =  exp(-R2s * TE_k)  for echo k in contrast w, else 0
  %
  % invCov(v) = J(v)' J(v) / sigma2_signal(v), stored in lower-triangular
  % row-major order matching KLdistsi (element order: (1,1),(2,1),(2,2),...).

  % Signal-space residual variance per voxel [n1 x n2 x n3].
  % SError.weighted{w} = rms of signal residuals for contrast w (per voxel).
  sigma2_NLS = zeros(n1, n2, n3);
  for w = 1:numcon
    nTE_w      = numel(weighted_data(w).TE);
    sigma2_NLS = sigma2_NLS + nTE_w * SError.weighted{w}.^2;
  end
  sigma2_NLS = sigma2_NLS / (nechos - nvec);

  % Optional spatial stabilisation of sigma2 (mirrors rsdx/rsdhat in R).
  if mscbw > 0 && exist('medianFilter3D', 'file')
    sigma2_s        = medianFilter3D(sigma2_NLS, mscbw, mask);
    sigma2_s(~mask) = mean(sigma2_s(mask));
  else
    sigma2_s = sigma2_NLS;
  end

  % Compute J'J blocks vectorised over masked voxels.
  nmask  = nnz(mask);
  R2s_v  = R2s(mask)';          % [1 x nmask]
  JtJ_11 = zeros(1, nmask);     % (R2s,  R2s ) block
  JtJ_w1 = zeros(numcon, nmask);% (A_w,  R2s ) blocks
  JtJ_ww = zeros(numcon, nmask);% (A_w,  A_w ) diagonal blocks

  for w = 1:numcon
    TEs_w  = weighted_data(w).TE(:);           % [nTE_w x 1]
    A_w_v  = extrapolated{w}(mask)';           % [1 x nmask]
    exp2   = exp(-2 * TEs_w .* R2s_v);        % [nTE_w x nmask]
    JtJ_11    = JtJ_11 + A_w_v.^2 .* sum(TEs_w.^2 .* exp2, 1);
    JtJ_w1(w,:) =      - A_w_v    .* sum(TEs_w    .* exp2, 1);
    JtJ_ww(w,:) =                    sum(             exp2, 1);
  end

  sigma2_v = sigma2_s(mask)';   % [1 x nmask]

  % Assemble invCov [nvd x nmask] in KLdistsi lower-triangular row-major order:
  % (1,1),(2,1),(2,2),(3,1),(3,2),(3,3),... Cross-contrast off-diagonals are 0.
  nvd    = nvec * (nvec + 1) / 2;
  invCov = zeros(nvd, nmask);
  idx    = 0;
  for row = 1:nvec
    for col = 1:row
      idx = idx + 1;
      if row == 1 && col == 1
        invCov(idx,:) = JtJ_11 ./ sigma2_v;
      elseif col == 1
        invCov(idx,:) = JtJ_w1(row-1,:) ./ sigma2_v;
      elseif row == col
        invCov(idx,:) = JtJ_ww(row-1,:) ./ sigma2_v;
      end
    end
  end

  % mpmData [nechos x nmask] — raw echo data restricted to mask.
  mpmData  = zeros(nechos, nnz(mask));
  echo_idx = 1;
  for w = 1:numcon
    nTE    = numel(weighted_data(w).TE);
    data_w = reshape(weighted_data(w).data, [], nTE);   % [nvoxel x nTE]
    mpmData(echo_idx:echo_idx+nTE-1, :) = data_w(mask(:), :)';
    echo_idx = echo_idx + nTE;
  end

  % PAWS smoothing — mirrors vpawscov2 call in smoothESTATICS.
  smoothedmpmData = hmri_calc_paws(modelCoeff, mpmData, invCov, mask, ...
                                   nechos, nvec, wghts, kstar, patchsize, ladjust, alpha);

  % ── Reconstruct denoised weighted_data ────────────────────────────────────
  denoised_weighted_data = weighted_data;
  echo_idx = 1;
  for w = 1:numcon
    nTE      = numel(weighted_data(w).TE);
    data_out = zeros(n1, n2, n3, nTE);
    for e = 1:nTE
      vol       = zeros(n1, n2, n3);
      vol(mask) = smoothedmpmData(echo_idx + e - 1, :);
      data_out(:,:,:,e) = vol;
    end
    denoised_weighted_data(w).data = data_out;
    echo_idx = echo_idx + nTE;
  end

end
