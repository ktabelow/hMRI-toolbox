function [denoised_weighted_data] = hmri_paws(weighted_data, params)
%main function that does PAWS

%%inputs
%weighted_data (cell): cell array of available contrasts
%params (str): 3 adaptive denoising parameters kstar, lambda, patchsize

%%outputs:
%denoised_weighted_data: adaptive denoised volumes

%call hmri_coreg before denoising

%call hmri_calc_R2s with 4 arguments (including design matrix D and SError)

%from R2s and extrapolates calculate residuals and variance (using design
%matrix)

%call PAWS with data, residuals, params (include here?)

%take out denoised_weighted_data
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end