function [denoised_weighted_data] = hmri_paws(weighted_data, params)
%main function that does PAWS

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
        res_var = res_var + (data[cccon].te -  extrapolated{ccon}*exp(-R2s*te))^2;
    end
    num_echo = num_echo + length(data[ccon].TE);
end
res_var = res_var / (num_echo - numconn + 1);

variance = res_var * inv(D'D);

%call PAWS with data, residuals, params (include here?)

ESTATICSmodel.extrapolated = extrapolated;
ESTATICSmodel.R2s = R2s;
ESTATICSmodel.nv = numconn + 1;
ESTATICSmodel.necho = num_echo;
mask = [];
% add ladjust = 1 to the list of defaults
hmri_calc_paws(ESTATICSmodel, dataToFit,mask,  kstar = 16, patchsize = 1, ladjust)

%take out denoised_weighted_data
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end