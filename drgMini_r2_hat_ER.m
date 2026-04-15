function [r2_hat_ER] = drgMini_r2_hat_ER(y_actual, y_predicted)
% Calculate rER², the fraction of the variance explained
%by the model and r_hat_ER2, the unbiased estimator of
%the variance
%
%Popsil and Bair, PLOS Computational Biology, 2021
%https://doi.org/10.1371/journal.pcbi.1009212

% r2er_n2m  Neuron-to-model approx. unbiased estimator of r2er
%
% Parameters
% ----------
% x : [1 x m] or [m x 1] array
%     Model predictions (m observations)
% y : [N x n x m] array
%     Data: N neurons x n repeats x m observations
%
% Returns
% -------
% hat_r2er : estimate of unbiased r^2 between model and expected value
% r2       : biased r^2 estimate
%

x = y_predicted(:)'; % Ensure column vectors
y = y_actual(:)';

% covariance_term = sum((nu - mean(nu)) .* (mu - mean(mu)));
% var_nu = sum((nu - mean(nu)).^2);
% var_mu = sum((mu - mean(mu)).^2);
% rER2 = (covariance_term^2) / (var_nu * var_mu);
%
% function [hat_r2er, r2] = r2er_n2m(x, y)


% Get dimensions: n repeats, m observations
dims = size(y);
n = dims(end-1);
m = dims(end);

% Estimate of trial-to-trial variability (sig2_hat)
sig2_hat = mean(var(y, 0, ndims(y)-1), ndims(y)); % unbiased variance (ddof=1)

% Mean-center the model predictions
x_ms = x - mean(x);

% Average across repeats
y_mean = mean(y, ndims(y)-1);

% Mean-center the data average
y_ms = y_mean - mean(y_mean, ndims(y_mean));

% Dot product of mean-centered model and data (biased numerator)
xy2 = sum(x_ms .* y_ms, ndims(y_ms)).^2;

% Individual variances of model and data
x2 = sum(x_ms.^2);
y2 = sum(y_ms.^2, ndims(y_ms));

% Biased denominator
x2y2 = x2 .* y2;

% Biased r^2
rER2 = xy2 ./ x2y2;

% Unbias numerator and denominator
ub_xy2 = xy2 - sig2_hat / n .* x2;
ub_x2y2 = x2y2 - (m - 1) * sig2_hat / n .* x2;

% Form ratio of individually unbiased estimates (unbiased r^2_er)
r2_hat_ER = ub_xy2 ./ ub_x2y2;

end
