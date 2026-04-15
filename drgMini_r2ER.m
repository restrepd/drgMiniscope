function [rER2] = drgMini_r2ER(y_actual, y_predicted)
    % Calculate rER², the fraction of the variance explained
    %by the model
    %
    %Popsil and Bair, PLOS Computational Biology, 2021
    %https://doi.org/10.1371/journal.pcbi.1009212
    %
    %Note that this is not r_carrot_ER2, the unbiased estimator of
    %the variance

    nu = y_predicted(:); % Ensure column vectors
    mu = y_actual(:);
    
    covariance_term = sum((nu - mean(nu)) .* (mu - mean(mu)));
    var_nu = sum((nu - mean(nu)).^2);
    var_mu = sum((mu - mean(mu)).^2);
    rER2 = (covariance_term^2) / (var_nu * var_mu);

    
end