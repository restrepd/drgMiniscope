function [rER2] = drgMini_r2ER(y_actual, y_predicted)
    % Calculate rER² (Normalized Signal Power Explained)
    %this is a method used to quantify goodness of fit
    %that, unlike R2, is not biased by trial to trial variability
    %Popsil and Bair, PLOS Computational Biology, 2021
    %https://doi.org/10.1371/journal.pcbi.1009212

    nu = y_predicted(:); % Ensure column vectors
    mu = y_actual(:);
    
    covariance_term = sum((nu - mean(nu)) .* (mu - mean(mu)));
    var_nu = sum((nu - mean(nu)).^2);
    var_mu = sum((mu - mean(mu)).^2);
    rER2 = (covariance_term^2) / (var_nu * var_mu);

    
end