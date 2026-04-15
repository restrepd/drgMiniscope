function [rER2] = drgMini_r2ER_XY(x_actual, y_actual, x_predicted, y_predicted)
    % Ensure column vectors
    x_actual = x_actual(:);
    y_actual = y_actual(:);
    x_predicted = x_predicted(:);
    y_predicted = y_predicted(:);

    % Stack actual and predicted data
    mu = [x_actual, y_actual];         % Actual data (N x 2)
    nu = [x_predicted, y_predicted];   % Predicted data (N x 2)

    % Calculate means
    mu_mean = mean(mu, 1);
    nu_mean = mean(nu, 1);

    % Center data
    mu_centered = mu - mu_mean;
    nu_centered = nu - nu_mean;

    % Covariance matrices
    covariance_matrix = (nu_centered' * mu_centered) / length(mu_centered);
    covariance_term = trace(covariance_matrix * covariance_matrix');
    
    % Variances
    var_nu = trace((nu_centered' * nu_centered) / length(nu_centered));
    var_mu = trace((mu_centered' * mu_centered) / length(mu_centered));
    
    % Calculate rER^2 (normalized signal power explained)
    rER2 = covariance_term^2 / (var_nu * var_mu);
end