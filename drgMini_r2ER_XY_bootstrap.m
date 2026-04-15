function [r2ER_observed, p_value] = drgMini_r2ER_XY_bootstrap(x_actual, y_actual, x_predicted, y_predicted, n_boot)
    % Calculate observed r2ER
    r2ER_observed = drgMini_r2ER_XY(x_actual,y_actual, x_predicted, y_predicted);
    
    % Generate bootstrap distribution under null hypothesis
    boot_stats = bootstrp(n_boot, @(x_a,y_a,x_p,y_p) drgMini_r2ER_XY(x_a,y_a,x_p(randperm(numel(x_p))),y_p(randperm(numel(y_p)))),...
                         x_actual, y_actual, x_predicted, y_predicted);
    
    % Calculate p-value (two-tailed)
    extreme_count = sum(abs(boot_stats) >= abs(r2ER_observed));
    p_value = extreme_count / n_boot;
end