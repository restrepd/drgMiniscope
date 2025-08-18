function [r2ER_observed, p_value] = drgMini_r2ER_bootstrap(y_actual, y_predicted, n_boot)
    % Calculate observed r2ER
    r2ER_observed = drgMini_r2ER(y_actual, y_predicted);
    
    % Generate bootstrap distribution under null hypothesis
    boot_stats = bootstrp(n_boot, @(y_a,y_p) drgMini_r2ER(y_a,y_p(randperm(numel(y_p)))),...
                         y_actual, y_predicted);
    
    % Calculate p-value (two-tailed)
    extreme_count = sum(abs(boot_stats) >= abs(r2ER_observed));
    p_value = extreme_count / n_boot;
end