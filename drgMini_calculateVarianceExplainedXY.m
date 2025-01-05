function percentExplained = drgMini_calculateVarianceExplainedXY(x,y,xfit,yfit)
% This method calculates the coefficient of determination (R-squared) and converts it to a percentage.
% The result tells you what percent of the variance in y is explained by the fit yfit

% Calculate the total sum of squares (SST):
SST = sum((x - mean(x)).^2)+sum((y - mean(y)).^2);

% Calculate the residual sum of squares (SSR)
SSR = sum((x - xfit).^2)+sum((y - yfit).^2);

% Calculate the coefficient of determination (R-squared)
R2 = 1 - (SSR / SST);

%Multiply by 100
percentExplained = R2;
end