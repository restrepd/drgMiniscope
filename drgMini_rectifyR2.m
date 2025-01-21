function R2rect = drgMini_rectifyR2(R2)
% This method provides a rectificaiton of R2 values that can become quite negative. 
% R2 calculates the coefficient of determination (R-squared).
% The result tells you what percent of the variance in y is explained by the fit yfit
% This is the same as R2 used by Glasser et al eNEURO 2020 
R2rect=R2;
R2rect(R2<0)=-(-R2(R2<0)/(0.5-R2(R2<0)));

end