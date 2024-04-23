function [R2] = drgGetR2(y,y_predicted)
%Calculates R2 as done in the Kording discriminant analysis
sum_difsq=0;
sum_difmeansq=0;
ymean=mean(y);
for ii=1:length(y)
    sum_difsq=sum_difsq+(y_predicted(ii)-y(ii))^2;
    sum_difmeansq=sum_difmeansq+(y(ii)-ymean)^2;
end
R2=1-(sum_difsq/sum_difmeansq);
end