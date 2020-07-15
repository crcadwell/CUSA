function [chi2stat,p] = ChiSquared(x)
total = sum(sum(x));
totalj = sum(x,1);
totali = sum(x,2);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        e(i,j) = totali(i)/total*totalj(j);
        chi2(i,j) = (x(i,j)-e(i,j))^2/e(i,j);
    end
end
chi2stat = sum(sum(chi2));
p = 1 - (chi2cdf(chi2stat,(size(x,1)-1)*(size(x,2)-1)));
end