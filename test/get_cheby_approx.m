function c = get_cheby_approx(fnc,intvl, n,m)

bma=0.5*(intvl(2)-intvl(1));
bpa=0.5*(intvl(2)+intvl(1));
f = zeros(n,1);

for k = 1:n
    y=cos(pi*(k-0.5)/n); 
    f(k)=fnc(y*bma+bpa);
end

fac=2.0/n;
c = zeros(1,m+1);
for j = 0:m 
    sum_temp = f.*cos(pi*j*([1:n]'-0.5)/n); 
    c(j+1)=fac*sum(sum_temp);
end

c(1) = 0.5*c(1);

end