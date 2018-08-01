function f_all = NMPC_dyn(f,x,n,N)

f_all = zeros((N+1)*n,1);
for j = 1:N+1
    f_all(1+(j-1)*n:j*n) = f(x(1+(j-1)*n:j*n));
end
    
    

end