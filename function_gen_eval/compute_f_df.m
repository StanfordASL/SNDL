function [f, dfdx_p] = compute_f_df(alpha,Zf,Jf,n,m,Nc)

%compute f(Xc) and dfdx_p(Xc) (i.e., first (n-m) rows of dfdx)

f = Zf*alpha;

dfdx_p = zeros(n-m,n,Nc);

for k = 1:Nc
    for j = 1:n-m
        dfdx_p(j,:,k) = alpha'*Jf(:,:,k,j);
    end
end

end