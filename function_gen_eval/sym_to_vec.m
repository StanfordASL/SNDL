function v = sym_to_vec(V,n)
%minimal vector encoding for symmetric matrix
%V: n x n symmetric matrix
%v: n(n+1)/2 length row vector

%pre-allocate vector
v = repmat(V(1,1),1,n*(n+1)/2);
v(1,1) = V(1,1);

%indices
i_end = cumsum(1:n);

for i = 2:n
    v(i_end(i-1)+1:i_end(i)) = V(i,1:i);
end


end