function V = vec_to_sym(v,n)
%symmetric matrix from minimal encoding vector
%V: n x n x N stack of symmetric matrices
%v: N x n(n+1)/2 stack of row vectors

%pre-allocate matrix
N = size(v,1);
V = repmat(v(1),n,n,N);

%reshape v
v = reshape(v', 1,size(v,2),N);
vp = permute(v,[2,1,3]);

%top-left corner
V(1,1,:) = v(1,1,:);

%indices
i_end = cumsum(1:n);

for j = 2:n
    V(j,1:j,:) = v(1,i_end(j-1)+1:i_end(j),:); %row
    V(1:j,j,:) = vp(i_end(j-1)+1:i_end(j),1,:); %column
end


end