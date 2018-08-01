function[Xc_i_up,constraint_viol,flag_done,max_slack] = update_constraint_set(eps_l,delta_wl,eps_wl,tol,W_p,W,partial_Wp,f,dfdx_p,lambda,N,Nc,n,m,Xc_i)

%sort points in Xc based on constraint violation
%eliminate points with eig(F), eig(w_lower*I - W) <= -tol
%generate replacements from X \ Xc

%% Check constraint violation

dWp_f = zeros((n-m),(n-m),N);
for j = 1:n-m
    dWp_f = dWp_f + partial_Wp(:,:,:,j).*reshape(f(j:n:end),1,1,N);
end

%there are two LMIs we care about
lmi_viol = zeros(N,2);

for k = 1:N
   
    F = -dWp_f(:,:,k) + dfdx_p(:,:,k)*W(:,1:n-m,k)+...
                W(1:n-m,:,k)*dfdx_p(:,:,k)' +2*(lambda+eps_l)*W_p(:,:,k);
    
    lmi_viol(k,1) = max(eig(F));  %want <0
    lmi_viol(k,2) = max(eig((delta_wl+eps_wl)*eye(n)-W(:,:,k))); %want <0
end

fprintf('CCM viol: %.3f, w_lower viol: %.3f, ', max(lmi_viol(:,1)),max(lmi_viol(:,2)));

max_slack = max(lmi_viol(:,1));

%take the max of the two violations for each datapoint
constraint_viol = max(lmi_viol,[],2);

%find top constraint violators
[sorted_all,idx_all] = sort(constraint_viol,'descend');

%idx_all(1:n_violating) are violations
n_violating = find(sorted_all<=1e-3,1)-1;
if isempty(n_violating)
    n_violating = N;
end
fprintf('n_viol:%d, ',n_violating);

%% Update Xc

%start:
[sorted,idx] = sort(constraint_viol(Xc_i),'descend');

%idx(1:n_keep) are above tolerance
n_keep = find(sorted <= -tol)-1;

n_disc = 0; 
if isempty(n_keep)
    %none below tolerance
    Xc_i_up = Xc_i;
    n_add_max = 50; %add 50 worst violators
else
    %discard some points
    Xc_i_up = Xc_i(idx(1:n_keep));
    n_disc = length(Xc_i)-length(Xc_i_up); 
    n_add_max = n_disc + 50; %replace the discarded and add up to 50
end
fprintf('n_disc:%d, ',n_disc);

%now add (at most n_disc) violating points in X \ Xc
n_added = 0; i = 1;
while (i <= n_violating) && (n_added < n_add_max)
    if (~ismember(idx_all(i),Xc_i))
        Xc_i_up = [Xc_i_up;idx_all(i)];
        n_added = n_added + 1;
    end
    i = i+1;
end

fprintf('Nc = %d \n', length(Xc_i_up));

%% set flag for all constraints satisfied

flag_done = n_violating==0;


end
