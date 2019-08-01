function[Xc_i_up,constraint_viol,flag_done,max_slack] = update_constraint_set(constants,W_p,W,partial_Wp,f,dfdx_p,lambda,Xc_i)

%sort points in Xc based on constraint violation
%eliminate points with eig(F), eig(w_lower*I - W) <= -tol
%generate replacements from X \ Xc

%% extract constants

n = constants.n;
m = constants.m;
eps_l = constants.eps_l;
delta_wl = constants.delta_wl;
eps_wl = constants.eps_wl;
tol = constants.tol;
eps_term = constants.eps_term;
N = constants.N;


%% Compute constraint violation

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

fprintf('Full Set: CCM viol: %.3f, w_lower viol: %.3f, ', max(lmi_viol(:,1)),max(lmi_viol(:,2)));

max_slack = max(lmi_viol(:,1));

%take the max of the two violations for each datapoint
constraint_viol = max(lmi_viol,[],2);

%find top constraint violators
[sorted_all,idx_all] = sort(constraint_viol,'descend');

%% Check for termination

%idx_all(1:n_term) are above termination tolerance
n_term = find(sorted_all<=eps_term,1)-1;
if isempty(n_term)
    n_term = N;
end
fprintf('n_viol:%d, ',n_term);

if n_term == 0
    Xc_i_up = [];
    flag_done = 1;
    return;
else
    flag_done = 0;
    %get actual number of violations
    %idx_all(1:n_violating) are above 0
    n_violating = find(sorted_all<=0,1)-1;
    if isempty(n_violating)
        n_violating = N;
    end
end    

%% Update Xc

%start:
[sorted,idx] = sort(constraint_viol(Xc_i),'descend');

%retain idx(1:n_keep) in Xc_i
n_keep = find(sorted <= -tol,1)-1;

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

%now add (at most n_add_max) violating points in X \ Xc
n_added = 0; i = 1;
while (i <= n_violating) && (n_added < n_add_max)
    if (~ismember(idx_all(i),Xc_i))
        Xc_i_up = [Xc_i_up;idx_all(i)];
        n_added = n_added + 1;
    end
    i = i+1;
end

fprintf('Nc = %d \n', length(Xc_i_up));

end