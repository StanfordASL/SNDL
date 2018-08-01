function [lambda_act,d] = check_lambda(n,m,X,X_dot,f_h,B_h,df_h,W_h,dW_h,E,E_prev,d_prev,lambda_prev,lambda_nom,dt,E_thresh)

%X: [x_1,...,x_Ke]
%X_dot: [dx_1, ..., dx_Ke]

Ke = size(X,2);

W = W_h(X')+repmat(0.1*eye(n),1,1,Ke);
dW = dW_h(X');
lambda_lim = inf*ones(Ke,1);

%% Evaluate constraint along geodesic

for k = 1:Ke
        
   eta_x = W(:,:,k)\X_dot(:,k);
   
   %check if eta_x in B_perp
   in_perp = norm(eta_x'*B_h(X(:,k)))<=0.001;
   
   if (~in_perp) %lambda constraint wont' be from here
       continue;
   else
       f = f_h(X(:,k));
       dfdx = df_h(X(:,k));
       
       dW_f = zeros(n);
       for i = 1:n
           dW_f = dW_f + dW(:,:,k,i)*f(i);
       end
       
       a = eta_x'*(-dW_f + dfdx*W(:,:,k) + W(:,:,k)*dfdx')*eta_x;
       b = eta_x'*W(:,:,k)*eta_x; % always >=0
       
       if (b>0.001)
           lambda_lim(k) = a/(-2*b);
       end
   end   
end

%% First update estimate of disturbance/mismatch

if (lambda_prev > 0)
    E_hat = E_prev*exp(-2*lambda_prev*dt) + (d_prev/(2*lambda_prev))*(1-exp(-2*lambda_prev*dt));
    d = d_prev +  2*lambda_prev*(E-E_hat)/(1-exp(-2*lambda_prev*dt));
else
    d = 0;
end

%% Now check for most limiting lambda along geodesic

lambda_allow = min(lambda_lim);

%pick lambda sensibly
if isinf(lambda_allow) %no limiting points here
    
    %smoothly choose desired lambda
    lambda_des = max(lambda_nom*(E/E_thresh),10);
    
    %E_des
    E_des = E*exp(-2*lambda_des*dt);
    
    fun = @(l) (E - (d/(2*l)))*exp(-2*l*dt) + (d/(2*l)) - E_des;
    [l,~,flag] = fzero(fun,max(lambda_des,lambda_prev));
    
    if (flag >0) && (l>0)
        lambda_act = l;
    else
        lambda_act = lambda_des;
    end    
    
%     fprintf('E: %.4f, d: %.4f, l_des: %.2f, l:%.2f, E_des: %.4f \n',E,d,lambda_des,lambda_max,E_des);

else
    lambda_act = lambda_allow;
end


end