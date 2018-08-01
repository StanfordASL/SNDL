function [f, B, df, dB] =  construct_dyn(om_f, alpha, om_b, beta,O,Lf, Lb,n,m)


phi_f_T = @(x) (1/sqrt(O))*( kron(kron(cos(x'*om_f'),[1,0])+kron(sin(x'*om_f'),[0,1]),Lf') );

phi_b_T = @(x) (1/sqrt(O))*( kron(kron(cos(x'*om_b'),[1,0])+kron(sin(x'*om_b'),[0,1]),Lb') );

f = @(x) phi_f_T(x)*alpha;
B = @(x) phi_b_T(x)*beta;

Jphi_j_f = @(x) (1/sqrt(O))*(kron(-diag(sin(om_f*x))*om_f,[1;0]) + kron(diag(cos(om_f*x))*om_f,[0;1]));
Jphi_j_b = @(x) (1/sqrt(O))*(kron(-diag(sin(om_b*x))*om_b,[1;0]) + kron(diag(cos(om_b*x))*om_b,[0;1]));

df = @Jacobian_f;

    function dfdx = Jacobian_f(x)
        dfdx = zeros(n,n);
        for i = 1:n
            dfdx(i,:) = alpha'*kron(Jphi_j_f(x),Lf(:,i));
        end
    end

dB = @Jacobian_B;

    function dBdx = Jacobian_B(x)
        dBdx = zeros(n,n,m); 
        for j = 1:m
            for i = 1:n
                dBdx(i,:,j) = beta(:,j)'*kron(Jphi_j_b(x),Lb(:,i));
            end
        end
    end

end