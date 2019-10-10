function [d_P,A_inv,Q_full] = get_dP(d_F,s_F,s_P,N_seg,T_seg,p_order,Q,pos_coef,vel_coef,acc_coef,jer_coef,snap_coef)
Q_full = Q(T_seg(1));
for ns = 2:N_seg
    Q_full = blkdiag(Q_full,Q(T_seg(ns)));
end
A = get_constraints(N_seg,T_seg,p_order,pos_coef,vel_coef,acc_coef,jer_coef,snap_coef);
A_inv = A\eye(p_order*N_seg);
R = A_inv'*Q_full*A_inv;
R_PP = R(s_F+1:end,s_F+1:end);
R_FP = R(1:s_F,s_F+1:end);
d_P = -R_PP\(R_FP'*d_F);

end