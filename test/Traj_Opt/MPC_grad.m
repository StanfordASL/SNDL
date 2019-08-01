function J = MPC_grad(xu, Prob, F, Tp)

xu_nom = Prob.user.xu_nom;

J = Tp*F*(xu-xu_nom);

end