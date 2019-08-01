function c = MPC_cost(xu, Prob, F, Tp)

xu_nom = Prob.user.xu_nom;

c = (Tp/2)*(xu-xu_nom)'*F*(xu-xu_nom);

end