function F = theta_solve(theta, betta, s, kappa, chi, eta, gama, mpl, b)
    LHS	= (1/betta-1+s)*kappa/chi*theta^(1-eta);
    RHS	= gama*(mpl-b) - (1-gama)*kappa*theta;
    F	= LHS-RHS;
end