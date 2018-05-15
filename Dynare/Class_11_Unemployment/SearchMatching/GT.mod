% Replication file for Gertler & Trigari (2009, JPE)

% 0 - flexible wages; 1 - lambda = 8/9; 2 - lambda = 11/12
@#define flex = 0

% 0 - nice irfs; 1 - nice moments
@#define moments = 0

close all;

var		y c i k n r w z a
		u v m p q x Lambda chi epsilon mu w_o ls theta;

varexo	epsilon_z;

parameters
% Structural parameters
	beta delta alpha rho_z sigma_z rho sigma eta sigma_m kappa b lambda
% Steady state values
	c_y_ss i_y_ss y_ss k_ss n_ss r_ss w_ss z_ss a_ss
	u_ss p_ss x_ss chi_ss epsilon_ss mu_ss ls_ss H_ss
% Derived parameters
	aleph aleph_a aleph_w varphi_a varphi_x varphi_p varphi_chi psi tau tau_1 tau_2 varsigma phi gamma_b gamma_o gamma_f
	;

% TABLE 1 Values of Parameters (p. 59)
	beta	= 0.997;	% Discount factor
	delta	= 0.008;	% Capital depreciation rate
	alpha	= 0.33;		% Production function parameter
	rho_z	= 0.983;	% Technology autoregressive parameter
	sigma_z	= 0.0075;	% Technology standard deviation
	rho		= 0.965;	% Survival rate
	sigma	= 0.5;		% Elasticity of matches to unemployment
	eta		= 0.5;		% Bargaining power parameter
	sigma_m	= 1;		% Matching function constant
	kappa	= 148.2;	% Adjustment cost parameter
	b		= 1.46;		% Unemployment flow value
	
@#if flex == 0
	lambda	= 0.001;	% Renegotiation frequency
@#endif
@#if flex == 1
	lambda	= 8/9;		% Renegotiation frequency
@#endif
@#if flex == 2
	lambda	= 11/12;	% Renegotiation frequency
@#endif
	
% Steady state values (ordering follows Appendix B in the 2006 Working Paper version)
	p_ss = 0.45;		% Job finding probability (p. 58; Shimer, 2005)

	z_ss = 1;
	n_ss = p_ss/(1-rho+p_ss);
	u_ss = 1-n_ss;
	x_ss = p_ss*u_ss/n_ss;
	mu_ss = 1/(1-lambda*beta);
	epsilon_ss = 1/(1-rho*lambda*beta);
	chi_ss = eta/(eta + (1-eta)*mu_ss/epsilon_ss);
	r_ss = 1/beta-(1-delta);
	k_y_ss = alpha/r_ss;
	i_y_ss = delta*k_y_ss;
	k_ss = (z_ss*k_y_ss)^(1/(1-alpha)) * n_ss;
	y_ss = z_ss * k_y_ss^alpha * n_ss^(1-alpha);
	% ac_ss = kappa/2*x_ss^2*n_ss;
	a_ss = (1-alpha)*z_ss*(k_ss/n_ss)^alpha;
	b_bar = b/(a_ss + kappa/2*x_ss^2);
	w_ss = chi_ss*(a_ss + kappa/2*x_ss^2 + kappa*p_ss*x_ss) + (1-chi_ss)*b;
	ls_ss = w_ss*n_ss/y_ss;
	c_y_ss = 1 - i_y_ss - kappa/2 * x_ss^2*n_ss/y_ss;
	
	H_ss = 4.09;

% Derived parameters
	aleph = beta/(kappa*x_ss);
	aleph_a = a_ss*aleph;
	aleph_w = w_ss*aleph;
	
	varphi_a = a_ss/w_ss;
	varphi_x = kappa*x_ss^2/w_ss;
	varphi_p = p_ss*beta*H_ss/w_ss;
	varphi_chi = chi_ss/(1-chi_ss)*kappa*x_ss/beta/w_ss;
	
	psi = chi_ss*mu_ss + (1-chi_ss)*epsilon_ss;
	tau = ((rho*beta*lambda)*psi) / (1 + (rho*beta*lambda)*psi);
	tau_1 = chi_ss*mu_ss*beta * ((rho-p_ss)*(x_ss*beta*lambda)*(lambda*mu_ss) + p_ss/eta)*(1-tau);
	tau_2 = chi_ss*mu_ss*(x_ss*beta*lambda)*(1-lambda*mu_ss)*(1-tau);
	varsigma = (1-lambda)*(1-tau)/lambda;
	phi = (1+tau_2) + varsigma + (tau/lambda - tau_1);
	
	gamma_b = (1+tau_2)/phi;
	gamma_o = varsigma/phi;
	gamma_f = (tau/lambda - tau_1)/phi;

model;
% Appendix B (pp. 82-84)
% B1-B20
	y = z + alpha*k(-1) + (1-alpha)*n(-1);
	y = c_y_ss*c + i_y_ss*i + (1-c_y_ss-i_y_ss)*(2*x+n(-1));
	m = sigma*u + (1-sigma)*v;
	n = n(-1) + (1-rho)*x;
	q = m - v;
	p = m - u;
	u = -(n_ss/u_ss) * n(-1);
	k = (1-delta)*k(-1) + delta*i;
	x = q + v - n(-1);
	0 = Lambda + (1-beta*(1-delta))*r(+1);
	Lambda = c - c(+1);
	x = aleph_a*a(+1) - aleph_w*w(+1) + Lambda + beta*x(+1);
	a = y - n(-1);
	r = y - k(-1);
	chi = -(1-chi_ss)*(mu-epsilon);
	epsilon = rho*lambda*beta*(Lambda + epsilon(+1));
	mu = (x_ss*lambda*beta)*x - (x_ss*lambda*beta)*(aleph_w*lambda*mu_ss)*(beta*lambda)*mu_ss*(w - w(+1)) 
		+ lambda*beta*(Lambda + mu(+1));
	w_o = chi_ss*varphi_a*a + (1-chi_ss)*varphi_p*p + (chi_ss*varphi_x + (1-chi_ss)*varphi_p)*x 
		+ varphi_chi*(chi - beta*(rho-p)*chi(+1));
	w = gamma_b*w(-1) + gamma_o*w_o + gamma_f*w(+1);
	z = rho_z*z(-1) + epsilon_z;
% Extras
	ls = w + n(-1) - y;
	theta = v - u;
end;

@#if moments == 0

	shocks;
		var epsilon_z; stderr 1;
	end;

	stoch_simul(order=1, noprint, irf=80, hp_filter=129600);

	@#if flex == 0
		irf0 = oo_.irfs;
		save irf_0 irf0;
	@#endif
	@#if flex == 1
		irf1 = oo_.irfs;
		save irf_1 irf1;
	@#endif
	
@#endif

@#if moments == 1

	shocks;
		var epsilon_z; stderr 100*sigma_z;
	end;

	stoch_simul(order=1, nofunctions, irf=0, hp_filter=129600);

	stoch_simul(order=1, periods=1000, irf=0, hp_filter=129600) u v;

	uu = oo_.endo_simul(strmatch('u', M_.endo_names, 'exact'),:)';
	vv = oo_.endo_simul(strmatch('v', M_.endo_names, 'exact'),:)';

	figure;
	scatter(uu, vv);
	title('Deviations from Hodrick-Prescott trend (%)');
	xlabel('Unemployment rate');
	ylabel('Vacancy rate');

	@#if flex == 0
		print('GT_BC_0','-depsc');
	@#endif
	@#if flex == 1
		print('GT_BC_1','-depsc');
	@#endif
	
@#endif