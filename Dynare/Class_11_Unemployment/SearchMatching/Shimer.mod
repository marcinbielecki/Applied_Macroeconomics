close all;

var			mpn w n u v theta m q J
			Output Employment Productivity Unemployment Wages Vacancies;
varexo		epsilon;
parameters	s chi eta kappa betta gama b rho sigma theta_ss;
	s		= 0.033;
	chi		= 0.45;
	eta		= 0.28;
	kappa	= 0.2098;
	betta	= 0.99;
	gama	= 0.28;
	b		= 0.4;
	rho		= 0.97;
	sigma	= 0.007;
	
	mpn_ss	= 1;
	// chi		= (1/betta-1+s)*kappa/(gama*(mpn_ss-b)-(1-gama)*kappa);
	
	x0		= 1;
	opt		= optimset('Display', 'off');
	sol		= fzero(@(x) theta_solve(x, betta, s, kappa, chi, eta, gama, mpn_ss, b), x0, opt);
	theta_ss= sol(1);
	
model;
	w		= gama*b + (1-gama)*(mpn + kappa*theta);
	J		= mpn - w + betta*(1-s)*J(+1);
	kappa/q	= betta * J(+1);
	u		= 1 - n;
	n		= (1-s)*n(-1) + m(-1);
	q		= chi * theta^(eta-1);
	m		= chi * v^eta * u^(1-eta);
	theta	= v / u;
	log(mpn)= rho*log(mpn(-1)) + epsilon;
	
	Output			= 100* log(mpn*n/(steady_state(mpn)*steady_state(n)));
	Employment		= 100* log(n/steady_state(n));
	Unemployment	= 100* log(u/steady_state(u));
	Productivity	= 100* log(mpn/steady_state(mpn));
	Wages			= 100* log(w/steady_state(w));
	Vacancies		= 100* log(v/steady_state(v));
end;

steady_state_model;
	mpn		= 1;
	theta	= theta_ss;
	q		= chi * theta^(eta-1);
	u		= s / (s + theta*q);
	v		= theta * u;
	n		= 1 - u;
	m		= s * n;
	w		= gama*b + (1-gama)*(mpn + kappa*theta);
	J		= kappa / (betta*q);
end;

steady;

shocks;
	var epsilon = sigma^2;
end;

stoch_simul;

stoch_simul(periods=1000, hp_filter=1600, nograph) Output Employment Productivity Unemployment Wages Vacancies;

uu = oo_.endo_simul(strmatch('Unemployment', M_.endo_names, 'exact'),:)';
vv = oo_.endo_simul(strmatch('Vacancies', M_.endo_names, 'exact'),:)';

figure;
scatter(uu, vv);
title('Deviations from Hodrick-Prescott trend (%)');
xlabel('Unemployment rate');
ylabel('Vacancy rate');
print('Shimer_BC','-depsc');