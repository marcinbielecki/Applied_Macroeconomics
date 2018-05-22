// close all;

var			y c inv k h r_k w z R Pi mc Delta p_tilde Num Den lambda d;
varexo		eps_z eps_d eps_R;
parameters	alffa betta delta eta mu phi theta sigmma
			rho_z rho_d
			gamma_R gamma_Pi gamma_y Pi_ss;
			
	alffa	= 0.33;		// capital share in output
	betta	= 0.99;		// quarterly discount factor
	delta	= 0.025;	// quarterly capital depreciation rate
	eta		= 2;		// inverse of Frisch elasticity of labor
	mu		= 1.33;		// average monopolistic competition markup
	theta	= 0.75;		// 1 - quarterly probability of resetting price
	sigmma	= 2;		// inverse of elasticity of intertemporal substitution
	
	rho_z	= 0.95;		// autoregressivity of supply shock
	rho_d	= 0.95;		// autoregressivity of demand shock
	
	gamma_R	= 0.8;		// autoregressivity of nominal interest rates
	gamma_Pi = 1.5;		// strength of central bank's reaction to inflation gap, should be > 1
	gamma_y	= 0.5;		// strength of central bank's reaction to output gap
	Pi_ss	= 1;		// (gross) quarterly inflation target - assumed 1 for simplicity
	
	// Set phi via formula so that h = 1/3 in steady state
	k_h_ss	= (alffa/mu/(1/betta-1+delta))^(1/(1-alffa));
	c_h_ss	= k_h_ss^alffa-delta*k_h_ss;
	h_ss	= 1/3;
	phi		= ((1-alffa)/mu*k_h_ss^alffa/(h_ss^(eta+sigmma)*c_h_ss^sigmma));

model;
	// (0) Marginal utility of consumption
	lambda	= d * c^(-sigmma);
	// (1) Euler equation
	1		= betta * lambda(+1)/lambda * (R/Pi(+1));
	// (2) Consumption-hours choice
	w		= phi * h^eta * c^sigmma;
	// (3) Asset markets
	// 1		= betta * lambda(+1)/lambda * (r_k(+1)+1-delta);
	inv		= delta * k;
	// (4) Real wages
	w		= mc * (1-alffa) * z * k(-1)^alffa * h^(-alffa);
	// (5) Capital rental rate
	r_k		= mc * alffa * z * k(-1)^(alffa-1) * h^(1-alffa);
	// (6) Production function
	y		= z * k(-1)^alffa * h^(1-alffa) / Delta;
	// (7) Price dispersion
	Delta	= theta*Delta(-1)*Pi^(mu/(mu-1)) + (1-theta)*p_tilde^(mu/(1-mu));
	// (8) Inflation dynamics
	Pi		= (theta / (1 - (1-theta) * p_tilde^(1/(1-mu))))^(1-mu);
	// (9) Optimal reset price
	p_tilde	= mu * Num / Den;
	// (10) Numerator
	Num		= lambda*mc*y + betta*theta*Pi^(mu/(mu-1))*Num(+1);
	// (11) Denominator
	Den		= lambda*y + betta*theta*Pi^(1/(mu-1))*Den(+1);
	// (12) Investment
	// k		= (1-delta)*k(-1) + inv;
	k		= steady_state(k);
	// (13) Output accounting
	y		= c + inv;
	// (14) TFP AR(1) process
	log(z)	= rho_z * log(z(-1)) + eps_z;
	// (15) Monetary policy rule
	R		= R(-1)^gamma_R * (steady_state(R) * (Pi/steady_state(Pi))^gamma_Pi * (y/steady_state(y))^gamma_y)^(1-gamma_R) * exp(eps_R);
	// (16) Demand shock AR(1) process
	log(d)	= rho_d * log(d(-1)) + eps_d;
end;

steady_state_model;
	z		= 1;
	d		= 1;
	Pi		= Pi_ss;
	R		= Pi / betta;
	r_k		= 1 / betta - (1-delta);
	p_tilde	= 1;
	Delta	= 1;
	mc		= p_tilde / mu;
	k_h		= (mc * alffa / r_k)^(1/(1-alffa));
	y_h		= k_h^alffa;
	inv_h	= delta * k_h;
	c_h		= y_h - inv_h;
	w		= mc * (1-alffa) * k_h^alffa;
	h		= (w / (phi * c_h^sigmma))^(1/(eta+sigmma));
	c		= c_h * h;
	inv		= inv_h * h;
	y		= y_h * h;
	k		= k_h * h;
	lambda	= c^(-sigmma);
	Num		= (lambda*mc*y)/(1-betta*theta);
	Den		= (lambda*y)/(1-betta*theta);
end;

shocks;
	var eps_z = 0.01^2;
	var eps_d = 0.01^2;
	var eps_R = 0.01^2;
end;

steady;
// check;

stoch_simul(order=1) y Pi R;