close all;

var			y c inv k h r w z y_h
			Output Consumption Investment
			Capital Hours Interest_rate
			Wages TFP Productivity Markup
			R mu;
varexo		eps_z eps_mu;
parameters	alffa betta delta rho_z sigma_z mu_ss rho_mu sigma_mu phi;
	alffa	= 0.33;
	betta	= 0.99;
	delta	= 0.025;
	
	rho_z	= 0.97;
	sigma_z = 0.0065;	//0.009;	//
	
	mu_ss	= 1.4;		//1;		//
	rho_mu	= 0.99;
	sigma_mu= 0.011;
	
	k_h_ss	= (alffa/mu_ss/(1/betta-1+delta))^(1/(1-alffa));
	h_ss	= 1/3;
	phi		= (1-h_ss)/h_ss*(1-alffa)/mu_ss*k_h_ss^alffa/(k_h_ss^alffa-delta*k_h_ss);

model;
	1		= betta * c/c(+1) * (1+r(+1));
	h		= 1 - phi * c/w;
	R		= r + delta;
	y		= z * k(-1)^alffa * h^(1-alffa);
	R		= 1/mu * alffa * y/k(-1);
	w		= 1/mu * (1-alffa) * y/h;
	inv		= k - (1-delta)*k(-1);
	y		= c + inv;
	log(z)	= rho_z*log(z(-1)) + eps_z;
	log(mu)	= (1-rho_mu)*log(mu_ss) + rho_mu*log(mu(-1)) + eps_mu;
	y_h		= y / h;

	Output			= 100* log(y/steady_state(y));
	Consumption		= 100* log(c/steady_state(c));
	Investment		= 100* log(inv/steady_state(inv));
	Capital			= 100* log(k/steady_state(k));
	Hours			= 100* log(h/steady_state(h));
	Interest_rate	= 400*    (r-steady_state(r));
	Wages			= 100* log(w/steady_state(w));
	TFP				= 100* log(z/steady_state(z));
	Productivity	= 100* log(y_h/steady_state(y_h));
	Markup			= 100* log(mu/steady_state(mu));
end;

steady_state_model;
	z	= 1;
	mu	= mu_ss;
	r	= 1/betta - 1;
	R	= r + delta;
	k_h	= (alffa/mu/(r+delta))^(1/(1-alffa));
	y_h	= k_h^alffa;
	w	= (1-alffa)/mu * y_h;
	c_h	= y_h - delta*k_h;
	h	= 1 / (1 + phi * c_h / w);
	k	= k_h * h;
	c	= c_h * h;
	y	= y_h * h;
	inv	= y - c;
end;

shocks;
	var eps_z = sigma_z^2;
	var eps_mu = sigma_mu^2;
end;

steady;
check;

stoch_simul(hp_filter=1600, nofunctions) Output Consumption Investment Hours Wages Productivity Markup TFP;