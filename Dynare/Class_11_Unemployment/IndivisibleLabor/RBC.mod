% Preamble: declare variable types
var			y c i k h w r z y_h							% Endogenous variables, 
			Output Consumption Investment Capital 		% separated by comma or whitespace
			Hours Wages InterestRate TFP Productivity;	% Semicolon is necessary to end the line
varexo		epsilon;									% Exogenous variables, usually shocks
parameters	alpha beta delta phi rho_z sigma_z;			% Parameter names

predetermined_variables k;								% Current k was decided in the past

% Declare parameter values
	alpha	= 0.33;		% Capital share of income
	beta	= 0.99;		% Households' discount factor
	delta	= 0.025;	% Quarterly depreciation rate
	phi		= 1.75;		% Disutility of labor
	rho_z	= 0.9622;	% Autocorrelation of deviation of TFP from trend
	sigma_z	= 0.00853;	% Standard deviation of shock to TFP

% Write down model equations: first order conditions and equilibrium conditions
model;
% 8 main equations
	1		= beta * c/c(+1) * (1+r(+1));	% Euler equation
	h		= 1 - phi*c/w;					% Consumption-hours choice
	y		= z * k^alpha * h^(1-alpha);	% Production function
	r		= alpha*y/k - delta;			% Real interest rate
	w		= (1-alpha)*y/h;				% Real hourly wage
	i		= k(+1) - (1-delta)*k;			% Investment
	y		= c + i;						% Output accounting
	log(z)	= rho_z*log(z(-1)) + epsilon;	% TFP AR(1) process

% Definition of labor productivity
	y_h		= y/h;

% Counterparts to our observables
	Output			= 100 * log(y/steady_state(y));
	Consumption		= 100 * log(c/steady_state(c));
	Investment		= 100 * log(i/steady_state(i));
	Capital			= 100 * log(k/steady_state(k));
	Hours			= 100 * log(h/steady_state(h));
	Wages			= 100 * log(w/steady_state(w));
	InterestRate	= 400 * (r-steady_state(r));
	TFP				= 100 * log(z/steady_state(z));
	Productivity	= 100 * log(y_h/steady_state(y_h));
end;

% Instructions for finding the steady state of the model
steady_state_model;
	z		= 1;
	r		= 1/beta - 1;
	k_h		= (alpha / (r+delta))^(1/(1-alpha));
	y_h		= k_h^alpha;
	w		= (1-alpha) * y_h;
	c_h		= y_h - delta*k_h;
	h		= 1 / (1 + phi * c_h / w);
	k		= k_h * h;
	c		= c_h * h;
	y		= y_h * h;
	i		= y - c;
end;

% Declare shock variance
shocks;
	% var epsilon = sigma_z^2;				% When calculating model statistical moments
	var epsilon	; stderr 0.01;						% When producing impulse response functions
end;

% Calculate steady state
steady;

% Check if model is stable
check;

% Simulate the model
stoch_simul(hp_filter=1600);