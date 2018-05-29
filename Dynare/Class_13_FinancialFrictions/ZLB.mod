close all;

var			y Pi R z d;
varexo		eps_z eps_d eps_R;
parameters	alffa betta eta sigmma theta
			gamma_R gamma_Pi gamma_y
			rho_z rho_d;
		
alffa	= 0.33;		// capital share in output
betta	= 0.99;		// quarterly discount factor
eta		= 2;		// parameter in the utility function, inverse of Frisch elasticity of labor
sigmma	= 2;		// parameter in the utility function, inverse of elasticity of intertemporal substitution
theta	= 0.75;		// price stickiness, average price duration = 1/(1-theta) = 4 quarters

gamma_R	= 0.85;		// autoregressivity of nominal interest rates
gamma_Pi = 1.5;		// strength of central bank's reaction to inflation gap, should be > 1
gamma_y	= 0.5;		// strength of central bank's reaction to output gap

rho_z	= 0.95;		// autoregressivity of supply shock
rho_d	= 0.95; 	// autoregressivity of demand shock

model;
	// (1) Dynamic IS curve
	y	= y(+1) - 1/sigmma * (R - Pi(+1)) + 1/sigmma * (1-rho_d) * d;
	// (2) Dynamic AS curve
	theta*Pi = betta*theta*Pi(+1) + ((1-betta*theta)*(1-theta)*(alffa+eta+sigmma*(1-alffa))/(1-alffa))*y  - ((1-betta*theta)*(1-theta)*(1+eta)/(1-alffa))*z;
	// (3) Taylor rule
	R	= gamma_R*R(-1) + (1-gamma_R)*(gamma_Pi*Pi + gamma_y*y) + eps_R;
	// R	= max(betta-1, gamma_R*R(-1) + (1-gamma_R)*(gamma_Pi*Pi + gamma_y*y) + eps_R);
	// (4-5) Shock processes
	z	= rho_z*z(-1) + eps_z;
	d	= rho_d*d(-1) + eps_d;
end;

shocks;
	var eps_d;
	periods 1:2;
	values -0.09;
end;

simul(periods = 200);

t=[1:40];
figure;
plotyy(t,400*(1/betta+R(1:40)-1),t,100*y(1:40),'plot');
legend('R','y');