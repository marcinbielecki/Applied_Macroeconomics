close all;
%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c i g ce n rk r q k x a h pi rn premium;
varexo e_rn e_g e_a;
parameters beta eta alph delt omeg eps G_Y C_Y I_Y Ce_Y Y_N X rho_a rho_g psi K_N R gam mu nu thet rho sig kap;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

C_Y		= 0.64;
Ce_Y	= 0.01;
I_Y		= 0.15;
G_Y		= 2;
K_N		= 2;
Y_N		= 0.282494996;
X		= 1.1;
beta    = 0.99;
R       = 1/beta;
alph    = 0.35; 
eta     = 3;
omeg    = 0.99;
delt    = 0.025;
rho_a   = 1;
rho_g   = 0.95;
psi     = 0.25;
Rk      = R + 0.02;
gam     = 1-0.0272;
mu      = 0.12;
thet    = 0.75;
rho     = 0.9;
sig     = 0.11;
kap     = ((1-thet)/thet)*(1-thet*beta);
eps     = (1-delt)/((1-delt) + ((alph/X)*(Y_N/K_N)));
nu      = 0.052092347;

%----------------------------------------------------------------
% 3. Model 
%----------------------------------------------------------------

model(linear);

%Aggregate Demand

y          = C_Y*c + I_Y*i + G_Y*g + Ce_Y*ce;                    //4.14
c          = -r + c(+1);                                         //4.15       
ce         = n;                                                  //4.16
rk(+1) - r = -nu*(n -(q + k));                                   //4.17
rk         =(1-eps)*(y - k(-1) - x) + eps*q - q(-1);             //4.18
q          = psi*(i - k(-1));                                    //4.19    

%Aggregate Supply 

y             = a + alph*k(-1) + (1-alph)*omeg*h;                //4.20
y - h - x - c = (eta^(-1))*h;                                    //4.21
pi            = kap*(-x) + beta*pi(+1);                          //4.22

%Evolution of State Variables

k = delt*i + (1-delt)*k(-1);                                     //4.23
n = gam*R*K_N*(rk - r(-1)) + r(-1) + n(-1);                      //4.24

% Monetary Policy Rules and Shocks

rn      = rho*rn(-1) + sig*pi(-1) - e_rn;                        //4.25 
g       = rho_g*g(-1) + e_g;                                     //4.26
a       = rho_a*a(-1) + e_a;                                     //4.27
rn      = r + pi(+1);                                            // Fisher Equation
premium = rk(+1) - r;                                        
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

check;
steady;

shocks;
var e_rn; stderr 0.25;  %in the paper a decline in the nom. int. rate by 25 basis points is considered
end;

nu_values = [0.052092347 0];

for j = 1:2;

    nu = nu_values(j);

    stoch_simul(irf=24, nograph);

    % save IRFs
    y_data(:,j)       = y_e_rn;
    i_data(:,j)       = i_e_rn;
    rn_data(:,j)      = rn_e_rn;
    premium_data(:,j) = premium_e_rn;
    n_data(:,j)       = n_e_rn;
    k_data(:,j)       = k_e_rn;
    q_data(:,j)       = q_e_rn;
    pi_data(:,j)      = pi_e_rn;
    rk_data(:,j)      = rk_e_rn;
    r_data(:,j)       = r_e_rn;
    c_data(:,j)       = c_e_rn;
    h_data(:,j)       = h_e_rn;

end;

%----------------------------------------------------------------
% 5. Plots
%----------------------------------------------------------------

length_irf = 12; % from the paper
tt = 0:length_irf;

figure
subplot(2,2,1);plot(tt,y_data(tt+1,1), tt,y_data(tt+1,2), '--');title('Output')
subplot(2,2,2);plot(tt,i_data(tt+1,1), tt,i_data(tt+1,2), '--');title('Investment')
subplot(2,2,3);plot(tt,rn_data(tt+1,1), tt,rn_data(tt+1,2), '--');title('Nominal Interest Rate')
legend('with Financial Accelerator' , 'without Financial Accelerator','Location', 'southeast');
subplot(2,2,4);plot(tt,premium_data(tt+1,1), tt,premium_data(tt+1,2), '--');title('Risk Premium')
