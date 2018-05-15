%% Prepare workspace
close all; clear; clc;

% Run two models to compare IRFs
dynare RBC nowarn
load('RBC_results.mat', 'oo_');
irf0 = oo_.irfs;
save irf_0.mat irf0;

dynare ind_lab nowarn
load('ind_lab_results.mat', 'oo_');
irf1 = oo_.irfs;
save irf_1.mat irf1;

%% Plot IRF comparison
close all; clear; clc;

load irf_0;
load irf_1;

horizon	= 1:1:40;
sub_width = 3;
sub_height = 3;
sub_tot = sub_width*sub_height;
shocks	= {'_epsilon'};
vars	= {'Output','Consumption','Investment','Capital','Hours','Wages','InterestRate','TFP','Productivity'};

for ii = 1:length(shocks)
    for jj = 1:length(vars)
        if mod(jj-1,sub_tot) == 0
            figure('units','normalized','outerposition',[0 0 1 1])
        end
        subplot(sub_height, sub_width, 1+mod(jj-1,sub_tot))
        hold on
        plot(horizon, eval(['irf0.' vars{jj}, shocks{ii}]), 'k-', 'LineWidth', 2)
        plot(horizon, eval(['irf1.' vars{jj}, shocks{ii}]), 'r--', 'LineWidth', 2)
        %plot(horizon, 0.*horizon, 'r-', 'LineWidth', 0.5) 
        title(vars{jj})
        hold off
    end
end

hL = legend({'divisible labor','indivisible labor'});
set(hL,...
    'Position',[0.326190479223922 0.00428571481789857 0.394047612980717 0.0499999989356313],...
    'Orientation','horizontal',...
    'FontSize',10);

print('ind_lab_irf','-depsc');