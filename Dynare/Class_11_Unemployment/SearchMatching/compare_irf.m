%% Plot IRF comparison
close all; clear; clc;

load irf_0;
load irf_1;

horizon	= 1:1:80;
sub_width = 3;
sub_height = 4;
sub_tot = sub_width*sub_height;
shocks	= {'_epsilon_z'};
vars	= {'y','w','ls','n','u','v','theta','p','x','k','i','z'};
   
for ii = 1:length(shocks)
    for jj = 1:length(vars)
        if mod(jj-1,sub_tot) == 0
            figure('units','normalized','outerposition',[0 0 1 1])
        end
        subplot(sub_height, sub_width, 1+mod(jj-1,sub_tot))
        hold on
        plot(horizon, eval(['irf0.' vars{jj}, shocks{ii}]), 'k-', 'LineWidth', 2)
        plot(horizon, eval(['irf1.' vars{jj}, shocks{ii}]), 'r-', 'LineWidth', 2)
        %plot(horizon, 0.*horizon, 'r-', 'LineWidth', 0.5) 
        title(vars{jj})
        hold off
    end
end

hL = legend({'flexible wages','staggered wages'});
set(hL,...
    'Position',[0.326190479223922 0.00428571481789857 0.394047612980717 0.0499999989356313],...
    'Orientation','horizontal',...
    'FontSize',10);

print('GT_irf','-depsc');