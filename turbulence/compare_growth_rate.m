close all; clear all;
format long

% Remove pre-existing var files
if exist("variables", 'dir')
    delete variables/*.mat;
end

%% INPUTS
% Set directories
set_directory(); load variables/directory.mat;
% Set multiscale parameters
set_multiscale(ncase); load variables/multiscale.mat;
% Set grids
set_space(ncase); load variables/space.mat;
% Set time steps
set_time(ncase); load variables/time.mat;
% Set line style
set_linestyle(ncase); load variables/linestyle.mat;
% Set options
set_options(ncase); load variables/options.mat;
% Set polydisperse bubbles
set_poly_bubbles(); load variables/poly_bubbles.mat;
% Set index
set_index(); load variables/index.mat;

% Load data
%% POST-PROCESS

% Integrated growth rate
f_growth_rate = figure("DefaultAxesFontSize", 18); 
plot([0 400],[0.014 0.014], 'k-.','LineWidth',2); hold on; grid on;
plot([0 400],[0.016 0.016], 'k--.','LineWidth',2);
plot([0 400],[0.022 0.022], 'k:','LineWidth',2);

for l = 1:ncase
    load(post_stat_dir(l)+"/momentum_thickness.mat");
    dmth = (mth(2:Nfiles(1)) - mth(1:Nfiles(1)-1)) ./ (time(2:Nfiles(1)) - time(1:Nfiles(1)-1));
    plot(time(1:Nfiles(1)-1),dmth,linestyle(l),'LineWidth',2);
end

axis([0 400 0 0.1]);
xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
ylabel('$\dot{\delta} / U_1$','interpreter','latex');
legend( "$0.014^{[1-4]}$","$0.016^{[5-6]}$","$0.022^{[1]}$", ...
        "Present, WENO3M",...
        "Present, WENO5M",...
        "Present, WENO7M",...
        'Interpreter','latex','location','northeast');
set(gca,'TickLabelInterpreter','latex');

saveas(f_growth_rate, output_dir+"/growth_rate.png"); 
close(f_growth_rate); 
disp("f_growth_rate saved");
