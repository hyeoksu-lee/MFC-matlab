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
f_tke_budget = figure("DefaultAxesFontSize", 18); 

for l = 1:ncase
    load(post_stat_dir(l)+"/tke_budget.mat");
    plot(y_norm,T,'ko-','LineWidth',2); hold on;
    plot(y_norm,P,'ro-','LineWidth',2);
    plot(y_norm,D,'bo-','LineWidth',2);
end

% axis([0 400 0 0.1]);
xlabel('$y / \delta_{\theta} $','interpreter','latex');
legend( "$T$","$P$","$\epsilon$",...
        'Interpreter','latex','location','northeast');
set(gca,'TickLabelInterpreter','latex');

saveas(f_tke_budget, output_dir+"/tke_budget.png"); 
close(f_tke_budget); 
disp("f_tke_budget saved");
