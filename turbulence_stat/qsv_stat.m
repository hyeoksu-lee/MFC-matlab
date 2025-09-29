close all; clear all;

%% Initialization
disp("Initialize run_turbulence ..."); tic;
% Read user inputs
set_user_inputs(); load variables/user_inputs.mat;
toc;

for i = 1:Nfiles
    % Read qsv_data
    filename = "../../qsv_data.dat"; disp(filename);
    fileID = fopen(filename,'rb');
    val = fread(fileID, [9, Inf], 'double');
    fclose(fileID);
    
    % Assign values
    qsv_info1 = val(1,:);
    qsv_info5 = val(2,:);
    pres = val(3,:);
    liutex_mag = val(4,:);
    omega1 = val(5,:);
    omega2 = val(6,:);
    omega3 = val(7,:);
    omega_xy = sqrt(omega1.^2 + omega2.^2);
    vs_proj = val(8,:);
    vs_res = val(9,:);
    
    % Vortex
    qsv_info5_v = qsv_info5(qsv_info1 == 1);
    pres_v = pres(qsv_info1 == 1);
    liutex_mag_v = liutex_mag(qsv_info1 == 1);
    omega_xy_v = omega_xy(qsv_info1 == 1);
    vs_proj_v = vs_proj(qsv_info1 == 1);
    vs_res_v = vs_res(qsv_info1 == 1);
    
    % QSV
    pres_qsv = pres(qsv_info5 == 1);
    liutex_mag_qsv = liutex_mag(qsv_info5 == 1);
    omega_xy_qsv = omega_xy(qsv_info5 == 1);
    vs_proj_qsv = vs_proj(qsv_info5 == 1);
    vs_res_qsv = vs_res(qsv_info5 == 1);
    
    % non-QSV
    pres_nqsv = pres_v(qsv_info5_v == 0);
    liutex_mag_nqsv = liutex_mag_v(qsv_info5_v == 0);
    omega_xy_nqsv = omega_xy_v(qsv_info5_v == 0);
    vs_proj_nqsv = vs_proj_v(qsv_info5_v == 0);
    vs_res_nqsv = vs_res_v(qsv_info5_v == 0);
    
    % Plot
    plot_pdf(pres_v,"$p$", [-3, 1, 2], [-3:0.002:1.5], "pdf_pres", timesteps(i));
    plot_pdf(omega_xy_v,"$\omega_{xy}$", [-3, 1, 2], [-3:0.002:1.5], "pdf_omegaxy", timesteps(i));
    
    plot_pdf(pres_qsv, "$p f_{QSV}$", [-3, 1, 2], [-3:0.002:1.5], "pdf_pres_qsv", timesteps(i));
    plot_pdf(omega_xy_qsv, "$\omega_{xy} f_{QSV}$", [-3, 1, 2], [-3:0.002:1.5], "pdf_omegaxy_qsv", timesteps(i));

    plot_pdf(pres_nqsv, "$p f_{nonQSV}$", [-3, 1, 2], [-3:0.002:1.5], "pdf_pres_nqsv", timesteps(i));
    plot_pdf(omega_xy_nqsv, "$\omega_{xy} f_{nonQSV}$", [-3, 1, 2], [-3:0.002:1.5], "pdf_omegaxy_nqsv", timesteps(i));
end
disp("End of program");

% plot_pdf
function plot_pdf(var, varname, varaxis, histbin, filename, timestep)

    % Create dir if not exist
    if ~exist("results/"+filename, "dir")
        mkdir("results/"+filename);
    end

    f1 = figure("DefaultAxesFontSize",18);
    histogram(reshape(var,[],1),histbin,'EdgeColor','k','LineWidth',1.5,'Normalization','pdf', 'DisplayStyle', 'stairs'); hold on; grid on;
    xlim([varaxis(1) varaxis(3)]); 
    set(gca, 'YScale', 'log');
    xlabel(varname,'interpreter','latex');
    ylabel('$PDF$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    saveas(f1,"results/"+filename+"/tstep_"+string(timestep),"png"); 
    close(f1);
end

% plot_jpdf
function plot_jpdf(var1, var1name, var1axis, ...
                   var2, var2name, var2axis, ...
                   filename, timestep)

    load variables/user_inputs.mat;
    
    % Create dir if not exist
    if ~exist("results/"+filename, "dir")
        mkdir("results/"+filename);
    end
    
    f1 = figure("DefaultAxesFontSize",18);

    x = reshape(var1,[],1);
    y = reshape(var2,[],1);
    [counts, xEdges, yEdges] = histcounts2(x, y, 100);

    % Convert histogram counts to probability density
    binWidthX = xEdges(2) - xEdges(1);
    binWidthY = yEdges(2) - yEdges(1);
    jointPDF = counts / (sum(counts(:)) * binWidthX * binWidthY);

    % Define bin centers
    xCenters = xEdges(1:end-1) + binWidthX/2;
    yCenters = yEdges(1:end-1) + binWidthY/2;

    % Plot joint PDF as a contour plot
    contourf(xCenters, yCenters, log(jointPDF'), 20, 'LineColor', 'none'); hold on;
    if (var1name == "$p$")
        plot([pv pv],[0 10],'r--','LineWidth',1.5);
    end
    % xlim([var1axis(1) var1axis(3)]); xticks([var1axis(1):var1axis(2):var1axis(3)]);
    % ylim([var2axis(1) var2axis(3)]); yticks([var2axis(1):var2axis(2):var2axis(3)]);
    colorbar; caxis([-10 6]);
    xlabel(var1name,'Interpreter','latex');
    ylabel(var2name,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % Save figure
    saveas(f1,"results/"+filename+"/tstep_"+string(timestep),"png"); 
    close(f1);
end
