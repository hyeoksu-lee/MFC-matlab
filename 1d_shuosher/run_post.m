close all; clear all;

Ncase = 7;
mfc_dir = [ "/p/global/hyeoksu/MFC/reynolds/case/1d_shuosher/N200/weno5m";
            "/p/global/hyeoksu/MFC/reynolds/case/1d_shuosher/N200/weno7m";
            "/p/global/hyeoksu/MFC/reynolds/case/1d_shuosher/N200/wcns6ld";
            "/p/global/hyeoksu/MFC/reynolds/case/1d_shuosher/N1000/weno5m";
            "/p/global/hyeoksu/MFC/reynolds/case/1d_shuosher/N1000/weno7m";
            "/p/global/hyeoksu/MFC/reynolds/case/1d_shuosher/N1000/wcns6ld";
            "/p/global/hyeoksu/MFC/reynolds/case/1d_shuosher/N3000/weno5m";];

% Space and time
m = [200 200 200 1000 1000 1000 3000];
timesteps = m*2;
linestyle = ["g<-","b>-","ro-","g<-","b>-","ro-","k-"];

% Figure
f1 = figure("DefaultAxesFontSize",18);
f1.Position = [100 100 1500 600];

% Run for each case
for i = 1:Ncase
    mp = m(i) + 1;
    dx = 10/m(i); x = 0:dx:10;

    % Read data
    filename = strcat(mfc_dir(i),"/restart_data/lustre_",int2str(timesteps(i)),".dat");
    disp(filename);
    fileID = fopen(filename,'r');
    A = fread(fileID,'double');
    fclose(fileID);

    % Data
    rho = A(1:mp); % Continuity

    % Plot
    if (i == 1 || i == 2 || i == 3)
        subplot(1,2,1);
        plot(x,rho,linestyle(i)); hold on;
    elseif(i == 4 || i == 5 || i == 6)
        subplot(1,2,2);
        plot(x,rho,linestyle(i)); hold on;
    else
        subplot(1,2,1);
        plot(x,rho,linestyle(i)); hold on;
        subplot(1,2,2);
        plot(x,rho,linestyle(i)); hold on;
    end
end

subplot(1,2,1);
legend("WENO5-M, N=200","WENO7-M, N=200","WCNS6-LD, N=200","WENO5-M, N=3000",'location','southwest');
subplot(1,2,2);
legend("WENO5-M, N=1000","WENO7-M, N=1000","WCNS6-LD, N=1000","WENO5-M, N=3000",'location','southwest');
saveas(f1,"results/results",'png');
close(f1);