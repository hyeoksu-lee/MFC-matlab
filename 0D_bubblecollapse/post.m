close all; clear all; clc;

Ncase = 2;
mfc_dir = [ "/Users/hyeoksu/MyWork/MFC-local/MFC/adap_dt_max/examples/0D_bubblecollapse_adap";
            "/Users/hyeoksu/MyWork/MFC-local/MFC/adap_dt_max/case";
            ];

% Space
m = 30; mp = m + 1;
xa = -200; xb = 200; x = linspace(xa,xb,mp); dx = x(2) - x(1);

% Time
tfinal = 0.05; dt = 0.0001; Nt = [501; 301]; t = 0:dt:tfinal;

linestyle = ["ro-","k<--","ro-","g<-","b>-","ro-","k-"];

% Figure
f1 = figure("DefaultAxesFontSize",18);
f1.Position = [100 100 800 600];

%
R = zeros(Ncase,Nt(1));

% Run for each case
for i = 1:Ncase
    for k = 1:Nt(i)
        % Read data
        filename = strcat(mfc_dir(i),"/restart_data/lustre_",int2str(k),".dat");
        disp(filename);
        fileID = fopen(filename,'r');
        A = fread(fileID,'double');
        fclose(fileID);
    
        % Data
        rho = A(       1:  mp); % Continuity
        mom = A(  mp + 1:2*mp); % Continuity
        E   = A(2*mp + 1:3*mp); % Continuity
        alf = A(3*mp + 1:4*mp); % Continuity
        nR  = A(4*mp + 1:5*mp); % Continuity
        nV  = A(5*mp + 1:6*mp); % Continuity
        n   = A(6*mp + 1:7*mp); % Continuity
        R(i,k) = nR(1)/n(1);
    end

    % Plot
    plot(t(1:Nt(i)),R(i,1:Nt(i)),linestyle(i),'LineWidth',2,'MarkerSize',10); hold on; grid on;
    xlabel("$t$",'Interpreter','latex');
    ylabel("$R/R0$",'Interpreter','latex');
end
legend("adap\_dt\_max\_iter = 200","adap\_dt\_max\_iter = 100 (Default)");
saveas(f1,"results",'png');
close(f1);