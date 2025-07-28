close all; clear all; clc;

Ncase = 2;
mfc_dir = [ "/Users/hyeoksu/MyWork/MFC-local/MFC/bubble/examples/0D_bubblecollapse_adap";
            "/Users/hyeoksu/MyWork/MFC-local/MFC/bubble/examples/0D_bubblecollapse_adap_multiscale";
            ];

% Space
m = 30; mp = m + 1;

% Time
tfinal = 0.05; dt = 0.0001; Nt = [501; 501]; t = 0:dt:tfinal;

% Refs
rho0 = 1000;
R0ref = 50.0e-06;
p0eq = 8236.0;
ub0 = sqrt(p0eq/rho0);
x0 = 4*R0ref;
p0 = 3*p0eq;
u0 = sqrt(p0/rho0);

% Rescale
rescale = [false true];

% Linestyle
linestyle = ["ro-","k<--","ro-","g<-","b>-","ro-","k-"];

% Figure
f1 = figure("DefaultAxesFontSize",18);
f1.Position = [100 100 1800 600];

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
        rho = A(     1:  mp); % Continuity
        mom = A(  mp+1:2*mp); % Continuity
        E   = A(2*mp+1:3*mp); % Continuity
        alf = A(3*mp+1:4*mp); % Continuity
        nR  = A(4*mp+1:5*mp); % Continuity
        nV  = A(5*mp+1:6*mp); % Continuity
        n   = A(6*mp+1:7*mp); % Continuity
        R(i,k) = nR(1)/n(1);

        if (rescale(i))
            R(i,k) = R(i,k)*(x0/R0ref);
        end
    end

    % Plot
    subplot(1,2,1);
    plot(t(1:Nt(i)),R(i,1:Nt(i)),linestyle(i),'LineWidth',2,'MarkerSize',10); hold on; grid on;
    xlabel("$t$",'Interpreter','latex');
    ylabel("$R/R_0$",'Interpreter','latex');
end
legend("single-scale","multi-scale");

subplot(1,2,2);
plot(t(1:min(Nt)),R(2,1:min(Nt)) - R(1,1:min(Nt)),linestyle(i),'LineWidth',2,'MarkerSize',10); hold on; grid on;
xlabel("$t$",'Interpreter','latex');
ylabel("$\epsilon (= R_{s}/R_0 - R_{m}/R_0)$",'Interpreter','latex');

saveas(f1,"results",'png');
close(f1);