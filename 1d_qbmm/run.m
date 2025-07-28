close all; clear all; clc;

% INPUT
Ncase = 2;

%
Nx = [100 100];
mp = Nx + 1;

%
vf0 = 1e-5;

% EQUATIONS 
num_fluids = 1;
num_dims   = 1;
model_eqns = 2;
nb         = 1;
bubbles    = true;
qbmm       = true;
adv_n      = false;

% TIME STEP 
Nfiles       = 20;
t_step_save  = [50 50];

% INDEXING 
contxb = 1;
contxe = num_fluids;
momxb  = contxe + 1;
momxe  = contxe + num_dims;
E_idx  = momxe + 1;
advxb  = E_idx + 1;
advxe  = E_idx + num_fluids;
sys_size = advxe;
if (model_eqns == 3)
    inteb = sys_size + 1;
    intee = sys_size + num_fluids;
    sys_size = intee;
end
if (bubbles)
    bubxb  = sys_size + 1;
    if (qbmm)
        bubxe  = sys_size + 6*nb;
    else
        bubxe  = sys_size + 2*nb;
    end
    sys_size = bubxe;
    if (adv_n)
        n_idx  = sys_size + 1;
        sys_size = n_idx;
    end
end

% INDEX TO READ
idx_read = [E_idx; advxe];

% CASE DIRECTORIES
mfc_dir = [ "/Users/hyeoksu/MyWork/MFC-local/MFC/bubble/case/1D_qbmm/multi";
            "/Users/hyeoksu/MyWork/MFC-local/MFC/MFC/examples/1D_qbmm";
            ];
output_dir = "./results";

linestyle1 = ["ko-", "r^-", "g-", "r-.", "go"];
linestyle2 = ["ko-", "b^-", "g-", "r-.", "go"];

rescale = [true false];
% ==============================================================================

% POST-PROCESS
for i = 1:Nfiles+1
    f1 = figure("DefaultAxesFontSize",20);
    % set(f1,'Position',[100 100 2500 2000]);
    for k = 1:Ncase
        var = zeros(mp(k),sys_size);

        % Read data from file
        s = sprintf('%06d',(i - 1)*t_step_save(k));
        for l = 1:length(idx_read)
            filename = strcat(mfc_dir(k),"/D/prim.",int2str(idx_read(l)),".00.",s,".dat");

            flag = 1; % Flag for file access
            if ~exist(filename,'file')
                disp(strcat(filename," does not exist."));
                flag = 0;
            end

            if (flag == 1)
                disp(filename);
                A = readmatrix(filename);
                x = A(:,1);
                var(:,idx_read(l)) = A(:,2);
            end
        end

        if (rescale(k))
            x = x/5;
            var(1:1:mp(k),E_idx) = var(1:1:mp(k),E_idx)*3;
        end

        yyaxis left;
        plot(x(1:mp(k)),var(1:mp(k),E_idx),linestyle1(k),'LineWidth',1,'MarkerSize',6); hold on;
        xlabel("$x/R_{0\mbox{ref}}$",'interpreter','latex');
        ylabel("$p/p_{eq}$",'interpreter','latex');
        xlim([-1000 1000]); 
        ylim([0 2]); yticks([0:0.5:2]);
        set(gca, 'YColor','r');
        set(gca,'TickLabelInterpreter','latex');

        yyaxis right;
        plot(x(1:mp(k)),var(1:mp(k),advxe)/vf0,linestyle2(k),'LineWidth',1,'MarkerSize',6); hold on;
        ylabel("$\alpha/\alpha_0$",'interpreter','latex');
        ylim([0 6]); yticks([0:2:6]);
        set(gca, 'YColor','b');
        set(gca,'TickLabelInterpreter','latex');

    end
    if (flag == 0) 
        break;
    end
    legend("multi-scale","single-scale",'location','northeast');
    saveas(f1,strcat(output_dir,"/t",int2str(i-1),".png"),'png');
    close(f1);
end
