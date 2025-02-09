function plot_Reynolds_stress(f_Reynolds_stress, y_norm, ruu, rvv, rww, ruv)

    load variables/linestyle.mat;
    start_idx = 23;
    
    figure(f_Reynolds_stress);

    % sqrt(ruu)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/ruu.dat");
    y_norm_ref1 = A(:,1); ruu_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/ruu.dat");
    y_norm_ref2 = A(:,1); ruu_ref2 = A(:,2);
    A = readmatrix("data/Wang_et_al_2022/ruu.dat");
    y_norm_ref3 = A(:,1); ruu_ref3 = A(:,2);
    subplot(2,2,1); hold on;
    % plot(y_norm_ref1,ruu_ref1,'g+','LineWidth',2);
    plot(y_norm_ref2,ruu_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,ruu_ref3,'r--','LineWidth',1);
    plot(y_norm(:,start_idx:end),ruu(:,start_idx:end),'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.2]);
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.2]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} u^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(rvv)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/rvv.dat");
    y_norm_ref1 = A(:,1); rvv_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/rvv.dat");
    y_norm_ref2 = A(:,1); rvv_ref2 = A(:,2);
    A = readmatrix("data/Wang_et_al_2022/rvv.dat");
    y_norm_ref3 = A(:,1); rvv_ref3 = A(:,2);
    subplot(2,2,2); hold on;
    % plot(y_norm_ref1,rvv_ref1,'r+');
    plot(y_norm_ref2,rvv_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,rvv_ref3,'r--','LineWidth',1);
    plot(y_norm(:,start_idx:end),rvv(:,start_idx:end),'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.2]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.2]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} v^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(rww)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/rww.dat");
    y_norm_ref1 = A(:,1); rww_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/rww.dat");
    y_norm_ref2 = A(:,1); rww_ref2 = A(:,2);    
    A = readmatrix("data/Wang_et_al_2022/rww.dat");
    y_norm_ref3 = A(:,1); rww_ref3 = A(:,2);
    subplot(2,2,3); hold on;
    % plot(y_norm_ref1,rww_ref1,'r+');
    plot(y_norm_ref2,rww_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,rww_ref3,'r--','LineWidth',1);
    plot(y_norm(:,start_idx:end),rww(:,start_idx:end),'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.2]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.2]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} w^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(-rvu)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/ruv.dat");
    y_norm_ref1 = A(:,1); ruv_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/ruv.dat");
    y_norm_ref2 = A(:,1); ruv_ref2 = A(:,2);   
    A = readmatrix("data/Wang_et_al_2022/ruv.dat");
    y_norm_ref3 = A(:,1); ruv_ref3 = A(:,2); 
    subplot(2,2,4); hold on;
    % plot(y_norm_ref1,ruv_ref1,'r+');
    plot(y_norm_ref2,ruv_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,ruv_ref3,'r--','LineWidth',1);
    plot(y_norm(:,start_idx:end),ruv(:,start_idx:end),'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.2]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.2]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{-\left< \overline{\rho} u^{\prime} v^{\prime}\right>}  / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
end