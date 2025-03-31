function set_linestyle(ncase)

    load variables/time.mat;

    linestyle = ["bo-", "g^--", "r+-.", "k-" ];

    linewidth = [1, 1, 1, 2];    

    % blue1 = [0.3569,0.8118,0.9569]; blue2 = [0.0196,0.0745,0.6706];
    % blueGrad = [linspace(blue1(1),blue2(1),max(Nfiles))', ...
    %             linspace(blue1(2),blue2(2),max(Nfiles))', ...
    %             linspace(blue1(3),blue2(3),max(Nfiles))'];

    save variables/linestyle.mat;
end