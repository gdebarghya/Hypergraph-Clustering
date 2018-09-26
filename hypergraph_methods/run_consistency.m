% This code compares the number of nodes misclustered by TTM, HOSVD and
% NH-Cut algorithms when uniform hypergraph has a planted bisection
% Two settings are studied (uncomment each set of parameters to run code)

clear all; clc;

trials = 50;         % no. of runs
k = 2;               % no. of clusters
s = 5:5:50;          % no. of nodes per cluster (varies in each plot
num = length(s);        


% Setting-1: varying m, fixed p,q
m = [2    3    4];          % order of hypergraph
p = [0.3  0.3  0.3];        % within cluster edge prob
q = [0.2  0.2  0.2];        % inter cluster edge prob


% % Setting-2: varying p, fixed m,q
% m = [3     3     3];          % order of hypergraph
% p = [0.3   0.25  0.225];      % within cluster edge prob
% q = [0.2   0.2   0.2];        % inter cluster edge prob


%%%%
n = k*s;
err = zeros(3,length(n),length(m),trials);
% Algo sequence (row wise): TTM, HOSVD, NH-Cut

tic
for t = 1:trials
    for i_exp = 1:length(m)
    disp(['Trial ' int2str(t) '/' int2str(trials) ' : Exp ' int2str(i_exp)])
        for i_n = 1:length(n)
            err(:,i_n,i_exp,t) = n(i_n)*planted_hypergraph(s(i_n),m(i_exp),k,p(i_exp),q(i_exp));
        end
    end
    toc
end

mean_err = mean(err,4);

set(gcf, 'Position', [1 1 1000 200]);
for i_exp = 1:length(m)
    subplot(1,1+length(m),i_exp)
    hold on
    plot(n,mean_err(1,:,i_exp),'-ks','LineWidth',1)
    plot(n,mean_err(2,:,i_exp),'-ko','LineWidth',1)
    plot(n,mean_err(3,:,i_exp),'-kx','LineWidth',1)
    hold off
    box on
%    ylim([0 0.5])
end
legend('TTM','HOSVD','NH-Cut','Location','EastOutside')
