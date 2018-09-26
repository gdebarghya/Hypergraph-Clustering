% Original code by Dohyung Park (Park et al., NIPS-2014)
% Modified by Debarghya Ghoshdastidar for the paper arXiv:1606.06516
% Modification: Added SCC, TTM, TeTrIS and SGC for comparison
%               Changed setting by adding noise to the generated points
%               Removed NSN+GSR and comparison based on NSE

clear all; close all; clc;

addpath otherSC
addpath SSC_ADMM_v1.1
addpath code2
addpath supp_material_RSCT/include
addpath SCC
addpath tetris

n_settings = 25;                  % # experimental settings

p_exp = 5*ones(1,n_settings);
d_exp = 3*ones(1,n_settings);
L_exp = 5*ones(1,n_settings);
n_exp = repmat(2:2:10,1,5).*d_exp;
o_exp = reshape(repmat((0:0.005:0.02),5,1),1,n_settings);
% p_exp = [ 5  5  5  5  5 10 10 10 10 10 20 20 20 20 20 35 35 35 35 35 50 50 50 50 50];      % ambient dimension
% d_exp = [ 3  3  3  3  3  6  6  6  6  6 12 12 12 12 12 21 21 21 21 21 30 30 30 30 30];      % subspace dimension
% L_exp = [ 5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5];    % # subspaces
% n_exp = [ 2  4  6  8 10  2  4  6  8 10  2  4  6  8 10  2  4  6  8 10  2  4  6  8 10] .* d_exp;

n_algo = 10;                 % # algorithms : SSC, GFS, LRR, TSC, NSN+Spectral, SCC, SGC, TeTrIS, TTM-uniform
n_trial = 50;               

E_exp = zeros(n_algo,n_settings,n_trial);     % Clustering error
%C_exp = zeros(n_algo,n_settings);     % Proportion of successful points
tic
for i_trial = 1:n_trial          % # trials
for i_exp = 1:n_settings
    
    disp(' ')
    disp(['Trial ' int2str(i_trial) '/' int2str(n_trial) ' : Setting ' int2str(i_exp) '/' int2str(n_settings)])
    
    %% Parameters
    p = p_exp(i_exp);                       % Ambient dimension
    L = L_exp(i_exp);                       % # subspaces
    d = d_exp(i_exp);                       % subspace dimension
    n = n_exp(i_exp)*ones(1,L);             % # sample points for each subspace
    N = sum(n);                             % Total # of samples
    o = o_exp(i_exp);
    %% True subspace generation   
    D0 = cell(1,L);
    for i=1:L
        D0{i} = orth(randn(p,d));
    end

    %% Data point generation
    A0 = zeros(1,N);           % True labels for sample points
    Y  = zeros(p,N);           % Sample points
    X0 = zeros(d,N);           % Weights for the sample points
    
    IDX = [1 cumsum(n)+1];
    for i=1:L
        X = normc(randn(d,n(i)));   % uniformly random unit vectors
    
        A0(:,IDX(i):IDX(i+1)-1) = i;
         Y(:,IDX(i):IDX(i+1)-1) = normc(D0{i}*X + o*randn(p,n(i)));
        X0(:,IDX(i):IDX(i+1)-1) = X;
    end
    
    %% SSC (l1-min)
    display('SSC..'); i_algo = 1;
    Z = admmLasso_mat_func(Y,false,20); Z_SSC = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
    
    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
   
    %% SSC-OMP (OMP)
    display('GFS..'); i_algo = 2;
    Z = OMPSC(Y,0); Z_GFS = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
    
    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
      
    %% LRR (nuclear norm minimization)
    display('LRR..'); i_algo = 3;
    Z = solve_lrr(Y,1e5); Z_LRR = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
    
    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
    
    %% TSC (NN)
    display('TSC..'); i_algo = 4; K = d-1;
    [A,Z,~] = TSC(Y,K,L); Z_TSC = Z;
        
    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);

    %% NSN + Spectral clustering
    display('NSN+Spectral..'); i_algo = 5; K = d-1; kmax = d;
    [Z,~] = NSN(Y,K); Z_NSN = Z;
    W = abs(Z) + abs(Z)';
    A = SpectralClusteringL(W,L);
    D = EstimateSubspace(Y,A,d,L);
        
    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0); 

    %% SCC
    display('SCC..'); i_algo = 6;
    opts.c = 250;
    [A,~] = scc(Y',d,L);

    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
        
    %% SGC
    display('SGC..'); i_algo = 7;
    [A,~] = tetris(Y',d,L);

    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
        
    % TETRIS
    display('TETRIS..'); i_algo = 8;
    [A,~,sigmaTetris] = tetris(Y',d,L);

    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
        
    % TTM-Uniform (same column)
    display('TTM-Uniform..'); i_algo = 9;
    opts.c = L*100;
    opts.sigma = sigmaTetris;
    [A,~] = ttm(Y',d,L,opts);

    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
      
    % TTM-Uniform (double column)
    display('TTM-Uniform2..'); i_algo = 10;
    opts.c = 2*L*100;
    opts.sigma = sigmaTetris;
    [A,~] = ttm(Y',d,L,opts);

    E_exp(i_algo,i_exp,i_trial) = computeCE(A,A0);
        
    %
%    mean(E_exp,3)
    toc
end
end

mE_exp = mean(E_exp,3)

close all;
n_algo = 10; algo = {'SSC', 'SSC-OMP', 'LRR', 'TSC', 'NSN+Spectral' 'SCC' 'SGC' 'TeTrIS' 'TTM-Uniform (100k)' 'TTM-Uniform (200k)'};
figure;
set(gcf, 'Position', [1 1 1000 400]);
for i_algo=1:n_algo
    subplot(2,n_algo/2,i_algo);
    pos = get(gca, 'Position');
%    set(gca, 'Position', [pos(1) pos(2) 0.1 0.2] );
    h = imagesc(flipud(reshape(mE_exp(i_algo,:),5,5)'),[0 1]); colormap(flipud(gray));
    set(gca, 'xTick', 1:5);
    set(gca, 'xTickLabel', 6:6:30);
    set(gca, 'yTick', 1:5);
    set(gca, 'yTickLabel', 0.2:-0.05:0);
    title(algo{i_algo},'FontSize',14);
    if i_algo == 1
        ylabel('Noise level (\sigma_a)','FontSize',16);
    end
    if i_algo == 8
        xlabel('Number of points for each subspace (n/k)','FontSize',16);
    end
end
colorbar;
