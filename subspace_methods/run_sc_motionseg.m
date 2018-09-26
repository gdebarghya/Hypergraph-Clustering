% Original code by Dohyung Park, Constantine Caramanis & Sujay Sanghavi (NIPS,2014)
% Modified by Debarghya Ghoshdastidar for the paper arXiv:1606.06516
% Modifications: Added TeTrIS and SGC for comparison, and average results
% over 20 independent runs (to account for randomness in k-means, k-flat, 
% SCC, SGC, TeTrIS etc.)

clear all; close all; clc;

addpath otherSC
addpath SSC_ADMM_v1.1
addpath code2
addpath supp_material_RSCT/include
addpath scc
addpath tetris
addpath sgc

NumTrial = 20;
NumAlgo = 10;
CE = zeros(NumAlgo,1,NumTrial);
ET = zeros(NumAlgo,1,NumTrial);

load Hopkins155_titles.mat; i_data = 0;
for i=1:length(Hopkins155_titles)
    fprintf('%d/%d %s \n',i,length(Hopkins155_titles),Hopkins155_titles{i});
    eval(['load Hopkins155/' Hopkins155_titles{i} '/' Hopkins155_titles{i} '_truth.mat' ]);
    
    L = max(s);

    if (L == 2)
        
        i_data = i_data + 1;
        
        N = size(x,2);
        F = size(x,3);
        p = 2*F;
        Y = reshape(permute(x(1:2,:,:),[1 3 2]),p,N);
        
        A0 = s; 
        [~,I] = sort(A0,1);
        A0 = A0(I);
        Y = Y(:,I);
        
        for i_trial = 1: NumTrial
            fprintf('%d/%d: Trial %d\n',i,length(Hopkins155_titles),i_trial);
%% K-means
        fprintf('Running K-means..\n'); i_algo = 1;
        
        tic;
        [A,~] = kmeans(Y',L,'emptyaction','singleton','replicates',10,'display','off');
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);
        
        
%% K-flats
        fprintf('Running K-flats..\n'); i_algo = 2;
        
        CE(i_algo,i_data,i_trial)  = 1;
        for i_repl = 1:10
            tic;
            A = Kflats(Y,3,L);
            ET(i_algo,i_data,i_trial)  = toc * 10;
            CE(i_algo,i_data,i_trial)  = min(computeCE(A,A0), CE(i_algo,i_data,i_trial));
        end

%% SSC
        fprintf('Running SSC..\n'); i_algo = 3;
        r = 0; affine = true; outlier = false; rho = 0.7; alpha = 800;
        tic
        [Z,A] = SSC(Y,r,affine,alpha,outlier,rho,s);
        
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);
       
%% LRR
        fprintf('Running LRR..\n'); i_algo = 4;
        % The following scripts are copied from the LRR source code.
        tic
        lambda = 4;
        %run lrr
        Z = solve_lrr(Y,lambda);
        %post processing
        [U,S,V] = svd(Z,'econ');
        S = diag(S);
        r = sum(S>1e-4*S(1));
        U = U(:,1:r);S = S(1:r);
        U = U*diag(sqrt(S));
        %U = normr(U);
        U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));
        LL = (U*U').^4;
        % spectral clustering
        D = diag(1./sqrt(sum(LL,2)));
        LL = D*LL*D;
        [U,S,V] = svd(LL);
        V = U(:,1:L);
        V = D*V;
        A = kmeans(V,L,'emptyaction','singleton','replicates',10,'display','off');
        
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);

%% SCC
        fprintf('Running SCC..\n'); i_algo = 5;
        tic

        opts.c = L*100;
        [A,~] = scc(Y',3,L,opts);
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);        


%% GFS / SSC-OMP
        fprintf('Running GFS..\n'); i_algo = 6;
        tic
        Z = OMPSC(Y,8);
        A = SpectralClustering(abs(Z)+abs(Z)',L);
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);
        
%% TSC
        fprintf('Running TSC..\n'); i_algo = 7;
        tic
        [A,Z] = TSC(Y,10,L);
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);
         
%% NSN+Spectral
        fprintf('Running NSN+Spectral..\n'); i_algo = 8;
        tic
        Y1 = Y./repmat(sqrt(sum(Y.^2,1)),p,1);
        Z = NSN(Y1,5,5,1e-4);
        A = SpectralClusteringL(Z+Z',L);
            
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);

%% SGC
        fprintf('Running SGC..\n'); i_algo = 9;
        tic
        
        opts.c = L*100;
        [A,~] = sgc(Y',3,L,opts);
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);        
        

%% TETRIS
        fprintf('Running TETRIS..\n'); i_algo = 10;
        tic
        
        opts.c = L*100;
        [A,~,~] = tetris(Y',3,L,opts);
        ET(i_algo,i_data,i_trial)  = toc;
        CE(i_algo,i_data,i_trial)  = computeCE(A,A0);        

        end
    end
end

disp('   k-means    k-flats    SSC       LRR        SCC    SSC-OMP      TSC       NSN      SGC      TETRIS')
meanError = (mean(mean(CE,3),2))'
medianError = (median(mean(CE,3),2))'
maxMeanError = (max(mean(CE,3),[],2))'
meanTime = (mean(mean(ET,3),2))'

