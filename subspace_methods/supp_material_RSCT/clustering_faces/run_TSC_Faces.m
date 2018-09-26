%--------------------------------------------------------------------------
% This is an adaption of the file run_SSC_Faces.m by Ehsan Elhamifar 
% from the paper
% `` Sparse subspace clustering: Algorithm, theory, and applications''
%
% ``avg_err'' and ``avg_err_whitening'' contain the results 
%
%
%
% Reinhard Heckel, 2013


clear all, close all

addpath '../include/';
addpath '../SSC_ADMM_v1.1/';

load YaleBCrop025.mat

%nSet = [2 3 4 5 6 7 8 9 10];

nSet = [3];

for i = 1:length(nSet)
      
    n = nSet(i);
    idx = Ind{n};   
    for j = 1:size(idx,1)
        X = [];
        Xlabels = [];
        for p = 1:n
            X = [X Y(:,:,idx(j,p))];         
            Xlabels = [Xlabels,ones(1, size( Y(:,:,idx(j,p)) , 2))*p];          
            %X = [X Xo{idx(j,p)}];
            %Xlabels = [Xlabels,ones(1, size( Xo{idx(j,p)} , 2))*p];
        end
        [D,N] = size(X);
        
        X = normc(X); 
        
        % without whitening
        q = 3;
        [labels,Z] = TSC(X,q,n);
        err{n}(j) = clustering_error(Xlabels,labels);
        clustering_error(Xlabels,labels)
        
        % with whitening
        
        W = 2;
        X = normc(X);
        [U,S,V] = svd(X);
        A = U(:,1:W);
        X = X - (A*pinv(A))*X;
        X = normc(X);      
        
        [labels,Z] = TSC(X,q,n);
        err_whitening{n}(j) = clustering_error(Xlabels,labels);
        
        
    end
    avg_err(n) = mean(err{n})
    avg_err_whitening(n) = mean(err_whitening{n});
end

avg_err

avg_err_whitening

