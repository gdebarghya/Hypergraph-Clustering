% 
% This function implements the TSC algorithm from the paper 
% ``Robust subspace clustering via thresholding'' by Reinhard Heckel and Helmut Bölcskei
% Reinhard Heckel, 2013
%
% X: m x N matrix of N data points
% q: input parameter of TSC
% L: number of clusters, optional. If not provided, L is estimated via the eigengap heuristic%
% labels: labels of the data points
% Z: adjacency matrix 
% nL: estimated number of clusters

function [labels,Z,nL] = TSC(X,q,L)

% normalize the data points
X = normc(X);

[m,N] = size(X);

Z = zeros(N,N);

for i=1:N
    corvec = abs(X'*X(:,i));
    corvec(i) = 0; % so TSC will not select it
    [el,order] = sort(corvec, 'descend');
    Z(i, order(1:q) ) = exp(-2*acos(el(1:q))); % better than squared arcsin
end

Z = Z + Z';

% (normalized) spectral clustering step
D = diag( 1./sqrt(sum(Z)+eps) );
Lap = speye(N) - D * Z * D;
[U,S,V] = svd(Lap);

%% estimate L, if not provided as input
if(nargin == 3)
    nL = L;
else    
    svals = diag(S);
    [ min_val , ind_min ] = min( diff( svals(1:end-1) ) ) ;
    nL = N - ind_min;
end
%%

V = V(:,N-nL+1:N);
V = normr(V);   % normalize rows

  
warning off;
maxiter = 1000;     % maximum number of iterations 
replicates = 200;   % number of replications
labels = kmeans(V,nL,'maxiter',maxiter,'replicates',replicates,'EmptyAction','singleton');



end


