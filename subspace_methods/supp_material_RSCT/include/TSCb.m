

function [labels,Z,nL] = TSCb(X,q,L)

% normalize
% X = normc(X);

[m,N] = size(X);

Z = zeros(N,N);

for i=1:N
    corvec = abs(X'*X(:,i));
    corvec(i) = 0; % so TSC will not select it
    [el,order] = sort(corvec, 'descend');
    %Z(i, order(1:q) ) = exp(-2*acos(el(1:q))); % better than squared arcsin
    Z(i, order(1:q) ) = abs(pinv(X(:, order(1:q) ))*X(:,i));
end

Z = Z + Z';

% (normalized) spectral clustering step






D = diag( 1./sqrt(sum(Z)+eps) );
Lap = speye(N) - D * Z * D;
[U,S,V] = svd(Lap);

%% estimate L
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
maxiter = 1000;     % Maximum number of iterations for KMeans 
replicates = 200;       % Number of replications for KMeans
labels = kmeans(V,nL,'maxiter',maxiter,'replicates',replicates,'EmptyAction','singleton');



end


