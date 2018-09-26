clear;

%%%%%%%%%%%%%%%%%%
% performance of TSC when subspaces intersect
%%%%%%%%%%%%%%%%%%

addpath '../../include'

m = 50;    % ambient dimension 

% original parameters
%L = 8;      % number of subspaces 
%nit = 100;
%dvec = 5:2:25;
%nl = 10:10:150;

% to obtain a figure fast 
L = 4;      % number of subspaces 
nit = 1;
dvec = 5:5:25;
nl = 10:50:150;

fde = zeros(length(dvec),length(nl));
ce = zeros(length(dvec),length(nl));
Le = zeros(length(dvec),length(nl));
Ls = zeros(length(dvec),length(nl));

for j=1:length(dvec)  % vary over d
  
    j
    
    for k=1:length(nl)% vary over nl
     k   
        
        
        N = L*nl(k);  % number of points overall

        for i=1:nit
            
            X = [];
            Xlabels = [];
            for l=1:L
                U = orth(randn(m,dvec(j)));
                A = normc(randn(dvec(j),nl(k))); 
                X = [X, U*A];  
                Xlabels = [Xlabels, ones(1,nl(k))*l];
            end            
            
            size(X)
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q = 2*max(3,ceil(nl(k)/20));
            [labels,A,Lest] = TSC(X,q);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % compute error 
            fde(j,k) = fde(j,k) + feature_detection_error(A,L)/nit;
            ce(j,k) = ce(j,k) + clustering_error(labels',Xlabels)/nit;
            
            if(L ~= Lest) 
                Le(j,k) = Le(j,k) + 1/nit;   
                if L < Lest
                    Ls(j,k) = Ls(j,k) + 1/nit;
                else
                    Ls(j,k) = Ls(j,k) - 1/nit;
                end
            end
        end
        
    end
end

% 
em2tikzf(Le,dvec,nl,'./LE_varyd.dat');
em2tikzf(Ls,dvec,nl,'./LS_varyd.dat');
em2tikzf(fde,dvec,nl,'./FDE_varyd.dat');
em2tikzf(ce,dvec,nl,'./CE_varyd.dat');

