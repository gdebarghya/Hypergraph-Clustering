%%%%%%%%%%%%%%%%%%
% performance of TSC when subspaces intersect
%%%%%%%%%%%%%%%%%%

clear;


addpath('../../include/')

% for reproducable results, seed the random number generator
s = RandStream('mcg16807','Seed',200);
RandStream.setGlobalStream(s);


L = 2;      % number of subspaces 
n = 200;    % ambient dimension 
d = 10;     % dimension of the subspaces
q = 15; 
pis = 20*d; % points in each subspace
N = L*pis;  % number of points overall
X = [];
nit = 2;


fde = zeros(d+1,1);
ce = zeros(d+1,1);
Le = zeros(d+1,1);
Ls = zeros(d+1,1);

for s=1:d
    
s
    for i=1:nit
        
        U1 = orth(randn(n,d)); % random orthonormal columns
        U2 = orth(randn(n,d-s)); % random orthonormal columns
        U2 = [U2, U1( : , 1:s )];

        % U2 and U1 are random subspaces, that intersect in at least s dimensions

        A1 = normc(randn(d,pis)); % random mtx, with columns uniformly on the unit sphere of reals^d, the coefficients of the points in the corresponding subspace
        A2 = normc(randn(d,pis)); % random mtx, with columns uniformly on the unit sphere of reals^d, the coefficients of the points in the corresponding subspace

        X = [U1*A1, U2*A2];  % without random mean
		X = normc(X);
        
		Xlabels = [ones(1,pis)*1, ones(1,pis)*2];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [labels,A,Lest] = TSC(X,8);   % clustering
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % compute error 
        fde(s+1) = fde(s+1) + feature_detection_error(A,L)/nit;
        ce(s+1) = ce(s+1) + clustering_error(Xlabels,labels')/nit;
		
        
        
        
        j = s+1;
        if(L ~= Lest) 
            Le(j) = Le(j) + 1/nit;   
            if L < Lest
                Ls(j) = Ls(j) + 1/nit;
            else
                Ls(j) = Ls(j) - 1/nit;
            end
        end
    end
end



t = (0:10)';

dlmwrite('./LS_SSintersect_TSC.dat', [t,Ls],'delimiter',' ');
dlmwrite('./LE_SSintersect_TSC.dat', [t,Le],'delimiter',' ');
dlmwrite('./FDE_SSintersect_TSC.dat', [t,fde],'delimiter',' ');
dlmwrite('./CE_SSintersect_TSC.dat', [t,ce],'delimiter',' ');

%%%%%%% 
