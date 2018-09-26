
%%%%%%%%%%%%%%%%%%
% performance of TSC when subspaces intersect
% and when a number of entries is errased 
%%%%%%%%%%%%%%%%%%

% m = 50;		% ambient dimension 
% L = 4;  		% number of subspaces 
% nit = 10;
% 
% del = 10; 	% number of deleltions
% frac = 1/3; 	% fraction of the subspaces which is shared..
% 
% dvec = 6:3:15;
% nl = 10:50:110;
% 
% filename = 'out.dat'
% 
% erasefrac(m,L,nit,del,frac,dvec,nl,filename)

function erasefrac(m,L,nit,del,frac,dvec,nl,filename)

addpath '../../include/'

fde = zeros(length(dvec),length(nl));
ce = zeros(length(dvec),length(nl));
Le = zeros(length(dvec),length(nl));

for j=1:length(dvec)  % vary over d
  j
    for k=1:length(nl)% vary over rho   
   k     
        pis = nl(k);
        
        N = L*pis;  % number of points overall
		
        for i=1:nit
    
            X = [];
            Xlabels = [];
            UC = orth(randn(m,dvec(j)*frac));  
            for l=1:L
                % orthonormal basis U consisting of the common part UC
                % and the individual part with 1-frac entries 
                U = orth([UC,randn(m,dvec(j)*(1-frac))] );
                A = normc(randn(dvec(j),pis)); 
                X = [X, U*A];  
                Xlabels = [Xlabels, ones(1,nl(k))*l];
            end

            % missing entries: We set randomly some to zero:
            for l = 1:L*pis
                p = randperm(m);
                X(p(1:del),l) = zeros(del,1);
            end

            %%%%%
            q = max(3,ceil(nl(k)/20));
            [labels,A,Lest] = TSC(X,q,L);
            %%%%%
            
            % compute error 
            ce(j,k) = ce(j,k) + clustering_error(labels',Xlabels)/nit;
        end
    end
end

em2tikzf(ce,dvec,nl,filename);

end
 







