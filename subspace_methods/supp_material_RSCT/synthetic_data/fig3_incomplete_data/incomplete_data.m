
% code to generate Figure 3

m = 50;    % ambient dimension 
L = 6;      % number of subspaces 
frac = 1/3; % fraction of the subspaces which is shared..


%%% original parameters 
% nit = 100;
% dvec = 3:3:21;
% nl = 10:10:150;

%%% to obtain a figure fast
nit = 1;
dvec = 3:9:21;
nl = 10:50:150;


filename = 'erasure_0.dat'
del = 0; % number of deleltions
erasefrac(m,L,nit,del,frac,dvec,nl,filename)

filename = 'erasure_5.dat'
del = 5; % number of deleltions
erasefrac(m,L,nit,del,frac,dvec,nl,filename)

filename = 'erasure_10.dat'
del = 10; % number of deleltions
erasefrac(m,L,nit,del,frac,dvec,nl,filename)

filename = 'erasure_15.dat'
del = 15; % number of deleltions
erasefrac(m,L,nit,del,frac,dvec,nl,filename)

filename = 'erasure_20.dat'
del = 20; % number of deleltions
erasefrac(m,L,nit,del,frac,dvec,nl,filename)


quit;

