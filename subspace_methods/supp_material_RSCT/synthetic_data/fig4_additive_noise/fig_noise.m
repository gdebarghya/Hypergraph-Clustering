
% Code to generate figure 4

m = 50;    % ambient dimension 


%%% original parameters
% L = 8;
% nit = 10;
% dvec = 5:2:25;
% nl = 10:10:150;

%%% parameters to produce a figure fast
L = 2;      % number of subspaces 
nit = 1;
dvec = 5:10:25;
nl = 10:50:150;



for sig=0:0.1:0.8
	sig
	filename = ['CE_sig', num2str(sig),'.dat'];
	run_TSC_noise(L,m,nit,dvec,nl,sig,filename);
end


quit;

