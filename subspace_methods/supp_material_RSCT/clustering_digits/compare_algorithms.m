%%%%%%%%%%%%%%
% This script compares the performance of TSCa and TSCb to SSC
% Reinhard Heckel, 2013
%%%%%%%%%%%%%%

clear; 

%% Problem parameters

% the number of digits in each subspace, n, over which is varied
%nSet = [25,50,75,100,125,150,175,200,225,250,275,300,325,350,375];
nSet = [25,100,175]; % to run faster

% the number of iterations 
numit = 10;

% the actual digits 
pl = [3,5,9]; % is actuall 2,4,8

addpath '../include/'
addpath '../SSC_ADMM_v1.1/'


L = length(pl);
%nSet = [150];

% for reproducable results, seed the random number generator
s = RandStream('mcg16807','Seed',200);
RandStream.setGlobalStream(s);

%%%%%% READ IMAGES
%% test images
images = loadMNISTImages('t10k-images.idx3-ubyte');
labels = loadMNISTLabels('t10k-labels.idx1-ubyte');

[labelssorted,IX] = sort(labels);
imgssorted = images(:,IX);

% beg, endd contain the indices of the begin (end) indices of the numbers
beg(1) = 1;
k = 1;
beg(k) = 1; 
for i =1:size(images,2) % for each point
    if labelssorted(i) == k-1;
    else
        endd(k) = i-1;
        k = k+1;
        beg(k) = i;
    end
end
endd(k) = size(images,2);
%%%%%% IMAGES READ

for i=1:length(nSet)
     n = nSet(i);
     for it=1:numit
         
         N = n*L;
         
		 % generate a problem instance
         T = [];
         Xlabels = [];
         X = [];
         for l=1:L
             U = beg(pl(l)):endd(pl(l));
             if(numit == 1)
                V = U(1:n); % take the first n images of each number
             else
                pn = randperm(length(U));
                V = U(pn(1:n)); % n random indices drawn without replacement from
             end
             
             Xlabels = [Xlabels, ones(1,n)*l];
             X = [X, imgssorted(:,V)];
         end
        X = normc(X);
       	
		q = max(3,ceil(n/20));
		% TSC 
		[labelTSC,A,nL] = TSC(X,q,L);
        ce_TSC(it) = clustering_error(Xlabels,labelTSC);
        
		% SSC
        r = 0; affine = false; outlier = false; rho = 1; alpha = 20;
        [ce_SSC(it),C] = SSC(X,r,affine,alpha,outlier,rho,Xlabels);

		% TSCb 
        [labelTSCb,A,nL] = TSCb(X,q,L);
        ce_TSCb(it) = clustering_error(Xlabels,labelTSCb);

     end
     
     ce_mean_TSC(i) = mean(ce_TSC);
     ce_mean_TSCb(i) = mean(ce_TSCb);
     ce_mean_SSC(i) = mean(ce_SSC);
    
     ce_std_TSCa(i) = std(ce_TSC);
	 ce_std_TSCb(i) = std(ce_TSCb);
	 ce_std_SSC(i) = std(ce_SSC);
       
end

dlmwrite('./TSCa.dat', [nSet',ce_mean_TSC',ce_std_TSCa'],'delimiter',' ');
dlmwrite('./TSCb.dat', [nSet',ce_mean_TSCb',ce_std_TSCb'],'delimiter',' ');
dlmwrite('./SSC.dat', [nSet',ce_mean_SSC',ce_std_SSC'],'delimiter',' ');
	
exit;
