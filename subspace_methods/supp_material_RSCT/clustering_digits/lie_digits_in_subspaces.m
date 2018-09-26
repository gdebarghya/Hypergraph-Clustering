%%%%%%%%%%%%%%
% Singular values of the matrices with columns consisting of vectorized images of a single digit
% Reinhard Heckel, 2013
%%%%%%%%%%%%%%

clear;

filename = 'sval.dat';

d = 20;

%%%%%% READ IMAGES
%% test images
images = loadMNISTImages('t10k-images.idx3-ubyte');
labels = loadMNISTLabels('t10k-labels.idx1-ubyte');
labels = labels + 1; % so the labels start from 1 and not from 0
[labelssorted,IX] = sort(labels);
imgssorted = images(:,IX);


% beg, endd contain the indices of the begin (end) indices of the numbers
beg(1) = 1;
k = 1;
beg(k) = 1; 
for i =1:length(labelssorted) % for each point
    if labelssorted(i) == k
    else
        endd(k) = i-1;
        k = k+1;
        beg(k) = i;
    end
end
endd(k) = length(labelssorted);
%%%%%% IMAGES READ


data(1:200,1) = 1:200;

for l=1:10
    l
    T = beg(l):endd(l);
    X = imgssorted(:,T);
    X = normc(X);
    [U{l},S,V] = svd(X);
    U{l} = U{l}(:,1:d);
    s = diag(S);
    data(1:200,l+1) = s(1:200);
    
end

dlmwrite(filename, data,'delimiter',' ');

exit;
