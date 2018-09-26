%%
% This script reads the faces from the CroppedYale database, resizes the
% images to 48x42 pixel, and stores the vectorized images in the array Xo.
% Xo{i} is a matrix, with the colums corresponding to the vectorized images
% of the ith person.
%
% Reinhard Heckel, 2013


clear;
% faces to read, T subseteq 1,...,39
T = [1:13 15:39];

Nfaces = length(T);
N = 1.

Xo = cell(length(T),1);
for i=1:length(T)

    if T(i) < 10
        direct = ['./CroppedYale/yaleB0' num2str(T(i)) '/'];
    else
        direct = ['./CroppedYale/yaleB' num2str(T(i)) '/'];
    end
    files = dir([direct '*.pgm']);
    num_files = numel(files);
    T(i)
    n(i) = num_files;
    Xtmp = [];
    for k = 1:n(i)
        files(k).name
        face = imread([direct files(k).name]);
        face = imresize(face, [48 42]);
        vface = im2double(reshape(face,1,[]));
        Xtmp = [Xtmp, vface' ];
    end
    Xo{i} = Xtmp;
end

save Xo

