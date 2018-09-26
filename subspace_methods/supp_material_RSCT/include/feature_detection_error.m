% This function computs the feature detection error
% A: N x N adjacency matrix, where N = Ln, and where the 
% rows 1+i*N:(i+1)*N correspond to a cluster
%
% Reinhard Heckel, 2013


function [f] = feature_detection_error(A,L)

N = length(A);
n = floor(N/L);

f = 0;

for l = 0:L-1
    for i = 1:n
        cp = i + l*n;
        f = f + norm( A(cp ,(1+l*n):(n*(1+l))  ))/norm(A(cp,:));
    end
end

f = 1-f/N;

end
