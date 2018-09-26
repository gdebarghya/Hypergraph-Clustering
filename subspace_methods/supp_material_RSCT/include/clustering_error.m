% This function computes the clustering error, defined as:
% 
% err(map) = 1 - sum_i delta(origlabel(i), map(label(i)))/N
%
% here map() is apermutation map for the cluster indices. The clustering
% error is defined as the minimum of err(map) over all possible permutation maps.
% The optimal matching is found via the Hungarian algorithm.
%
% the entries in origlabel and label must be in 1:L
%
% Reinhard Heckel, 2013

function [err] = clustering_error(origlabel,label)

% construct the contingency table between the two labels, which contains
% the number of cluster label co-occurences 

N = length(label);

[labelssorted, ind] = sort(origlabel);

origlabel = origlabel(ind); % now the original labels are sorted
label = label(ind);

%% beg, endd will contain the indices of the begin (end) indices of the respective cluster
beg(1) = 1;
k = 1;
beg(k) = 1; 

for i =1:N % for each point
    if origlabel(i) == k
    else
        endd(k) = i-1;
        k = k+1;
        beg(k) = i;
    end
end
endd(k) = N;

%%
L = length(beg);
Q = length(unique(label));

A = zeros(L,Q);

for l = 1:L
    T = beg(l):endd(l); % the indices of coeff within the same cluster
    
    % construct a row of K, 
    S = label(T);
    unv = unique(S);
    his = histc(S,unv); 
    
    % labelssorted(T(1)) is the current row
    A(labelssorted(T(1)), unv) = his;
    
end

A = -A;
[Matching,Cost] = Hungarian(A);  % find the maximal weighted matching of the bipartite Graph defined by A;
                                 % the rows correspond to the nodes on one
                                 % side, the columns to the nodes on the
                                 % other side
                                 
% ``Matching'' contains the optimal maping 
% construct the vector containing the matching from ``Matching''
for l=1:Q
    if(find(Matching(:,l) == 1))
        matchv(l) = find(Matching(:,l)==1);
    else
        matchv(l) = l;
    end
end

for i = 1:N
    labelreordered(i) = matchv(label(i));
end

err = sum(origlabel ~= labelreordered)/N;

end
