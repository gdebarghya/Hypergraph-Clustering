function [A] = LabelPoints(Y,D)
    % For each column of Y, find the nearest subspace among the column
    % spaces of D{1}, ..., D{L}.
    %
    % Y : (p*N matrix) Data matrix whose columns are the data points
    % D : (1*L cells of p*d matrix) Orthogonal bases of L subspaces
    % A : (1*N vector) Output label

    L = numel(D);
    [p,N] = size(Y);
    
    proj = zeros(L,N);
    for i=1:L
        proj(i,:) = sum((D{i}'*Y).^2,1);
    end
    
    [~,A] = max(proj,[],1);
end