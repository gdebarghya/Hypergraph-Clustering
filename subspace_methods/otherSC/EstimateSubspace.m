function [D] = EstimateSubspace(Y,A,d,L)
    % Estimate L d-dimensional subspaces representing the data matrix
    % For each label, find the d-dim approximation of the columns of Y with
    % the label
    %
    % Y : (p*N matrix) Data matrix whose columns are the data points
    % A : (1*N vector) Labels of data points
    % D : (1*L cells of p*d matrix) Orthogonal bases of L subspaces
    
    D = cell(1,L);
    for i=1:L
        M = Y(:,find(A == i));
        [U,~,~] = svds(M,d);
        D{i} = U;
        if size(U,2) < d
            D{i} = [D{i} zeros(size(U,1),d-size(U,2))];
        end
    end
end

