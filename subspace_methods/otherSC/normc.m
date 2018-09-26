function Xnorm = normc(X)
    Xnorm = X./repmat(sqrt(sum(X.^2,1)),size(X,1),1);
end