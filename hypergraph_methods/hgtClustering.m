%% 
% This function clusters a weighted hypergraph H
% returning a number of clusters not exceeding max_k
% optionally restricting the number of iterations for
% single cluster estraction.
% 
% H: hypergraph structure
%    H.nVertices: number of vertices (N)
%    H.nEdges:    number of edges (M)
%    H.edges:     MxK matrix of indices to vertices representing the edges
%    H.w:         Mx1 vector edge weights
%
% max_k: maximum number of desired clusters
%
% max_iter: maximum number of iterations per cluster extraction

function [clustering,clusters]=hgtClustering(H,max_k,max_iter)
  if(~exist('max_iter','var'))
    max_iter=Inf;
  end
  
  %H=H.clone();
  n=H.nVertices;
  mask=zeros(n,1)>0;
  clustering=zeros(n,1);
  clusters=zeros(max_k,n);
  clustered=0;
  k=1;
  while sum(mask)<n && k<=max_k
    x=rand(n,1)+1000;
    x(mask)=0;
    x=x/sum(x);
    [x,Fx,niter]=hgtExtractCluster(H,x,max_iter);
    %disp(['Cluster ' num2str(k) 'extracted, niter=' num2str(niter)]);
    clusters(k,:)=x;
    den=x'*Fx;
    if den>0
      cluster=(Fx/den)>0.98;
      %cluster(x<(1E-4/(n-sum(mask))))=0;
    else
      cluster=~mask;
    end
    clustered=clustered+sum(cluster);
    clustering(cluster)=k;
    mask(cluster)=true;
    k=k+1;
    %REMOVE UNUSEFUL EDGES
    m=zeros(H.nEdges,1)>1;
    for i=1:H.nEdges
      m(i)=max(clustering(H.edges(i,:)))==0;
    end
%    H.removeEdges(m);
  H.edges=H.edges(m,:);
  H.nEdges=sum(m);
  H.w=H.w(m);
  end
  %clusters=clusters(1:k-1,:);
  %if(min(clusters)==0)
  %  clusters=clusters+1;
  %end
end
