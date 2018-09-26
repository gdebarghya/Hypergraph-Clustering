% Hypergraph clustering
% 
% H: hypergraph structure
%    H.nVertices: number of vertices (N)
%    H.nEdges:    number of edges (M)
%    H.edges:     MxK matrix of indices to vertices representing the edges
%    H.w:         Mx1 vector edge weights
%
% x: Nx1 starting point in the simplex
%
% maxiter: maximum number of iterations

function [x,Fx,niter]=hgtExtractCluster(H,x,maxiter)
  if ~exist('maxiter','var')
    maxiter=1000;
  end
  niter=0;
  toll=1e-6;
  error=2*toll;
  old_x=x;
  while niter<maxiter && error>toll

   
    Fx=zeros(H.nVertices,1);
    for i=1:H.nEdges
      edge=H.edges(i,:);
      tmp=prod(x(edge))*H.w(i);
      if tmp>0
        Fx(edge)=Fx(edge)+tmp./x(edge);
      end 
    end
    x=x.*Fx;
    xFx=sum(x);
    if xFx==0
      return;
    end
    x=x/xFx;
    
    error=norm(x-old_x);
    old_x=x;
    
    niter=niter+1;
    %fprintf('iter=%d error=%f\n',niter,error);
  end
end

