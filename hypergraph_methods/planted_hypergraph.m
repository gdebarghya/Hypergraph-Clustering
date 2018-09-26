function err = planted_hypergraph(s,m,k,p,q)

n = s*k;
true_labels = zeros(n,1);
for i = 1:k
    true_labels((i-1)*s + (1:s)) = i;
end

A = zeros(n*ones(1,m));
H = zeros(n,nchoosek(n,m));
Hcounter = 0;

for i1 = 1:n
    for i2 = i1:n
        % Case: m=2
        if (m==2)
            if (min(true_labels([i1 i2]))==max(true_labels([i1 i2])))
                prob = p;
            else
                prob = q;
            end
            if (rand(1)>prob)
                continue
            end
            Hcounter = Hcounter + 1;
            H([i1 i2],Hcounter) = 1;

            A(i1,i2) = 1;
            A(i2,i1) = 1;
            continue
        end
        
        for i3 = i2:n
            %Case: m=3
            if (m==3)
                if (min(true_labels([i1 i2 i3]))==max(true_labels([i1 i2 i3])))
                    prob = p;
                else
                    prob = q;
                end
                if (rand(1)>prob)
                    continue
                end
                Hcounter = Hcounter + 1;
                H([i1 i2 i3],Hcounter) = 1;

                A(i1,i2,i3) = 1;
                A(i1,i3,i2) = 1;
                A(i3,i1,i2) = 1;
                A(i3,i2,i1) = 1;
                A(i2,i3,i1) = 1;
                A(i2,i1,i3) = 1;
                continue
            end
            
            %Case: m=4
            for i4 = i3:n
                if (min(true_labels([i1 i2 i3 i4]))==max(true_labels([i1 i2 i3 i4])))
                    prob = p;
                else
                    prob = q;
                end
                if (rand(1)>prob)
                    continue
                end
                Hcounter = Hcounter + 1;
                H([i1 i2 i3 i4],Hcounter) = 1;

                A(i1,i2,i3,i4) = 1;
                A(i1,i2,i4,i3) = 1;
                A(i1,i3,i2,i4) = 1;
                A(i1,i3,i4,i2) = 1;
                A(i1,i4,i2,i3) = 1;
                A(i1,i4,i3,i2) = 1;
                A(i2,i1,i3,i4) = 1;
                A(i2,i1,i4,i3) = 1;
                A(i2,i3,i1,i4) = 1;
                A(i2,i3,i4,i1) = 1;
                A(i2,i4,i1,i3) = 1;
                A(i2,i4,i3,i1) = 1;
                A(i3,i1,i2,i4) = 1;
                A(i3,i1,i4,i2) = 1;
                A(i3,i2,i1,i4) = 1;
                A(i3,i2,i4,i1) = 1;
                A(i3,i4,i1,i2) = 1;
                A(i3,i4,i2,i1) = 1;
                A(i4,i1,i2,i3) = 1;
                A(i4,i1,i3,i2) = 1;
                A(i4,i2,i1,i3) = 1;
                A(i4,i2,i3,i1) = 1;
                A(i4,i3,i1,i2) = 1;
                A(i4,i3,i2,i1) = 1;
                
            end
        end
    end
end
H = H(:,1:Hcounter);

%%%%%% TTM (G&D, arxiv:1606.06516)
W = A;
for i = 1:(m-2)
    W = sum(W,(m+1-i));
end
D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
L = W.*dd;
[vec,val] = eig(L,'nobalance');
temp = sortrows([diag(val) vec'],-1);
evecs = temp(1:k,2:end)';
for i = 1:n
    if (norm(evecs(i,:))>0)
        evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
    end
end
idx = kmeans(evecs,k,'emptyaction','singleton');
err(1,1) = computeCE (idx,true_labels);


%%%%%% HOSVD (Govindu,CVPR-2005; G&D,NIPS-2014)
A = reshape(A,n,n^(m-1));           
W = A*A';
D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
L = W.*dd;
[vec,val] = eig(L,'nobalance');
temp = sortrows([diag(val) vec'],-1);
evecs = temp(1:k,2:end)';
for i = 1:n
    if (norm(evecs(i,:))>0)
        evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
    end
end
idx = kmeans(evecs,k,'emptyaction','singleton');
err(2,1) = computeCE (idx,true_labels);


%%%%%% NHCut (Zhou et la,NIPS-2007; G&D, AnnStat-2016)
W = H*H'/m;
D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
L = W.*dd;
[vec,val] = eig(L,'nobalance');
temp = sortrows([diag(val) vec'],-1);
evecs = temp(1:k,2:end)';
for i = 1:n
    if (norm(evecs(i,:))>0)
        evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
    end
end
idx = kmeans(evecs,k,'emptyaction','singleton');
err(3,1) = computeCE (idx,true_labels);
