function [] = em2tikzf(A,dvec,rho,filename)
%

fid=fopen(filename, 'w');


for i=1:length(dvec)
    for j=1:length(rho)
        str = sprintf('%d %d %d',rho(j),dvec(i), A(i,j));
        fprintf(fid,'%s\n',str);
        % disp(str);
    end
end


 fclose(fid);

end

