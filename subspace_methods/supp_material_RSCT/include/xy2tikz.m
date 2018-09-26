function [] = xy2tikzf(a,b,filename)
fid=fopen(filename, 'w');
        
		str = sprintf('%d %d', a(i), b(i));
        disp(str);
        fprintf(fid,'%s\n',str);

for i=1:length(a)
        str = sprintf('(%d, %d)', a(i), b(i));
        disp(str);
        fprintf(fid,'%s\n',str);
end
fclose(fid);
end
