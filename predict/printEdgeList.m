function printEdgeList( cd_adjmat )
%PRINTEDGELIST 此处显示有关此函数的摘要
%   此处显示详细说明



[rows, cols]=size(cd_adjmat);

fid=fopen (  'cd_edgelist', 'w');

for i=1:cols
    
    fprintf(fid, '%d', i);
   
    for j=1:rows
        
       if cd_adjmat(j,i)==1
          
           fprintf(fid, ' %d', j+cols);
           
       end
        
    end
    fprintf(fid, '\n');
   
end

for i=1:rows
    
    fprintf(fid, '%d', i+cols);
   
    for j=1:cols
        
       if cd_adjmat(i,j)==1
          
           fprintf(fid, ' %d', j);
           
       end
        
    end
    fprintf(fid, '\n');
   
end



fclose(fid);


end

