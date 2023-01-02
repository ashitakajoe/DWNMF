function [MD_mat_new] = IWKNN( MD_mat, MM_mat, DD_mat, K1,K2, r )

[rows,cols]=size(MD_mat);
% [rows,~]=size(MD_mat);
% rows=sum(MD_mat(:));
% save rows rows;
% cols=2;
y_m=zeros(rows,cols);  
y_d=zeros(rows,cols);  

% [nc,nd]=size(MD_mat);
% interaction1=zeros(rows,2);
% x=1;
% for i=1:nc
%     for j=1:nd
%         if MD_mat(i,j)==1
%             interaction1(x,1)=i;
%             interaction1(x,2)=j;
%             x=x+1;
%         end
%     end
% end
% save interaction1 interaction1;


knn_network_m = KNN( MM_mat, K1 );  %for circRNA
for i = 1 : rows   
         w=zeros(1,K1);
        [sort_m,idx_m]=sort(knn_network_m(i,:),2,'descend'); 
        sum_m=sum(sort_m(1,1:K1));   
        for j = 1 : K1
            w(1,j)=r^(j-1)*sort_m(1,j); 
            y_m(i,:) =  y_m(i,:)+ w(1,j)* MD_mat(idx_m(1,j),:); 
        end                      
            y_m(i,:)=y_m(i,:)/sum_m;              
end

knn_network_d = KNN( DD_mat , K2 );  %for disease
for i = 1 : cols   
        w=zeros(1,K2);
        [sort_d,idx_d]=sort(knn_network_d(i,:),2,'descend');
        sum_d=sum(sort_d(1,1:K2));
        for j = 1 : K2
            w(1,j)=r^(j-1)*sort_d(1,j);
            y_d(:,i) =  y_d(:,i)+ w(1,j)* MD_mat(:,idx_d(1,j)); 
        end                      
            y_d(:,i)=y_d(:,i)/sum_d;               
end

a1=1;
a2=1;
y_md=(y_m*a1+y_d*a2)/(a1+a2);  

 for i = 1 : rows
     for j = 1 : cols
         MD_mat_new(i,j)=max(MD_mat(i,j),y_md(i,j));
     end    
 end

end

function [ knn_network ] = KNN( network , k )
    [rows, cols] = size( network );
    %v=diag(X,k),XΪ����vΪ����,ȡ����X�ĵ�K���Խ���Ԫ��Ϊ����v
    network= network-diag(diag(network)); 
    knn_network = zeros(rows, cols);
    %sort(A,dim,mode)��modeΪ'ascend'ʱ��������������modeΪ'descend'ʱ�����н�������.
    %�Ծ���ָ���ķ�������������򣬲����������ľ���
    %��dim=1ʱ���Ծ����ÿһ�����򣨼�����һά�����������ţ�����dim=2ʱ���Ծ����ÿһ�����򣨼����ڶ�ά�����������ţ���
    %Y�е�Ԫ�ض�Ӧ��ԭʼ����X�е���һ��Ԫ�ء��������ǿ�����������[Y,I] = sort(X,DIM,MODE)
    %I�����������У�����ʾY�е�Ԫ����X��Ԫ�صĶ�Ӧ
    [sort_network,idx]=sort(network,2,'descend');
    for i = 1 : rows
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);
    end
end


