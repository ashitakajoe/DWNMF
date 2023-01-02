%circRNA
%*******************************************************************************************************************
clear;

%*********************************************************************************************************************
load dataset.mat %585*88


newData1=cd_adjmat;%585*88    dataset中
Y = newData1; 

clear newData1
interaction=cd_adjmat;%585*88  
% nc:the number of circRNA
% nd:the number of disease
% pp:the number of known diseae-circRNA associations

[nc,nd]=size(interaction);%585*88  circRNA*disease
%interaction: adjacency matrix for the disease-circrna association network
%*******************************************************************************************************************
%Calculate the network topology similarity
[rows, cols]=size(cd_adjmat);          
 printEdgeList(cd_adjmat); 
        
dos('python __main__.py --input cd_edgelist --output cd_embeddings --workers 1 --representation-size 256');
        
cd_embeddings=textread('cd_embeddings','','headerlines',1);
cd_embeddings=sortrows(cd_embeddings,1);
%save cd_embeddings cd_embeddings;
dis_vectors=cd_embeddings(1:cols,2:257);
circ_vectors=cd_embeddings(cols+1:cols+rows,2:257);
deep_dis_sim=cosSim(dis_vectors);
deep_circ_sim=cosSim(circ_vectors);
%save deep_dis_sim deep_dis_sim;
%save deep_circ_sim deep_circ_sim; 

SD=deep_dis_sim;

SC=deep_circ_sim;

Sd=zeros(nd,nd);
 for i=1:nd
    for j=1:nd
       Sd(i,j) = a1*SD(i,j)+(1-a1)*disease_sim(i,j);
    end
 end 

%circRNA function similarity
circSim=miRNASS( interaction, disease_sim );



 Sc=zeros(nc,nc);
for i=1:nc
    for j=1:nc
       Sc(i,j) = a2*SC(i,j)+(1-a2)*circSim(i,j); 
    end
end 

%*******************************************************************************************************************


% 对 Y 进行预处理
new_Y_IWKNN = IWKNN( interaction, Sc, Sd, K1,K2, r );  

% Laplacian Matrices   
Dm = diag(sum(Sc));
Dd = diag(sum(Sd));
Lm = Dm - Sc;
Lm = (Dm^(-0.5))*Lm*(Dm^(-0.5));
Ld = Dd - Sd;
Ld = (Dd^(-0.5))*Ld*(Dd^(-0.5));

% initialize A & B
[u,s,v] = svds(new_Y_IWKNN,K);   
A = u*(s^0.5);
B = v*(s^0.5);
BB=sqrt(sum(B.*B,1) + 1e-6);  
BBA=1./BB;  
Dij = diag (BBA); 

AA=sqrt(sum(A.*A,1) + 1e-6);  
AAB=1./AA;  
Dij1 = diag (AAB); 

lambda_m_Lm = lambda_m*Lm;
lambda_d_Ld = lambda_d*Ld;
lambda_l_eye_K = lambda_l*eye(K);
lambda_l_eye_K1 = lambda_l*Dij*eye(K);
lambda_l_eye_K2 = lambda_l*Dij1*eye(K);
clear Dij BBA BB Dij1 AAB


for z=1:num_iter
     A = (new_Y_IWKNN*B  - lambda_m_Lm*A) / (B'*B +lambda_l*(B'*B)+lambda_l_eye_K2 );
     B = (new_Y_IWKNN'*A - lambda_d_Ld*B) / (A'*A +lambda_l*(A'*A)+lambda_l_eye_K1);
end 

% compute prediction matrix
newY = A*B';
%save newY newY;


% obtain the final prediction result, 1st column is the score, 2nd columan is the corresponding disease ID, 3rd is circRNA ID
finalprediction=[];
for i=1:nd
    finalprediction=[finalprediction;newY(:,i)];
end

finalpredictiondisease=[];
for i=1:nd
    finalpredictiondisease=[finalpredictiondisease,1:nc];
end
finalprediction(:,2)=finalpredictiondisease';

finalpredictionmiRNA=[];
for i=1:nd
    finalpredictionmiRNA=[finalpredictionmiRNA;i.*ones(nc,1)];
end
finalprediction(:,3)=finalpredictionmiRNA;

discard=finalprediction(find(finalprediction(:,1)==-10000),:);
finalprediction=setdiff(finalprediction,discard,'rows');
%sortrows mean ranking according to given column,here 1st column, -1 means
result=sortrows(finalprediction,-1);


%load the matrix for the correspondence relation between disease ID and names:

[a1,~,a]=xlsread('D:\MATLAB\DWNMF- master\diseasename-list.xlsx');
%load the matrix for the correspondence relation between circRNA ID and names: 
[b1,~,b]=xlsread('D:\MATLAB\DWNMF- master\circRNA-name-list.xlsx');

m=length(result);
%define two variables to store the correspondence relation of disease and circRNA, respectively:
c=zeros(m,1);
f=zeros(m,1);
%define two cells to store the result after the replacement:
d={};
e={};
%obtain the correspondence relation for circRNA ID in the third column of result matrix and their names:
for i=1:m
    if(~isempty(find(b1==result(i,2))))
        c(i,1)=find(b1==result(i,2));
        d{i,1}=b{c(i,1),2};
    else
        d{i,1}=result(i,2);
    end
end

%obtain the correspondence relation for disease in the third column of result matrix and their names: 
for i=1:m
    if(~isempty(find(a1==result(i,3))))
        f(i,1)=find(a1==result(i,3));
        e{i,1}=a{f(i,1),2};
    else
        e{i,1}=result(i,3);
    end
end

xlswrite('D:\MATLAB\DWNMF- master\final prediction result.xlsx',result((1:m),1),'A1:A50830');
xlswrite('D:\MATLAB\DWNMF- master\final prediction result.xlsx',d,'B1:B50830');
xlswrite('D:\MATLAB\DWNMF- master\final prediction result.xlsx',e,'C1:C50830');







