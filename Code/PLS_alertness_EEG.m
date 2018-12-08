%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               PLS ANALYSIS ON AEC FUNCTIONAL CONNECTIVITY
%       
%       Uses output AEC from Brainstorm to test contrast between groups
%       Conditions are all within participants 2x2 study design
%       CUED vs NOT CUED and Subjective report SEEN vs NOTSEEN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read data from files
fNOTSEENCUE= dir("~/Desktop/AEC_alertness/NOTSEENCUE/0*");
fNOTSEENNOCUE= dir("~/Desktop/AEC_alertness/NOTSEENNOCUE/0*");
fSEENCUE= dir("~/Desktop/AEC_alertness/SEENCUE/0*");
fSEENNOCUE= dir("~/Desktop/AEC_alertness/SEENNOCUE/0*");
%% Set parameters for PLS
con=[fNOTSEENCUE, fNOTSEENNOCUE, fSEENCUE, fSEENNOCUE];
conditions = {'NOTSEENCUE','NOTSEENOCUE','SEENCUE','SEENNOCUE'};
ncond = numel(conditions);
option.method = 1;
option.num_perm = 5000;
option.num_boot = 5000;
%option.stacked_designdata = [-1 1 1 -1]';

count=1;
%set number of participants, features, and freq band of choice
nsub=30;
nsubj = [nsub];
nfeat=2080;
data=zeros(nsub*ncond, nfeat);
freqband= 1; % 1= delta; 2=theta; 3=alpha; 4=beta; 5=gamma;

for j=1:4
   
    for i=1:nsub
        load(strcat(con(i,j).folder, "/", con(i,j).name));
        disp(con(i,j).name)
        data(count,:)= squeeze(TF(:,1,freqband))';
        count=count+1;
        
    end
    
end

datamat{1} = data;
result = pls_analysis(datamat,nsubj,ncond,option);

result.perm_result.sprob
figure; bar(result.v(:,1));
%figure; hist(result.u(:,1));
%% percent covariance explained 
(result.s(1)^2)/sum(result.s.^2)

%%bootstrap ratio
ratio=result.boot_result.compare_u(:,1)./result.boot_result.u_se(:,1);


d=data(1:4:120,:);
d=mean(d,1);
n=64;
v=d;
A = triu(ones(n));
A(A~=0) = 1:n*(n+1)/2;
A = A + triu(A,1).';
A1 = reshape(v(A),n,n);
d=data(2:4:120,:);
d=mean(d,1);
n=64;
v=d;
A = triu(ones(n));
A(A~=0) = 1:n*(n+1)/2;
A = A + triu(A,1).';
A2 = reshape(v(A),n,n);
d=data(3:4:120,:);
d=mean(d,1);
n=64;
v=d;
A = triu(ones(n));
A(A~=0) = 1:n*(n+1)/2;
A = A + triu(A,1).';
A3 = reshape(v(A),n,n);
d=data(4:4:120,:);
d=mean(d,1);
n=64;
v=d;
A = triu(ones(n));
A(A~=0) = 1:n*(n+1)/2;
A = A + triu(A,1).';
A4 = reshape(v(A),n,n);
figure;subplot(2,2,1);imagesc(A1);subplot(2,2,2);imagesc(A2);subplot(2,2,3);imagesc(A3);subplot(2,2,4);imagesc(A4);

n=64;
v=result.u(:,1);
A = triu(ones(n));
A(A~=0) = 1:n*(n+1)/2;
A = A + triu(A,1).';
A = reshape(v(A),n,n);
figure;imagesc(A);