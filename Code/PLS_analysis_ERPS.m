%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               PLS ANALYSIS ON AEC FUNCTIONAL CONNECTIVITY
%       
%       Uses output AEC from Brainstorm to test contrast between groups
%       Conditions are all within participants 2x2 study design
%       CUED vs NOT CUED and Subjective report SEEN vs NOTSEEN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read data from files
fNOTSEENCUE= dir("~/Desktop/export_data_brainstorm/NOTSEENCUE/0*");
fNOTSEENNOCUE= dir("~/Desktop/export_data_brainstorm/NOTSEENNOCUE/0*");
fSEENCUE= dir("~/Desktop/export_data_brainstorm/SEENCUE/0*");
fSEENNOCUE= dir("~/Desktop/export_data_brainstorm/SEENNOCUE/0*");
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
nfeat=19264;
data=zeros(nsub*ncond, nfeat);
eegdata=zeros(64,301,nsub*4);

for j=1:4
   
    for i=1:nsub
        load(strcat(con(i,j).folder, "/", con(i,j).name));
        disp(con(i,j).name)
        data(count,:)= reshape(eeg.F(1:64,:).',1,[]);
        eegdata(:,:,count)=eeg.F(1:64,:);
        count=count+1;
        
    end
    
end

datamat{1} = data;
result = pls_analysis(datamat,nsubj,ncond,option);

result.perm_result.sprob
figure; bar(result.v(:,1));
figure; hist(result.u(:,1));
%% percent covariance explained 
(result.s(1)^2)/sum(result.s.^2)

%%bootstrap ratio
ratio=result.boot_result.compare_u(:,1)./result.boot_result.u_se(:,1);

% d=mean(data(1:4:120,:),1);
% dd=reshape(d,301,64)';
% 
% d2=mean(data(3:4:120,:),1);
% dd2=reshape(d2,301,64)';
% 
% cue=dd2-dd;
% 
% plot(cue(17,:))
% gg=mean(eegdata(:,:,1:4:120),3);
% gg2=mean(eegdata(:,:,3:4:120),3);
% cued=gg2-gg;
% figure;plot(cued(17,:))
% 
% d=mean(data(1:4:120,:),1);
% dd=reshape(d,301,64)';

load("~/Desktop/McGill/EEG_alertness/channel_mat.");
channelmat=table2cell(channelmat);
d2=result.u(:,1);
dd2=reshape(d2,301,64)';
figure; imagesc(dd2);
set(gca,'Ytick',1:64)
set(gca,'YtickLabel',channelmat)
set(gca,'Xtick',1:10:301)
set(gca,'XtickLabel',eeg.Time(1:10:301))


d2=result.u(:,2);
dd2=reshape(d2,301,64)';
figure; imagesc(dd2); 
set(gca,'Ytick',1:64)
set(gca,'YtickLabel',channelmat)
set(gca,'Xtick',1:10:301)
set(gca,'XtickLabel',eeg.Time(1:10:301))
