%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               PLS ANALYSIS ON AEC FUNCTIONAL CONNECTIVITY
%       
%       Uses output AEC from Brainstorm to test contrast between groups
%       Conditions are all within participants 2x2 study design
%       CUED vs NOT CUED and Subjective report SEEN vs NOTSEEN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read data from files
fCUE=dir("/Users/jasondsc/Desktop/new_CNV/cue/0*");
fNOCUE=dir("/Users/jasondsc/Desktop/new_CNV/nocue/0*");
%% Set parameters for PLS
con=[fCUE, fNOCUE];
conditions = {'CUE','NOCUE'};
ncond = numel(conditions);
option.method = 1;
option.num_perm = 10000;
option.num_boot = 10000;
%option.stacked_designdata = [-1 -0.5 1 0.5]';

count=1;
%set number of participants, features, and freq band of choice
nsub=30;
nsubj = [nsub];
nfeat=19264;
data=zeros(nsub*ncond, nfeat);
eegdata=zeros(64,301,nsub*4);

for j=1:2
   
    for i=1:nsub
        load(strcat(con(i,j).folder, "/", con(i,j).name));
        disp(con(i,j).name)
        data(count,:)= reshape(avg(1:64,:).',1,[]);
        eegdata(:,:,count)=avg(1:64,:);
        count=count+1;
        
    end
    
end

datamat{1} = data;
result = pls_analysis(datamat,nsubj,ncond,option);

result.perm_result.sprob
figure; bar(result.v(:,1));
figure; bar(result.boot_result.orig_usc(:,1)); hold on; errorbar(result.boot_result.orig_usc(:,1),[result.boot_result.ulusc(:,1)-result.boot_result.llusc(:,1) ], '.');
figure; bar(result.boot_result.orig_usc(:,2)); hold on; errorbar(result.boot_result.orig_usc(:,2),[result.boot_result.ulusc(:,2)-result.boot_result.llusc(:,2) ], '.');

%figure; hist(result.u(:,1));
%% percent covariance explained 
(result.s(1)^2)/sum(result.s.^2)

load("~/Desktop/McGill/EEG_alertness/channel_mat.");
channelmat=table2cell(channelmat);
d2=result.boot_result.compare_u(:,1);
d2(abs(d2)<2.58)=0;
d2(d2>=2.58)=1;
d2(d2<=-2.58)=-1;
dd2=reshape(d2,301,64)';
figure; imagesc(dd2); 
set(gca,'Ytick',1:64)
set(gca,'YtickLabel',channelmat)
set(gca,'Xtick',1:10:301)
set(gca,'XtickLabel',eeg.Time(1:10:301))


d2=result.boot_result.compare_u(:,1);
dd2=reshape(d2,301,64)';
f=figure; imagesc(dd2); 
set(gca,'Ytick',1:64)
set(gca,'YtickLabel',channelmat)
set(gca,'Xtick',1:20:301)
set(gca,'XtickLabel',eeg.Time(1:20:301))





notseencue=reshape(mean(data(1:30,:),1)',301,64)';
notseennocue=reshape(mean(data(31:60,:),1)',301,64)';
seencue=reshape(mean(data(61:90,:),1)',301,64)';
seennocue=reshape(mean(data(91:120,:),1)',301,64)';

%VAN
plot(-1*mean(notseencue([14,15, 16,17,18,20,46,47,50,51,52],:))); hold on;
plot(-1*mean(seencue([14,15, 16,17,18,20,46,47,50,51,52],:))); hold on;
figure
plot(-1*mean(notseennocue([14,15, 16,17,18,20,46,47,50,51,52],:))); hold on;
plot(-1*mean(seennocue([14,15, 16,17,18,20,46,47,50,51,52],:))); hold on;
set(gca,'Xtick',1:20:301)
set(gca,'XtickLabel',eeg.Time(1:20:301))

plot(-1*(mean(seencue([14,15, 16,17,18,20,46,47,50,51,52],:))-mean(notseencue([14,15, 16,17,18,20,46,47,50,51,52],:))));
hold on;
plot(-1*(mean(seennocue([14,15, 16,17,18,20,46,47,50,51,52],:))-mean(notseennocue([14,15, 16,17,18,20,46,47,50,51,52],:))));

%% P100 posterior
figure; subplot(1,2,1)
plot(-1*mean(notseencue([18, 19, 20,53,51,52],:))); hold on;
plot(-1*mean(seencue([18, 19, 20,53,51,52],:))); hold on;
set(gca,'Xtick',1:20:301)
set(gca,'XtickLabel',eeg.Time(1:20:301))
subplot(1,2,2)
plot(-1*mean(notseennocue([18, 19, 20,53,51,52],:))); hold on;
plot(-1*mean(seennocue([18, 19, 20,53,51,52],:))); hold on;
set(gca,'Xtick',1:20:301)
set(gca,'XtickLabel',eeg.Time(1:20:301))

plot(-1*(mean(seencue([18, 19, 20,53,51,52],:))-mean(notseencue([18, 19, 20,53,51,52],:))));
hold on;
plot(-1*(mean(seennocue([18, 19, 20,53,51,52],:))-mean(notseennocue([18, 19, 20,53,51,52],:))));


load("~/Desktop/McGill/EEG_alertness/channel_mat.");
channelmat=table2cell(channelmat);

%% plot  for all electrodes 
load('~/Desktop/pls_analysis_eeg_alertness/channel_locs.mat');
cue=reshape(mean(data(1:30,:),1)',301,64)';
nocue=reshape(mean(data(31:60,:),1)',301,64)';
figure;
ch=reshape(channellocs', 1, []);
for j=1:numel(channellocs);
    
    if (ch(j) ~=0)
        i=ch(j);
        dp1=-1*cue(i,:);
        dp2=-1*nocue(i,:);
        subplot(9,11,j); plot(dp1(:,101:301),'Color', [0.9, 0.1, 1]); hold on; plot(dp2(:,101:301),'b');
        %axis([0 200 -2.5*10^-6 2.5*10^-6])
        title(channelmat{i})
        set(gca,'Xtick',[])
        set(gca,'Ytick',[])
        
    end
end


%% plot  for all electrodes 
load('~/Desktop/pls_analysis_eeg_alertness/channel_locs.mat');
notseencue=reshape(mean(data(1:30,:),1)',301,64)';
notseennocue=reshape(mean(data(31:60,:),1)',301,64)';
seencue=reshape(mean(data(61:90,:),1)',301,64)';
seennocue=reshape(mean(data(91:120,:),1)',301,64)';
figure;
ch=reshape(channellocs', 1, []);
for j=1:numel(channellocs);
    
    if (ch(j) ~=0)
        i=ch(j);
        dp1=-1*notseencue(i,:);
        dp2=-1*seencue(i,:);
        dp3=-1*notseennocue(i,:);
        dp4=-1*seennocue(i,:);
        subplot(9,11,j); plot((dp2(:,51:251)+dp4(:,51:251))/2,'Color', [0.9, 0.1, 1]); hold on; plot((dp1(:,51:251)+dp3(:,51:251))/2,'b');
        %axis([0 200 -2.5*10^-6 2.5*10^-6])
        title(channelmat{i})
        set(gca,'Xtick',[])
        set(gca,'Ytick',[])
        
    end
end
