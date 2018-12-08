%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               PLS ANALYSIS ON SPECTRAL POWER TO CUE
%       
%       Uses output AEC from Brainstorm to test contrast between groups
%       Conditions are all within participants 2x2 study design
%       CUED vs NOT CUED and Subjective report SEEN vs NOTSEEN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read data from files
fCUE= dir("/Users/jasondsc/Desktop/McGill/EEG_alertness/export_raw_spect/cue/0*");
fNOCUE= dir("/Users/jasondsc/Desktop/McGill/EEG_alertness/export_raw_spect/nocue/0*");
%% Set parameters for PLS
con=[fCUE, fNOCUE];
conditions = {'fCUE','fNOCUE'};
ncond = numel(conditions);
option.method = 1;
option.num_perm = 1000;
option.num_boot = 1000;
%option.stacked_designdata = [-1 -0.5 1 0.5]';

count=1;
%set number of participants, features, and freq band of choice
nsub=30;
nsubj = [nsub];
nfeat=3203200;
data=zeros(nsub*ncond, nfeat);
%eegdata=zeros(64,301,nsub*4);

for j=1:2
   
    for i=1:nsub
        load(strcat(con(i,j).folder, "/", con(i,j).name));
        disp(con(i,j).name)
        data(count,:)= reshape(TF(:,:,:),1,[]);
        %eegdata(:,:,count)=eeg.F(1:64,:);
        count=count+1;
        
    end
    
end

datamat{1} = data;
result = pls_analysis(datamat,nsubj,ncond,option);

result.perm_result.sprob
figure; bar(-1*result.v(:,1));
figure; bar([1],-1*result.v([1],1)); hold on;bar([2],-1*result.v([2],1)); legend( "cued", "uncued")
figure; bar(result.boot_result.orig_usc(:,1)); hold on; errorbar(1:2,result.boot_result.orig_usc(:,1)',result.boot_result.ulusc(:,1),result.boot_result.llusc(:,1), '.');

load("~/Desktop/McGill/EEG_alertness/channel_mat.");
channelmat=table2cell(channelmat);


%% plot unthreshold for all electrodes 
load('~/Desktop/pls_analysis_eeg_alertness/channel_locs.mat');
figure;
d1=result.boot_result.compare_u(:,1);
dd1=reshape(d1,64,1001,50);
ch=reshape(channellocs', 1, []);
for j=1:numel(channellocs);
    
    if (ch(j) ~=0)
        i=ch(j);
        dp1=squeeze((dd1(i,:,:)))';
        subplot(9,11,j); imagesc(-1*dp1(:,300:900));
        title(channelmat{i})
        set(gca,'Ydir','Normal')
        set(gca,'Xtick',[])
        set(gca,'Ytick',[])
        caxis([-4,4])
    end
end

%% plot threshold for all electrodes
figure;
d1=result.boot_result.compare_u(:,1);
dd1=reshape(d1,64,1001,50);
ch=reshape(channellocs', 1, []);
for j=1:numel(channellocs);
    
    if (ch(j) ~=0)
        i=ch(j);
        dp1=squeeze((dd1(i,:,:)))';
        dp1(abs(dp1)<2.58)=0;
        dp1(dp1>=2.58)=1;
        dp1(dp1<=-2.58)=-1;
        subplot(9,11,j); imagesc(-1*dp1(:,300:900));
         title(channelmat{i})
        set(gca,'Ydir','Normal')
        set(gca,'Xtick',[])
        set(gca,'Ytick',[])
        caxis([-1,1])
    end
end

%% plot thresholded maps frontal central posterior
dd1=reshape(d1,64,1001,50);
dp1=squeeze(mean(dd1([2,3, 7, 39, 29, 35, 36, 59 ,62 ,40, 30, 24, 25, 58],:,:)))';
dp1(abs(dp1)<2.58)=0;
dp1(dp1>=2.58)=1;
dp1(dp1<=-2.58)=-1;
figure; subplot(1,3,1); imagesc(-1*dp1(:,400:900)); hold on;
line([102,102], [0,51], 'Color', 'k','LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r','LineWidth',1.5);
title('Frontal Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:500)
set(gca,'XtickLabel',Time(400:80:900))

dp1=squeeze(mean(dd1([13,14,15,45,46,53,11,12,22,23,44,54,55],:,:)))';
dp1(abs(dp1)<2.58)=0;
dp1(dp1>=2.58)=1;
dp1(dp1<=-2.58)=-1;
subplot(1,3,2); imagesc(-1*dp1(:,400:900)); hold on;
line([102,102], [0,51], 'Color', 'k','LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r','LineWidth',1.5);
title('Central Parietal Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:500)
set(gca,'XtickLabel',Time(400:80:900))


dp1=squeeze(mean(dd1([16, 17, 18, 47 ,48 ,49, 50, 51, 19, 20, 52],:,:)))';
dp1(abs(dp1)<2.58)=0;
dp1(dp1>=2.58)=1;
dp1(dp1<=-2.58)=-1;
subplot(1,3,3); imagesc(-1*dp1(:,400:900)); hold on;
line([102,102], [0,51], 'Color', 'k', 'LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r', 'LineWidth',1.5);
title('Posterior Occipital Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:500)
set(gca,'XtickLabel',Time(400:80:900))


%% plot unthresholded maps frontal central posterior
dd1=reshape(d1,64,1001,50);
dp1=squeeze(mean(dd1([2,3, 7, 39, 29, 35, 36, 59 ,62 ,40, 30, 24, 25, 58],:,:)))';
figure; subplot(1,3,1); imagesc(-1*dp1(:,400:900)); hold on;
line([102,102], [0,51], 'Color', 'k','LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r','LineWidth',1.5);
title('Frontal Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:500)
set(gca,'XtickLabel',Time(400:80:900))

dp1=squeeze(mean(dd1([13,14,15,45,46,53,11,12,22,23,44,54,55],:,:)))';
subplot(1,3,2); imagesc(-1*dp1(:,400:900)); hold on;
line([102,102], [0,51], 'Color', 'k','LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r','LineWidth',1.5);
title('Central Parietal Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:500)
set(gca,'XtickLabel',Time(400:80:900))


dp1=squeeze(mean(dd1([16, 17, 18, 47 ,48 ,49, 50, 51, 19, 20, 52],:,:)))';
subplot(1,3,3); imagesc(-1*dp1(:,400:900)); hold on;
line([102,102], [0,51], 'Color', 'k', 'LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r', 'LineWidth',1.5);
title('Posterior Occipital Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:500)
set(gca,'XtickLabel',Time(400:80:900))


%% plot thresholded maps frontal central posterior
dd1=reshape(data,60,64,1001,50);
dd1=squeeze(mean(dd1(1:30,:,:,:),1));
dd2=reshape(data,60,64,1001,50);
dd2=squeeze(mean(dd2(31:60,:,:,:),1));
dd2=dd2-mean(dd2(:,376:500,:),2);
dd1=dd1-mean(dd1(:,376:500,:),2);

dp1=squeeze(mean(dd1([2,3, 7, 39, 29, 35, 36, 59 ,62 ,40, 30, 24, 25, 58],:,:)))';
dp2=squeeze(mean(dd2([2,3, 7, 39, 29, 35, 36, 59 ,62 ,40, 30, 24, 25, 58],:,:)))';
dp1=(dp1-dp2)*10^12;
figure; subplot(1,3,1); imagesc(dp1(2:50,400:700)); hold on;
line([102,102], [0,51], 'Color', 'k','LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r','LineWidth',1.5);
title('Frontal Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:300)
set(gca,'XtickLabel',Time(400:80:700))


dp1=squeeze(mean(dd1([13,14,15,45,46,53,11,12,22,23,44,54,55],:,:)))';
dp2=squeeze(mean(dd2([13,14,15,45,46,53,11,12,22,23,44,54,55],:,:)))';
dp1=(dp1-dp2)*10^12;
subplot(1,3,2); imagesc(dp1(2:50,400:700)); hold on;
line([102,102], [0,51], 'Color', 'k','LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r','LineWidth',1.5);
title('Central Parietal Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:300)
set(gca,'XtickLabel',Time(400:80:700))


dp1=squeeze(mean(dd1([16, 17, 18, 47 ,48 ,49, 50, 51, 19, 20, 52],:,:)))';
dp2=squeeze(mean(dd2([16, 17, 18, 47 ,48 ,49, 50, 51, 19, 20, 52],:,:)))';
dp1=(dp1-dp2)*10^12;
subplot(1,3,3); imagesc(dp1(2:50,400:700)); hold on;
line([102,102], [0,51], 'Color', 'k', 'LineWidth',1.5);
line([282,282], [0,51], 'Color', 'r', 'LineWidth',1.5);
title('Posterior Occipital Bundle');
set(gca,'Ydir','Normal')
set(gca,'Xtick',0:80:300)
set(gca,'XtickLabel',Time(400:80:700))

