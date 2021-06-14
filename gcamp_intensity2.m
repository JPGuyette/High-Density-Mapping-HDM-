function gcamp_intensity2(dataset,pos,frames)
chkdir
    top_dir=[dataset,'\'];
     % Get absolute image locations in order of location priority
    mkimageloc(dataset);
    load([top_dir,'files\fileloc_tif.mat'])

% Load image stack
%curim=[0;0;0;0];

curim=imread(fileloc.paths{1});
getpos=0;
if ~exist('pos','var')
    getpos=1;
elseif isempty(pos)
    getpos=1;
end
if getpos == 1;
    fig=figure;
    imshow(curim);
    h=impoly;
    pos=wait(h);
    close(fig)
end


[h,w,z]=size(curim);z;

% Make target folders
mkdir([top_dir,'intensitymaps']);

% get images and average intensities
mask=poly2mask(pos(:,1),pos(:,2),h,w);

n = length(fileloc.paths);
imstack = zeros(h,w,z,'uint8');
tic;
for i = 1:n
    curim=imread(fileloc.paths{i});
    curim(~mask)=NaN;
    imstack(:,:,i)=curim;
end
toc;

% Get peak values and select cycles
avgint=squeeze(nanmean(nanmean(imstack)));
[maxi,mini]=peakfind1(avgint,1);
[f,x]=pfft(avgint,1);
maxlist=find(maxi);
minlist=find(mini);
maxlist2=maxlist;
minlist2=minlist;

if minlist2(1)>maxlist2(1)
    maxlist2(1)=[];
end
if minlist2(end)>maxlist2(end)
   minlist2(end)=[]; 
end 



% calculate regional minima and minima framemaps
qtcyc=round(mean(maxlist2-minlist2)/2); % Length of a quarter cycle
ncyc=length(minlist2);
minstack=zeros(size(imstack,1),size(imstack,2),ncyc);
maxstack=minstack;
minframes=uint16(minstack);
maxframes=minframes;

for i = 1:ncyc
    tempstack=double(imstack(:,:,...
        max(minlist2(i)-qtcyc,1):min(minlist2(i)+qtcyc,n)));
    [minstack(:,:,i),minframes(:,:,i)]=min(tempstack,[],3);
    tempstack=double(imstack(:,:,...
        max(maxlist2(i)-qtcyc,1):min(maxlist2(i)+qtcyc,n)));
   [maxstack(:,:,i),maxframes(:,:,i)]=max(tempstack,[],3);
    
end
clearvars tempstack



% Calculate intensity heatmaps
perstack=zeros(size(maxstack));
for i= 1:ncyc
   perinc=(maxstack(:,:,i)-minstack(:,:,i)); %<<---Calculation for heatmaps
   perinc(~mask)=NaN;
   perinc=perinc/255;
   perstack(:,:,i)=perinc;
%    tofile=[ind2rgb(uint8(perstack(:,:,i)*255),jet),...
%         getcolorbar([0,1],size(curim,1),'jet')];
    q=ind2rgb(uint8(perstack(:,:,i)*255),jet);
    w=getcolorbar([0,1],size(curim,1),'jet');
    tofile=[q,w];
   imwrite(tofile,[top_dir,...
       'intensitymaps\cycle_',num2str(i),'.tif'])
end
    tofile=[ind2rgb(uint8(mean(perstack(:,:,i),3)*255),jet),...
        getcolorbar([0,1],size(curim,1),'jet')];
   imwrite(tofile,[top_dir,...
       'intensitymaps\cycle_AVG.tif'])
       


% Calculate average intensity curves
avgint=zeros(n,1);
int=zeros(n,1);
j=1;
curmin=zeros(size(curim),'uint8');
for i =1:n
    curim=imstack(:,:,i);
    if i==minlist2(j);
        curmin=curim;
        j=j+1;
    end
    
    if j>ncyc
        break
    end   
    avgint(i)=mean(mean((curim),1),2);%<<---Calculation for intensity curve (averaged over the selected region)
    int(i)=mean(mean(curim,1),2);
end
xcel=cell(n+1,2);
xcel(1,:)={'Frame','Average intensity change from baseline'};
xcel(2:end,1)=num2cell((1:n)');
xcel(2:end,2)=num2cell(avgint);
xlswrite([top_dir,'intensity.xlsx'],xcel);
save([top_dir,'intensitylocation.mat'],'dataset','mask','pos')


% peakvals = max(imstack,[],3);
% minvals = min(imstack,[],3);
% avgvals= mean(imstack,[],3);
% F=(avgvals)./minvals;
% imshow(avgvals)
% 
% Fstack = imstack-repmat(minvals,[1 1 n]);
% y=Fstack(:,:,frames);%to pick a certain amount of frames
% [~,Ftime] = max(y,[],3);
% imagesc(Ftime)
% K = filter2(fspecial('average',9),Ftime);
% imagesc(K)
%plot mean
% plot(squeeze(nanmean(nanmean(Fstack))))
% avg=squeeze(nanmean(nanmean(Fstack)));
% xlswrite('fluo.xlsx',[(1:n)',(avg)])%writes to excel


% make a figure of time mapC
% figure
% imagesc(Ftime)
% colormap(jet)

%save time map image
% outdir=[dataset,'\Ftime'];
% mkdir(outdir);
% Fscale=Ftime/n; %scales F time to the number of frames in set
% Fscale8=uint8(Fscale*255); %make Fscale into a 8 bit image
% Frgb = ind2rgb(Fscale8,jet(255)); %turns Fscale8 into color
% cbar=getcolorbar([0 n],size(Fscale8,1),'jet'); %makes color bar and scales it to Fscale8 height
% imwrite([Frgb,cbar],[outdir,'\timemap.tif']) %writes and saves 