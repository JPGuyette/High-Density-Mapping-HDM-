function bar=getcolorbar(scale,height,varargin)
% type=hsv;
% fig=figure(100);
% bar=repmat(permute(jet(2000),[1,3,2]),[1,100,1]);
% imshow(bar)
% numlines=1+(scale(2)-scale(1))/increment;
% [x1,y1]=dataxy2figxy(100,0);
% [x2 y2]=dataxy2figxy(100,110);
% annotation('line',[x1,y1],[x2 y2],'LineWidth',3)
% text(120,0,num2str(scale(2)),'FontSize',16,'Color','k')
% close all
% fig=figure('Position',[0 0 1200 1200]);
% imagesc(scale);
% h=colorbar;
% set(h,'FontSize',18)
% pause(1)
% x=screencapture(fig);
% bar=x(138:1001,1035:1150,:);
% bar=zeros(100,100,3);
% %close(fig)
% bar=double(bar)/255;
% bar=imresize(bar,height/size(bar,1));

if nargin==3
    mapname=varargin{1};
else
    mapname='jet';
end
fullout=ones(2040,400,3);
bar=repmat(permute(jtfheatmap(mapname,2000),[1,3,2]),[1,100,1]);
fullout(21:2020,1:100,:)=bar;
inc=(scale(2)-scale(1))/10;
for i=1:11
    fullout((18:24)+((i-1)*200),1:110,:)=0;
    curnum=scale(1)+(i-1)*inc;
    txtnum=num2str(curnum);
    if curnum>=0;
        txtnum=[' ',txtnum]; %#ok<AGROW>
    end
    texter=text2im(txtnum);
    texter=repmat(texter,[1 1 3]);
    if size(texter,2)>285
        texter=texter(:,1:185,:);
    end
    texter=imresize(texter,2);
    fullout((1:40)+((i-1)*200),115:(115+size(texter,2)-1),:)=texter;
end
bar=imresize(fullout,height/2040);
bar(bar>1)=1;
bar(bar<0)=0;
end

