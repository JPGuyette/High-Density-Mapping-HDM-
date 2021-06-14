function split_tiff_stack_sa(prefix,rawdatafile,firstfilenum,brighten)
% Compiled by ADMIN\kjhansen on 01/07/2017 13:55:41
% split_tiff_stack Version 2015.04.0002
% 
% split_tiff_stack Splits a stack of tiff into the more universally
% accessible sequence of tiffs. Can be run without inputs.
%
% INPUTS:
%   prefix: prefix for output files
%   rawdatafile: Path to tiss file containing stack. Leave empty to ask
%       user for file
%   firstfilenum: number for first file

% Process rawdatafile input

if ~exist('rawdatafile','var')
    rawdatafile=[];
end
if isempty(rawdatafile)
    [rawdatafile,dirin]=uigetfile({'*.tif','*.tiff'},'Select tiff stack',...
        '','MultiSelect','on');
    if ischar(rawdatafile)
        rawdatafile=[dirin,rawdatafile];
    else
        for i=1:length(rawdatafile)
            rawdatafile{i}=[dirin,rawdatafile{i}];
        end
    end
end

if iscell(rawdatafile)
    N=length(rawdatafile);
    rdf=rawdatafile;
else
    N=1;
    rdf{1}=rawdatafile;
end
% Get directory
tmp=rdf{1};
d=strfind(tmp,'\');
d=d(end);
D=tmp(1:d);
% Process firstfilenum input
if ~exist('firstfilenum','var');
    firstfilenum=1;
end
if isempty(firstfilenum);
    firstfilenum=1;
end

cnt=zeros(N,1);
for i=1:N
    t=imfinfo(rdf{i});
    cnt(i)=length(t);
end

nums=firstfilenum:1:firstfilenum+sum(cnt)-1;
% Write individual tiffs
tstart=now();
bar=waitbar(0,'Calculating time remaining...');
tot=sum(cnt);
if ~exist('brighten','var')
    brighten=[];
end
if ischar(brighten) 
    if strcmpi(brighten,'imagej')
       x=t(1).ImageDescription;
       x1=strfind(x,'min=');
       x2=strfind(x,'max=');
       brighten=zeros(2,1);
       brighten(1)=str2double(x(x1+4:x2-1));
       brighten(2)=str2double(x(x2+4:end));
       bd=t(1).BitDepth;
       brighten=brighten./(2^bd);
    end
    
end
    
for i=1:N
    for j=1:cnt(i)
        if i>1
            k=sum(cnt(1:i-1))+j;
        else
            k=j;
        end
        curim=imread(rdf{i},j);
        curim=doubleim(curim);
        if ~isempty(brighten)
            curim=imadjust(curim,brighten,[0 1]);            
        end
        imwrite(curim,[D,prefix,'_',numpad(nums(k),4),'.tif']);
        trem=(k/(tstart-now())*(tot-k));
        waitbar(k/tot,bar,['Time remaining: ',...
            datestr(trem,'hh:mm:ss')])
    end
end

% Move source stacks
if ~exist([D,'srcstack'],'dir')
    mkdir([D,'srcstack'])
end
for i=1:N
    movefile(rdf{i},[D,'srcstack\',prefix,...
        '_stack_',numpad(i,2),'.tif'])
end  
end

function outim=doubleim(inputim)
%doubleim Converts an integer class image to a double class image

    if isa(inputim,'integer')
       cur=class(inputim);
       x=strfind(cur,'int');
       n=str2double(cur(x+3:end));
       outim=double(inputim)/(2^n);
    else

        switch class(inputim)
            case {'logical','single'}
                outim=double(inputim);
            case 'double'
                outim=inputim;                
            otherwise
                error('The specified variable is not a recognized image class')
        end
        %Scale image if necessary
        if min(outim(:))<0
            outim=outim-min(outim(:));
            warning('hdm2:imbelow0',['The specified image has',...
                ' values below 0...scaling image'])
        end
        if max(outim(:))>1
            outim=outim/max(outim(:));
            warning('hdm2:imabove1',['The specified image has',...
                ' values above 1...scaling image'])
        end
            
    end
end

function num = numpad(num,len)
    if isnumeric(num)
        num = num2str(num);
    end
    
    if length(num) < len
        num = numpad(['0',num],len);
    end
end

