function varargout=mkfileloc(dataset,ext,varargin)
if nargin==3
    tardir=varargin{1};
end
chkdir
    
    
top_dir=[dataset,'\'];
    if exist([top_dir,'files\fileloc_',ext,'.mat'],'file')
        x=load([top_dir,'files\fileloc_',ext,'.mat']);
        tardir=x.fileloc.tardir;
        % Priority 1: imageloc.mat pointer to the raw data folder
        % DO NOTHING! load([top_dir,'files\imageloc.mat'])
    elseif ~isempty(dir([top_dir,'files\*.',ext]))
        % Priority 2: files stored in the projects folder
        tardir=[top_dir,'files'];
        list=dir([top_dir,'files\*',ext]);
        mkfilelist(dataset,tardir,list,ext);
    elseif ~isempty(dir([top_dir,'*',ext]))
        % Priority 3: files in the ds folder but not sorted into the
        % files directory yet
        list=dir([top_dir,'*.',ext]);
        tardir=[top_dir,'files'];
        movefile([top_dir,'*.',ext],[top_dir,'files\']);
        mkfilelist(dataset,tardir,list,ext);
    elseif exist('tardir','var')
        %tardir=[top_dir,'files'];
        list=dir([tardir,'\*',ext]);
        mkfilelist(dataset,tardir,list,ext);
    else
        % Request that the user chose the location of the raw files
        if nargin==4
            [tarfile,tardir]=uigetfile(...
                ['\\research.wpi.edu\gaudettelab\Raw Data\*.',ext],...
            'Select the target file');
            list=dir([tardir,tarfile]);
        else
            tardir=uigetdir('\\research.wpi.edu\gaudettelab\Raw Data\',...
                'Select directory with the target files');
            list=dir([tardir,'\*',ext]);
        end
        mkfilelist(dataset,tardir,list,ext);
        if nargin==4
            xx=load([dataset,'\files\fileloc_',ext,'.mat']);
            fileloc=xx.fileloc;
            num=find(~cellfun(@isempty,strfind(xx.fileloc.paths,tarfile)));
            fileloc.paths=fileloc.paths(num);
            fileloc.sha256=fileloc.sha256(num);
            fileloc.number={1};
            save([dataset,'\files\fileloc_',ext,'.mat'],'fileloc');
        end
    end
    
    if nargout == 1
        varargout{1}=tardir;
    end
    
end