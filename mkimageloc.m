function mkimageloc(dataset)
chkdir
top_dir=[dataset,'\'];
    if exist([top_dir,'files\fileloc_tif.mat'],'file')
        % Priority 1: imageloc.mat pointer to the raw data folder
        % DO NOTHING! load([top_dir,'images\imageloc.mat'])
    elseif exist([top_dir,'images\imageloc.mat'],'file')
        copyfile([top_dir,'images\fileloc_tif.mat'],...
            [top_dir,'images\imageloc.mat'])
        % Priority 1: imageloc.mat pointer to the raw data folder
        % DO NOTHING! load([top_dir,'images\imageloc.mat'])
    elseif ~isempty(dir([top_dir,'images\*.tif']))
        % Priority 2: images stored in the projects folder
        list=dir([top_dir,'images\*.tif']);
        mkfilelist(dataset,[dataset,'\images'],list,'tif');
    elseif ~isempty(dir([top_dir,'*.tif']))
        % Priority 3: images in the ds folder but not sorted into the
        % images directory yet
        list=dir([top_dir,'*.tif']);
        movefile([top_dir,'*.tif'],[top_dir,'images\']);
       mkfilelist(dataset,[dataset,'\images'],list,'tif');
    else
        % Request that the user chose the location of the raw images
        tardir=uigetdir('\\research.wpi.edu\gaudettelab\Raw Data\',...
            'Select directory with the target tif files');
        list=dir([tardir,'\*.tif']);
        mkfilelist(dataset,tardir,list,'tif');
    end
end