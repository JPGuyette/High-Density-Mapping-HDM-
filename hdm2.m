function varargout=hdm2_sa(dataset,framerate,frames,subsize,subshift,varargin)
% Compiled by BIOMED-NB03\jfavreau on 02/20/2015 14:02:13
% hdm2 Version 2015.01.0007
% 
% HDM2 mimics the C# version of HDM used in the lab previously. Outputs a
% file named 'rawhdm2.mat' in the dataset\matlab_data folder. This file can
% be used by other functions to calculate parameters from the HDM data.
%
%   Required Inputs:
%   dataset: Dataset to be analyzed (Defines the folder to be used for this
%       analysis
%   frames: frame numbers to be analyzed
%   subsize: Size of subimages to be used in the form [y,x]
%   subshift: Distance to shift between subimages in the form [y,x].
%       subsize must be evenly divisible by subshift
%   
%   Optional Inputs:
%   roi: Nx2 array of pixel points representing the location of the ROI
%       (user is required to select this if it is not defined)
%   subpix: Subpixel algorithm to be used. Default:foroosh
%   rigidbodymotion: logical defining whether to use rigid body motion or
%       not. Default:RBM is on
%   importrbm: Numeric defining whether to calculate rigid body shifting(0)
%       , to import RBM from a rawhdm.mat file (1) or to have the user
%       select a different region for RBS tracking. Has no effect on
%       analysis if rigidbodymotion is set to 0. Default is 0.
%   roitype: string specifying what type of roi to request from the user.
%       Default: 'impoly'
%   comparison:
%       'seq'(default): sequential comparison of 1->2, 2->3. 3->4...
%       'each': comparison to first image e.g. 1->2, 1->3, 1->4...
%   subpix: Subpixel algorithm to be used. Options are 'Foroosh' (default),
%       'none','esinc','favreau1'
%   showprogress:
%       1 (default) show progress bars
%       0 Suppress progress bars
%   checkframerate:
%       1 (default) Check with user if framerate doesnt match source
%       0 Ignore framerate mismatch
%
%   Version 2.0 fixes, enhancements over C# version
%   - Non rectangular ROIs can be selected and tracked
%   - Rigid body shifting uses an adaptive sized rectangle to maximize
%       accuracy of pixel level motion
%   - ROI selection forces all pixels used in HDM analysis to be WITHIN the
%       selected ROI
%
%   (C)2014 John T Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute
%   --


        
%% Process inputs
    p=inputParser;
    p.addRequired('dataset',@ischar);
    p.addRequired('framerate',@isnumeric);
    p.addRequired('frames',@isnumeric);
    p.addRequired('subsize',@isnumeric)
    p.addRequired('subshift',@isnumeric)
    p.addParamValue('roi',0,@(x)(isnumeric(x) && size(x,2)==2))
    p.addParamValue('subpix','foroosh',...
        @(x) any(strcmp(x,{'none','esinc','foroosh','favreau1'})));
    p.addParamValue('rigidbodymotion',1,@(x) x==0||x==1)
    p.addParamValue('roitype','impoly',...
          @(x) any(strcmp(x,{'impoly','imrect','imsquare','imcircle',...
          'imellipse'})));
    p.addParamValue('comparison','seq',...
        @(x) any(strcmp(x,{'seq','each'})))
    p.addParamValue('showprogress',1,@isnumeric)
    p.addParamValue('importrbm',0,@(x) x==0||x==1||x==2)
    p.addParamValue('checkframerate',1,@(x) any(x==[0 1]))
    p.parse(dataset,framerate,frames,subsize,subshift,varargin{:});
    comp=p.Results.comparison;
    prog=p.Results.showprogress;
    cfr=p.Results.checkframerate;
    
    
   

    %% Gaudette projects folder protection
    chkdir
    clearvars compname curuser buildin buildyear sessdate x y
    
    %% Create folders as necessary
    if ~exist(dataset,'dir')
        mkdir(dataset)
        y=1;
    else
        y=0;
    end
    top_dir=[dataset,'\'];
    subdirs={'matlab_data','images','files'};
    for i =1:length(subdirs)
        if~exist([top_dir,subdirs{i}],'dir')
            mkdir([top_dir,subdirs{i}]);
        end
    end
    
    rbm=p.Results.rigidbodymotion;
    roi=p.Results.roi;
    roitype=p.Results.roitype;
    subpix=p.Results.subpix;
    N=length(frames);
    framerate=p.Results.framerate;
    if p.Results.importrbm==1
        if ~exist([top_dir,'matlab_data\rbmdata.mat'],'file')
            [matfile,matpath]=uigetfile('*.mat',['Select the mat file ',...
                'containing rigid body motion data...']);
            rbmdata=load([matpath,matfile]);
            save([top_dir,'matlab_data\rbmdata.mat'],'rbmdata')
        else
            load([top_dir,'matlab_data\rbmdata.mat'])
        end
        rbmshift=rbmdata.rawhdm.xxyy(2:end,:)-...
            rbmdata.rawhdm.xxyy(1:end-1,:);
        rbm=2;
    elseif p.Results.importrbm==2
       error('RBM shifting by choosing a different region is not yet supported')
    else
        rbmshift=zeros(N,4);
    end
    
    if numel(subsize)==1
        subsize=[subsize subsize];
    end
    if numel(subshift)==1
        subshift=[subshift subshift];
    end
    
    %% Get absolute image locations in order of location priority
    dataset=p.Results.dataset;
    loc=mkfileloc(dataset,'tif');
    mkfileloc(dataset,'exp',loc);
    %% Check framerate input if possible
    exp=load([top_dir,'files\fileloc_exp.mat']);
    if ~isempty(exp.fileloc.paths)
        f=fopen(exp.fileloc.paths{1});
        w=textscan(f,'%s');
        x=strfind(w{1},'Framerate=');
        y=~cellfun(@isempty,x);
        w=w{1};
        w=w(y);
        w=w{1};
        z=strfind(w,'=');
        fr=str2double(w(z+1:end));
        if fr~=(framerate*abs(frames(2)-frames(1))) && cfr
            x=questdlg(['The specified framerate(',num2str(framerate),...
                ' fps) does not match the record framerate (',...
                num2str(fr/abs(frames(2)-frames(1))),' fps',...
                ' \n Which framerate would you like to use?'],...
                ['Please verify your selected framerate. ds=',dataset],...
                ['Specified: ',num2str(framerate),' fps'],...
                ['Actual: ',num2str(fr/abs(frames(2)-frames(1))),' fps'],...
                ['Actual: ',num2str(fr/abs(frames(2)-frames(1))),' fps']);
            if strcmp(x,['Actual: ',num2str(fr/abs(frames(2)-frames(1))),...
                    ' fps'])
                framerate=fr/abs(frames(2)-frames(1));
            end
        end
    end
    load([top_dir,'files\fileloc_tif.mat'])
    
    % Ensure that all target frames are found
    if ~all(ismember(frames,fileloc.number))
        error('hdm2:missingframes', ...
            'Not all target frames are available in the dataset')
    end
    
%% Process image sequence in order
    % Get the frames that we want to analyze
    [~,tarind]=ismember(frames,fileloc.number);
    images=fileloc.paths(tarind);
    
    % Ensure that the images are in numerical order
%     [~,ind]=sort(fileloc.number(tarind));
%     images=images(ind);
    
    % If an ROI is not specified, request the user for an ROI in the first
    % frame
    im1 = imread(images{1});
    if size(im1,3)>1;im1=rgb2gray(im1);end;
    if roi==0
        [roi, BB,shape, roiim]=hdmgetroi(im1,subsize,subshift,roitype);
        imwrite(roiim,[top_dir,'matlab_data\hdmroi.tif']);
    else
        %[~,~,BB,shape]=getsubimages(roi,im1,subsize,subshift,roitype);
        im1=imread(images{1});
        if size(im1,3)~=3
            im1show=repmat(im1(:,:,1),[1 1 3]);
        else
            im1show=im1;
        end
        if isa(im1show,'uint8')
            im1show=double(im1show)/255;
        elseif isa(im1show,'uint16')
            im1show=double(im1show)/(2^16);
        end
        [BB,shape,mask]=disptrueroi(im1,subsize,subshift,roi);
        imwrite(mask.*im1show,[top_dir,'matlab_data\hdmroi.tif']);
    end
    tstart=now();
    if prog==1,bar=waitbar(0,'Processing image sequence...');end
    
    % Preallocate data output
    hdmtimer=tic;
    comptime=0;
    subimcount=0;
    
    udata=NaN([shape,N-1]);
    vdata=NaN([shape,N-1]);
    xxyy=NaN([N-1,4]);
    pixX=NaN([shape,N-1]);
    pixY=NaN([shape,N-1]);
    targetroi=roi;
    im1=imread(images{1});
    if size(im1,3)>1;im1=rgb2gray(im1);end;
    
    roilist=NaN([size(roi),N]);
    roilist(:,:,1)=targetroi;
    im2=im1;
    if prog==1,waitbar(1/N,bar);end
    for k=1:N-1
        % Write values to output matrices
        
        switch comp
            case 'seq'
                im1=im2;
                if k==1
                    xxyy(k,:)=[BB(1),BB(1)+BB(3)-1,BB(2),BB(2)+BB(4)-1];
                else
                    xxyy(k,:)=xxyy(k-1,:)+temphdm.xxyyshift;
                end
            case 'each'
                if k==1
                    xxyy(k,:)=[BB(1),BB(1)+BB(3)-1,BB(2),BB(2)+BB(4)-1];
                else
                    xxyy(k,:)=xxyy(1,:)+temphdm.xxyyshift;
                end
        end
        im2=imread(images{k+1});
        if size(im2,3)~=1;im2=rgb2gray(im2);end
        [temphdm,subim]=roicorr(im1,im2,subsize,subshift,targetroi,...
            'subpix',subpix,'rigidbodymotion',rbm,'roitype',roitype,...
            'importrbm',rbmshift(k,:));
        subimcount=subimcount+subim.count;
        comptime=comptime+subim.comptime;
        % Write displacement data to output matrices
        pixX(:,:,k)=subim.pixX;
        pixY(:,:,k)=subim.pixY;
        udata(:,:,k)=temphdm.udata;
        vdata(:,:,k)=temphdm.vdata; % Positive direction indicates downward motion
        switch comp
            case 'seq'
                targetroi=targetroi+repmat(temphdm.xxyyshift([1 3]),...
                    [size(targetroi,1),1]); %Shift ROI!
        end
        trem=(now()-tstart)*(N-k+1)/(k+1);
        if prog==1,waitbar((k+1)/N,bar,...
                ['Time remaining: ',datestr(trem,'HH:MM:SS')]);end
        roilist(:,:,k+1)=targetroi; %% REMOVED SHIFTING HERE
    end
    if prog==1, delete(bar);end
    use=subim.use; 
    rawhdm.xxyy=[xxyy;xxyy(k-1,:)+temphdm.xxyyshift;];
    rawhdm.udata=udata;
    rawhdm.vdata=vdata;% Positive shifts are displacements in the downward direction
    rawhdm.subsize=subsize;
    rawhdm.framerate=framerate;
    rawhdm.subshift=subshift; 
    rawhdm.roilist=roilist;
    switch nargout
        case 0
            
            
            save([top_dir,'matlab_data\rawhdm2.mat'],...
                'rawhdm','frames','pixX','pixY','use','subimcount','comptime',...
                'roi','comp');
        case 1
            output.rawhdm=rawhdm;
            output.frames=frames;
            output.pixX=pixX;
            output.pixY=pixY;
            output.use=use;
            output.subimcount=subimcount;
            output.comptime=comptime;
            output.roi=roi;
            output.comp=comp;
            varargout{1}=output;
            save([top_dir,'matlab_data\rawhdm2.mat'],...
                'rawhdm','frames','pixX','pixY','use','subimcount','comptime',...
                'roi','comp');
    end
    toc(hdmtimer);

end

        
        

function [varargout] = GCE_curbuild(i) %#ok<INUSD>
    GCEloc='';

    % Get newest build year
        files=dir([GCEloc,'~GCE\~builds']);
        buildyear=files(end).name;
        
    % Get newest build
        files=dir([GCEloc,'~GCE\~builds\',buildyear]);
        curbuild=files(end).name;
        
    % Handle output
        switch nargout 
            case {1,0}
                varargout(1)={[buildyear,'.',curbuild]};
            case 2
                varargout={curbuild,buildyear};
            case 3
                varargout={curbuild(end-3:end),curbuild(1:2),buildyear};
            otherwise
                error('Too many output arguments');
        end
        if exist('i','var')
            build.version=[buildyear,'.',curbuild];
            build.matlabver=version;
            build.accessdate=datestr(now(),'yyyy-mm-dd');
            build.computer=lower(deblank(getenv('COMPUTERNAME')));
            logs=dir([GCEloc,'~GCE\~log\*',build.computer,'.txt']);
            build.logfile=[GCEloc,'~GCE\~log\',logs(end).name];
            
            varargout{1}=build;
        end
            
end

function Hash = DataHash(Data, Opt)
% DATAHASH - Checksum for Matlab array of any type
% This function creates a hash value for an input of any type. The type and
% dimensions of the input are considered as default, such that UINT8([0,0]) and
% UINT16(0) have different hash values. Nested STRUCTs and CELLs are parsed
% recursively.
%
% Hash = DataHash(Data, Opt)
% INPUT:
%   Data: Array of these built-in types:
%           (U)INT8/16/32/64, SINGLE, DOUBLE, (real or complex)
%           CHAR, LOGICAL, CELL (nested), STRUCT (scalar or array, nested),
%           function_handle.
%   Opt:  Struct to specify the hashing algorithm and the output format.
%         Opt and all its fields are optional.
%         Opt.Method: String, known methods for Java 1.6 (Matlab 2009a):
%              'SHA-1', 'SHA-256', 'SHA-384', 'SHA-512', 'MD2', 'MD5'.
%            Known methods for Java 1.3 (Matlab 6.5):
%              'MD5', 'SHA-1'.
%            Default: 'MD5'.
%         Opt.Format: String specifying the output format:
%            'hex', 'HEX':      Lower/uppercase hexadecimal string.
%            'double', 'uint8': Numerical vector.
%            'base64':          Base64 encoded string, only printable
%                               ASCII characters, 33% shorter than 'hex'.
%            Default: 'hex'.
%         Opt.Input: Type of the input as string, not case-sensitive:
%             'array': The contents, type and size of the input [Data] are
%                      considered  for the creation of the hash. Nested CELLs
%                      and STRUCT arrays are parsed recursively. Empty arrays of
%                      different type reply different hashs.
%             'file':  [Data] is treated as file name and the hash is calculated
%                      for the files contents.
%             'bin':   [Data] is a numerical, LOGICAL or CHAR array. Only the
%                      binary contents of the array is considered, such that
%                      e.g. empty arrays of different type reply the same hash.
%             Default: 'array'.
%
% OUTPUT:
%   Hash: String, DOUBLE or UINT8 vector. The length depends on the hashing
%         method.
%
% EXAMPLES:
% % Default: MD5, hex:
%   DataHash([])                % 7de5637fd217d0e44e0082f4d79b3e73
% % MD5, Base64:
%   Opt.Format = 'base64';
%   Opt.Method = 'MD5';
%   DataHash(int32(1:10), Opt)  % bKdecqzUpOrL4oxzk+cfyg
% % SHA-1, Base64:
%   S.a = uint8([]);
%   S.b = {{1:10}, struct('q', uint64(415))};
%   Opt.Method = 'SHA-1';
%   DataHash(S, Opt)            % ZMe4eUAp0G9TDrvSW0/Qc0gQ9/A
% % SHA-1 of binary values:
%   Opt.Method = 'SHA-1';
%   Opt.Input  = 'bin';
%   DataHash(1:8, Opt)          % 826cf9d3a5d74bbe415e97d4cecf03f445f69225
%
% NOTE:
%   Function handles and user-defined objects cannot be converted uniquely:
%   - The subfunction ConvertFuncHandle uses the built-in function FUNCTIONS,
%     but the replied struct can depend on the Matlab version.
%   - It is tried to convert objects to UINT8 streams in the subfunction
%     ConvertObject. A conversion by STRUCT() might be more appropriate.
%   Adjust these subfunctions on demand.
%
%   MATLAB CHARs have 16 bits! In consequence the string 'hello' is treated as
%   UINT16('hello') for the binary input method.
%
%   DataHash uses James Tursa's smart and fast TYPECASTX, if it is installed:
%     http://www.mathworks.com/matlabcentral/fileexchange/17476
%   As fallback the built-in TYPECAST is used automatically, but for large
%   inputs this can be more than 100 times slower.
%   For Matlab 6.5 installing typecastx is obligatory to run DataHash.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2011-2012 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also: TYPECAST, CAST.
% FEX:
% Michael Kleder, "Compute Hash", no structs and cells:
%   http://www.mathworks.com/matlabcentral/fileexchange/8944
% Tim, "Serialize/Deserialize", converts structs and cells to a byte stream:
%   http://www.mathworks.com/matlabcentral/fileexchange/29457
% Jan Simon, "CalcMD5", MD5 only, faster C-mex, no structs and cells:
%   http://www.mathworks.com/matlabcentral/fileexchange/25921

% $JRev: R-k V:011 Sum:kZG25iszfKbg Date:28-May-2012 12:48:06 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLFile\DataHash.m $
% History:
% 001: 01-May-2011 21:52, First version.
% 007: 10-Jun-2011 10:38, [Opt.Input], binary data, complex values considered.
% 011: 26-May-2012 15:57, Fails for binary input and empty data.

% Main function: ===============================================================
% Java is needed:
if ~usejava('jvm')
   error(['JSimon:', mfilename, ':NoJava'], ...
      '*** %s: Java is required.', mfilename);
end

% typecastx creates a shared data copy instead of the deep copy as Matlab's
% TYPECAST - for a [1000x1000] DOUBLE array this is 100 times faster!
persistent usetypecastx
if isempty(usetypecastx)
   usetypecastx = ~isempty(which('typecastx'));  % Run the slow WHICH once only
end

% Default options: -------------------------------------------------------------
Method    = 'MD5';
OutFormat = 'hex';
isFile    = false;
isBin     = false;

% Check number and type of inputs: ---------------------------------------------
nArg = nargin;
if nArg == 2
   if isa(Opt, 'struct') == 0   % Bad type of 2nd input:
      error(['JSimon:', mfilename, ':BadInput2'], ...
         '*** %s: 2nd input [Opt] must be a struct.', mfilename);
   end
   
   % Specify hash algorithm:
   if isfield(Opt, 'Method')
      Method = upper(Opt.Method);
   end
   
   % Specify output format:
   if isfield(Opt, 'Format')
      OutFormat = Opt.Format;
   end
   
   % Check if the Input type is specified - default: 'array':
   if isfield(Opt, 'Input')
      if strcmpi(Opt.Input, 'File')
         isFile = true;
         if ischar(Data) == 0
            error(['JSimon:', mfilename, ':CannotOpen'], ...
               '*** %s: 1st input is not a file name', mfilename);
         end
         
         if exist(Data, 'file') ~= 2
            error(['JSimon:', mfilename, ':FileNotFound'], ...
               '*** %s: File not found: %s.', mfilename, Data);
         end
         
      elseif strncmpi(Opt.Input, 'bin', 3)  % Accept 'binary'
         isBin = true;
         if (isnumeric(Data) || ischar(Data) || islogical(Data)) == 0
            error(['JSimon:', mfilename, ':BadDataType'], ...
               '*** %s: 1st input is not numeric, CHAR or LOGICAL.', mfilename);
         end
      end
   end
   
elseif nArg ~= 1  % Bad number of arguments:
   error(['JSimon:', mfilename, ':BadNInput'], ...
      '*** %s: 1 or 2 inputs required.', mfilename);
end

% Create the engine: -----------------------------------------------------------
try
   Engine = java.security.MessageDigest.getInstance(Method);
catch
   error(['JSimon:', mfilename, ':BadInput2'], ...
      '*** %s: Invalid algorithm: [%s].', mfilename, Method);
end

% Create the hash value: -------------------------------------------------------
if isFile
   % Read the file and calculate the hash:
   FID = fopen(Data, 'r');
   if FID < 0
      error(['JSimon:', mfilename, ':CannotOpen'], ...
         '*** %s: Cannot open file: %s.', mfilename, Data);
   end
   Data = fread(FID, Inf, '*uint8');
   fclose(FID);
   
   Engine.update(Data);
   if usetypecastx
      Hash = typecastx(Engine.digest, 'uint8');
   else
      Hash = typecast(Engine.digest, 'uint8');
   end

elseif isBin             % Contents of an elementary array:
   if isempty(Data)      % Nothing to do, Engine.update fails for empty input!
      Hash = typecastx(Engine.digest, 'uint8');
   elseif usetypecastx   % Faster typecastx:
      if isreal(Data)
         Engine.update(typecastx(Data(:), 'uint8'));
      else
         Engine.update(typecastx(real(Data(:)), 'uint8'));
         Engine.update(typecastx(imag(Data(:)), 'uint8'));
      end
      Hash = typecastx(Engine.digest, 'uint8');
      
   else                  % Matlab's TYPECAST is less elegant:
      if isnumeric(Data)
         if isreal(Data)
            Engine.update(typecast(Data(:), 'uint8'));
         else
            Engine.update(typecast(real(Data(:)), 'uint8'));
            Engine.update(typecast(imag(Data(:)), 'uint8'));
         end
      elseif islogical(Data)               % TYPECAST cannot handle LOGICAL
         Engine.update(typecast(uint8(Data(:)), 'uint8'));
      elseif ischar(Data)                  % TYPECAST cannot handle CHAR
         Engine.update(typecast(uint16(Data(:)), 'uint8'));
         Engine.update(typecast(Data(:), 'uint8'));
      end
      Hash = typecast(Engine.digest, 'uint8');
   end
   
elseif usetypecastx  % Faster typecastx:
   Engine = CoreHash_(Data, Engine);
   Hash   = typecastx(Engine.digest, 'uint8');
   
else                 % Slower built-in TYPECAST:
   Engine = CoreHash(Data, Engine);
   Hash   = typecast(Engine.digest, 'uint8');
end

% Convert hash specific output format: -----------------------------------------
switch OutFormat
   case 'hex'
      Hash = sprintf('%.2x', double(Hash));
   case 'HEX'
      Hash = sprintf('%.2X', double(Hash));
   case 'double'
      Hash = double(reshape(Hash, 1, []));
   case 'uint8'
      Hash = reshape(Hash, 1, []);
   case 'base64'
      Hash = fBase64_enc(double(Hash));
   otherwise
      error(['JSimon:', mfilename, ':BadOutFormat'], ...
         '*** %s: [Opt.Format] must be: HEX, hex, uint8, double, base64.', ...
         mfilename);
end

end

% ******************************************************************************
function Engine = CoreHash_(Data, Engine)
% This mothod uses the faster typecastx version.

% Consider the type and dimensions of the array to distinguish arrays with the
% same data, but different shape: [0 x 0] and [0 x 1], [1,2] and [1;2],
% DOUBLE(0) and SINGLE([0,0]):
Engine.update([uint8(class(Data)), typecastx(size(Data), 'uint8')]);

if isstruct(Data)                    % Hash for all array elements and fields:
   F      = sort(fieldnames(Data));  % Ignore order of fields
   Engine = CoreHash_(F, Engine);    % Catch the fieldnames
   
   for iS = 1:numel(Data)            % Loop over elements of struct array
      for iField = 1:length(F)       % Loop over fields
         Engine = CoreHash_(Data(iS).(F{iField}), Engine);
      end
   end
   
elseif iscell(Data)                  % Get hash for all cell elements:
   for iS = 1:numel(Data)
      Engine = CoreHash_(Data{iS}, Engine);
   end
      
elseif isnumeric(Data) || islogical(Data) || ischar(Data)
   if isempty(Data) == 0
      if isreal(Data)                % TRUE for LOGICAL and CHAR also:
         Engine.update(typecastx(Data(:), 'uint8'));
      else                           % typecastx accepts complex input:
         Engine.update(typecastx(real(Data(:)), 'uint8'));
         Engine.update(typecastx(imag(Data(:)), 'uint8'));
      end
   end
   
elseif isa(Data, 'function_handle')
   Engine = CoreHash(ConvertFuncHandle(Data), Engine);
   
else  % Most likely this is a user-defined object:
   try
      Engine = CoreHash(ConvertObject(Data), Engine);
   catch
      warning(['JSimon:', mfilename, ':BadDataType'], ...
         ['Type of variable not considered: ', class(Data)]);
   end
end

end

% ******************************************************************************
function Engine = CoreHash(Data, Engine)
% This methods uses the slower TYPECAST of Matlab
% See CoreHash_ for comments.

Engine.update([uint8(class(Data)), typecast(size(Data), 'uint8')]);

if isstruct(Data)                    % Hash for all array elements and fields:
   F      = sort(fieldnames(Data));  % Ignore order of fields
   Engine = CoreHash(F, Engine);     % Catch the fieldnames
   for iS = 1:numel(Data)            % Loop over elements of struct array
      for iField = 1:length(F)       % Loop over fields
         Engine = CoreHash(Data(iS).(F{iField}), Engine);
      end
   end
elseif iscell(Data)                  % Get hash for all cell elements:
   for iS = 1:numel(Data)
      Engine = CoreHash(Data{iS}, Engine);
   end
elseif isempty(Data)
elseif isnumeric(Data)
   if isreal(Data)
      Engine.update(typecast(Data(:), 'uint8'));
   else
      Engine.update(typecast(real(Data(:)), 'uint8'));
      Engine.update(typecast(imag(Data(:)), 'uint8'));
   end
elseif islogical(Data)               % TYPECAST cannot handle LOGICAL
   Engine.update(typecast(uint8(Data(:)), 'uint8'));
elseif ischar(Data)                  % TYPECAST cannot handle CHAR
   Engine.update(typecast(uint16(Data(:)), 'uint8'));
elseif isa(Data, 'function_handle')
   Engine = CoreHash(ConvertFuncHandle(Data), Engine);
else  % Most likely a user-defined object:
   try
      Engine = CoreHash(ConvertObject(Data), Engine);
   catch
      warning(['JSimon:', mfilename, ':BadDataType'], ...
         ['Type of variable not considered: ', class(Data)]);
   end
end

end

% ******************************************************************************
function FuncKey = ConvertFuncHandle(FuncH)
%   The subfunction ConvertFuncHandle converts function_handles to a struct
%   using the Matlab function FUNCTIONS. The output of this function changes
%   with the Matlab version, such that DataHash(@sin) replies different hashes
%   under Matlab 6.5 and 2009a.
%   An alternative is using the function name and name of the file for
%   function_handles, but this is not unique for nested or anonymous functions.
%   If the MATLABROOT is removed from the file's path, at least the hash of
%   Matlab's toolbox functions is (usually!) not influenced by the version.
%   Finally I'm in doubt if there is a unique method to hash function handles.
%   Please adjust the subfunction ConvertFuncHandles to your needs.

% The Matlab version influences the conversion by FUNCTIONS:
% 1. The format of the struct replied FUNCTIONS is not fixed,
% 2. The full paths of toolbox function e.g. for @mean differ.
FuncKey = functions(FuncH);

% ALTERNATIVE: Use name and path. The <matlabroot> part of the toolbox functions
% is replaced such that the hash for @mean does not depend on the Matlab
% version.
% Drawbacks: Anonymous functions, nested functions...
% funcStruct = functions(FuncH);
% funcfile   = strrep(funcStruct.file, matlabroot, '<MATLAB>');
% FuncKey    = uint8([funcStruct.function, ' ', funcfile]);

% Finally I'm afraid there is no unique method to get a hash for a function
% handle. Please adjust this conversion to your needs.

end

% ******************************************************************************
function DataBin = ConvertObject(DataObj)
% Convert a user-defined object to a binary stream. There cannot be a unique
% solution, so this part is left for the user...

% Perhaps a direct conversion is implemented:
DataBin = uint8(DataObj);

% Or perhaps this is better:
% DataBin = struct(DataObj);

end

% ******************************************************************************
function Out = fBase64_enc(In)
% Encode numeric vector of UINT8 values to base64 string.

Pool = [65:90, 97:122, 48:57, 43, 47];  % [0:9, a:z, A:Z, +, /]
v8   = [128; 64; 32; 16; 8; 4; 2; 1];
v6   = [32, 16, 8, 4, 2, 1];

In  = reshape(In, 1, []);
X   = rem(floor(In(ones(8, 1), :) ./ v8(:, ones(length(In), 1))), 2);
Y   = reshape([X(:); zeros(6 - rem(numel(X), 6), 1)], 6, []);
Out = char(Pool(1 + v6 * Y));

end


function chkdir
err=0;
    if strcmp(pwd,'')
       x=questdlg(['MatLab is currently working in the Gaudette Lab Projects folder.',...
           '  Typically analysis is performed in a subfolder of Projects.',...
           ' Would you like to continue?'],...
           'Screw up data management for the entire lab?',...
           'Continue ruining data management for the entire lab',...
           'Cancel and manually change the working directory',...
           'Cancel and manually change the working directory');
       if strcmp(x,'Continue ruining data management for the entire lab')
           x=questdlg(['Seriously?.',...
           '  You really want to F*** everything up for the whole lab.',...
           ' because your too !#%$@%#$^ lazy to change the directory?'],...
           'Screw up data management for the entire lab?',...
           'Yes, F*** it all',...
           'Fine, I''ll change directories',...
           'Fine, I''ll change directories');
            if strcmp(x,'Yes, F*** it all')
                compname = getenv('COMPUTERNAME');
                compname = lower(deblank(compname));
                % Get the Last name of the current user
                curuser=[getenv('userdomain'),'\',getenv('username')];
                % Determine target build ID
                if ~exist('buildin','var')
                    [buildin,buildyear]=GCE_curbuild;
                end
                buildid=[buildyear,'.',buildin]; 
               % Create session log file
                sessdate=datestr(now(),'yyyy_mm_dd_HH:MM:SS');
                %logs=dir([GCEloc,'~GCE\~log\',sessdate,'*']);
                fpath='f---info.txt';
                logfile=fopen(fpath,'at+');
                fprintf(logfile,'\n%s','############');
                fprintf(logfile,'\n%s',sessdate);
                fprintf(logfile,'\n%s',['Build: ',buildid]);
                fprintf(logfile,'\n%s',['User: ',curuser]);
                fprintf(logfile,'\n%s',['Workstation: ',compname]);
                fprintf(logfile,'\n%s',['Matlab Version: ',version]);
                fclose('all');
                x=0;
                bar=waitbar(0,'F***ing up the lab. Please wait...');
                while x<30
                    pause(1)
                    x=x+1;
                    waitbar(x/30,bar);
                end
            else
                err=1;
            end
       else
           err=1 ;
       end
       if err==1
          error('User decided not to F*** up the lab'); 
       end
    end
    clearvars compname curuser buildin buildyear sessdate x y
end

        function [BB,shape,mask]=disptrueroi(im1,subsize,subshift,pos)
            if size(im1,3)~=3
                im1show=repmat(im1(:,:,1),[1 1 3]);
            else
                im1show=im1;
            end
            % Get the subimage positions and draw them on the plot
            [U,C,BB,shape]=getsubimages(pos,im1,subsize,subshift);
            % Create mask with overall position
            mask=ones(size(im1show));
            mask(:,:,2)=mask(:,:,2).*(U*0.6);
            mask(:,:,3)=mask(:,:,2);
            % Add sample subimage at centroid
            mask(C(2):C(2)-1+subsize(1),C(1):C(1)-1+subsize(2),2)=1;
            mask(C(2):C(2)-1+subsize(1),C(1):C(1)-1+subsize(2),3)=0.6;
            mask(C(2):C(2)-1+subsize(1),C(1):C(1)-1+subsize(2),1)=0.6;
            % Add overlapping subimage in each direction to show shift
            mask(C(2)+subsize(1):C(2)-1+subsize(1)+subshift(1),C(1):C(1)-1+subsize(2),2)=0.6;
            mask(C(2)+subsize(1):C(2)-1+subsize(1)+subshift(1),C(1):C(1)-1+subsize(2),1)=0.6;
            mask(C(2)+subsize(1):C(2)-1+subsize(1)+subshift(1),C(1):C(1)-1+subsize(2),3)=1;

            mask(C(2):C(2)-1+subsize(1),C(1)+subsize(2):C(1)-1+subsize(2)+subshift(2),2)=0.6;
            mask(C(2):C(2)-1+subsize(1),C(1)+subsize(2):C(1)-1+subsize(2)+subshift(2),1)=0.6;
            mask(C(2):C(2)-1+subsize(1),C(1)+subsize(2):C(1)-1+subsize(2)+subshift(2),3)=1;

            mask(mask==0)=1;
            
    %         ann=annotation('rectangle',[C(1)+subshift(1),C(2),subsize(1),subsize(2)],'FaceColor','b','FaceAlpha',0.5);
    %         rectangle('Position',[C(1)-subshift(1),C(2),subsize(1),subsize(2)],'FaceColor','g');
    %         rectangle('Position',[C(1),C(2)+subshift(2),subsize(1),subsize(2)],'FaceColor','k');
    %         rectangle('Position',[C(1),C(2)-subshift(2),subsize(1),subsize(2)],'FaceColor','y');
            %setPosition(h,roi);
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

function out=esinc(args,x)
    A=args(1);
    B=args(2);
    C=args(3);
    if x==C
        out=A;
    else
        out=A*exp(-(B*(x-C)).^2).*((sin(pi()*B*(x-C)))./(pi()*B*(x-C)));
    end
end

function sse=esinc_fit(params,c)
    sse=sum((c-esinc(params,[-1 0 1])).^2);
end

 function [varargout]=getsubimages(pos,im1,subsize,subshift)
        %% Performs 2 separate functions
        %1 given 2 outputs, outputs a logical pixel mask
        %2 given 1 output, outputs a logical struct mask
        
        % Get properties of selected region
        mask=poly2mask(pos(:,1),pos(:,2),size(im1,1),size(im1,2));
        R=regionprops(mask,'Centroid','BoundingBox');
        %Get top left corner of center subimage
        xc=floor(R.Centroid(1,1))-(subsize(2)/2);
        yc=floor(R.Centroid(1,2))-(subsize(1)/2);
        % Get top left pixel of top left subimage in bounding box
        R.BoundingBox(1:2)=floor(R.BoundingBox(1:2));
        xs=min(xc:-subshift(2):R.BoundingBox(1));
        ys=min(yc:-subshift(1):R.BoundingBox(2));
        % Get right and bottom boundaries for subimages
        xbound=R.BoundingBox(1)+R.BoundingBox(3)-subsize(2);
        ybound=R.BoundingBox(2)+R.BoundingBox(4)-subsize(1);
        % Create location grids (top left corner of each region)
        xarray=xs:subshift(2):xbound;
        yarray=ys:subshift(1):ybound;
        if nargout==1;
            [X,Y]=meshgrid(xarray,yarray);
            X=X+(subsize(2)-1)/2;
            Y=Y+(subsize(1)-1)/2;
        end
        if nargout==4
            U=false(size(mask));
        else
            U=false(length(yarray),length(xarray));
        end
        % Determine if subimages fall within the region
        for I =1:length(xarray)
            for J=1:length(yarray)
                i=xarray(I);
                j=yarray(J);
                cursub=mask(j:j+subsize(1)-1,i:i+subsize(2)-1);
                if all(cursub)
                    if nargout==4
                        U(j:j+subsize(1)-1,i:i+subsize(2)-1)=true;
                    else
                        U(J,I)=true;
                    end
                end
            end
        end
        switch nargout
            case 4
                varargout={U,[xc,yc],R.BoundingBox,...
                    [length(yarray),length(xarray)]};
            case 1
                getsubimout.pixX=X;
                getsubimout.pixY=Y;
                getsubimout.use=U;
                getsubimout.centroid=[xc yc];
                getsubimout.shape=size(U);
                getsubimout.pixelmask=mask;
                getsubimout.boundingbox=R.BoundingBox;
                varargout{1}=getsubimout;
        end
    end

    function [roiL,varargout]=hdmgetroi(im1L,subsizeL,subshiftL,roitypeL)
    global BB
    global shape
    global im1
    global subsize
    global subshift
    global im1show
    global mask
    global roi
    global hh
    global roitype
    global h
    im1=im1L;
    subsize=subsizeL;
    subshift=subshiftL;
    roitype=roitypeL;
        % Show image
        im1=doubleim(im1);
        if size(im1,3)~=3
            im1show=repmat(im1(:,:,1),[1 1 3]);
        else
            im1show=im1;
        end
        close all
        fig=figure;
        hh=imshow(im1show);
        hold on;
        x=[];
        switch roitype
            case 'imrect'
                h=imrect;
            case 'imsquare'
                h=imrect(gca,'PositionConstraintFcn',@forcecircle);
            case 'impoly'
                h=impoly;
            case 'imellipse'
                h=imellipse;
            case 'imcircle'
                h=imellipse(gca,'PositionConstraintFcn',@forcecircle);
%                 setConstrainedPosition(h,forcecircle(getPosition(h)))
%                 setPositionConstraintFcn(h,@forcecircle);
        end
        addNewPositionCallback(h,@roicallback);
        roi=getPosition(h);

        switch roitype
            case {'imrect','imsquare'}
                roi=imrect2impoly(roi);
            case {'imellipse','imcircle'}
                roi=getVertices(h);
        end
        [BB,shape,mask]=disptrueroi(im1,subsize,subshift,roi);
        hold on;
        delete(hh);
        hh=imshow(mask.*im1show);
        uistack(hh,'bottom')
        uiwait(fig);
        
        if nargout>1
            varargout{1}=BB;
        end
        if nargout>2
            varargout{2}=shape;
        end
        if nargout>3
            varargout{3}=mask.*im1show;
        end
        roiL=roi;
    

    end
    
        function roicallback(pos)
            global BB
            global shape
            global im1
            global subsize
            global subshift
            global im1show
            global mask
            global roi
            global hh
            global roitype
            global h
            
            switch roitype
                case {'imellipse','imcircle'}
                    pos=getVertices(h);
                case {'imrect','imsquare'}
                    pos=imrect2impoly(pos);
            end
            [BB,shape,mask]=disptrueroi(im1,subsize,subshift,pos);
                        hold on;
            delete(hh);
            hh=imshow(mask.*im1show);
            uistack(hh,'bottom')
            roi=pos;
    %         ann=annotation('rectangle',[C(1)+subshift(1),C(2),subsize(1),subsize(2)],'FaceColor','b','FaceAlpha',0.5);
    %         rectangle('Position',[C(1)-subshift(1),C(2),subsize(1),subsize(2)],'FaceColor','g');
    %         rectangle('Position',[C(1),C(2)+subshift(2),subsize(1),subsize(2)],'FaceColor','k');
    %         rectangle('Position',[C(1),C(2)-subshift(2),subsize(1),subsize(2)],'FaceColor','y');
            %setPosition(h,roi);
        end
function setpos=forcecircle(pos)
    dd=mean(pos(3:4));
    dx=dd-pos(3);
    dy=dd-pos(4);
    setpos=[pos(1:2),0,0]+[-dx/2,-dy/2,dd,dd];
end

function roi=imrect2impoly(I)
% Converts the [xmin,ymin, width, height] vector output by imrect to the
% more standard 2 x N vector output of impoly,imroi etc

roi=[I(1),I(2);...
    I(1),I(2)+I(4);...
    I(1)+I(3),I(2)+I(4);...
    I(1)+I(3),I(2)];
end

function mkfilelist(dataset,tardir,~,ext)
% mkimagelist: Makes a list of files and calculates sha512 hashes of each
% image. Stores hashes, image numbering and absolute paths in a structure
% in the projects folder. Calculates hashes using: http://www.mathworks.com
% /matlabcentral/fileexchange/31272-datahash
%
%   Inputs:
%   dataset: folder in which to putput image list variable
%   tardir: the directory in which the images are stored. Can be absolute
%       or relative. use '' to specify current working directory
%   filelist: List of file information as generated by the 'dir' command
%
tic
%% Create output struct and destination folders
    top_dir=[dataset,'\'];
    outdir=[top_dir,'files\'];
    if ~exist(outdir,'dir')
        mkdir(outdir)
    end
    if ~strcmp(tardir(end),'\')
        tardir=[tardir,'\'];
    end
    filelist=dir([tardir,'*',ext]);
    if strcmp(tardir,'')
        tardir=pwd;
    end
    N=length(filelist);
    paths=cell(N,1);
    sha256=cell(N,1);
    number=NaN([N,1]);
    % Prepare hashing structure
    opt.Method='SHA-256';
    opt.format='HEX';

% Check if the folder has already been hashed. If so, check the dirhash
% if exist([tardir,'srchash_',ext,'.mat'],'file');
%     load([tardir,'srchash_',ext,'.mat']);
%     % Drop srchash.mat from hash value
%     curdir=dir([tardir,'*',ext]);
%     for j=1:length(curdir)
%         if strcmp(curdir(j).name,['srchash_',ext,'.mat'])>0
%             curdir(j)=[];
%             break;
%         end
%     end
%     curdirhash=DataHash(curdir,opt);
%     if strcmp(srchash.dirhash,curdirhash) 
%         fileloc=srchash;
%         save([outdir,'fileloc_',ext,'.mat'],'fileloc');
%     else
%         warning('HDM:mkfilehash:hashchange',...
%             'The source directory has changed')
%         %keyboard;
%         calchash;
%     end
% else
    curdir=dir([tardir,'*',ext]);
    calchash;
% end
    toc;
    

    

    
    
    function calchash
        %% Loop through files and fill in imageloc structure
        bar=waitbar(0,'Analyzing directory...');
        tstart=now();
        for i=1:N
            paths{i}=[tardir,filelist(i).name];
            fid=fopen(paths{i});
            %curfile=fread(fid);
            fclose('all');
            if length(filelist(i).name)>7
                if ~isnan(str2double(filelist(i).name(end-7:end-4)))
                    number(i)=str2double(filelist(i).name(end-7:end-4));
                end
            end
            sha='temp';
            %sha=DataHash(curfile,opt);
            sha256{i}=sha;
            tleft= (N-i)*(now()-tstart)/i;
            waitbar(i/N,bar,['Analyzing directory. Time remaining: '...
                ,datestr(tleft,'HH:MM:SS')]);
        end
        delete(bar)
        fileloc.paths=paths;
        fileloc.sha256=sha256;
        fileloc.dirhash=DataHash(curdir,opt);
        fileloc.number=number;
        fileloc.dataset=dataset;
        fileloc.tardir=tardir;
    
    %% Output the imageloc structure to the source folder
        srchash=fileloc;
        %save([tardir,'srchash_',ext,'.mat'],'srchash');
        save([outdir,'fileloc_',ext,'.mat'],'fileloc');
    end
end



function [varargout]=pixcorr(subim1,subim2,subpix)
% [disp]=pixcorr(subim1,subim2,subpix)
% [xdisp,ydisp]=pixcorr(subim1,subim2,subpix)
% This method performs a 2D phase correlation on a pair of subimages
% (subim1 and subim2).
%
%   Inputs:
%   subim1: First subimage (matrix of values ranging from 0 to 1
%   subim2: Second subimage. Same size as subim1
%   subpix: Optional. A string defining the subpixel algorithm to be used.
%       Current options are: 'foroosh'(default), 'esinc', 'hdm', and 'none'

    global hamm
    % Parse input
    if ~exist('subpix','var')
        subpix='foroosh';
    end
    if size(subim1)~=size(subim2)
        error('hdm:pixcorr:sizemismatch','Subimages must be the same size')
    end
    subim1=doubleim(subim1);
    subim2=doubleim(subim2);
   
    % Create hamming window
    winsize=size(subim1);
    if size(hamm)~=size(subim1)
        hamm1=repmat(hamming(winsize(1)),[1 winsize(2)]);
        hamm2=repmat(hamming(winsize(2))',[winsize(1),1]);
        hamm=hamm1.*hamm2;
    end
%     % Remove DC component and normalize regions
%     subim1=(subim1-mean(mean(subim1)))/max(max(abs(subim1)));
%     subim2=(subim1-mean(mean(subim2)))/max(max(abs(subim2)));
    
    % Apply hamming window
    subim1=subim1.*hamm;
    subim2=subim2.*hamm;
     
    % Apply phase correlation to the subimages
       
    ft1=fft2(subim1);
    ft2=fft2(subim2);
    xcorr=conj(ft1).*ft2;
    ab=abs(xcorr);
    ab(ab==0)=1;
    xcorr=xcorr./ab;
    %xcorr(1,1)=0;%??????
    ifft_corr=ifft2(xcorr,'symmetric');
    icorr=fftshift(ifft_corr);
    
    % Find the maximum peak
    [maxy, maxpos]=max(sqrt(hamm).*icorr); % 15-Oct-2014 Added .*sqrt(hamm) This should focus on peaks near the middle.
    [Cmax, x0]=max(maxy);
    y0=maxpos(x0);
    halfrange=floor((winsize/2)+1);

    
    switch subpix
        case 'none'
        % Return pixel level displacement only
            displacement=[x0-halfrange(2),halfrange(1)-y0];
            
        case {'foroosh2','foroosh'}
             if x0==1 || x0==winsize(2) || y0==1 || y0==winsize(1)
%                  evalin('caller','roifail=roifail+1;')
%                  roifail=evalin('caller','roifail');
%                 warning('hdm:pixcorr:peakatedge','Peak is at edge of window')
%                 subim=evalin('caller','subim');
%                 subsize=evalin('caller','subsize');
%                 L=evalin('caller','L');
%                 NN=numel(subim.use);
%                 curim=evalin('caller','im1');
%                 curim=repmat(curim,[1 1 3]);
%                 temp.X=reshape(subim.pixX,[NN 1]);
%                 temp.Y=reshape(subim.pixY,[NN 1]);
%                 mask=zeros(size(curim(:,:,1)));
%                 mask(...
%                    (temp.Y(L)-((subsize(1)-1)/2):temp.Y(L)+((subsize(1)-1)/2)),...
%                    (temp.X(L)-((subsize(2)-1)/2):temp.X(L)+((subsize(2)-1)/2)))...
%                    =1;
%                mask=mask*.5;
%                mask(mask==0)=1;
%                curim(:,:,1)=curim(:,:,1).*mask;
%                curim(:,:,3)=curim(:,:,3).*mask;
%                close all;figure(10);
%                subplot(2,2,1),hold on;
%                imshow(subim1)
%                subplot(2,2,2)
%                imshow(subim2)
%                subplot(2,2,3)
%                imshow(curim)
%                subplot(2,2,4)
%                surf(icorr)
%                 
%                 % Peak is on edge of subregion
% 
%                 % Have user select main peak
%                 dcm=datacursormode;
%                 input('Select the main peak then press enter\n');
%                 
%                 info=getCursorInfo(dcm);
%                 close all;
%                 x0=info.Position(1);
%                 y0=info.Position(2);
              end
            % In the x direction:
            % Find maximum side peak
            if x0==1
                cx=icorr(y0,x0+1);
                xsp=2;
            elseif x0==winsize(2)
                cx=icorr(y0,x0-1);
                xsp=winsize(2)-1;
            elseif icorr(y0,x0-1)>icorr(y0,x0+1)
                cx=icorr(y0,x0-1);
                xsp=x0-1;
            else
                cx=icorr(y0,x0+1);
                xsp=x0+1;
            end
            
            % Calculate subpixel displacements
%                 dx1=cx/(cx+Cmax);
%                 dx2=cx/(cx-Cmax);
%             if abs(dx1)<1
%                 dx=abs(dx1)*(xsp-x0);
%             elseif abs(dx2)<1
%                 dx=abs(dx2)*(xsp-x0);
%             elseif abs(dx1)<1 && abs(dx2)<1
%                 disp('PRoblem with foroosh: both dx values meet criteria');
%                 keyboard
%             end
            dx=(cx/(cx+Cmax))*(xsp-x0);
            % In the y direction:
            % Find maximum side peak
            
            if y0==1
                cy=icorr(y0+1,x0);
                ysp=2;
            elseif y0==winsize(1)
                cy=icorr(y0-1,x0);
                ysp=winsize(1)-1;
            elseif icorr(y0-1,x0)>icorr(y0+1,x0)
                cy=icorr(y0-1,x0);
                ysp=y0-1;
            else
                cy=icorr(y0+1,x0);
                ysp=y0+1;
            end
            
            % Calculate subpixel displacements
%                 dx1=cx/(cx+Cmax);
%                 dx2=cx/(cx-Cmax);
%             if abs(dx1)<1
%                 dx=abs(dx1)*(xsp-x0);
%             elseif abs(dx2)<1
%                 dx=abs(dx2)*(xsp-x0);
%             elseif abs(dx1)<1 && abs(dx2)<1
%                 disp('PRoblem with foroosh: both dx values meet criteria');
%                 keyboard
%             end
            dy=(cy/(cy+Cmax))*(ysp-y0);
            displacement=[x0-halfrange(2)+dx,-(y0-halfrange(1)+dy)];
        case 'forooshold'
            if x0==1 || x0==winsize(2) || y0==1 || y0==winsize(2)
                % Peak is on edge of subregion
                close all;
                surf(icorr)
                % Have user select main peak
                dcm=datacursormode;
                input('Select the main peak then press enter\n');
                
                info=getCursorInfo(dcm);
                close all;
                x0=info.Position(2);
                y0=info.Position(1);
            end
            % For the x-direction
            % find the maximum of the two side peaks
            
            if x0==1
                cx=icorr(y0,x0+1);
                xdir=-1;
            elseif x0==winsize(2)
                cx=icorr(y0,x0-1);
                xdir=1;
            else
                if icorr(y0,x0+1)>icorr(y0,x0-1)
                    cx=icorr(y0,x0+1);
                    xdir=1;
                else
                    cx=icorr(y0,x0-1);
                    xdir=-1;
                end
            end
            % Get subpixel displacement
            dx1=cx/(cx+Cmax);
            dx2=cx/(cx-Cmax);
            if xdir==sign(dx1) && xdir==sign(dx2)
                if dx1<1 && dx1>-1
                    dx=dx1;
                elseif dx2<1 && dx2>-1
                    dx=dx2;
                else
                    error('hdm:pixcorr:ambigsubpix',...
                        'Ambiguous foroosh subpixel displacement')
                end
            elseif xdir==sign(dx1)
                dx=dx1;
            elseif xdir==sign(dx2)
                dx=dx2;
            elseif dx1==0 && dx2==0
                dx=0;
            end
            
            % For the y-direction
            % find the maximum of the two side peaks
            if y0==1
                cy=icorr(y0+1,x0);
                ydir=1;
            elseif y0==winsize(1)
                cy=icorr(y0-1,x0);
                ydir=-1;
            else
                if icorr(y0+1,x0)>icorr(y0-1,x0)
                    cy=icorr(y0+1,x0);
                    ydir=1;
                else
                    cy=icorr(y0-1,x0);
                    ydir=-1;
                end
            end
            % Get subpixel displacement
            dy1=cy/(cy+Cmax);
            dy2=cy/(cy-Cmax);
            if ydir==sign(dy1) && ydir==sign(dy2)
                if dy1<1 && dy1>-1
                    dy=dy1;
                elseif dy2<1 && dy2>-1
                    dy=dy2;
                else
                    error('hdm:pixcorr:ambigsubpix',...
                        'Ambiguous foroosh subpixel displacement')
                end
            elseif ydir==sign(dy1)
                dy=dy1;
            elseif ydir==sign(dy2)
                dy=dy2;
            elseif dy1==0 && dy2==0
                dy=0;
            end
            displacement=[(x0-halfrange(2))+dx,-((y0-halfrange(1))+dy)];
                
        case 'hdm'
        % Determine the subpixel level displacement using the method 
        % described by kelly et al
            % Find larger of the two side peaks in x
            if x0 ~= 1 && x0 ~=winsize(2)
                if icorr(y0,x0+1)>icorr(y0,x0-1)
                    cx=icorr(y0,x0+1);
                    x=x0+1;
                    xpos=1;
                else
                    cx=icorr(y0,x0-1);
                    x=x0-1;
                    xpos=0;
                end
            elseif x0==1
                cx=icorr(y0,x0+1);
                x=x0+1;
                xpos=1;
            else
                cx=icorr(y0,x0-1);
                x=x0-1;
                xpos=0;
            end
            
            % Find larger of the two side peaks in y
            if y0 ~= 1 && y0 ~=winsize(1)
                if icorr(y0+1,x0)>icorr(y0-1,x0)
                    cy=icorr(y0+1,x0);
                    y=y0+1;
                    ypos=1;
                else
                    cy=icorr(y0-1,x0);
                    y=y0-1;
                    ypos=0;
                end
            elseif y0==1
                cy=icorr(y0+1,x0);
                y=y0+1;
                ypos=1;
            else
                cy=icorr(y0-1,x0);
                y=y0-1;
                ypos=0;
            end
            
            % Calculate two solutions in x:
            tempx1= (x*cx-x0*Cmax)/(cx-Cmax);
            tempx2= (x*cx+x0*Cmax)/(cx+Cmax);
            if(xpos==1)
                if tempx1<(x0+1) && tempx1>x0
                    disx=tempx1;
                end
                if tempx2<(x0+1) && tempx2>x0
                    disx=tempx2;
                end 
            end
            if(xpos==0)
                if tempx1<(x0) && tempx1>x0-1
                    disx=tempx1;
                end
                if tempx2<(x0) && tempx2>x0-1
                    disx=tempx2;
                end 
            end
            
            
             % Calculate two solutions in y:
            tempy1= (y*cy-y0*Cmax)/(cy-Cmax);
            tempy2= (y*cy+y0*Cmax)/(cy+Cmax);
            if(ypos==1)
                if tempy1<(y0+1) && tempy1>y0
                    disy=tempy1;
                end
                if tempy2<(y0+1) && tempy2>y0
                    disy=tempy2;
                end 
            end
            if(ypos==0)
                if tempy1<(y0) && tempy1>y0-1
                    disy=tempy1;
                end
                if tempy2<(y0) && tempy2>y0-1
                    disy=tempy2;
                end 
            end
            
            displacement=[disx-size(subim1,2)/2+1,disy-size(subim1,1)/2+1];            
        case 'esinc'
            % Calculate subpixel displacements using methods described by
            % Argyriou et al
            
            % Set up fitting options
            lb=[0 0 -.8];
            ub=[6 1 .8];
            start=[.6 .6 .001];
            
                % Start with the default options
                options = optimset;
                % Modify options setting
                options = optimset(options,'Display', 'off');
                options = optimset(options,'Algorithm', 'interior-point');

            
            
            % Determine x direction displacement
            vals=x0-1:x0+1;
            c=icorr(y0,vals);
            params = ...
            fmincon(@(x)esinc_fit(x,c),start,[],[],[],[],lb,ub,[],options);
            dx=params(3);
            % Graph X outcome
%             close all;
%             xvals=-1:.01:1;
%             figure
%             subplot(1,2,1)
%             plot([-1 0 1],c,'r*')
%             hold on;
%             plot(xvals,esinc(params,xvals),'k-')
%             plot(params(3),esinc(params,params(3)),'b+')
            % Determine y direction displacement
            vals=y0-1:y0+1;
            c=icorr(vals,x0)';
            params = ...
            fmincon(@(x)esinc_fit(x,c),start,[],[],[],[],lb,ub,[],options);
            dy=params(3);
           % Graph Y outcome
%             subplot(1,2,2)
%             plot([-1 0 1],c,'r*')
%             hold on;
%             plot(xvals,esinc(params,xvals),'k-')
%             plot(params(3),esinc(params,params(3)),'b+')
            % output displacement
            displacement=[(x0-halfrange(2))+dx,-((y0-halfrange(1))+dy)];
        case 'favreau1'
            % Isolate 3x3 region around the main peak
            tarreg=icorr(y0-1:y0+1,x0-1:x0+1);
            % Normalize the region around the peak
            tarreg=tarreg/tarreg(2,2);
            % Find the dominant quadrant
            tarline=reshape(tarreg,[1 9]);
            quads=[sum(tarline([1 2 4])),sum(tarline([4 7 8])),...
                sum(tarline([6 9 8])),sum(tarline([2 3 6]))];
            [~,q]=max(quads);
            
            switch q
                case 1
                    P=tarline([2 1 4 5]);
                case 2
                    P=tarline([5 4 7 8]);
                case 3
                    P=tarline([6 5 8 9]);
                case 4
                    P=tarline([3 2 5 6]);
            end
            dx(1)=(P(3)+P(4))/(sum(P));
            dx(2)=(P(2)+P(3))/(sum(P));
            switch q
                case 1
                    dx(1)=dx(1)-1;
                case 3
                    dx(2)=dx(2)-1;
                case 4
                    dx(1)=dx(1)-1;
                    dx(2)=dx(2)-1;
            end
            if any(abs(dx)>0.5)
                %keyboard
            end
            displacement=[(x0-halfrange(2))+dx(1),-((y0-halfrange(1))-dx(2))];        

            
    end
    switch nargout
        case 1
            varargout{1}=[1 -1].*displacement; % Downward direction is positive
        case 2
            varargout{1}=displacement(1);
            varargout{2}=-1*displacement(2);% Downward direction is positive
    end
    
end

function [varargout]=roicorr(im1,im2,subsize,subshift,varargin)
% roicorr: Function to correlate subimages in a region of interest. First 
%   divides a region of interest into all subimages that can fit inside.
%   Next, finds the rigid body shift of the region. Finally, calculates
%   pixel or subpixel level shifts for each region and outputs the
%   displacement data.
%
%   REQUIRED INPUTS:
%       im1: integer or double image to be analyzed
%       im2: integer or double image to be compared to im1
%       subsize: Size of subimage to analyze. Can be a single number or an
%           array of two numbers. Single number indicates a square subimage
%           while an array indicates [ysize,xsize]
%       subshift: Pixel shift between subimages. follows same convention as
%           subsize.
%
%   OPTIONAL INPUTS:
%       roi: Region of interest to analyze. Can eiter be an integer 0
%           (requests user input) or the output of the function hdmgetroi.
%       subpix: Method to calculate subpixel displacements. Can be foroosh
%           (default), none, esinc, or favreau1
%       roitype: The shape of the roi to be defined. Only used if roi=0.
%           Can be: 'impoly'(default), 'imrect' or 'imfreehand'
%       rigidbodymotion: Numeric input defining if rigid body motion should
%           be calculated(1,default), ignored(0), or is provided(2)
%       importrbm: Only used if rigidbodymotion=2. Uses the specified U and
%           V shifts as the rigid body motion.
%
%   (C)2013 John T Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute
%   --
%% Check function input for accuracy
    tic;
    %roifail=evalin('caller','roifail');
    % Parse input
    p=inputParser;
    imvalid=@(x) (isa(x,'double')...
        && max(max(x))<=1 && min(min(x))>= 0) || (isa(x,'integer'));
    p.addRequired('im1',imvalid);
    p.addRequired('im2',imvalid);
    p.addRequired('subsize',@isnumeric)
    p.addRequired('subshift',@isnumeric)
    p.addOptional('roi',0)
    p.addParamValue('subpix','foroosh',...
        @(x) any(strcmp(x,{'none','esinc','foroosh','favreau1'})));
    p.addParamValue('roitype','impoly',...
          @(x) any(strcmp(x,{'impoly','imrect','imsquare','imcircle',...
          'imellipse'})));
    p.addParamValue('rigidbodymotion',1,@(x) x==0||x==1||x==2)
    p.addParamValue('importrbm',[0 0 0 0],@(x) length(x)==4)
    p.parse(im1,im2,subsize,subshift,varargin{:});
    
    im1=p.Results.im1;
    im2=p.Results.im2;
    im1=doubleim(im1);
    im2=doubleim(im2);
    
    subsize=p.Results.subsize;
    if length(subsize)==1
        subsize=[subsize,subsize];
    end
    subshift=p.Results.subshift;
    if length(subshift)==1;
        subshift=[subshift,subshift];
    end
    if any(mod(subsize,subshift))
        error('hdm:roicorr:invalidshift',...
            'Subimage size must be evenly divisible by shift size')
    end
    subpix=p.Results.subpix;
    roitype=p.Results.roitype;
    rbm=p.Results.rigidbodymotion;
    importrbm=p.Results.importrbm;
    
%% Have user select ROI
    % Process ROI, have user select it if it is not specified
    if p.Results.roi==0
        roi=hdmgetroi(im1,subsize,subshift,roitype);
        disp(roi)
    else
        roi=p.Results.roi;
    end
    if size(im1)~=size(im2)
        error('hdm:roicorr:sizemismatch','Subimages must be the same size')
    end
    
%% Populate subimages in the ROI
    
    subim=getsubimages(roi,im1,subsize,subshift);
    
    
%% Compute rigid body motion
switch rbm
    case 1
        % Get rectangle within ROI and centered on centroid
        c=subim.centroid;
        incX=1;
        incY=1;
        % Increase square size until a boundary is exceeded
        while all(subim.pixelmask(c(2)-incY:c(2)+incY,c(1)-incX:c(1)+incX))
            incX=incX+1;
            incY=incY+1;
        end
        incX=incX-1;
        incY=incY-1;

        % Increase vertically until boundary is exceeded
        while all(subim.pixelmask(c(2)-incY:c(2)+incY,c(1)-incX:c(1)+incX))
            incY=incY+1;
        end
        incY=incY-1;

        % Increase horizontally until a boundary is exceeded
        while all(subim.pixelmask(c(2)-incY:c(2)+incY,c(1)-incX:c(1)+incX))
            incX=incX+1;
        end
        incX=incX-1;

        % Compute rigid body displacements to the nearest pixel
        [RU,RV]=pixcorr(im1(c(2)-incY:c(2)+incY,c(1)-incX:c(1)+incX),...
            im2(c(2)-incY:c(2)+incY,c(1)-incX:c(1)+incX),'none');
        rawhdm.xxyyshift=[RU RU RV RV]; % Positive Y direction indicates downward motion
    case 0
        rawhdm.xxyyshift=[0 0 0 0];
        RU=0;RV=0;
    case 2
        rawhdm.xxyyshift=importrbm;
        RU=importrbm(1);
        RV=importrbm(3);
    otherwise
        error('HDM2:roicorr:badrbm',['The rigid body motion value spec'...
            ,'ified is invalid'])
end
    
%% Get frame rois for analysis
    Ns=numel(subim.use);
    curcorr=getframes(RU,RV);
    curcorr.udata=zeros([Ns,1]);
    curcorr.vdata=curcorr.udata;
    
%% Compute subimage displacements
    % Loop through subimages and determine displacements
    N=size(curcorr.frame1,3);
    tic;
    tempuse=reshape(subim.use,[Ns, 1]);
    frame1=curcorr.frame1;
    frame2=curcorr.frame2;
    udata=curcorr.udata;
    vdata=curcorr.vdata;
    m=0;
    for L=1:Ns
       if tempuse(L)
           m=m+1;
           [udata(L),vdata(L)]=pixcorr(frame1(:,:,m),frame2(:,:,m),subpix);
%            %TEMP
%            tempu=reshape(udata,subim.shape);
%            surf(tempu)
%            keyboard
       else
           udata(L)=NaN;
           vdata(L)=NaN;
       end
    end
    udata=udata+RU;
    vdata=vdata+RV;
    rawhdm.udata=reshape(udata,subim.shape);
    rawhdm.vdata=reshape(vdata,subim.shape);

    subim.comptime=toc;
    subim.count=Ns;
    varargout{1}=rawhdm;
    varargout{2}=subim;
    if nargout > 2
        varargout{3}=roi;
    end
    
    
%% HELPER FUNCTIONS (EMBEDDED)

    function getframeout=getframes(RU,RV)
       % Uses the information in subim struct to determine the target ROIs
       N=numel(subim.use);
       temp.use=reshape(subim.use,[N 1]);
       temp.X=reshape(subim.pixX,[N 1]);
       temp.Y=reshape(subim.pixY,[N 1]);
       getframeout.frame1=zeros([subsize,numel(subim.use(subim.use==1))]);
       getframeout.frame2=getframeout.frame1;
       j=0;
       for i = 1:N
           if temp.use(i)==1
               j=j+1;
               getframeout.frame1(:,:,j)=im1(...
                   (temp.Y(i)-((subsize(1)-1)/2):temp.Y(i)+((subsize(1)-1)/2)),...
                   (temp.X(i)-((subsize(2)-1)/2):temp.X(i)+((subsize(2)-1)/2)));
               getframeout.frame2(:,:,j)=im2(...
                   RV+(temp.Y(i)-((subsize(1)-1)/2):temp.Y(i)+((subsize(1)-1)/2)),...
                   RU+(temp.X(i)-((subsize(2)-1)/2):temp.X(i)+((subsize(2)-1)/2)));
               %^^^ im2 has the rigid body motion added in
           end
       end
               
    end
end



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
        % Request that the user choose the location of the raw files
        if nargin==4
            [tarfile,tardir]=uigetfile(...
                ['',ext],...
            'Select the target file');
            list=dir([tardir,tarfile]);
        else
            tardir=uigetdir('Z:\Megan OBrien\Raw Data','Select the folder');
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

