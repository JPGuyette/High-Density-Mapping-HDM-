function gcamp_intensity2_sa(dataset,pos,frames)
% Compiled by ADMIN\kjhansen on 01/07/2017 13:56:01
% gcamp_intensity2 Version 2015.04.0002
% 
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
%    tofile=[ind2rgb(uint8(perstack(:,:,i)*255),parula),...
%         getcolorbar([0,1],size(curim,1),'parula')];
    q=ind2rgb(uint8(perstack(:,:,i)*255),parula);
    w=getcolorbar([0,1],size(curim,1),'parula');
    tofile=[q,w];
   imwrite(tofile,[top_dir,...
       'intensitymaps\cycle_',num2str(i),'.tif'])
end
    tofile=[ind2rgb(uint8(mean(perstack(:,:,i),3)*255),parula),...
        getcolorbar([0,1],size(curim,1),'parula')];
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


peakvals = max(imstack,[],3);
minvals = min(imstack,[],3);
avgvals= mean(imstack,[],3);
F=(avgvals)./minvals;
imshow(avgvals)

Fstack = imstack-repmat(minvals,[1 1 n]);
y=Fstack(:,:,frames);%to pick a certain amount of frames
[~,Ftime] = max(y,[],3);
imagesc(Ftime)
K = filter2(fspecial('average',9),Ftime);
imagesc(K)
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





function mapout=jtfheatmap(mapname,size)
    if strcmp(mapname(1),'-')
        if strcmp(mapname(2),'-')
            mn=mapname(3:end);
            mapname=mapname(3:end);
        else
            mn=mapname(2:end);
        end
    else
        mn=mapname;
    end
    
    mapout=eval([mn,'(',num2str(size),');']);
    
    if strcmp(mapname(1),'-')
        mapout=flipud(mapout);
    end
end

function varargout=linreg(x,y)
% [m, b] = linreg(u,x)
% [m, b, R2] = linreg(u,x)
%   Outputs the slope of x w.r.t. y using a computationally efficient
%   matrix algebr approach.
%
%   INPUT:
%       x: x values
%       y: y values
%
%   OUTPUT:
%       m: Slope of line
%       b: y-intercept of line
%       R2: (Optional) R-squared value for goodness of fit
%
%

% Check input
    x=squeezecol(x);
    y=squeezecol(y);
    if numel(x)>length(x) || numel(y)>length(y)
        error('Inputs must be arrays')
    end
    
    if length(x) ~= length(y)
        error('Input vectors must be the same length')
    end
    
    
% Calculate necessary values
    n=numel(x);
    m=(n*sum(x.*y)-sum(x).*sum(y))/(n*sum(x.^2)-(sum(x)).^2);
    b=(1/n)*(sum(y)-m*sum(x));
    SSerr=sum((y-(m*x+b)).^2);
    SStot=sum((y-mean(y)).^2);
    R2=1-(SSerr/SStot);
    
    switch nargout
        case {0,1}
            out.m=m;
            out.b=b;
            out.R2=R2;
            varargout={out};
        case 2
            varargout={m,b};
        case 3
            varargout={m,b,R2};
        otherwise
            error('Function has too many outputs');
    end
end

function out=squeezecol(in)
% Forces 1d variable to be a column array
    in=squeeze(in);
    if size(in,2)>1
        in=in';
    end
    if size(in,2)>1
        error('Input must only have data in one dimension')
    end
    out=in;
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
        tardir=uigetdir('',...
            'Select directory with the target tif files');
        list=dir([tardir,'\*.tif']);
        mkfilelist(dataset,tardir,list,'tif');
    end
end

function map = parula(n)
%PARULA Blue-green-orange-yellow color map
%   PARULA(M) returns an M-by-3 matrix containing a colormap. 
%   The colors begin with dark purplish-blue and blue, range
%   through green and orange, and end with bright yellow. PARULA is named
%   after a bird, the tropical parula, which has these colors.
%
%   PARULA returns a colormap with the same number of colors as the current
%   figure's colormap. If no figure exists, MATLAB creates one.
%
%   EXAMPLE
%
%   This example shows how to reset the colormap of the current figure.
%
%       colormap(parula)
%
%   See also AUTUMN, BONE, COLORCUBE, COOL, COPPER, FLAG, GRAY, HOT, HSV,
%   JET, LINES, PINK, PRISM, SPRING, SUMMER, WHITE, WINTER, COLORMAP,
%   RGBPLOT.

%   Copyright 2013 The MathWorks, Inc.

if nargin < 1
   n = size(get(gcf,'Colormap'),1);
end

values = [
   0.2081 0.1663 0.5292
   0.2091 0.1721 0.5411
   0.2101 0.1779 0.5530
   0.2109 0.1837 0.5650
   0.2116 0.1895 0.5771
   0.2121 0.1954 0.5892
   0.2124 0.2013 0.6013
   0.2125 0.2072 0.6135
   0.2123 0.2132 0.6258
   0.2118 0.2192 0.6381
   0.2111 0.2253 0.6505
   0.2099 0.2315 0.6629
   0.2084 0.2377 0.6753
   0.2063 0.2440 0.6878
   0.2038 0.2503 0.7003
   0.2006 0.2568 0.7129
   0.1968 0.2632 0.7255
   0.1921 0.2698 0.7381
   0.1867 0.2764 0.7507
   0.1802 0.2832 0.7634
   0.1728 0.2902 0.7762
   0.1641 0.2975 0.7890
   0.1541 0.3052 0.8017
   0.1427 0.3132 0.8145
   0.1295 0.3217 0.8269
   0.1147 0.3306 0.8387
   0.0986 0.3397 0.8495
   0.0816 0.3486 0.8588
   0.0646 0.3572 0.8664
   0.0482 0.3651 0.8722
   0.0329 0.3724 0.8765
   0.0213 0.3792 0.8796
   0.0136 0.3853 0.8815
   0.0086 0.3911 0.8827
   0.0060 0.3965 0.8833
   0.0051 0.4017 0.8834
   0.0054 0.4066 0.8831
   0.0067 0.4113 0.8825
   0.0089 0.4159 0.8816
   0.0116 0.4203 0.8805
   0.0148 0.4246 0.8793
   0.0184 0.4288 0.8779
   0.0223 0.4329 0.8763
   0.0264 0.4370 0.8747
   0.0306 0.4410 0.8729
   0.0349 0.4449 0.8711
   0.0394 0.4488 0.8692
   0.0437 0.4526 0.8672
   0.0477 0.4564 0.8652
   0.0514 0.4602 0.8632
   0.0549 0.4640 0.8611
   0.0582 0.4677 0.8589
   0.0612 0.4714 0.8568
   0.0640 0.4751 0.8546
   0.0666 0.4788 0.8525
   0.0689 0.4825 0.8503
   0.0710 0.4862 0.8481
   0.0729 0.4899 0.8460
   0.0746 0.4937 0.8439
   0.0761 0.4974 0.8418
   0.0773 0.5012 0.8398
   0.0782 0.5051 0.8378
   0.0789 0.5089 0.8359
   0.0794 0.5129 0.8341
   0.0795 0.5169 0.8324
   0.0793 0.5210 0.8308
   0.0788 0.5251 0.8293
   0.0778 0.5295 0.8280
   0.0764 0.5339 0.8270
   0.0746 0.5384 0.8261
   0.0724 0.5431 0.8253
   0.0698 0.5479 0.8247
   0.0668 0.5527 0.8243
   0.0636 0.5577 0.8239
   0.0600 0.5627 0.8237
   0.0562 0.5677 0.8234
   0.0523 0.5727 0.8231
   0.0484 0.5777 0.8228
   0.0445 0.5826 0.8223
   0.0408 0.5874 0.8217
   0.0372 0.5922 0.8209
   0.0342 0.5968 0.8198
   0.0317 0.6012 0.8186
   0.0296 0.6055 0.8171
   0.0279 0.6097 0.8154
   0.0265 0.6137 0.8135
   0.0255 0.6176 0.8114
   0.0248 0.6214 0.8091
   0.0243 0.6250 0.8066
   0.0239 0.6285 0.8039
   0.0237 0.6319 0.8010
   0.0235 0.6352 0.7980
   0.0233 0.6384 0.7948
   0.0231 0.6415 0.7916
   0.0230 0.6445 0.7881
   0.0229 0.6474 0.7846
   0.0227 0.6503 0.7810
   0.0227 0.6531 0.7773
   0.0232 0.6558 0.7735
   0.0238 0.6585 0.7696
   0.0246 0.6611 0.7656
   0.0263 0.6637 0.7615
   0.0282 0.6663 0.7574
   0.0306 0.6688 0.7532
   0.0338 0.6712 0.7490
   0.0373 0.6737 0.7446
   0.0418 0.6761 0.7402
   0.0467 0.6784 0.7358
   0.0516 0.6808 0.7313
   0.0574 0.6831 0.7267
   0.0629 0.6854 0.7221
   0.0692 0.6877 0.7173
   0.0755 0.6899 0.7126
   0.0820 0.6921 0.7078
   0.0889 0.6943 0.7029
   0.0956 0.6965 0.6979
   0.1031 0.6986 0.6929
   0.1104 0.7007 0.6878
   0.1180 0.7028 0.6827
   0.1258 0.7049 0.6775
   0.1335 0.7069 0.6723
   0.1418 0.7089 0.6669
   0.1499 0.7109 0.6616
   0.1585 0.7129 0.6561
   0.1671 0.7148 0.6507
   0.1758 0.7168 0.6451
   0.1849 0.7186 0.6395
   0.1938 0.7205 0.6338
   0.2033 0.7223 0.6281
   0.2128 0.7241 0.6223
   0.2224 0.7259 0.6165
   0.2324 0.7275 0.6107
   0.2423 0.7292 0.6048
   0.2527 0.7308 0.5988
   0.2631 0.7324 0.5929
   0.2735 0.7339 0.5869
   0.2845 0.7354 0.5809
   0.2953 0.7368 0.5749
   0.3064 0.7381 0.5689
   0.3177 0.7394 0.5630
   0.3289 0.7406 0.5570
   0.3405 0.7417 0.5512
   0.3520 0.7428 0.5453
   0.3635 0.7438 0.5396
   0.3753 0.7446 0.5339
   0.3869 0.7454 0.5283
   0.3986 0.7461 0.5229
   0.4103 0.7467 0.5175
   0.4218 0.7473 0.5123
   0.4334 0.7477 0.5072
   0.4447 0.7482 0.5021
   0.4561 0.7485 0.4972
   0.4672 0.7487 0.4924
   0.4783 0.7489 0.4877
   0.4892 0.7491 0.4831
   0.5000 0.7491 0.4786
   0.5106 0.7492 0.4741
   0.5212 0.7492 0.4698
   0.5315 0.7491 0.4655
   0.5418 0.7490 0.4613
   0.5519 0.7489 0.4571
   0.5619 0.7487 0.4531
   0.5718 0.7485 0.4490
   0.5816 0.7482 0.4451
   0.5913 0.7479 0.4412
   0.6009 0.7476 0.4374
   0.6103 0.7473 0.4335
   0.6197 0.7469 0.4298
   0.6290 0.7465 0.4261
   0.6382 0.7460 0.4224
   0.6473 0.7456 0.4188
   0.6564 0.7451 0.4152
   0.6653 0.7446 0.4116
   0.6742 0.7441 0.4081
   0.6830 0.7435 0.4046
   0.6918 0.7430 0.4011
   0.7004 0.7424 0.3976
   0.7091 0.7418 0.3942
   0.7176 0.7412 0.3908
   0.7261 0.7405 0.3874
   0.7346 0.7399 0.3840
   0.7430 0.7392 0.3806
   0.7513 0.7385 0.3773
   0.7596 0.7378 0.3739
   0.7679 0.7372 0.3706
   0.7761 0.7364 0.3673
   0.7843 0.7357 0.3639
   0.7924 0.7350 0.3606
   0.8005 0.7343 0.3573
   0.8085 0.7336 0.3539
   0.8166 0.7329 0.3506
   0.8246 0.7322 0.3472
   0.8325 0.7315 0.3438
   0.8405 0.7308 0.3404
   0.8484 0.7301 0.3370
   0.8563 0.7294 0.3336
   0.8642 0.7288 0.3300
   0.8720 0.7282 0.3265
   0.8798 0.7276 0.3229
   0.8877 0.7271 0.3193
   0.8954 0.7266 0.3156
   0.9032 0.7262 0.3117
   0.9110 0.7259 0.3078
   0.9187 0.7256 0.3038
   0.9264 0.7256 0.2996
   0.9341 0.7256 0.2953
   0.9417 0.7259 0.2907
   0.9493 0.7264 0.2859
   0.9567 0.7273 0.2808
   0.9639 0.7285 0.2754
   0.9708 0.7303 0.2696
   0.9773 0.7326 0.2634
   0.9831 0.7355 0.2570
   0.9882 0.7390 0.2504
   0.9922 0.7431 0.2437
   0.9952 0.7476 0.2373
   0.9973 0.7524 0.2310
   0.9986 0.7573 0.2251
   0.9991 0.7624 0.2195
   0.9990 0.7675 0.2141
   0.9985 0.7726 0.2090
   0.9976 0.7778 0.2042
   0.9964 0.7829 0.1995
   0.9950 0.7880 0.1949
   0.9933 0.7931 0.1905
   0.9914 0.7981 0.1863
   0.9894 0.8032 0.1821
   0.9873 0.8083 0.1780
   0.9851 0.8133 0.1740
   0.9828 0.8184 0.1700
   0.9805 0.8235 0.1661
   0.9782 0.8286 0.1622
   0.9759 0.8337 0.1583
   0.9736 0.8389 0.1544
   0.9713 0.8441 0.1505
   0.9692 0.8494 0.1465
   0.9672 0.8548 0.1425
   0.9654 0.8603 0.1385
   0.9638 0.8659 0.1343
   0.9623 0.8716 0.1301
   0.9611 0.8774 0.1258
   0.9600 0.8834 0.1215
   0.9593 0.8895 0.1171
   0.9588 0.8958 0.1126
   0.9586 0.9022 0.1082
   0.9587 0.9088 0.1036
   0.9591 0.9155 0.0990
   0.9599 0.9225 0.0944
   0.9610 0.9296 0.0897
   0.9624 0.9368 0.0850
   0.9641 0.9443 0.0802
   0.9662 0.9518 0.0753
   0.9685 0.9595 0.0703
   0.9710 0.9673 0.0651
   0.9736 0.9752 0.0597
   0.9763 0.9831 0.0538
   ];

P = size(values,1);
map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');


function varargout=peakfind1(signal,Fs,passrange)
%   peakfind1 Finds peaks in a highly repetitive time sequence. This 
%   analysis is based on a fourier transform yielding the dominant
%   frequency in the time series. 
%
%   2013 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%^^^%     end
    if ~exist('passrange','var')

%     if ~exist('target','var')
%         target='both';
        passrange=[0 Inf];
    end
    % Remove linear trend from signal
    x=(1:length(signal))';
    [m,b]=linreg(x,signal);
    signal=signal-(m*x+b);
    
    % Remove DC component of signal
    signal=signal-mean(signal);
    
    % Take fourier transform and extract the dominant frequency
    [P,freq]=pfft(signal,Fs);
    P(freq<passrange(1)|freq>passrange(2))=0;
    maxfreq=freq(P==max(P));
    %disp(maxfreq)
    N=length(signal);
    win=2*round(Fs/(2*maxfreq));
    w=ceil(0.9*(win/2));
    
    % Zero pad signal
    signal2=[zeros(w,1);signal;zeros(w,1)];
    N2=length(signal2);
    
    
    maxindex=false(N,1);
    minindex=false(N,1);
    s=0;
    for i=1+w:1:N2-w
        j=i-w;
        if signal2(i)==max(signal2(i-w:i+w)) && s~=1
            maxindex(j)=true;
            s=1;
        end
        if signal2(i)==min(signal2(i-w:i+w)) && s~=-1
            minindex(j)=true;
            s=-1;
        end
    end
    peakindex=maxindex | minindex;
    switch nargout
        case {1,0}
            varargout{1}=peakindex;
        case 2
            varargout{1}=maxindex;
            varargout{2}=minindex;
        case 3
            varargout{1}=peakindex;
            varargout{2}=maxindex;
            varargout{3}=minindex;
    end
end

function [X,freq]=pfft(x,Fs)
%   pfft: Returns the positive real portion of the fft of a 1D signal
%   [X,freq]=pfft(x,Fs)
%   x: initial signal
%   Fs: Signal sampling frequency
%   X: Fourier domain of series
%   Corresponding frequency spectrum of fourier domain series
%
%   2012 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%^^^
    % Detrend signal
    x=x-mean(x);
    N=length(x); %get the number of points
    freq=((0:N/2-1)/(N/2)*(Fs/2))';
    X=fft(x,N);
    X=abs(X(1:length(freq))).^2;
    
end



function imtext=text2im(text,hgt)
if ~exist('hgt','var')
    hgt=20;
end
% text2im - generates an image, containing the input text

text=text+0;                    % converting string into Ascii-number array

laenge=length(text);
imtext=zeros(20,18*laenge);     % Preparing the resulting image

for i=1:laenge
    
    code=text(i);               % Ascii code of i.th letter

    if code==0
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;     %   .
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   .
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   .
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   I       Complete
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   I       List of all
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           % \ I /     256 ASCII
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %  \I/      Letters
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   V
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   .
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   .
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;           %   .
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==1
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==2
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==3
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==4
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==5
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==6
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==7
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==8
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==9
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==10
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==11
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==12
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==13
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==14
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==15
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==16
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==17
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==18
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==19
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==20
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==21
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==22
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==23
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==24
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==25
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==26
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==27
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==28
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==29
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==30
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==31
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==32
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==33
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==34
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==35
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==36
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==37
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==38
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==39
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==40
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==41
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==42
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==43
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==44
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1];
    elseif code==45
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==46
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==47
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==48
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==49
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==50
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==51
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==52
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==53
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==54
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==55
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==56
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==57
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==58
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==59
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1];
    elseif code==60
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==61
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==62
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==63
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==64
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1;
     0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==65
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==66
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==67
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==68
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==69
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==70
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==71
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==72
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==73
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==74
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==75
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==76
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==77
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==78
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==79
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==80
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==81
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1];
    elseif code==82
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==83
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==84
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==85
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==86
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==87
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==88
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==89
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==90
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==91
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==92
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==93
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==94
    TxtIm=[1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==95
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1];
    elseif code==96
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==97
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==98
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==99
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==100
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==101
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==102
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==103
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1];
    elseif code==104
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==105
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==106
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1];
    elseif code==107
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==108
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==109
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==110
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==111
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==112
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==113
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1];
    elseif code==114
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==115
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==116
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==117
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==118
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==119
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==120
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==121
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==122
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==123
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1];
    elseif code==124
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==125
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1];
    elseif code==126
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==127
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==128
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1];
    elseif code==129
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==130
    TxtIm=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==131
    TxtIm=[1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==132
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==133
    TxtIm=[1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==134
    TxtIm=[1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==135
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1];
    elseif code==136
    TxtIm=[1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==137
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==138
    TxtIm=[1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==139
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==140
    TxtIm=[1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==141
    TxtIm=[1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==142
    TxtIm=[0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==143
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==144
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==145
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==146
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==147
    TxtIm=[1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==148
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==149
    TxtIm=[1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==150
    TxtIm=[1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==151
    TxtIm=[1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==152
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==153
    TxtIm=[1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==154
    TxtIm=[1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==155
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1];
    elseif code==156
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     0,0,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==157
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==158
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==159
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==160
    TxtIm=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==161
    TxtIm=[1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==162
    TxtIm=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==163
    TxtIm=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==164
    TxtIm=[1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==165
    TxtIm=[1,1,1,1,0,0,0,0,0,0,1,1,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==166
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==167
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==168
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==169
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==170
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==171
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==172
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==173
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==174
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==175
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==248
    TxtIm=[1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1;
     0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1];
    elseif code==177
    TxtIm=[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1;
     0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1];
    elseif code==178
    TxtIm=[0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1];
    elseif code==179
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==180
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==181
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==182
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==183
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==184
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==185
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==186
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==187
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==188
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==189
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==190
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==191
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==192
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==193
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==194
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==195
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==196
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==197
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==198
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==199
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==200
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==201
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==202
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==203
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==204
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==205
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==206
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==207
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==208
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==209
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==210
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==211
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==212
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==213
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==214
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==215
    TxtIm=[1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1];
    elseif code==216
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==217
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==218
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==219
    TxtIm=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    elseif code==220
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    elseif code==221
    TxtIm=[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1];
    elseif code==222
    TxtIm=[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0];
    elseif code==223
    TxtIm=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==224
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==225
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==226
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==227
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==228
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==229
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==230
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==231
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==232
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==233
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==234
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==235
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==236
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==237
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==238
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==239
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==240
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==241
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==242
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==243
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==244
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1];
    elseif code==245
    TxtIm=[1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==246
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==247
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==176
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==249
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==250
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==251
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==252
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==253
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==254
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    elseif code==255
    TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    end

    imtext(1:20,((i-1)*18+1):i*18)=TxtIm;
end

if hgt ~= 20
    scl=hgt/20;
    imtext=round(imresize(imtext,scl));
end

