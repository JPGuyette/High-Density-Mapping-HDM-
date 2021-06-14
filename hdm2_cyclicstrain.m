function varargout = hdm2_cyclicstrain_sa(dataset,field,method,varargin)
% Compiled by BIOMED-NB03\jfavreau on 02/20/2015 14:03:33
% hdm2_cyclicstrain Version 2015.01.0007
% 

    if nargin==4
        freqpass=varargin{1};
    else
        freqpass=[0, Inf];
    end
%% Load data
    if ischar(dataset)
    top_dir=[dataset,'\'];
        if exist([top_dir,'matlab_data\rawhdm2.mat'],'file')
            load([top_dir,'matlab_data\rawhdm2.mat']);
        elseif exist([top_dir,'matlab_data\rawhdm.mat'],'file')
            rawhdm.subshift=rawhdm.shift; %#ok<NODEF>
        else
            error('hdm or hdm2 data do not exist for this dataset')
        end
    else
        error('hdm_greenstrain_cyclic is not configured to analyze this dataset')
    end
    if ~exist([top_dir,'matlab_data\avg_strain_tensor.mat'],'file')
        hdm2_avgstrain(dataset,0)
    end
    load([top_dir,'matlab_data\avg_strain_tensor.mat'])
    
%% Process inputs
    %[y,x,time]=size(rawhdm.udata);
%     y=1:y;
%     x=1:x;
    switch lower(field)
        case {'eyy','ey','e22','y'}
            fielddata=squeeze(strain.green(2,2,:));
        case {'exx','ex','e11','x'}
            fielddata=squeeze(strain.green(1,1,:));
        case {'exy','eyx','e12','e21','xy'}
            fielddata=squeeze(strain.green(2,2,:));
        case {'area','a'};
            fielddata=strain.area;
        case {'e1'}
            Exx=squeeze(strain.green(1,1,:));
            Eyy=squeeze(strain.green(2,2,:));
            Exy=squeeze(strain.green(1,2,:));
            fielddata=((Exx+Eyy)/2)+sqrt(((Exx-Eyy)/2).^2+Exy.^2);
        case {'e2'}
            Exx=squeeze(strain.green(1,1,:));
            Eyy=squeeze(strain.green(2,2,:));
            Exy=squeeze(strain.green(1,2,:));
            fielddata=((Exx+Eyy)/2)-sqrt(((Exx-Eyy)/2).^2+Exy.^2);
        otherwise
            error('Invalid field input')
    end
    switch method
        case 'userselect'
            error('This input has not been implemented')
        case 'userfreq'
            error('This input has not been implemented')
        case 'maxfreqmins'
            % Calculate the maximum frequency and then find minima
            [pkpoints,zeropoints]=peakfind1(fielddata,rawhdm.framerate,freqpass);
        case 'maxfreqmaxs'
            [zeropoints,pkpoints]=peakfind1(fielddata,rawhdm.framerate,freqpass);
        otherwise
            error('This input has not been implemented')
    end
    x=find(pkpoints);
    y=find(zeropoints);
    if x(1)<y(1)
        pkpoints(x(1))=0;
    end
    if x(end)>y(end)
        pkpoints(x(end))=0;
    end
        

%% Preallocate variables
   
    Ncyc=numel(find(zeropoints))-1;
    strainsum=cell(Ncyc,7);
    pk=find(zeropoints);
    cycstrain=zeros(2,2,pk(Ncyc+1)-pk(1)+1);
    cycarea=zeros(pk(Ncyc+1)-pk(1)+1,1);
    % loop through data
    for c= 1:Ncyc
        
        for f=pk(c):pk(c+1)
            F=[1 0;0 1;];
            A=1;
            for ff=(pk(c)+1):f
                F=F*incgrad.defgrad(:,:,ff); % iterate deformation gradients
                A=A*incgrad.astretch(ff,1);% iterate area stretch ratios
            end
            cycstrain(:,:,f)= 0.5*(F'*F-eye(2)); % Calculate green's strain tensor
            cycarea(f,1)=A-1;
        end
        Exx=cycstrain(1,1,pk(c):pk(c+1));
        Eyy=cycstrain(2,2,pk(c):pk(c+1));
        Exy=cycstrain(1,2,pk(c):pk(c+1));
        E1=((Exx+Eyy)/2)+sqrt(((Exx-Eyy)/2).^2+Exy.^2);
        E2=((Exx+Eyy)/2)-sqrt(((Exx-Eyy)/2).^2+Exy.^2);
        strainsum(c,:)={c,max(cycstrain(1,1,pk(c):pk(c+1)))-...
            min(cycstrain(1,1,pk(c):pk(c+1))),...
            max(cycstrain(2,2,pk(c):pk(c+1)))-...
            min(cycstrain(2,2,pk(c):pk(c+1))),...
            max(cycstrain(1,2,pk(c):pk(c+1)))-...
            min(cycstrain(1,2,pk(c):pk(c+1))),...
            max(E1)-min(E1),max(E2)-min(E2),...
            max(cycarea(pk(c):pk(c+1)))-min(cycarea(pk(c):pk(c+1)))};
    end
    cycstrain(:,:,1:pk(1)-1)=[];
    cycarea(1:pk(1)-1)=[];
    Exx=squeeze(cycstrain(1,1,:));
    Eyy=squeeze(cycstrain(2,2,:));
    Exy=squeeze(cycstrain(1,2,:));
    E1=((Exx+Eyy)/2)+sqrt(((Exx-Eyy)/2).^2+Exy.^2);
    E2=((Exx+Eyy)/2)-sqrt(((Exx-Eyy)/2).^2+Exy.^2);
    
    cycout.Frame=num2cell((pk(1):pk(Ncyc+1))');
    cycout.Exx=Exx;
    cycout.Eyy=Eyy;
    cycout.Exy=Exy;
    cycout.E1=E1;
    cycout.E2=E2;
    cycout.area=cycarea;
    cycout.resetframes=pk;
    cycout.peakframes=find(pkpoints);
    cycout.mode=lower(field);
    
    if nargout==0
        save([top_dir,'matlab_data\cyclicstrain.mat'],'dataset','cycout')
        % Compile output
        xcel=cell(size(cycstrain,3)+1,7);
        xcel(1,:)={'Frame','Exx','Eyy','Exy','E1','E2','Area';};
        xcel(2:end,:)=num2cell([(pk(1):pk(Ncyc+1))',...
            squeeze(cycstrain(1,1,:)),squeeze(cycstrain(2,2,:)),...
            squeeze(cycstrain(1,2,:)),E1,E2,cycarea]);
        xcelsum=cell(Ncyc+1,7);
        xcelsum(1,:)={'Cycle','Exx','Eyy','Exy','E1','E2','Area'};
        xcelsum(2:end,:)=strainsum;
        
        try
            xlswrite([top_dir,'greenstrain_',dataset,'.xlsx'],xcel,...
                'Cyclic strain')
            xlswrite([top_dir,'greenstrain_',dataset,'.xlsx'],xcelsum,...
                'Cycle Summary')
            y=0;
        catch
            warning('hdm2:wkbkreopen',['The specified workbook'...
                ,' was already open. The file will be closed then reopened'])
            
            wbkname = ['greenstrain_',dataset,'.xlsx'];
            h = actxGetRunningServer('Excel.Application');
            h.WorkBooks.Item(wbkname).Save;
            h.WorkBooks.Item(wbkname).Close;
            xlswrite([top_dir,'greenstrain_',dataset,'.xlsx'],xcel,...
                'Cyclic strain')
            xlswrite([top_dir,'greenstrain_',dataset,'.xlsx'],xcelsum,...
                'Cycle Summary')
            y=1;
        end
        if y==1
            winopen([top_dir,'greenstrain_',dataset,'.xlsx'])
        end
    else
        
        varargout{1}=cycout;
    end
    
    
    
end

function F=dispgrad(U,V,X,Y)
% DISPGRAD: Using the known initial position of a field of points as well
% as their displacements in the X and Y directions, this function returns
% the average deformation gradient tensor of the field.
% F=dispgrad(U,V,X,Y)
%
%   2014 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute
%
%^^^
    [dudx,dudy,~]=fitplane(X,Y,U);
    [dvdx,dvdy,~]=fitplane(X,Y,V);
    F=[dudx,dudy;dvdx,dvdy;]+eye(2);
end

function [varargout] = fitplane(x,y,z)
% FITPLANE: Returns the equation for a plane fit to x y and z data

Np=numel(z);
% Check input
if ~(Np==numel(y) && Np==numel(z))
    error('fitplane:badinput', 'All three inputs must have the same number of elements')
end
keep=~isnan(z);
xx = x(keep);
yy = y(keep);
zz = z(keep);
N=numel(zz);
c=[xx yy ones(N,1)];
coef = c\zz;
if nargout==1
    varargout{1}=coef;
else
    for i=3:-1:1
        varargout{i}=coef(i);
    end
end
end

function varargout = hdm2_avgstrain( dataset,xl )
% hdm2_avgstrain Calculates strains from compiled HDM2 data
%   Calculates average green and cauchy strains accross the HDM region of
%   interest. This data is compiled and output in two excel spreadsheets.
%
%   dataset: folder of dataset that has been compiled
%   xl: binary defining whether or not to output to excel file
%
%   2014 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%^^^

tic;

%% Import data
    if ischar(dataset)
    top_dir=[dataset,'\'];
        if exist([top_dir,'matlab_data\rawhdm2.mat'],'file')
            load([top_dir,'matlab_data\rawhdm2.mat']);
        elseif exist([top_dir,'matlab_data\rawhdm.mat'],'file')
            rawhdm.subshift=rawhdm.shift;
            
        else
            error('hdm or hdm2 data do not exist for this dataset')
        end
    else
        xl=0;
        rawhdm=dataset;
    end
%% Process inputs
    [y,x,time]=size(rawhdm.udata);
%     y=1:y;
%     x=1:x;

%% Preallocate variables
    incgrad.defgrad=zeros(2,2,time);
    incgrad.astretch=zeros(time,1);
    incgrad.astretchmap=zeros(y-1,x-1,time);
    incstrain.cauchy = zeros(2,2,time);
    incstrain.green = zeros(2,2,time);
    strain.cauchy = zeros(2,2,time);
    strain.green = zeros(2,2,time);
    strain.area=zeros(time,1);
    [X,Y]=meshgrid(0:rawhdm.subshift:x*rawhdm.subshift-rawhdm.subshift,...
        0:rawhdm.subshift:y*rawhdm.subshift-rawhdm.subshift);
    XX=(X-repmat(mean(X,2),[1,x]));
    YY=(Y-repmat(mean(Y,1),[y,1]));
        
%% Loop through data and calculate incremental strains

    for f = 1:time
        % Calculate slopes using matrix algebra
        U=rawhdm.udata(:,:,f);
        V=rawhdm.vdata(:,:,f);
        % Calculate displacement gradients and
        % Store the incremental deformation gradient
        incgrad.defgrad(:,:,f)=dispgrad(U,V,XX,YY);
        
        % Calculate area strains
        s1=double(rawhdm.subshift(1));
        s2=double(rawhdm.subshift(2));
        s=s1;
        U1=U(1:end-1,1:end-1);
        U2=U(1:end-1,2:end);
        U3=U(2:end,1:end-1);
        U4=U(2:end,2:end);
        V1=V(1:end-1,1:end-1);
        V2=V(1:end-1,2:end);
        V3=V(2:end,1:end-1);
        V4=V(2:end,2:end);
        
        adata2=.5*(2*s1*s2+s1*(-U1+U2-U3+U4)+s2*(-V1-V2+V3+V4)+U1.*V2-U2.*V1-U1.*V3+U3.*V1+U2.*V4-U4.*V2-U3.*V4+U4.*V3);
        incgrad.astretch(f,1)=((nanmean(adata2(:))-s^2)/(s^2))+1;
        incgrad.astretchmap(:,:,f)=((adata2-s^2)/(s^2))+1;
        F=incgrad.defgrad(:,:,1);
        for ff=2:f
            F=F*incgrad.defgrad(:,:,ff);
        end
        incstrain.green(:,:,f)=0.5*...
            (incgrad.defgrad(:,:,f)'*incgrad.defgrad(:,:,f)-eye(2));
        incstrain.cauchy(:,:,f)=0.5*...
            (incgrad.defgrad(:,:,f)+incgrad.defgrad(:,:,f)')-eye(2);
        strain.green(:,:,f)= 0.5*(F'*F-eye(2));
        strain.cauchy(:,:,f)= 0.5*(F+F')-eye(2);
        strain.area(f,1)=(prod(incgrad.astretch(1:f,1)))-1;
    end

%% Calculate frequencies from incremental strains
    [powerspec.E11,powerspec.freq]=...
        pfft(squeeze(incstrain.green(1,1,:)),rawhdm.framerate);
    [powerspec.E22,~]=...
        pfft(squeeze(incstrain.green(2,2,:)),rawhdm.framerate);
    [powerspec.E12,~]=...
        pfft(squeeze(incstrain.green(1,2,:)),rawhdm.framerate);
    [powerspec.e11,~]=...
        pfft(squeeze(incstrain.cauchy(1,1,:)),rawhdm.framerate);
    [powerspec.e22,~]=...
        pfft(squeeze(incstrain.cauchy(2,2,:)),rawhdm.framerate);
    [powerspec.e12,~]=...
        pfft(squeeze(incstrain.cauchy(1,2,:)),rawhdm.framerate);
    [powerspec.area,~]=...
        pfft(squeeze(incgrad.astretch),rawhdm.framerate);
    timeseries=(0:(1/rawhdm.framerate):((time-1)/rawhdm.framerate))';
    
    Exx=squeeze(strain.green(1,1,:));
    Eyy=squeeze(strain.green(2,2,:));
    Exy=squeeze(strain.green(1,2,:));
    E1=((Exx+Eyy)/2)+sqrt(((Exx-Eyy)/2).^2+Exy.^2);
    E2=((Exx+Eyy)/2)-sqrt(((Exx-Eyy)/2).^2+Exy.^2);
    if exist('atan2','builtin')
        Q=0.5*atan2(2*(Exy),Exx-Eyy);
        Q=Q*180/pi();
        %Q(Q<0)=Q(Q<0)+360;
    else
        Q=zeros(size(Exx));
        warning('Function atan2 does not exist')
    end
%% Prepare data for Excel output
    % For Green's Strain
    disp('>> Creating output matrices');
        % Create a frame matrix
        frame=cell(time,1);
        for i=1:time
            frame{i,1}=['''(',num2str(i-1),',',num2str(i),')'];
        end

    % Create output matrices
        xcel1=cell(time+1,9);
        xcel2=cell(length(powerspec.freq)+1,5);
        xcel3=cell(time+1,6);

        % Create headers
            xcel1(1,:)={'Frame','Time(s)','E11','E22','E12','E1','E2','Angle','area'};
            xcel2(1,:)={'Frequency(Hz)','11 Power','22 Power','12 Power','area power'};
            xcel3(1,:)={'Frame','Time(s)','inc_E11','inc_E22','inc_E12','inc_area'};
        % Place data in output matrix
            xcel1(2:(time+1),:)=horzcat(frame,...
                cellstr(num2str(timeseries))...
                ,cellstr(num2str(Exx))...
                ,cellstr(num2str(Eyy))...
                ,cellstr(num2str(Exy))...
                ,cellstr(num2str(E1))...
                ,cellstr(num2str(E2))...
                ,cellstr(num2str(Q))...
                ,cellstr(num2str(squeeze(strain.area))));
            xcel2(2:(length(powerspec.freq)+1),:)=...
                horzcat(cellstr(num2str(powerspec.freq))...
                ,cellstr(num2str(powerspec.E11))...
                ,cellstr(num2str(powerspec.E22))...
                ,cellstr(num2str(powerspec.E12))...
                ,cellstr(num2str(powerspec.area)));
            xcel3(2:(time+1),:)=horzcat(frame,...
                cellstr(num2str(timeseries))...
                ,cellstr(num2str(squeeze(incstrain.green(1,1,:))))...
                ,cellstr(num2str(squeeze(incstrain.green(2,2,:))))...
                ,cellstr(num2str(squeeze(incstrain.green(1,2,:))))...
                ,cellstr(num2str(squeeze(incgrad.astretch))));
%     % For Cauchy strain
% 
%     % Create output matrices
%         xcel4=cell(time+1,5);
%         xcel5=cell(length(powerspec.freq)+1,3);
%         xcel6=cell(time+1,6);
% 
%         % Create headers
%             xcel4(1,:)={'Frame','Time(s)','E11','E22','E12'};
%             xcel5(1,:)={'Frequency(Hz)','11 Power','22 Power','12 Power'};
%             xcel6(1,:)={'Frame','Time(s)','inc_E11','inc_E22','inc_E12'};
%         % Place data in output matrix
%             xcel4(2:(time+1),:)=horzcat(frame,...
%                 cellstr(num2str(timeseries))...
%                 ,cellstr(num2str(squeeze(strain.cauchy(1,1,:))))...
%                 ,cellstr(num2str(squeeze(strain.cauchy(2,2,:))))...
%                 ,cellstr(num2str(squeeze(strain.cauchy(1,2,:)))));
%             xcel5(2:(length(powerspec.freq)+1),:)=...
%                 horzcat(cellstr(num2str(powerspec.freq))...
%                 ,cellstr(num2str(powerspec.e11))...
%                 ,cellstr(num2str(powerspec.e22)));
%             xcel6(2:(time+1),:)=horzcat(frame,...
%                 cellstr(num2str(timeseries))...
%                 ,cellstr(num2str(squeeze(incstrain.cauchy(1,1,:))))...
%                 ,cellstr(num2str(squeeze(incstrain.cauchy(2,2,:))))...
%                 ,cellstr(num2str(squeeze(incstrain.cauchy(1,2,:)))));
            
%% Output data to excel files
    
    if xl==1
        warning('off','MATLAB:xlswrite:AddSheet')
        ds = strrep(dataset,'\','-');
        for i = 1:3
            try
                % Greens strain sheet
                disp('>> Attempting to write Excel Data file');
                xlswrite([top_dir,'\greenstrain_',ds,'.xlsx']...
                    ,xcel1,'Strain Data');
                xlswrite([top_dir,'\greenstrain_',ds,'.xlsx']...
                    ,xcel2,'Frequency Data');
                xlswrite([top_dir,'\greenstrain_',ds,'.xlsx']...
                    ,xcel3,'Incremental Strains');
    %             % Cauchy strain sheet
    %             xlswrite([top_dir,'\cauchystrain_',ds,'.xlsx']...
    %                 ,xcel4,'Strain Data');
    %             xlswrite([top_dir,'\cauchystrain_',ds,'.xlsx']...
    %                 ,xcel5,'Frequency Data');
    %             xlswrite([top_dir,'\cauchystrain_',ds,'.xlsx']...
    %                 ,xcel6,'Incremental Strains');
                disp('>> Excel data successfully written');
                break;
            catch err
                disp(['>> Waiting for user to close Excel File (attempt #',num2str(i),')'])
                disp('>> Press any key to try again');
                err.message,err.stack
                pause;
            end
        end
        warning('on','MATLAB:xlswrite:AddSheet')
    end
    tfinal=toc;

    if nargout==1
        out.time=time;
        out.strain=strain;
        out.powerspec=powerspec;
        out.incstrain=incstrain;
        out.incgrad=incgrad;
        out.tfinal=tfinal;
        varargout{1}=out;        
    elseif nargout == 5
        varargout= {time,strain,powerspec,incstrain,incgrad};
    else
        %% Save the data to a matlab file and clear the memory
        clearvars -except time strain top_dir Pf1 Pf2 freq incstrain...
            tfinal incgrad adata2
        save(strcat(top_dir,'\matlab_data\avg_strain_tensor.mat'))

        disp(['Average Strain calculations completed in: ',...
            datestr(datenum(0,0,0,0,0,tfinal),'HH:MM:SS')])
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



