function varargout = hdm_greenstrain_cyclic_standalone(dataset,field,method,varargin)
% Compiled by HORUS\John T Favreau on 02/24/2014 14:52:03
% hdm_greenstrain_cyclic Version 2014.02.0006
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
            rawhdm.subshift=rawhdm.shift;
        else
            error('hdm or hdm2 data do not exist for this dataset')
        end
    else
        error('hdm_greenstrain_cyclic is not configured to analyze this dataset')
    end
    if ~exist([top_dir,'matlab_data\avg_strain_tensor.mat'],'file')
        avg_strain_green(dataset,0)
    end
    load([top_dir,'matlab_data\avg_strain_tensor.mat'])
    
%% Process inputs
    [y,x,time]=size(rawhdm.udata);
%     y=1:y;
%     x=1:x;
    switch lower(field)
        case {'eyy','ey','e22','y'}
            fielddata=squeeze(strain.green(2,2,:));
        case {'exx','ex','e11','x'}
            fielddata=squeeze(strain.green(1,1,:));
        case {'exy','eyx','e12','e21','xy'}
            fielddata=squeeze(strain.green(2,2,:));
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
            [~,zeropoints]=peakfind1(fielddata,rawhdm.framerate,freqpass);
        case 'maxfreqmaxs'
            [zeropoints,~]=peakfind1(fielddata,rawhdm.framerate,freqpass);
        otherwise
            error('This input has not been implemented')
    end

%% Preallocate variables
   
    Ncyc=numel(find(zeropoints))-1;
    strainsum=cell(Ncyc,4);
    pk=find(zeropoints);
    cycstrain=zeros(2,2,pk(Ncyc+1)-pk(1)+1);
    % loop through data
    for c= 1:Ncyc
        
        for f=pk(c):pk(c+1)
            F=[1 0;0 1;];
            for ff=(pk(c)+1):f
                F=F*incgrad.defgrad(:,:,ff);
            end
            cycstrain(:,:,f)= 0.5*(F'*F-eye(2));
        end
        strainsum(c,:)={c,max(cycstrain(1,1,pk(c):pk(c+1)))-...
            min(cycstrain(1,1,pk(c):pk(c+1))),...
            max(cycstrain(2,2,pk(c):pk(c+1)))-...
            min(cycstrain(2,2,pk(c):pk(c+1))),...
            max(cycstrain(1,2,pk(c):pk(c+1)))-...
            min(cycstrain(1,2,pk(c):pk(c+1)))};
    end
    cycstrain(:,:,1:pk(1)-1)=[];
    % Compile output
    xcel=cell(size(cycstrain,3)+1,4);
    xcel(1,:)={'Frame','Exx','Eyy','Exy';};
    xcel(2:end,:)=num2cell([(pk(1):pk(Ncyc+1))',...
        squeeze(cycstrain(1,1,:)),squeeze(cycstrain(2,2,:)),...
        squeeze(cycstrain(1,2,:))]);
    xcelsum=cell(Ncyc+1,4);
    xcelsum(1,:)={'Cycle','Exx','Eyy','Exy'};
    xcelsum(2:end,:)=strainsum;
    xlswrite([top_dir,'greenstrain_',dataset,'.xlsx'],xcel,...
        'Cyclic strain')
    xlswrite([top_dir,'greenstrain_',dataset,'.xlsx'],xcelsum,...
        'Cycle Summary')
 
    
    
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
    x=squeeze(x);
    y=squeeze(y);
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



function varargout = avg_strain_green( dataset,xl )
% avg_strain_green2 Calculates strains from compiled HDM2 data
%   Calculates average green and cauchy strains accross the HDM region of
%   interest. This data is compiled and output in two excel spreadsheets.
%
%   dataset: folder of dataset that has been compiled
%   xl: binary defining whether or not to output to excel file
%
%   2013 John Favreau
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
        UX=(U-repmat(nanmean(U,2),[1,x]));
        UY=(U-repmat(nanmean(U,1),[y,1]));
        VX=(V-repmat(nanmean(V,2),[1,x]));
        VY=(V-repmat(nanmean(V,1),[y,1]));
        % Calculate displacement gradients
        dudx=nanmean(nansum(UX.*XX)/sum(XX.^2));
        dudy=nanmean(nansum(UY.*YY)/sum(YY.^2));
        dvdx=nanmean(nansum(VX.*XX)/sum(XX.^2));
        dvdy=nanmean(nansum(VY.*YY)/sum(YY.^2));
        % Store the incremental deformation gradient
        incgrad.defgrad(:,:,f)=[dudx,dudy;dvdx,dvdy;]+eye(2);
        
        % Calculate area strains
        s=double(rawhdm.subshift(1));
        U1=U(1:end-1,1:end-1);
        U2=U(1:end-1,2:end);
        U3=U(2:end,1:end-1);
        U4=U(2:end,2:end);
        V1=V(1:end-1,1:end-1);
        V2=V(1:end-1,2:end);
        V3=V(2:end,1:end-1);
        V4=V(2:end,2:end);
        adata=0.5*(U3.*V4 - (U4+s).*V3 + (U4+s).*(V2+s)...
                    -(U2+s).*V4 + (U2+s).*(V1+s) - U1.*(V2+s) ...
                    + U1.*V3 - U3.*(V1+s));
        incgrad.astretch(f,1)=((nanmean(adata(:))-s^2)/s^2)+1;
        incgrad.astretchmap(:,:,f)=((adata-s^2)/s^2)+1;
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
%% Prepare data for Excel output
    % For Green's Strain
    disp('>> Creating output matrices');
        % Create a frame matrix
        frame=cell(time,1);
        for i=1:time
            frame{i,1}=['''(',num2str(i-1),',',num2str(i),')'];
        end

    % Create output matrices
        xcel1=cell(time+1,6);
        xcel2=cell(length(powerspec.freq)+1,5);
        xcel3=cell(time+1,6);

        % Create headers
            xcel1(1,:)={'Frame','Time(s)','E11','E22','E12','area'};
            xcel2(1,:)={'Frequency(Hz)','11 Power','22 Power','12 Power','area power'};
            xcel3(1,:)={'Frame','Time(s)','inc_E11','inc_E22','inc_E12','inc_area'};
        % Place data in output matrix
            xcel1(2:(time+1),:)=horzcat(frame,...
                cellstr(num2str(timeseries))...
                ,cellstr(num2str(squeeze(strain.green(1,1,:))))...
                ,cellstr(num2str(squeeze(strain.green(2,2,:))))...
                ,cellstr(num2str(squeeze(strain.green(1,2,:))))...
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
    else
        %% Save the data to a matlab file and clear the memory
        clearvars -except time strain top_dir Pf1 Pf2 freq incstrain...
            tfinal incgrad adata
        save(strcat(top_dir,'\matlab_data\avg_strain_tensor.mat'))

        disp(['Average Strain calculations completed in: ',...
            datestr(datenum(0,0,0,0,0,tfinal),'HH:MM:SS')])
    end
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

