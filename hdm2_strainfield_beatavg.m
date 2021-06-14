function hdm2_strainfield_beatavg_sa(dataset,res,frameranges,fieldrange,mapname)
% Compiled by BIOMED-NB03\jfavreau on 02/26/2015 16:07:00
% hdm2_strainfield_beatavg Version 2015.02.0002
% 
    top_dir=[dataset,'\'];
    load([top_dir,'matlab_data\rawhdm2.mat'])
    if ischar(frameranges)
        
            switch frameranges
                case 'expand'
                    strain=hdm2_cyclicstrain(dataset,'E1','maxfreqmins');
                case 'contract'
                    strain=hdm2_cyclicstrain(dataset,'E2','maxfreqmaxs');
            end
        
        frameranges=cell(length(strain.peakframes),1);
        for k=1:length(frameranges)
            frameranges{k}=strain.resetframes(k):strain.peakframes(k);
        end
    end
    J=length(frameranges);
    for j=1:J
        SF=hdm2_strainfield(dataset,res,frameranges{j},fieldrange,mapname);
        if j==1
            straindata=SF;
            im1=frameranges{j};
            im1=im1(1);
        else 
            straindata.E1=SF.E1+straindata.E1;
            straindata.E2=SF.E2+straindata.E2;
            straindata.Exx=SF.Exx+straindata.Exx;
            straindata.Eyy=SF.Eyy+straindata.Eyy;
            straindata.Exy=SF.Exy+straindata.Exy;
            straindata.area=SF.area+straindata.area;
        end
    end
    straindata.E1=straindata.E1/J;
    straindata.E2=straindata.E2/J;
    straindata.Exx=straindata.Exx/J;
    straindata.Eyy=straindata.Eyy/J;
    straindata.Exy=straindata.Exy/J;
    straindata.area=straindata.area/J; 
    outdir=[dataset,'\heatmaps\beatavg\'];
    
    
    
    
    mkdir(outdir);
    save([outdir,'strainfields.mat'],'dataset','straindata')
     % Create overlay images
    imlist=load([dataset,'\files\fileloc_tif.mat']);
    curim=imread(imlist.fileloc.paths{im1});
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.E1,...
        0.5,rawhdm.roilist(:,:,im1),fieldrange,mapname);
    imwrite(outim,[outdir,'E1overlay.tif'])
    imwrite(pix,[outdir,'0_E1strain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.E2,...
        0.5,rawhdm.roilist(:,:,im1),fieldrange,mapname);
    imwrite(outim,[outdir,'E2overlay.tif'])
    imwrite(pix,[outdir,'0_E2strain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.Exx,...
        0.5,rawhdm.roilist(:,:,im1),fieldrange,mapname);
    imwrite(outim,[outdir,'Exxoverlay.tif'])
    imwrite(pix,[outdir,'0_Exxstrain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.Eyy,...
        0.5,rawhdm.roilist(:,:,im1),fieldrange,mapname);
    imwrite(outim,[outdir,'Eyyoverlay.tif'])
    imwrite(pix,[outdir,'0_Eyystrain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.Exy,...
        0.5,rawhdm.roilist(:,:,im1),fieldrange,mapname);
    imwrite(outim,[outdir,'Exyoverlay.tif'])
    imwrite(pix,[outdir,'0_Exystrain.tif']);
%     [~,~,outim]=transover(curim,straindata.xxyy,straindata.angle*(180/pi()),...
%         0.5,rawhdm.roilist(:,:,im1),[0 360],mapname);
%     imwrite(outim,[outdir,'Angle.tif'])
%     [~,~,outim]=transover(curim,straindata.xxyy,udisp,...
%         0.5,rawhdm.roilist(:,:,im1),[0 50],mapname);
%     imwrite(outim,[outdir,'Udispoverlay.tif'])
%     [~,~,outim]=transover(curim,straindata.xxyy,vdisp,...
%         0.5,rawhdm.roilist(:,:,im1),[0 50],mapname);
%     imwrite(outim,[outdir,'Vdispoverlay.tif'])
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.area,...
        0.5,rawhdm.roilist(:,:,im1),fieldrange,mapname);
    imwrite(outim,[outdir,'areaoverlay.tif'])
    imwrite(pix,[outdir,'0_areastrain.tif']);
%     [~,~,outim]=transover(curim,straindata.xxyy,straindata.angle*(180/pi()),...
%         0.5,rawhdm.roilist(:,:,f),[0 360],mapname);
%     imwrite(outim,[outdir,'angleoverlay.tif'])
        
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
            rawhdm.subshift=rawhdm.shift; %#ok<NODEF>
            
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

function varargout = hdm2_cyclicstrain(dataset,field,method,varargin)

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

function [varargout]=hdm2_strainfield(dataset,res,framerange,fieldrange,mapname)
% HDM2_strainfield Calculates strains from a displacement field
% Calculates strains on a displacement field with the given shift and  
% resolution parameters. Returns a structured array with Exx, Eyy, Exy,
% E1, E2, and rotation angle.
%   
%   dataset: dataset that has been analyzed with hdm2
%   res: number of points to calculate slopes across for strain
%       calculations'MarkerSize',10)
%   framerange: array of frames over which to sum displacements
%
%   2013 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%--------------------------------------------------------------------------
if isempty(mapname)
    mapname='jet';
end
    if ischar(dataset)
        top_dir=[dataset,'\'];
        load([top_dir,'matlab_data\rawhdm2.mat']);
    else
        rawhdm=dataset;
        dataset=rawhdm.dataset;
        top_dir=[dataset,'\'];
    end
    % Load data
    
    if ~exist([top_dir,'heatmaps'],'dir')
        mkdir([top_dir,'heatmaps'])
    end
    outdir=[top_dir,'heatmaps\strain_',...
            num2str(framerange(1)),'_',num2str(framerange(end)),'\'];
    if~exist(outdir,'dir')
        mkdir(outdir)
    end
    
    udisp=sum(rawhdm.udata(:,:,framerange),3);
    vdisp=sum(rawhdm.vdata(:,:,framerange),3);
    shift=rawhdm.subshift;    
    buf=(res-1)/2;
    [y,x]=size(udisp);
    xvector=0:shift:shift*(res-1);
    yvector=0:shift:shift*(res-1);
    [XX,YY]=meshgrid(xvector,yvector);
    % Calculate strains
    straintensor=nan(2,2,y,x);
    ang=zeros(y,x);
    s=rawhdm.subshift;
    tic
%     tstart=now();
%     bar=waitbar(0,'Processing displacement data...');
    % Calculate using linear regression formulas
    for i=buf+1:y-buf
        for j=buf+1:x-buf
            if ~any(isnan([udisp(i-buf:i+buf,j-buf:j+buf),...
                    vdisp(i-buf:i+buf,j-buf:j+buf)]))
            
                F=dispgrad(udisp(i-buf:i+buf,j-buf:j+buf),...
                    vdisp(i-buf:i+buf,j-buf:j+buf),XX,YY);
                % Green's Strain tensor
                straintensor(:,:,i,j)=0.5*(F'*F-eye(2));
                % Deformation gradient
                if any(any(isnan(F)))
                    ang(i,j)=NaN;
                else
                    R=F/(sqrtm(F'*F));
                    ang(i,j)=acos(R(1,1));
                end
            end
%             trem=(now()-tstart)*(N-k+1)/(k+1);
%             waitbar((k+1)/N,bar,['Time remaining: ',datestr(trem,'HH:MM:SS')]);
        end
    end
    %% Calculate area strain field
    if ~exist([top_dir,'matlab_data\avg_strain_tensor.mat'],'file')
        hdm2_avgstrain(dataset,0);
    end
    load([top_dir,'matlab_data\avg_strain_tensor.mat'])
    h=fspecial('gaussian',res,0.5);
    temparea=prod(incgrad.astretchmap(:,:,framerange),3)-1;
    straindata.area=filter2(h,temparea);
    f=framerange(1);
    %s=rawhdm.subshift;
    straindata.tensor=straintensor;
    straintensor=permute(straintensor,[3 4 1 2]);
    straindata.Exx=straintensor(:,:,1,1);
    straindata.Eyy=straintensor(:,:,2,2);
    straindata.Exy=straintensor(:,:,1,2);
    straindata.angle=(180/pi())*(1/2)*...
    atan2(2*straindata.Exy,(straindata.Exx-straindata.Eyy));
    straindata.E1=((straindata.Exx+straindata.Eyy)/2)+...
        sqrt(((straindata.Exx-straindata.Eyy)/2).^2+straindata.Exy.^2);
    straindata.E2=((straindata.Exx+straindata.Eyy)/2)-...
        sqrt(((straindata.Exx-straindata.Eyy)/2).^2+straindata.Exy.^2);
    % Strip edge effects from data
%     straindata.E1=straindata.E1(buf+1:end-buf,buf+1:end-buf);
%     straindata.E2=straindata.E2(buf+1:end-buf,buf+1:end-buf);
%     straindata.angle=straindata.angle(buf+1:end-buf,buf+1:end-buf);
%     straindata.Exx=straindata.Exx(buf+1:end-buf,buf+1:end-buf);
%     straindata.Exy=straindata.Exy(buf+1:end-buf,buf+1:end-buf);
%     straindata.Eyy=straindata.Eyy(buf+1:end-buf,buf+1:end-buf);
     straindata.xxyy=rawhdm.xxyy(f,:)+[buf*s(2),-buf*s(2),buf*s(1),-buf*s(1)];
    % Create overlay images
    imlist=load([dataset,'\files\fileloc_tif.mat']);
    
    curim=imread(imlist.fileloc.paths{f});
    if size(curim,3)~=1
        curim=curim(:,:,1);
    end
    
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.E1,...
        0.5,rawhdm.roilist(:,:,f),fieldrange,mapname);
    imwrite(outim,[outdir,'E1overlay.tif'])
    imwrite(pix,[outdir,'0_E1strain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.E2,...
        0.5,rawhdm.roilist(:,:,f),fieldrange,mapname);
    imwrite(outim,[outdir,'E2overlay.tif'])
    imwrite(pix,[outdir,'0_E2strain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.Exx,...
        0.5,rawhdm.roilist(:,:,f),fieldrange,mapname);
    imwrite(outim,[outdir,'Exxoverlay.tif'])
    imwrite(pix,[outdir,'0_Exxstrain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.Eyy,...
        0.5,rawhdm.roilist(:,:,f),fieldrange,mapname);
    imwrite(outim,[outdir,'Eyyoverlay.tif'])
    imwrite(pix,[outdir,'0_Eyystrain.tif']);
    [~,pix,outim]=transover(curim,straindata.xxyy,straindata.Exy,...
        0.5,rawhdm.roilist(:,:,f),fieldrange,mapname);
    imwrite(outim,[outdir,'Exyoverlay.tif'])
    imwrite(pix,[outdir,'0_Exystrain.tif']);
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.angle*(180/pi()),...
        0.5,rawhdm.roilist(:,:,f),[0 360],mapname);
    imwrite(outim,[outdir,'Angle.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,udisp,...
        0.5,rawhdm.roilist(:,:,f),[0 50],mapname);
    imwrite(outim,[outdir,'Udispoverlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,vdisp,...
        0.5,rawhdm.roilist(:,:,f),[0 50],mapname);
    imwrite(outim,[outdir,'Vdispoverlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.area,...
        0.5,rawhdm.roilist(:,:,f),fieldrange,mapname);
    imwrite(outim,[outdir,'areaoverlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.angle,...
        0.5,rawhdm.roilist(:,:,f),[-180 180],mapname);
    imwrite(outim,[outdir,'angleoverlay.tif'])

%     % Output SAC surface
%     [SAC1,SAC2,~]=size(beatdata.average.astretchmap);
%     beats(beatgroup(b)).strain.SAC = (ones(SAC1,SAC2)-...
%         beatdata.average.astretchmap(:,:,beatdata.average.ES))*100;
%     beats(beatgroup(b)).strain.SW=trapz(beatdata.average.pressure,beatdata.average.astretchmap,3);
%     filt2=fspecial(filttype,[rawhdm.res rawhdm.res]);
%     beats(beatgroup(b)).strain.SAC_s=filter2(filt2,padarray(beats(beatgroup(b)).strain.SAC,[buf buf],'replicate'),'valid');
%     beats(beatgroup(b)).strain.SW_s=filter2(filt2,padarray(beats(beatgroup(b)).strain.SW,[buf buf],'replicate'),'valid');
%     beats(beatgroup(b)).strain.Exx_s=filter2(filt2,padarray(beats(beatgroup(b)).strain.Exx,[buf buf],'replicate'),'valid');
%     beats(beatgroup(b)).strain.Eyy_s=filter2(filt2,padarray(beats(beatgroup(b)).strain.Eyy,[buf buf],'replicate'),'valid');
%     beats(beatgroup(b)).strain.angle=ang;
%     beats(beatgroup(b)).strain.angle2=(1/2)*atan2((2*beats(beatgroup(b)).strain.Exy),(beats(beatgroup(b)).strain.Exx-beats(beatgroup(b)).strain.Eyy));
    straindata.res=res;
    switch nargout
        case  1
            varargout{1}=straindata;
            save([outdir,'strainfields.mat'],'straindata')
        case 0
            save([outdir,'strainfields.mat'],'straindata')
    end
    
    close all;
    fclose('all');
    toc;
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
end

function [raw,pix,im]=transover(ol_image,ol_roi,overlay,trans,~,range,varargin)
% transover Creates overlays a matrix on an image given the location of the
% matric
%
%   ol_image: Matlab variable containing the image to overlay
%   ol_roi: xxyy position of the region of interest on the image (pixels)
%   overlay: Matrix to overlay in the ROI
%   trans: Opacity of the overlay
%   pos: NO LONGER USED
%   range: The strain range to show
%   colormap: (optional, default is now parula)
%
%   2013 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%^^^
if nargin==7
    mapname=varargin{1};
else
    mapname='parula';
end
% Create and resize a pixel-sized interpolation of the matrix image
        raw = overlay;
        pixyx=[1+ol_roi(4)-ol_roi(3),1+ol_roi(2)-ol_roi(1)];
        matyx=size(raw);
        maxratio=max(pixyx/matyx);
        n=nextpow2(maxratio);
        pix = interp2(raw,n);
        pix=imresize(pix,pixyx);
        kill=isnan(pix);
        
    % Create an overlay for just the roi
%         ol_mask=poly2mask(pos(:,1)-ol_roi(1)+1,pos(:,2)-ol_roi(3)+1,...
%             1+ol_roi(4)-ol_roi(3),1+ol_roi(2)-ol_roi(1));
%         
%         pix=pix.*ol_mask;
        im=pix;
    
    % Calculate a pixel-sized rgb image for the overlay
        
        if ~exist('range','var');
            minval = min(min(im));
            maxval=max(max(im));
        else
            minval=range(1);
            maxval=range(2);
        end
        im=im-minval;
        im=im/(maxval-minval);
        im=im*255;
        im=im+1;
        im=round(im);
        im=ind2rgb(im,jtfheatmap(mapname,256));
        pix=im;
        pix(repmat(kill,[1 1 3]))=NaN;
        close gcf

    % Construct strain overlay map and place in the ROI
        
        ol_image=doubleim(ol_image);
        ol_image=repmat(ol_image,[1 1 3]);
        im_roi=ol_image(ol_roi(3):ol_roi(4),ol_roi(1):ol_roi(2),:);
        c_a=trans*im;
        C=c_a+(1-trans)*im_roi;
        ignore=repmat(kill,[1 1 3]);
        C(ignore)=im_roi(ignore);
        [imheight,imwidth,~]=size(ol_image);
%         fullmask=poly2mask(pos(:,1),pos(:,2),imheight,imwidth);
% %         mask=repmat(fullmask(ol_roi(3):ol_roi(4),ol_roi(1):ol_roi(2),:),[1 1 3]);
% %         C(~mask)=im_roi(~mask);
        ol_image(ol_roi(3):ol_roi(4),ol_roi(1):ol_roi(2),:)=C;
        cbar=getcolorbar(range,size(ol_image,1),mapname);
        im=zeros(imheight,imwidth+size(cbar,2),3);
        im(:,1:imwidth,:)=ol_image;
        im(:,imwidth+1:end,:)=cbar;
        cbar2=getcolorbar(range,size(pix,1),mapname);
        pix=[pix,cbar2];
end

function scale=getcbar(min,max)

end

