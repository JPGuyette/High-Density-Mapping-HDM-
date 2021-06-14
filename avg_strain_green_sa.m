function varargout = avg_strain_green_sa( dataset,xl )
% Compiled by HORUS\John T Favreau on 03/13/2014 08:49:55
% avg_strain_green Version 2014.02.0008
% 
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
        adata=0.5*(U3.*V4 - (U4+s).*V3 + (U4+s).*(V2+s)...
                    -(U2+s).*V4 + (U2+s).*(V1+s) - U1.*(V2+s) ...
                    + U1.*V3 - U3.*(V1+s));
        adata2=.5*(2*s1*s2+s1*(-U1+U2-U3+U4)+s2*(-V1-V2+V3+V4)+U1.*V2-U2.*V1-U1.*V3+U3.*V1+U2.*V4-U4.*V2-U3.*V4+U4.*V3);
        incgrad.astretch(f,1)=((nanmean(adata2(:))-s^2)/s^2)+1;
        incgrad.astretchmap(:,:,f)=((adata2-s^2)/s^2)+1;
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
        Q(Q<0)=Q(Q<0)+360;
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



