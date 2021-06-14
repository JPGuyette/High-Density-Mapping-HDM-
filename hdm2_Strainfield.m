function [varargout]=hdm2_strainfield(dataset,res,framerange,fieldrange)
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
    % Load data
    top_dir=[dataset,'\'];
    if ~exist([top_dir,'heatmaps'],'dir')
        mkdir([top_dir,'heatmaps'])
    end
    outdir=[top_dir,'heatmaps\strain_',...
            num2str(framerange(1)),'_',num2str(framerange(end)),'\'];
    if~exist(outdir,'dir')
        mkdir(outdir)
    end
    load([top_dir,'matlab_data\rawhdm2.mat'])
    udisp=sum(rawhdm.udata(:,:,framerange),3);
    vdisp=sum(rawhdm.vdata(:,:,framerange),3);
    shift=rawhdm.subshift;    
    buf=(res-1)/2;
    [y,x]=size(udisp);
    xvector=0:shift:shift*(res-1);
    yvector=0:-shift:-shift*(res-1);
    % Calculate strains
    straintensor=zeros(2,2,y,x);
    ang=zeros(y,x);
    tic
%     tstart=now();
%     bar=waitbar(0,'Processing displacement data...');
    % Calculate using linear regression formulas
    for i=1:y
        for j=1:x
            x2=sum((xvector-mean(xvector)).^2);
            if j<1+buf
                % Calculate dudx at each point
                    u=udisp(i,1:j+buf);
                    ux2=sum((u-mean(u)).*(xvector(1:length(u))-mean(xvector(1:length(u)))));
                    dudx1=ux2/x2;
                % Calculate dvdx at each point
                    u=vdisp(i,1:j+buf);
                    ux2=sum((u-mean(u)).*(xvector(1:length(u))-mean(xvector(1:length(u)))));
                    dvdx1=ux2/x2;
            elseif j>x-buf
                % Calculate dudx at each point
                    u=udisp(i,j-buf:end);
                    ux2=sum((u-mean(u)).*(xvector(1:length(u))-mean(xvector(1:length(u)))));
                    dudx1=ux2/x2;
                % Calculate dvdx at each point
                    u=vdisp(i,j-buf:end);
                    ux2=sum((u-mean(u)).*(xvector(1:length(u))-mean(xvector(1:length(u)))));
                    dvdx1=ux2/x2;
            else
                % Calculate dudx at each point
                    u=udisp(i,j-buf:j+buf);
                    ux2=sum((u-mean(u)).*(xvector-mean(xvector)));
                    dudx1=ux2/x2;
                % Calculate dvdx at each point
                    u=vdisp(i,j-buf:j+buf);
                    ux2=sum((u-mean(u)).*(xvector-mean(xvector)));
                    dvdx1=ux2/x2;
            end
            if i<1+buf%--------------------------xvector below changed to yvector!!!!-Strain issue!
                % Calculate dudy at each point
                    u=udisp(1:i+buf,j);
                    ux2=sum((u-mean(u))'.*(yvector(1:length(u))-mean(yvector(1:length(u)))));
                    dudy1=ux2/x2;
                % Calculate dvdx at each point
                    u=vdisp(1:i+buf,j);
                    ux2=sum((u-mean(u))'.*(yvector(1:length(u))-mean(yvector(1:length(u)))));
                    dvdy1=ux2/x2;
            elseif i>y-buf
                % Calculate dudy at each point
                    u=udisp(i-buf:end,j);
                    ux2=sum((u-mean(u))'.*(yvector(1:length(u))-mean(yvector(1:length(u)))));
                    dudy1=ux2/x2;
                % Calculate dvdx at each point
                    u=vdisp(i-buf:end,j);
                    ux2=sum((u-mean(u))'.*(yvector(1:length(u))-mean(yvector(1:length(u)))));
                    dvdy1=ux2/x2;
            else
                % Calculate dudy at each point
                    u=udisp(i-buf:i+buf,j);
                    ux2=sum((u-mean(u))'.*(yvector-mean(yvector)));
                    dudy1=ux2/x2;
                % Calculate dvdx at each point
                    u=vdisp(i-buf:i+buf,j);
                    ux2=sum((u-mean(u))'.*(yvector-mean(yvector)));
                    dvdy1=ux2/x2;
            end
            F=[dudx1,dudy1;dvdx1,dvdy1;]+eye(2);
            % Green's Strain tensor
            straintensor(:,:,i,j)=0.5*(F'*F-eye(2));
            % Deformation gradient
            if any(any(isnan(F)))
                ang(i,j)=NaN;
            else
                R=F/(sqrtm(F'*F));
                ang(i,j)=acos(R(1,1));
            end
%             trem=(now()-tstart)*(N-k+1)/(k+1);
%             waitbar((k+1)/N,bar,['Time remaining: ',datestr(trem,'HH:MM:SS')]);
        end
    end
    f=framerange(1);
    s=rawhdm.subshift;
    straindata.tensor=straintensor;
    straintensor=permute(straintensor,[3 4 1 2]);
    straindata.Exx=straintensor(:,:,1,1);
    straindata.Eyy=straintensor(:,:,2,2);
    straindata.Exy=straintensor(:,:,1,2);
    straindata.angle=ang;
    straindata.E1=straindata.Exx.*... 
            cos(straindata.angle).^2+...
            straindata.Eyy.*...
            sin(straindata.angle).^2+...
            2*straindata.Exy.*...
            sin(straindata.angle).*...
            cos(straindata.angle); 
    straindata.E2=straindata.Exx.*...
            sin(straindata.angle).^2+...
            straindata.Eyy.*...
            cos(straindata.angle).^2-...
            2*straindata.Exy.*...
            sin(straindata.angle).*...
            cos(straindata.angle);
    % Strip edge effects from data
    straindata.E1=straindata.E1(buf:end-buf,buf:end-buf);
    straindata.E2=straindata.E2(buf:end-buf,buf:end-buf);
    straindata.angle=straindata.angle(buf:end-buf,buf:end-buf);
    straindata.Exx=straindata.Exx(buf:end-buf,buf:end-buf);
    straindata.Exy=straindata.Exy(buf:end-buf,buf:end-buf);
    straindata.Eyy=straindata.Eyy(buf:end-buf,buf:end-buf);
    straindata.xxyy=rawhdm.xxyy(f,:)+[buf*s(2),-buf*s(2),buf*s(1),-buf*s(1)];
    % Create overlay images
    imlist=load([dataset,'\files\fileloc_tif.mat']);
    
    curim=imread(imlist.fileloc.paths{f});
    if size(curim,3)~=1
        curim=curim(:,:,1);
    end
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.E1,...
        0.7,rawhdm.roilist(:,:,f),fieldrange);
    imwrite(outim,[outdir,'E1overlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.E2,...
        0.7,rawhdm.roilist(:,:,f),fieldrange);
    imwrite(outim,[outdir,'E2overlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.Exx,...
        0.7,rawhdm.roilist(:,:,f),fieldrange);
    imwrite(outim,[outdir,'Exxoverlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.Eyy,...
        0.7,rawhdm.roilist(:,:,f),fieldrange);
    imwrite(outim,[outdir,'Eyyoverlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.Exy,...
        0.7,rawhdm.roilist(:,:,f),fieldrange);
    imwrite(outim,[outdir,'Exyoverlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,straindata.angle*(180/pi()),...
        0.7,rawhdm.roilist(:,:,f),[0 360]);
    imwrite(outim,[outdir,'Angle.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,udisp,...
        0.7,rawhdm.roilist(:,:,f),[0 50]);
    imwrite(outim,[outdir,'Udispoverlay.tif'])
    [~,~,outim]=transover(curim,straindata.xxyy,vdisp,...
        0.7,rawhdm.roilist(:,:,f),[0 50]);
    imwrite(outim,[outdir,'Vdispoverlay.tif'])

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
    switch nargout
        case  1
            varargout{1}=straindata;
    end
    straindata.res=res;
    save([outdir,'strainfields.mat'],'straindata')
    close all;
    fclose('all');
    toc;
end
