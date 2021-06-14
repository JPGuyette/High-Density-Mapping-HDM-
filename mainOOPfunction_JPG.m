%% Orientational Order Parameter (OOP) Calculator
% Beta version 2.0
% date: 06/13/2013
%
% The following code is part of the development of COBRA: the Cytoskeletal
% Organization Benchmark and Representation Assay.
%
% This file contain a complete version of the main software that does not
% uses a GUI. This is helpful for experienced users, for debgging and
% upgrading the software. The code features several functions that
% communicate updating a single 'handlesLikeStructure' structure. This should make them
% easily compatible with the GUIDE-generated GUIs.

% written by: Francesco S. Pasqualini
%
% DISEASE BIOPHYSICS GROUP CONFIDENTIAL
% Unpublished Copyright (c) 2013,
% DISEASE BIOPHYSICS GROUP, All Rights Reserved.
%
% NOTICE:  All information contained herein is, and remains the property of
% the DISEASE BIOPHYSICS GROUP. The intellectual and technical concepts
% contained herein are proprietary to the DISEASE BIOPHYSICS GROUP and may
% be covered by U.S. and Foreign Patents, patents in process, and are
% protected by trade secret or copyright law.
% Dissemination of this information or reproduction of this material is
% strictly forbidden unless prior written permission is obtained
% from the DISEASE BIOPHYSICS GROUP.  Access to the source code contained
% herein is hereby forbidden to anyone except current
% DISEASE BIOPHYSICS GROUP employees, managers or contractors who have
% executed Confidentiality and Non-disclosure agreements explicitly
% covering such access.
%
% The copyright notice above does not evidence any actual or intended
% publication or disclosure  of  this source code, which includes
% information that is confidential and/or proprietary, and is a trade
% secret, of  the DISEASE BIOPHYSICS GROUP.
% ANY REPRODUCTION, MODIFICATION, DISTRIBUTION, PUBLIC  PERFORMANCE,
% OR PUBLIC DISPLAY OF OR THROUGH USE  OF THIS  SOURCE CODE  WITHOUT  THE
% EXPRESS WRITTEN CONSENT OF COMPANY IS STRICTLY PROHIBITED,
% AND IN VIOLATION OF APPLICABLE LAWS AND INTERNATIONAL TREATIES.
% THE RECEIPT OR POSSESSION OF  THIS SOURCE CODE AND/OR RELATED INFORMATION
% DOES NOT CONVEY OR IMPLY ANY RIGHTS  TO REPRODUCE, DISCLOSE OR
% DISTRIBUTE ITS CONTENTS, OR TO MANUFACTURE, USE, OR SELL ANYTHING THAT
% IT  MAY DESCRIBE, IN WHOLE OR IN PART.




function handlesLikeStructure = mainOOPfunction()
%% INITIALIZE
[PATHSTR,NAME,EXT] = fileparts(pwd);
addpath(genpath([PATHSTR,...
    '\Dependencies']));

clc
clear all
close all


    handlesLikeStructure.customSettingsFlag = 1;
    handlesLikeStructure = changeFitSettings(handlesLikeStructure);

handlesLikeStructure.typeOfAnalysis = menu('mtitle','image','batch');


if handlesLikeStructure.typeOfAnalysis == 1

    handlesLikeStructure.fileIndex = 0;
    % the name of the fields saved in in the handlesLikeStructure structure saved by
    % each function are commented after its usage
   
    handlesLikeStructure = selectImage(handlesLikeStructure);
    % orglImage, folderName and fileName
    handlesLikeStructure = thresSelection(handlesLikeStructure);
% thValues, thCoherency
    handlesLikeStructure = analyzeImageFile(handlesLikeStructure);
    % fltImage, ortSet, ortDensity, qltMetricCoverage,
    % qltMetricIntensity, qltMetricCoherency
    
    handlesLikeStructure = calculateOOPs(handlesLikeStructure);
    % metricsVector
     
elseif handlesLikeStructure.typeOfAnalysis == 2

    handlesLikeStructure = selectFolder(handlesLikeStructure);
    % folderName, fileList
    numFile = length(handlesLikeStructure.fileList);
    handlesLikeStructure = thresSelection(handlesLikeStructure);
% thValues, thCoherency
    for (fileIndex = 1:numFile)
        handlesLikeStructure.fileIndex = fileIndex; % fileIndex
        
        handlesLikeStructure = selectImage(handlesLikeStructure);
        % orglImage and fileName

        handlesLikeStructure = analyzeImageFile(handlesLikeStructure);
        % fltImage, ortSet, ortDensity, qltMetricCoverage,
        % qltMetricIntensity, qltMetricCoherency
                
        handlesLikeStructure = calculateOOPs(handlesLikeStructure);
        % metricsVector
        pause(1)
    end
else
    errordlg('please restart and select one option')
end

handlesLikeStructure = datasetFormation(handlesLikeStructure);

[saveFileName saveFilePath] = uiputfile;

saveFileString = [saveFilePath, saveFileName,'.csv'];
export(handlesLikeStructure.Dataset,'file', saveFileString,'Delimiter',',')

function handlesLikeStructure = selectFolder(handlesLikeStructure) % folderName, fileList

handlesLikeStructure.folderName = uigetdir(getuserdir);

% handlesLikeStructure.fileList = dir(fullfile(...
%     handlesLikeStructure.folderName, '*.tif*'));

handlesLikeStructure.fileList = dir(fullfile(...
    handlesLikeStructure.folderName, '*.tif*'));

function handlesLikeStructure = selectImage(handlesLikeStructure)


if handlesLikeStructure.fileIndex == 0
    % opens user-selected image
    
    [handlesLikeStructure.fileName, handlesLikeStructure.folderName] = ...
        uigetfile({'*.bmp;*.jpg;*.tif;*.tif*;*.png;','All Image Files';...
        '*.*','All Files' },...
        'Pick an image',...
        getuserdir);
    
    currentFilename = [handlesLikeStructure.folderName,...
        handlesLikeStructure.fileName];
    handlesLikeStructure.orglImage = imread(currentFilename);
%     figure(1)
%     imshow(handlesLikeStructure.orglImage,[])
    
elseif  handlesLikeStructure.fileIndex > 0
    % or opens the next image in the queue
    
    currentFilename = [handlesLikeStructure.folderName,'\',...
        handlesLikeStructure.fileList(handlesLikeStructure.fileIndex).name];
    handlesLikeStructure.orglImage = imread(currentFilename);
%     figure(1)
%     imshow(handlesLikeStructure.orglImage,[])
else
    % alternatively returns an error
    errordlg(['for developer: somehow the counter for',...
        ' the fileList has a negative value'])
end

function handlesLikeStructure = thresSelection(handlesLikeStructure) % thValues, thCoherency

% prompt the user to imput numerical values for threshold
prompt = {'Intensity:','Coherency:'};
dlg_title = 'Specify orientation inclusion threshold';
num_lines = 1;
def = {'0.2','0.4'};
UserSelections = inputdlg(prompt,dlg_title,num_lines,def);
UserSelections = str2double(UserSelections);

% pass the thresholds to a structure that will be used later in the
% computation
handlesLikeStructure.thValues = UserSelections(1);
handlesLikeStructure.thCoherency = UserSelections(2);

function handlesLikeStructure = analyzeImageFile(handlesLikeStructure)

im_HSV = rgb2hsv(handlesLikeStructure.orglImage);
handlesLikeStructure.fltImage = zeros(size(im_HSV));
% extract Orientations
Orientation = im_HSV(:,:,1);
% extract Coherency
Coherency = im_HSV(:,:,2);
% extract Gray Values for the original image
Values = im_HSV(:,:,3);

% apply defined thresholds
NonBackgroundOrientations = ...
    Orientation(Values>handlesLikeStructure.thValues);
HighConfidenceEdgeOrientations = ...
    Orientation(Values>handlesLikeStructure.thValues |...
    Coherency>=(handlesLikeStructure.thCoherency*max(Coherency(:))));

handlesLikeStructure.fltImage(:,:,1) = Orientation;
handlesLikeStructure.fltImage(:,:,2) = Coherency.*(Coherency>=...
    (handlesLikeStructure.thCoherency*max(Coherency(:))));
handlesLikeStructure.fltImage(:,:,3) = Values.*...
    (Values>=handlesLikeStructure.thValues);

handlesLikeStructure.fltImage = hsv2rgb(handlesLikeStructure.fltImage);
handlesLikeStructure.ortSet = HighConfidenceEdgeOrientations;


%% CONTROL THIS STEP
% change the expression in the hsv format [0 1] to radians in the -pi, pi
% range
intervalRad = linspace(-pi,pi,360);

fsp_hue_360rad = @(x) mod(x*2*pi,2*pi)-pi;
y = handlesLikeStructure.ortSet;
if(max(y(:))<=1 & min(y(:))>=0)
    y = fsp_hue_360rad(y);
end
[counts_angles] = hist(y,intervalRad);
handlesLikeStructure.ortDensity= counts_angles/trapz(intervalRad,counts_angles);
% handlesLikeStructure.ortDensity= handlesLikeStructure.ortDensity(:);
handlesLikeStructure.ortDensity = handlesLikeStructure.ortDensity;
handlesLikeStructure.ortSet = y;
handlesLikeStructure.intervalRad = intervalRad;

% TISSUE QUALITY METRICS
% the following 2 numbers are to be considered for further implementation.
% The idea is that good tissues required minimal thresholding so, for a
% given threshold, the code will remove more pixels from a bad tissue than
% from a good one. To quanitify this we can easily do as
% 1-%ofPixelRemovedAtGivenThreshold.

handlesLikeStructure.qltMetricCoverage = ...
    numel(Values(Values>0.05))/numel(Values);

handlesLikeStructure.qltMetricIntensity = 1-...
    numel(NonBackgroundOrientations)/numel(Orientation);
handlesLikeStructure.qltMetricCoherency = ...
    handlesLikeStructure.qltMetricIntensity-...
    ((numel(HighConfidenceEdgeOrientations)-...
    numel(NonBackgroundOrientations))/...
    numel(Orientation));

figure(2)
subplot(1,3,1)
imshow(handlesLikeStructure.fltImage,[])


function handlesLikeStructure = calculateOOPs(handlesLikeStructure)
% metricsVector
intervalRad = handlesLikeStructure.intervalRad;
intervalRadPlot = rad2deg(intervalRad/2); 
% return the axial data in the range [-pi, pi]
OOPr = circ_r(handlesLikeStructure.ortSet(:));

try
    
    [cf_, gof_] = Trimodal_VonMisesFit(handlesLikeStructure);
    cfvalues = coeffvalues(cf_);
    fhat = cfvalues(1);
    khat1 = cfvalues(2);
    khat2 = abs(cfvalues(3));
    
    mhat1 = cfvalues(4);
    mhat2 = ((-1)^(cfvalues(3)<0).*cfvalues(5));
    
    ghat = cfvalues(6);
    OOP1 = circ_r(circ_vmrnd(mhat1,khat1,1e3));
    OOP2 = circ_r(circ_vmrnd(mhat2,khat2,1e3));
    gof = gof_.rsquare;
    
    gamma = 1 - fhat - ghat;
    
    if (handlesLikeStructure.fileIndex == 0)
        handlesLikeStructure.metricsVector = ...
            [OOPr, OOP1, OOP2,gamma,gof,...
            handlesLikeStructure.thValues...
            handlesLikeStructure.thCoherency,...
            handlesLikeStructure.qltMetricCoverage,...
            handlesLikeStructure.qltMetricIntensity,...
            handlesLikeStructure.qltMetricCoherency,...
            ];
    else
        handlesLikeStructure.metricsVector(...
            handlesLikeStructure.fileIndex,:) = ...
            [OOPr, OOP1, OOP2,gamma,gof,...
            handlesLikeStructure.thValues...
            handlesLikeStructure.thCoherency,...
            handlesLikeStructure.qltMetricCoverage,...
            handlesLikeStructure.qltMetricIntensity,...
            handlesLikeStructure.qltMetricCoherency,...
            ];
    end
    % plot the values of all OOP indicators, as well as the graph of OOP vs
    % spread paraneter k and circular standard deviation
    
    bestFit_VM_1 = 0.5*ghat*circ_vmpdf(intervalRad,0,0)+...
        fhat*circ_vmpdf(intervalRad,mhat1,khat1);
    bestFit_VM_2 = 0.5*ghat*circ_vmpdf(intervalRad,0,0)+...
        (1-fhat-ghat)*circ_vmpdf(intervalRad,mhat2,khat2);
    AreaFit = trapz(intervalRad,bestFit_VM_1)+...
        trapz(intervalRad,bestFit_VM_2);
    
    
    figure(2)
    subplot(1,3,2:3)
    bar(intervalRadPlot,handlesLikeStructure.ortDensity,...
        'EdgeColor','g','FaceColor','g')
    hold on
    plot(intervalRadPlot,bestFit_VM_1,'k','LineWidth',2)
    plot(intervalRadPlot,bestFit_VM_2,'r','LineWidth',2)
    hold off
    legend({'OOP', 'OOP 1', 'OOP 2'},'FontSize',14,'FontWeight','bold')
    legend boxoff
    
catch err
    if (handlesLikeStructure.fileIndex == 0)
        handlesLikeStructure.metricsVector =...
            [OOPr,zeros(1,4),...
            handlesLikeStructure.thValues,...
            handlesLikeStructure.thCoherency,...
            handlesLikeStructure.qltMetricCoverage,...
            handlesLikeStructure.qltMetricIntensity,...
            handlesLikeStructure.qltMetricCoherency,...
            ];
    else
        handlesLikeStructure.metricsVector(...
            handlesLikeStructure.fileIndex,:) = ...
            [OOPr,zeros(1,4),...
            handlesLikeStructure.thValues,...
            handlesLikeStructure.thCoherency,...
            handlesLikeStructure.qltMetricCoverage,...
            handlesLikeStructure.qltMetricIntensity,...
            handlesLikeStructure.qltMetricCoherency,...
            ];
    end
    display(['could not perform fitting:']);
    getReport(err, 'extended')
    
    
    figure(2)
    subplot(1,3,2:3)
    bar(intervalRadPlot,handlesLikeStructure.ortDensity,'EdgeColor','g','FaceColor','g')
    % plot(intervalRad,handlesLikeStructure.ortDensity,'g')
    
    legend({'OOP', 'OOP 1', 'OOP 2'},'FontSize',14,'FontWeight','bold')
    legend boxoff
    
    handlesLikeStructure.intervalRadPlot = intervalRadPlot;
end


function handlesLikeStructure = datasetFormation(handlesLikeStructure)

VariableNames = ...
    {'OOP','OOP1','OOP2','gamma','rsq',...
    'thIntensity','thCoherency','qltMetricCoverage',...
    'qltMetricIntensity','qltMetricCoherency'};

handlesLikeStructure.Dataset = ...
    dataset({handlesLikeStructure.metricsVector, VariableNames{:}});

if (handlesLikeStructure.fileIndex == 0)
 handlesLikeStructure.Dataset.FileNames = ...
     handlesLikeStructure.fileName;
else
    handlesLikeStructure.Dataset.FileNames = ...
        {handlesLikeStructure.fileList.name}';
end
handlesLikeStructure.Dataset = ...
    [handlesLikeStructure.Dataset(:,end),...
    handlesLikeStructure.Dataset(:,1:end-1)];

function handlesLikeStructure = changeFitSettings(handlesLikeStructure)

if(handlesLikeStructure.customSettingsFlag == 1) 
    handlesLikeStructure.userDefinedLower = [0 2 1 deg2rad(-45) deg2rad(45) 0];
    handlesLikeStructure.userDefinedUpper = [1 100 100 deg2rad(45) deg2rad(315) 1];
    handlesLikeStructure.userDefinedGuess = [0.1, 10,10,deg2rad(0),deg2rad(180),.1];
else
    handlesLikeStructure.userDefinedLower = [0 1 1 deg2rad(-45) deg2rad(135) 0];
    handlesLikeStructure.userDefinedUpper = [1 100 100 deg2rad(45) deg2rad(225) 1];
    handlesLikeStructure.userDefinedGuess = [0.1, 5,5,deg2rad(0),deg2rad(90),.1];
end

function [cf_, gof_] = ...
    Trimodal_VonMisesFit(handlesLikeStructure)


% %% default for now (in the GUI are manually select)
% handlesLikeStructure.userDefinedLower = [0 1 1 deg2rad(-90) deg2rad(91) 0];
% handlesLikeStructure.userDefinedUpper = [1 100 100 deg2rad(90) deg2rad(180) 1];
% handlesLikeStructure.userDefinedGuess = [0.1, 5,5,deg2rad(0),deg2rad(90),.1];
%
intervalRad = handlesLikeStructure.intervalRad(:);
handlesLikeStructure.ortDensity = handlesLikeStructure.ortDensity(:);

ok_ = isfinite(intervalRad) & isfinite(handlesLikeStructure.ortDensity);

if(handlesLikeStructure.customSettingsFlag == 1)
    opt = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',handlesLikeStructure.userDefinedLower,...
        'Upper',handlesLikeStructure.userDefinedUpper,...
        'StartPoint',handlesLikeStructure.userDefinedGuess);
else
    opt = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',handlesLikeStructure.defaultLower ,...
        'Upper',handlesLikeStructure.defaultUpper,...
        'StartPoint',handlesLikeStructure.defaultGuess);
end

ft_ = fittype(['0*circ_vmpdf(x,0,0) + f*circ_vmpdf(x,mhat1,khat1)+'...
    '(1-f-g)*circ_vmpdf(x,mhat2,khat2)'],...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'f', 'khat1', 'khat2', 'mhat1', 'mhat2','g'},...
    'options',opt);


[cf_, gof_] = fit(intervalRad(ok_),handlesLikeStructure.ortDensity(ok_),ft_);




%% SUPPORT ROUTINES
% % In the final version unccomment this line so that the code starts from
% % the user home directory (in Windows)

function userDir = getuserdir

if ispc
    userDir = winqueryreg('HKEY_CURRENT_USER',...
        ['Software\Microsoft\Windows\CurrentVersion\' ...
         'Explorer\Shell Folders'],'Personal');
else
    userDir = char(java.lang.System.getProperty('user.home'));
end



%function userDir = getuserdir
% During development let's use this to simplify things
% 
% [PATHSTR,NAME,EXT] = fileparts(pwd);
% userDir = [PATHSTR,...
%     '\DemoImages\OOP\pCM'];
% 
% userDir = 'C:\Users\Francesco\Documents\COBRA_Apr2013Material\aaaaaa\Rotated\Analysis - MosaicISO\ImageJpreProcessed'









