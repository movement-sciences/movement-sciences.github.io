%ComparativeAvianScaling_Allbirds_25

% ** First **:  Browse to the directory containing the dataset and make it your
%current directory in Matlab

rootDir = cd;

studyInfo.loadDir = rootDir;
studyInfo.saveDir = rootDir;

%% Read in the literature metadata from the spreadsheet.
[num,~,raw] =  xlsread([ studyInfo.loadDir filesep 'LiteratureGaitDataSummary.xlsx'],'SortedByGroup','','basic'); %#ok<XLSRD>

scalingData = num(:,6:end);
scalingInfo_Hdrs = raw(1,7:end);

studyInfoData = raw(2:end,1:6);
studyInfo_Hdrs = raw(1,1:6);

t_TypeID = studyInfoData(:,(strcmpi('Classification',studyInfo_Hdrs)==1));
studyInfo.studyIndex = (1:length(t_TypeID))';

studyInfo.TypeID = studyInfoData(studyInfo.studyIndex,(strcmpi('Type',studyInfo_Hdrs)==1));
studyInfo.AuthorList = studyInfoData(studyInfo.studyIndex,(strcmpi('Author',studyInfo_Hdrs)==1));
studyInfo.YearList = studyInfoData(studyInfo.studyIndex,(strcmpi('Year',studyInfo_Hdrs)==1));
studyInfo.SpeciesID = studyInfoData(studyInfo.studyIndex,(strcmpi('Species',studyInfo_Hdrs)==1));
studyInfo.ClassID = studyInfoData(studyInfo.studyIndex,(strcmpi('Classification',studyInfo_Hdrs)==1));

bodyMass = scalingData(studyInfo.studyIndex,(strcmpi('bodyMass',scalingInfo_Hdrs)==1));
hipHeight = scalingData(studyInfo.studyIndex,(strcmpi('HipHeight',scalingInfo_Hdrs)==1));
segSum = scalingData(studyInfo.studyIndex,(strcmpi('SegSum',scalingInfo_Hdrs)==1));
Lnorm = scalingData(studyInfo.studyIndex,(strcmpi('Lnorm',scalingInfo_Hdrs)==1)); %Reference leg length based on body mass scaling
postureIdx = hipHeight./segSum; %#ok<NASGU>
lengthIdx = segSum./Lnorm; %#ok<NASGU>

nStudies = size(studyInfo.studyIndex,1);

%% Load in the gait data:

dataFileName = 'Compiled_Data_SP_SL.xlsx';
refLength = Lnorm; %Select option for reference leg length: segSum hipHeight Lnorm

%Currently we are testing just one reference length,  
% but the analysis in the paper analyzed alternative possible definitions
% of reference length and how this influenced the observed scaling trends
nRefLengths = size(refLength,2);

for j = 1:nRefLengths
    %Load data and normalise based on reference length
    [normGait(j),SI_Gait(j)] = GetGaitData(dataFileName,studyInfo,refLength(:,j)); %#ok<AGROW>
end
%% Plot the raw data in SI units

cmap = colormap(jet(nStudies));
close all;

fontProperties = {'FontName','Arial','FontSize',14};

SIGaitFig_H = figure;
set(SIGaitFig_H,'Color',[1 1 1]);
for i = 1:nStudies
    hold on;
    subplot(2,1,1)
    hold on;
    h1 = plot(SI_Gait(1).x_speed{i,:},SI_Gait(1).strideFreq{i,:},'.');
    set(h1,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
    legend off;
    if i == nStudies
        ylabel('Stride freq. (Hz)')
        legend([studyInfo.SpeciesID(1:nStudies)],'Box','off','Location','northeast');
    end
    set(gca,'Box', 'off',fontProperties{:})

    subplot(2,1,2)
    hold on;
    h2 = plot(SI_Gait(1).x_speed{i,:},SI_Gait(1).strideLength{i,:}, '.');
    set(h2,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
    legend off;
    if i == nStudies
        ylabel('stride length (m)')
        %legend([studyInfo.SpeciesID(1:nStudies)],'Box','off','Location','northeast');
    end
end
xlabel('Speed (m/s)')
    set(gca,'Box', 'off',fontProperties{:})

set(SIGaitFig_H,'InvertHardcopy','off');
figFilename = [rootDir filesep 'BirdGaitData_SI.pdf'];
figure(SIGaitFig_H)
print(figFilename,'-dpdf','-bestfit','-vector');


%% Plot data in in SI units side by side with dimensionless quantities to compare

cmap = colormap(jet(nStudies));

NormGaitFig_H = figure;
set(NormGaitFig_H,'Color',[1 1 1]);

for i = 1:nStudies
    hold on;
    subplot(2,2,1)
    hold on;
    h1 = plot(SI_Gait(1).x_speed{i,:},SI_Gait(1).strideFreq{i,:},'.');
    set(h1,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
    legend off;
    if i == nStudies
        ylabel('Stride freq. (Hz)')
        %legend([studyInfo.SpeciesID(1:nStudies)],'Box','off','Location','northeast');
    end
    set(gca,'Box', 'off',fontProperties{:})

    subplot(2,2,3)
    hold on;
    h2 = plot(SI_Gait(1).x_speed{i,:},SI_Gait(1).strideLength{i,:}, '.');
    set(h2,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
    legend off;
    if i == nStudies
        ylabel('stride length (m)')
    end
    xlabel('Speed (m/s)')
    set(gca,'Box', 'off',fontProperties{:})

  subplot(2,2,2)
    hold on;
    h1 = plot(normGait(1).x_speed{i,:},normGait(1).strideFreq{i,:},'.');
    set(h1,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
    legend off;
    if i == nStudies
        ylabel('relative stride freq.')
        %legend([studyInfo.SpeciesID(1:nStudies)],'Box','off','Location','northeast');
    end
    set(gca,'Box', 'off',fontProperties{:})

    subplot(2,2,4)
    hold on;
    h2 = plot(normGait(1).x_speed{i,:},normGait(1).strideLength{i,:}, '.');
    set(h2,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
    legend off;
    if i == nStudies
        ylabel('relative stride length')
    end
    xlabel('dimensionless speed')
    set(gca,'Box', 'off',fontProperties{:})

end

set(NormGaitFig_H,'InvertHardcopy','off');
figFilename = [rootDir filesep 'BirdGaitData_SI_v_Norm.pdf'];
figure(NormGaitFig_H)
print(figFilename,'-dpdf','-bestfit','-vector');


%% Curve fits
%Uncomment if you want to play with fitting curves to the data

% SF_fitSlope =  nan(nStudies,nRefLengths);
% SF_fitConstant =  nan(nStudies,nRefLengths);
% SL_fitSlope =  nan(nStudies,nRefLengths);
% SL_fitConstant =  nan(nStudies,nRefLengths);
% 
% % Create curve fits
%     for i = 1:nStudies
%         fitType = 'poly1';
%         myfittype = fittype('a + b*log(x)',...
%             'dependent',{'y'},'independent',{'x'},...
%             'coefficients',{'a','b'});
%         mode = 1;
%         %Fit 2nd order polynomial to gait data
%         [p2Fit_SF{i,j}, p2_SF_gof(i),p2_SF_output(i)] = createPolyFit(normGait(j).x_speed{i,:}, normGait(j).strideFreq{i,:},myfittype,mode); %#ok<*AGROW>
%         [p2Fit_SL{i,j}, p2_SL_gof(i),p2_SL_output(i)] = createPolyFit(normGait(j).x_speed{i,:}, normGait(j).strideLength{i,:},fitType,mode);
% 
%         t_cfs = coeffvalues(p2Fit_SL{i,j});
%         SL_fitSlope(i,j) = t_cfs(1,1);
%         SL_fitConstant(i,j) = t_cfs(1,2);
%         clear t_cfs;
%         %Note coeff and constant terms in reverse order for SF fit because
%         % it uses a differently specified fit model
%         t_cfs = coeffvalues(p2Fit_SF{i,j});
%         SF_fitSlope(i,j) = t_cfs(1,2);
%         SF_fitConstant(i,j) = t_cfs(1,1);
%         clear t_cfs;
%     end
% 
%     %% plot data with resulting curve fits
%     cmap = colormap(jet(nStudies));
% 
%     spdFitFH =figure;
%     hold on;
%     for i = 1:nStudies
%         subplot(2,1,1)
%         hold on;
%         h1 = plot(p2Fit_SF{i,j},normGait(j).x_speed{i,:},normGait(j).strideFreq{i,:},'.','PredFunc');
%         set(h1,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
%         legend off;
%         if i == nStudies
%             ylabel('rel. stride freq.')
%         end
% 
%         subplot(2,1,2)
%         hold on;
%         h2 = plot(p2Fit_SL{i,j},normGait(j).x_speed{i,:},normGait(j).strideLength{i,:}, '.','PredFunc');
%         set(h2,'MarkerFaceColor',cmap(i,:),'Color',cmap(i,:));
%         legend off;
%         if i == nStudies
%             ylabel('rel. stride length')
%         end
%     end
%     xlabel('rel. speed')
%     [~] = SaveFigAsPDF(spdFitFH,studyInfo.saveDir,'SpeedFits');

% %% Plot dimensionless data as a function of log body mass,
% % comparing species at a single reference dimensionless speed
% 
% %First get the SF and SL values at the reference speed for each species
% refSpeed = 0.75;
% refSpdData_SF = nan(nStudies,nRefLengths);
% refSpdData_SL = nan(nStudies,nRefLengths);
% 
% for j = 1:nRefLengths
%     for i = 1:nStudies
%         refSpdData_SF(i,j) = feval(p2Fit_SF{i,j},refSpeed)';
%         refSpdData_SL(i,j) = feval(p2Fit_SL{i,j},refSpeed)';
%     end
% end
% 
% %Plot the data as a function of body mass at the reference speed.
% classId = categorical(studyInfo.ClassID);
% logBodyMass = log10(bodyMass);
% 
% scalingFig = figure;
% hold on;
% subplot(2,1,1)
% gscatter(logBodyMass,refSpdData_SF,classId);
% ylabel('relative stride length');  xlabel('log body mass (kg)')
% set(gca,'XTick',log10([0.1 1 10 100]), 'XTickLabel',{'0.1' '1' '10' '100'});
% subplot(2,1,2)
% gscatter(logBodyMass,refSpdData_SF,classId);
% ylabel('relative stride length');  xlabel('log body mass (kg)')
% set(gca,'XTick',log10([0.1 1 10 100]), 'XTickLabel',{'0.1' '1' '10' '100'});
% 
% [~] = SaveFigAsPDF(scalingFig,studyInfo.saveDir,'SL_SF_v_BodyMass_Scatter');

%% support function below:

function [normGait,SI_Gait] = GetGaitData(dataFileName,studyInfo, refLength)

nStudies = size(studyInfo.studyIndex,1);

normGait.x_speed = cell(nStudies,1);
normGait.stridePeriod = cell(nStudies,1);
normGait.strideFreq = cell(nStudies,1);
normGait.strideLength = cell(nStudies,1);
normGait.vMax = [];

SI_Gait.x_speed = cell(nStudies,1);
SI_Gait.stridePeriod = cell(nStudies,1);
SI_Gait.strideFreq = cell(nStudies,1);
SI_Gait.strideLength = cell(nStudies,1);
SI_Gait.vMax = [];

% normGait(j).swingFreq = cell(nStudies,1);
% normGait(j).stanceFreq = cell(nStudies,1);
% normGait(j).dutyFactor = cell(nStudies,1);

for i = 1:nStudies
    %Load gait data for current study from Excel spreadsheet
    t_sheetName = [studyInfo.AuthorList{i} '_' num2str(studyInfo.YearList{i}) '_' studyInfo.SpeciesID{i}];
    [gaitData,gaitHdrs,~] =  xlsread([studyInfo.loadDir filesep dataFileName],t_sheetName,'','basic');
    %{'Speed_mps','SwingPeriod','StancePeriod','strideLength_m','DF','StridePeriod_s'}

    %organize gait data
    x_speed = gaitData(:,(strcmpi('Speed_mps',gaitHdrs)==1));
    stridePeriod = gaitData(:,(strcmpi('StridePeriod_s',gaitHdrs)==1));
    strideLength = gaitData(:,(strcmpi('strideLength_m',gaitHdrs)==1));

    SI_Gait.x_speed{i,:} = x_speed;
    SI_Gait.stridePeriod{i,:} = stridePeriod;
    SI_Gait.strideFreq{i,:} = stridePeriod.^-1;
    SI_Gait.strideLength{i,:} = strideLength;

    %Calculate normalised variables:
    normGait.vNorm = ((9.81.*refLength(i)).^0.5);
    normGait.fNorm = ((9.81./refLength(i)).^0.5);
    normGait.tNorm = (((9.81./refLength(i)).^0.5)).^-1;
    normGait.lNorm = refLength(i);

    normGait.x_speed{i,:} = x_speed./normGait.vNorm;
    normGait.stridePeriod{i,:} = stridePeriod./normGait.tNorm;
    normGait.strideFreq{i,:} = (stridePeriod.^-1)./normGait.fNorm;
    normGait.strideLength{i,:} = strideLength./normGait.lNorm;

    normGait.vMax(i) = max(normGait.x_speed{i,:});
    SI_Gait.vMax(i) = max(SI_Gait.x_speed{i,:});

end

end


function [fitresult, gof,output] = createPolyFit(xData, yData,fitName,mode)
%CREATEFIT(C_W_V,C_W_SL)
%  Create a fit.%
%  Data for 'untitled fit 1' fit:
%      xData
%      yData
%  fitName: e.g., 'poly1' or 'poly2'
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xData, yData );

% Set up fittype
ft = fittype( fitName );

% % Set up options.
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% if mode == 2
% opts.Robust = 'Bisquare';
% end

% Fit model to data.
%[fitresult, gof,output] = fit( xData, yData, ft, opts );

[fitresult, gof,output] = fit( xData, yData, ft );

end

