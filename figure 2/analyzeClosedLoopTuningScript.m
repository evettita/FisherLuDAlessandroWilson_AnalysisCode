%% analyzeClosedLoopTuningScript
% Analysis code for plotting data obtained using barRandLocON stimulus
% displays all trials from an expeirment in order as well as
% averaged tuning curves for the full expeirment
% 
% Written to look at data for Figure 2 of the manuscript, modified form
%   
% Yvette Fisher 11/8/18
%% %% Load in closed loop data set 
clear all; close all;
ephysSettings;

% INPUT wanted trials here
possibleTrials_1barBefore = 1:7;

% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = true;
trialFilesListClosedLoop = extractTrialsWithCertainStimulusName( 'closedLoop_vertStripeON_270world', includeIfNameContainsString, possibleTrials_1barBefore);

%% Plot all selected trials 
plotSelectedTrials( trialFilesListClosedLoop )

% At this stage in the script exclude any trials from analysis that do not meet quality
% criterea:
%-ClosedLoop_vertStripe270 stimulus trials are only taken from those collected before any stimuli with inverted CL or 2-3 bars, especially 2 bar closed loop stimuli was shown to the fly
%-Ave Vm for each trial is within 15 mv of first trial 
%-Ave Vm for each trial is below -20mV
%-Spike amplitude within 50% of initial size =  not dramatic evidence of decrease in Ra
%-Remove trials with large hyper pol/ depot events that look like that crazy thing
% Fly visits all bar locations in 4 min trial for the trial to be included
% Must have at least 2 ClosedLoop_vertStripe270 trials to be in the data
% set

%% loop over all the files and extract the medianFiltered version of barPosition data
MEDIAN_FILTER_WIDTH_SEC = 0.04; % sec
ORDER_MEDFILTER = MEDIAN_FILTER_WIDTH_SEC * settings.sampRate;

for fileNum = 1 : length ( trialFilesListClosedLoop )
    
    cd( trialFilesListClosedLoop( fileNum ).folder );
    % load current file for current trial
    load( trialFilesListClosedLoop(fileNum).name ); 
    xPanelPos{ fileNum } =  data.xPanelPos ;
    medFilteredVoltage { fileNum } = medfilt1( data.voltage, ORDER_MEDFILTER , 'truncate' ); % Median filtering of the trace

    %load/store trial number
    trialNums( fileNum ) = trialMeta.trialNum;
    
    trial_stimPatternNum ( fileNum ) = stimulus.panelParams.patternNum;
    
    %load/store stimulus name
    trialStimulusName{ fileNum } = stimulus.name;
end


TOTAL_DEGREES_IN_PATTERN = 360;
DEGREE_PER_LED_SLOT = 360 / 96;  % 
MIDLINE_POSITION = 34;%
EDGE_OF_SCREEN_POSITION = 72; % last LED slot in postion units
EDGE_OF_SCREEN_DEG = (EDGE_OF_SCREEN_POSITION -  MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

POSSIBLE_BAR_LOCATIONS = 2:2:71;
barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

 MINBARPOS = 1;
 MAXBARPOS = 72;
 
% makes barPositions match the edges of the interval used for open loop
% presentation
binIntervalSize_degrees = mean(  barPositionDegreesFromMidline(2 : end) - barPositionDegreesFromMidline(1 : end-1) );
barPositions = [ barPositionDegreesFromMidline - (binIntervalSize_degrees / 2), barPositionDegreesFromMidline(end) + (binIntervalSize_degrees / 2) ];

% extract data for heat map
meanVoltageByPosAcrossTrials=[];
numberOfSamplesByPosAcrossTrials=[];

for fileNum = 1 : length ( trialFilesListClosedLoop )
    % 
    patternOffsetPosition = decodePatternOffset( trial_stimPatternNum( fileNum) );
    
    xBarPostion =  xPanelPos{fileNum} + patternOffsetPosition;
    xBarPostion = mod( xBarPostion , MAXBARPOS );
    
    barPositionDegreesVsMidline = (  xBarPostion - MIDLINE_POSITION) * DEGREE_PER_LED_SLOT;
    
    meanVoltageByPos = [];
    semVoltageByPos = [];
    currLogical = [];
    numberOfSamplesInduced = [];
    secondsSpendInPosition = [];
    
    for i = 1: (length ( barPositions ) -1 )
        currLogical = ( barPositions(i) < barPositionDegreesVsMidline) & (barPositionDegreesVsMidline < barPositions (i + 1)) ;
        
        % store mean voltage values, empty = NaN
        meanVoltageByPos(i) = mean( medFilteredVoltage{fileNum}( currLogical ) );
        semVoltageByPos (i) = std( medFilteredVoltage{fileNum}( currLogical ) ) / sqrt( length (medFilteredVoltage{fileNum})) ;
        
        % store number of samples that exist for current bin
        numberOfSamplesInduced(i) = sum( currLogical );
        secondsSpendInPosition(i) = numberOfSamplesInduced(i) / settings.sampRate ; % convert samples to seconds
    end
    
    % store voltage values
    meanVoltageByPosAcrossTrials( fileNum, : ) = meanVoltageByPos; 
    numberOfSamplesByPosAcrossTrials ( fileNum, : ) = numberOfSamplesInduced;
end
  
middleOfBarPositionValues = round( mean( [barPositions( 1:end - 1) ; barPositions( 2:end ) ]) );

% Plot Heatmap
figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');

imagesc( meanVoltageByPosAcrossTrials ,'AlphaData', ~isnan(meanVoltageByPosAcrossTrials)); hold on;
set(gca,'color',0*[1 1 1]); 

c = colorbar;
c.Label.String = ' Vm  (mV)';


set(gca, 'xtick', 1: length(middleOfBarPositionValues) );
set(gca, 'xticklabel', middleOfBarPositionValues )
xlabel(' Pattern position (deg)')

set(gca, 'ytick', 1:length(trialNums) )
% build tick labels
tickLabels = [];
for i = 1: length (trialNums)
    tickLabels{i} = [ 'trial ' num2str( trialNums(i)) ',' trialStimulusName{i} ];
end
set(gca, 'yticklabel', tickLabels )
niceaxes

title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ] )

% Create a closed loop tuning curve by using the information about samples in each bin to make sure the
% average is correctly weighted based on the data amount for each bin
% 1)take response from each trial the matrix multiple it by the number of samples in
% that bin
% 2) sum the Vm signals for each bar position,
% 3) divide the number of the total number of samples for the bar position
% in the data set
voltageWeightedBySampleNum = meanVoltageByPosAcrossTrials .* numberOfSamplesByPosAcrossTrials;
sumVoltageByBarPos =  nansum(voltageWeightedBySampleNum);
sumSampleByBarPos = nansum(numberOfSamplesByPosAcrossTrials);

voltageByBarPosition = sumVoltageByBarPos ./ sumSampleByBarPos;

figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');
imagesc( voltageByBarPosition );
c = colorbar;
c.Label.String = ' Vm  (mV)';


set(gca, 'xtick', 1: length(middleOfBarPositionValues) );
set(gca, 'xticklabel', middleOfBarPositionValues )
xlabel(' Pattern position (deg)')

set(gca, 'ytick', 1:length(trialNums) )
% build tick labels
tickLabels = [];
for i = 1: length (trialNums)
    tickLabels{i} = [ 'trial ' num2str( trialNums(i)) ',' trialStimulusName{i} ];
end
set(gca, 'yticklabel', tickLabels )
niceaxes;

title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ] )


%% create plot for supplemental figure, Plot 2 EPG cells from the same fly:
TOTAL_DEGREES_IN_PATTERN = 360;
DEGREE_PER_LED_SLOT = 360 / 96; 
MIDLINE_POSITION = 34;
EDGE_OF_SCREEN_POSITION = 72; % last LED slot in postion units
EDGE_OF_SCREEN_DEG = (EDGE_OF_SCREEN_POSITION -  MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

POSSIBLE_BAR_LOCATIONS = 2:2:71;
barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

% closed loop data 4/3/18 fly 113 2 fills
twoRFs = [-36.0261827000000,-36.5012248900000,-36.4153547600000,-35.9661848300000,-37.5929094500000,-36.4996730700000,-39.2506758600000,-38.1202957100000,-38.5558632100000,-40.0056465600000,-40.6671602500000,-39.6547487900000,-41.2346608200000,-40.7858873700000,-41.4825391700000,-38.7292944800000,-35.3811916100000,-35.3661376600000,-33.5831814000000,-33.3607651700000,-33.2677294600000,-32.5916925200000,-32.8526500100000,-34.6333442000000,-36.0963436500000,-35.4182440700000,-35.0301238800000,-35.4748362600000,-33.7960856000000,-35.9654613100000,-37.0080588700000,-36.5090180000000,-39.8503516000000,-38.1774360900000,-36.8306262100000;-39.3625295100000,-39.5551826400000,-39.8609980200000,-40.7053161700000,-39.9523413000000,-40.8327924200000,-39.6940546800000,-40.5158011800000,-40.9846813100000,-40.4161667000000,-39.1655837200000,-38.7031752500000,-34.9429731300000,-41.4424115800000,-43.1507576300000,-41.0354601900000,-39.9665824400000,-36.4790737400000,-37.3108664500000,-34.2678216000000,-36.6918558900000,-35.0959459300000,-31.1766547700000,-31.0051825000000,-32.7382163400000,-33.6792377000000,-33.3999530600000,-33.6347358700000,-34.2088637300000,-33.0789394600000,-33.4910933100000,-33.1769654700000,-37.3927192000000,-37.0544063200000,-39.1644789200000];
% open loop fly 113
%twoRFs = [0.376362840000000,1.07270135300000,-0.623160552000000,-2.31035009500000,-2.06939229300000,-1.97138053400000,-0.845985682000000,-1.31160489200000,-0.519532785000000,-2.55616170300000,-2.60169263300000,-3.49347374400000,-1.96180256500000,-2.22315754800000,-2.33074791500000,-1.52437522400000,-1.10797164000000,0.559284946000000,2.40233722900000,2.33174510400000,1.83682456700000,2.71983070000000,2.46654431700000,2.40099126600000,-0.265080213000000,1.83364309000000,0.747280084000000,-0.0629772210000000,1.46412385700000,-0.637148314000000,0.152201025000000,-0.283898980000000,-1.56413479000000,-2.02194767500000,-2.55219002200000;-2.67730324900000,-2.26321198600000,-4.19545741600000,-4.82090191600000,-4.64421568500000,-4.21277224200000,-3.60416903400000,-3.05096199000000,-3.61545548800000,-2.91573140500000,-4.40956447500000,-4.59536392600000,-3.66596057500000,-1.76768763200000,-0.172912466000000,-1.71478934700000,-2.47906137700000,0.397644168000000,0.671504199000000,2.26425443900000,1.15906374200000,2.35714597500000,0.543872484000000,0.574566409000000,1.21770652200000,0.111790856000000,0.787858347000000,-0.124986525000000,0.819751325000000,-0.569442342000000,1.76171040900000,1.71558940000000,-0.948208647000000,-1.43893349800000,-1.16245587200000];

% Smooth the data with median filter if user wants to:
% order for median filter smooth
ORDER_MEDFILTER =2; %    
filteredRFs = medfilt1( twoRFs , ORDER_MEDFILTER, 'truncate', 2); % Median filtering of the trace
 
figure('Position',[50, 50, 1000, 800]); set(gcf, 'Color', 'w');
plot( barPositionDegreesFromMidline, filteredRFs ,'Color', [ 0 , 0, 0] );        
niceaxes; box off
xlim([ -120 135])

%%
function [] = plotSelectedTrials( trialFilesList )
% PLOTSELECTEDTRIALS
% helper function for OL/CL analysis
% Loop over all the files and extract the data for all the trials selected
ephysSettings;

dataWholeExp = struct();

for fileNum = 1 : length ( trialFilesList )
    
    cd( trialFilesList( fileNum ).folder );
    % load current file for current trial
    load( trialFilesList(fileNum).name );
    
    
    if( fileNum == 1) % first trial set up dataWholeExp variable with subfields that match data
        dataWholeExp = data;
    else
        % concatenate all the fields in data.___ into a new large struct will
        % whole data set in it still organized by trial order
        dataWholeExp = [ dataWholeExp , data ];
        
    end
end

% Concatinate variables from whole single experiment
voltage = cat(1, dataWholeExp(:).voltage);
xPanelPos = cat(1, dataWholeExp(:).xPanelPos);
ficTracAngularPosition = cat(1, dataWholeExp(:).ficTracAngularPosition);

% extract and concatinate Integrated X if exists
if(isfield( dataWholeExp ,'ficTracIntx') )
    ficTracIntx = cat(1, dataWholeExp(:).ficTracIntx);
end

% Plot data from whole experiment:
FigHand = figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');

timeArray = (1  :  length(voltage) ) / settings.sampRate; % seconds

ax(1) = subplot(2,1,1);
plot( timeArray,  voltage , 'DisplayName' , 'membrane Voltage' ); hold on;
ax(2) = subplot(2,1,2);
plot( timeArray,  xPanelPos , 'DisplayName' , 'panel position' ); hold on;
linkaxes(ax,'x');
legend('show')

end