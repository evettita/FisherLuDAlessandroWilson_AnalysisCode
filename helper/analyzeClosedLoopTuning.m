function [ voltageByBarPosition ] = analyzeClosedLoopTuning( trialFilesListClosedLoop , varargin)
%ANALYZECLOSEDLOOPTUNING extracts data from a series of trials for a
%particular recording, plots the heat maps for the user to look at and
% outputs the mean heading tuning curve for that whole batch of trials
% measured in closed loop
%
% INPUT  trialFilesList - directory and info for files the contain trials
% to be analzyzed
%       trialRemovalInfo - incase there is a portion of a trial that
%       needs to be removed this struct will contain the details of which 'trial'
%       and what 'period' need to be removed from the analysis
%
% OUTPUT  voltageByBarPosition - heading tuning curve for whole batch
%       also plots the data for the user to check
%
% Yvette Fisher 11/17/18
ephysSettings;

if( nargin > 1 )
    trialRemovalInfo = varargin{1};
    REMOVE_SECTION_OF_TRACE = true;
else
    REMOVE_SECTION_OF_TRACE = false;
end

% loop over all the files and extract the medianFiltered version of barPosition data
MEDIAN_FILTER_WIDTH_SEC = 0.04; % sec
ORDER_MEDFILTER = MEDIAN_FILTER_WIDTH_SEC * settings.sampRate;

for fileNum = 1 : length ( trialFilesListClosedLoop )
    
    cd( trialFilesListClosedLoop( fileNum ).folder );
    % load current file for current trial
    load( trialFilesListClosedLoop(fileNum).name );
    
    xPanelPos{ fileNum } =  data.xPanelPos ;
    medFilteredVoltage { fileNum } = medfilt1( data.voltage, ORDER_MEDFILTER , 'truncate' ); % Median filtering of the trace
    
    % check if need to remove section of the trial that user tagged as "bad"
    if( REMOVE_SECTION_OF_TRACE &&  sum( trialMeta.trialNum == [trialRemovalInfo(:).trial] ) > 0 )
        % find the correct index
        indexForRemoval = find( trialMeta.trialNum == [trialRemovalInfo(:).trial]);
        
        periodToRemove = trialRemovalInfo(indexForRemoval).period;% seconds
        startFrameToIgnore = periodToRemove(1)  * settings.sampRate;% start period
        endFrameToIgnore = periodToRemove(2)  * settings.sampRate;% end period
        
        % replace filteredVoltage and xPanelPos with only data from correct
        % frames
        medFilteredVoltage { fileNum } = medFilteredVoltage { fileNum }( [1 : startFrameToIgnore endFrameToIgnore: end]);
        xPanelPos{ fileNum } =  data.xPanelPos( [1 : startFrameToIgnore endFrameToIgnore: end]);
    end
    
    %load/store trial number
    trialNums( fileNum ) = trialMeta.trialNum;
    
    trial_stimPatternNum ( fileNum ) = stimulus.panelParams.patternNum;
    
    %load/store stimulus name
    trialStimulusName{ fileNum } = stimulus.name;
end

% Ajust into degrees 270 degree panels
TOTAL_DEGREES_IN_PATTERN = 360;
DEGREE_PER_LED_SLOT = 360 / 96;  %
MIDLINE_POSITION = 34;%
EDGE_OF_SCREEN_POSITION = 72; % last LED slot in postion units
EDGE_OF_SCREEN_DEG = (EDGE_OF_SCREEN_POSITION -  MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

POSSIBLE_BAR_LOCATIONS = 2:2:71;
barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

MINBARPOS = 1;
MAXBARPOS = 72;
% makes barPositions match the edges of the interval used for open-loop
% presentation
binIntervalSize_degrees = mean(  barPositionDegreesFromMidline(2 : end) - barPositionDegreesFromMidline(1 : end-1) );
barPositions = [ barPositionDegreesFromMidline - (binIntervalSize_degrees / 2), barPositionDegreesFromMidline(end) + (binIntervalSize_degrees / 2) ];

% extract data for heat map
meanVoltageByPosAcrossTrials=[];
numberOfSamplesByPosAcrossTrials=[];

for fileNum = 1 : length ( trialFilesListClosedLoop )
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

% used the information about samples in each bin to make sure the
% average is correctly weighted

% 1)take response from each trial and matrix multiple it by the number of samples in
% that bin
% 2) sum the Vm signals for each bar position,
% 3) divide the number of the total number of samples for the bar position
% in the data set
voltageWeightedBySampleNum = meanVoltageByPosAcrossTrials .* numberOfSamplesByPosAcrossTrials;
sumVoltageByBarPos =  nansum(voltageWeightedBySampleNum, 1);
sumSampleByBarPos = nansum(numberOfSamplesByPosAcrossTrials, 1);

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

end

