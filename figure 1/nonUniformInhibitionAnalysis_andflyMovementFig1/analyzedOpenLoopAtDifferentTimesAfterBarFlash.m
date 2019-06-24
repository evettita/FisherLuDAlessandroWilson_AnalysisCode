%% analyzedOpenLoopAtDifferentTimesAfterBarFlash
% Analysis code for plotting data obtained using barRandLocON stimulus
% displays both all trials from an expeirment in order as well as
% averaged tuning curves for the full expeirment
%
%  added to analysis to look at Vm in addition to delta Vm and
% to look at persistance of effect of a bar flash on Vm....
% 
%
% Yvette Fisher 3/5/19
%% 1) Load in open loop data set 
clear all
close all
ephysSettings;

% INPUT wanted trials here
possibleTrials_1barBefore = [1,3,5,7];


% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = false;
trialFilesList = extractTrialsWithCertainStimulusName( 'barRandLocON()', includeIfNameContainsString, possibleTrials_1barBefore);

%% extract delta Vm and Absolute Vm 
ephysSettings;

for fileNum = 1 : length ( trialFilesList )
    
    cd( trialFilesList( fileNum ).folder );
    % load current file for current trial
    load( trialFilesList(fileNum).name ); 
   
    % save TrialTime
    processedData(fileNum).trialStartTime = trialMeta.trialStartTime; 
    % same Trial number 
    processedData(fileNum).trialNum = trialMeta.trialNum; 
    processedData(fileNum).exptInfo = exptInfo;
    processedData(fileNum).stimulus = stimulus;
    
    % possible X-Positions
    minXPos = 1;
    maxXPos = 71;
    DURATION_TO_PLOT_BEFORE_FLASH = 0.5;
    DURATION_TO_PLOT_PER_EPOCH = 1; % sec
    
    % set how long the stimulus step (data.xPanelPos) needs to be constant
    % after change to be included in data set
    SEARCH_LENGTH = 0.3; % seconds
    seqFramesToLookFor = SEARCH_LENGTH * settings.sampRate;
    
    FINAL_STIMULUS_PERIOD =  4; %seconds
    numOfFramesToIgnoreAtEnd = FINAL_STIMULUS_PERIOD * settings.sampRate;
    
    for j = ( minXPos : maxXPos)  % loop over all possible X-pos values and pull out traces for each that has data
        % find index where the pos value has changed
        change = find( diff (data.xPanelPos) ~= 0 ) + 1;
        
        % find which of those indexes the panels was at the currect position j
        indofJ = change (data.xPanelPos (change)  == j) ;
        
        % remove any indexs too close to the end of the trace to cause errors
        indofJ  = indofJ (indofJ < ( numel( data.xPanelPos )  - numOfFramesToIgnoreAtEnd));
        
        epochStartInds = [];
        
        % only take those where there are not change in position values for 100 ms after the start
        for ii = 1: length( indofJ )
            averagePosDecode = mean ( data.xPanelPos (indofJ(ii) + 1 : indofJ(ii) + seqFramesToLookFor ) );
            
            if(averagePosDecode == j) % if not other changes then take this index             
                epochStartInds = [epochStartInds indofJ(ii) ];
            end
        end

        % check that this epoch step was not empty for some reason
        if( ~ isempty(epochStartInds) )
            
            stimulusStep= [];
            voltageByPosition = [];
            for i = 1 : ( numel(epochStartInds) )
            
            currEpochStart = epochStartInds(i) - (DURATION_TO_PLOT_BEFORE_FLASH * settings.sampRate);
            currEpochEnd = epochStartInds(i) + ( DURATION_TO_PLOT_PER_EPOCH * settings.sampRate);
            
            % store current Voltage trace
            voltageByPosition(i, :) = data.scaledVoltage( currEpochStart : currEpochEnd);
           
            % same stimulus step
            stimulusStep(i, :)  = data.xPanelPos( currEpochStart: currEpochEnd );
            
            end
             % save the data traces for this bar position:
            barPositionResp{j} = voltageByPosition; 
            stimulusStepByPosition{j} = stimulusStep;
        end
             
    end
    % same traces & stimulus steps for all bar positions
    processedData(fileNum).barPositionResp =  barPositionResp;
    processedData(fileNum).stimulusStep = stimulusStepByPosition;
end

% PLOTTING over trials across the experiment in order
clearvars -except processedData exptInfo stimulus settings trialMeta trialFilesList
% plot average RF using the whole data set:

% store traces, and average Vm for each time period "pre" "post" and during stimulus
POSSIBLE_BAR_LOCATIONS = 2:2:71;

for kk = POSSIBLE_BAR_LOCATIONS 
    % Plotting responses over the time course of the whole experiment:
    BAR_POSITION_TO_PLOT = kk;
    
    % start time of the first trial
    firstTrialStartTime = processedData(1).trialStartTime;

    for i = 1 : length( processedData )
        
        currTraces = processedData(i).barPositionResp{BAR_POSITION_TO_PLOT};
        currTimeArray = (1  :  length(currTraces) ) / settings.sampRate; % seconds

        %remove median filter step:
        currTracesFiltered = currTraces;
        
        % Extract average from each epoch period
        PRE_BASELINE_START_TIME = 0.25;
        PRE_BASELINE_END_TIME = 0.5;
        
        POST_BASELINE_START_TIME = 1.25;
        POST_BASELINE_END_TIME = 1.5;
        
        RESP_START_TIME = 0.75;
        RESP_END_TIME = 1;
        
        % save average values of Vm for summary plot
        preBaselinesIndexes = PRE_BASELINE_START_TIME < currTimeArray & currTimeArray < PRE_BASELINE_END_TIME ;
        postBaselineIndexes = POST_BASELINE_START_TIME < currTimeArray & currTimeArray < POST_BASELINE_END_TIME ;
        responsePeriodIndexes = RESP_START_TIME < currTimeArray & currTimeArray < RESP_END_TIME ;
        
        % save pre baseline mean and std
        preBaselineVoltage = mean ( currTracesFiltered (:, preBaselinesIndexes) , 2 );
        ave_preBaselineVoltage(i) = mean ( preBaselineVoltage );
        sem_preBaselineVoltage(i) =  std( preBaselineVoltage) / sqrt( numel( preBaselineVoltage));
        
        % same post baseline mean & std
        postBaselineVoltage = mean ( currTracesFiltered (:, postBaselineIndexes) , 2 );
        ave_postBaselineVoltage(i) = mean ( postBaselineVoltage );
        sem_postBaselineVoltage(i) =  std( postBaselineVoltage) / sqrt( numel( postBaselineVoltage));
        
        % same responses mean & std
        voltage = mean ( currTracesFiltered (:, responsePeriodIndexes) , 2 );
        ave_Voltage(i) = mean ( voltage );
        sem_Voltage(i) =  std( voltage) / sqrt( numel( voltage));
        
        % calculate time since first trial of this type
        timeElapsedSinceFirstExpTrial(i) =  caldiff( [datetime( firstTrialStartTime ), datetime( processedData(i).trialStartTime )] ) ;
        trialNum(i) = processedData(i).trialNum;
    end

    trialTimeInMinutes = minutes( time(timeElapsedSinceFirstExpTrial) );
    
    % save voltage difference from each time point to build heat map
    voltageDiff(BAR_POSITION_TO_PLOT / 2, :) = ave_Voltage - ave_preBaselineVoltage;
    % save voltage at that time point
    membraneVoltageResp(BAR_POSITION_TO_PLOT / 2, :) = ave_Voltage;
    membraneVoltagePost(BAR_POSITION_TO_PLOT / 2, :) =  ave_postBaselineVoltage;
    
end

%% % Plot Delta Vm
 deltaVm = true;
 plotHeatMapRFs( voltageDiff,  stimulus, exptInfo, trialNum , deltaVm )

%% % Plot Vm bar response
 deltaVm = false;
plotHeatMapRFs( membraneVoltageResp,  stimulus, exptInfo, trialNum , deltaVm );

%% Plot Vm 250ms post flash
 deltaVm = false;
meanResponse = plotHeatMapRFs( membraneVoltagePost,  stimulus, exptInfo, trialNum , deltaVm );
% 

%%
function [ meanResponse ] = plotHeatMapRFs( response, stimulus, exptInfo, trialNum , deltaVm)
figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');

POSSIBLE_BAR_LOCATIONS = 2:2:71;
DEGREE_PER_LED_SLOT = 360 / 96;  %
MIDLINE_POSITION = 34; % LED position where the fly is aligned to for all 230 deg EPG recordings 11/2017  - present

barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;
imagesc(response'); box off; hold on

c = colorbar;

if (deltaVm)
    colormap('bluewhitered'); % zero set to white in this plotting
    c.Label.String = ' Delta Vm  (mV) ( response - baseline)';
else
    c.Label.String = ' Vm  (mV)';
end

set(gca, 'xtick', POSSIBLE_BAR_LOCATIONS / 2 )
set(gca, 'xticklabel', round(barPositionDegreesFromMidline))
set(gca, 'ytick', 1:length(trialNum) )

% build tick labels
tickLabels = [];
for i = 1: length (trialNum)
    tickLabels{i} = [ 'trial ' num2str( trialNum(i))  ];
end

set(gca, 'yticklabel', tickLabels )
xlabel( ' bar center (deg) ' );
title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ' stim: ' num2str(stimulus.name) ] )
niceaxes

% AVERAGE across all Open LOOP trials
figure('Position',[50, 50, 1000, 400]);
set(gcf, 'Color', 'w');

meanResponse = mean(response',1);
imagesc( meanResponse ); box off

c = colorbar;
if (deltaVm)
    colormap('bluewhitered'); % zero set to white in this plotting
    c.Label.String = ' Delta Vm  (mV) ( response - baseline)';
else
    c.Label.String = ' Vm  (mV)';
end

set(gca, 'xtick', POSSIBLE_BAR_LOCATIONS / 2 )
set(gca, 'xticklabel', round(barPositionDegreesFromMidline))

set(gca, 'ytick', 1:length(trialNum) )

% build tick labels
tickLabels = [];
for i = 1: length (trialNum)
    tickLabels{i} = [ 'trial ' num2str( trialNum(i))  ];
end

title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ] )

% AVERAGE across all Open LOOP trials
figure('Position',[50, 50, 1000, 400]);
set(gcf, 'Color', 'w');
plot(barPositionDegreesFromMidline, meanResponse); box off
ylabel('Vm');
niceaxes
end

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

ax(1) = subplot(3,1,1);
plot( timeArray,  voltage , 'DisplayName' , 'membrane Voltage' ); hold on;
ax(2) = subplot(3,1,2);
plot( timeArray,  xPanelPos , 'DisplayName' , 'panel position' ); hold on;
ax(3) = subplot(3,1,3);
plot( timeArray, ficTracAngularPosition , 'DisplayName' , ' ball position (Yaw)' ); hold on;

linkaxes(ax,'x');
legend('show')

end