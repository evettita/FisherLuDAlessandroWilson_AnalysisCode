% plotBarRandLoc
% Builds a receptive field plot using the bar @ random location stimulus
% (barRandLoc)
% 
% Yvette Fisher 11/2018
%% Plot data to look at from aquireTrial code
ephysSettings
FigHand = figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');

timeArray = (1  :  length(data.current) ) / settings.sampRate; % seconds

% check if scaled current exists
if(isfield(data,'scaledCurrent'))
    current = data.scaledCurrent; % for voltage clamp traces
else
    current = data.current;
end

ax(1) = subplot(3,1,1);
% plot current trace
plot(timeArray, current); hold on;

low_pass_cutoff= 1e3; % 1 kHz
filteredCurrent = lowPassFilter( current, low_pass_cutoff, settings.sampRate );
hold on; plot(timeArray,filteredCurrent,'k')

if( isfield( stimulus, 'visualTriggerCommand' ) )
    plot(timeArray, stimulus.visualTriggerCommand); hold on;
end

% FIND PEAKS using signal processing function
 MIN_PEAK_HEIGHT = 13; % pA
 WIDTHS_REQUIRED = 10; 
 
% spikes detection for Ring neuron recordings
[~ , spikeIndex] = findpeaks((-1)*filteredCurrent, 'MinPeakHeight', MIN_PEAK_HEIGHT, 'MinPeakWidth',WIDTHS_REQUIRED);

hold on;
scatter( timeArray(spikeIndex) , filteredCurrent(spikeIndex) );

% converts spike times from sample rate into seconds...
spikeTimes = timeArray(spikeIndex);

% solve for an array interspike interval between each spike
interSpikeIntervals = diff(spikeTimes);

title('Current & stimulus Command time course');
xlabel('time(s)')
ylabel('pA');

% check if scaled voltage exists
if(isfield(data,'scaledVoltage'))
    voltage = data.scaledVoltage; % for current clamp traces
else
    voltage = data.voltage;
end

% plot voltage trace
ax(2) = subplot(3,1,2);
plot( timeArray, voltage); hold on;
scatter( timeArray(spikeIndex) , voltage(spikeIndex) );
title('voltage');
xlabel('time(s)')
ylabel('mV');

% plot panel data if it was aquired
if(isfield( data, 'xPanelPos') )
    % plot x and y panel pos
    ax(3) = subplot(3,1,3);
    plot( timeArray, data.xPanelPos, 'DisplayName','panel x'); hold on;
    plot(timeArray, data.yPanelPos, 'DisplayName','panel y');
    title('voltage');
    xlabel('time(s)')
    ylabel('Panel X and Y values (V)');
    legend('show')
end

linkaxes(ax,'x');
suptitle( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ' trial#: ' num2str(trialMeta.trialNum) ' stim: ' num2str(stimulus.name) ] )

%% raw VOLTAGE or CURRENT PLOT  Extract traces from each different X-pos 

% 180 degree panels 
 DEGREE_PER_LED_SLOT = 360 / 96;  
 MIDLINE_POSITION = 20; % approx. LED position where the fly is aligned to for all 180 deg EPG recordings 8/2017 - 11/2017 

% check if scaled voltage exists
if(isfield(data,'scaledVoltage'))
    dataTrace = data.scaledVoltage; % for current clamp traces
    OFFSET_FOR_PLOT = 0; % pA but this is abitrary
else
    dataTrace = data.scaledCurrent;
    OFFSET_FOR_PLOT = 0; % mV but this is abitrary
end

OFFSET_FOR_aveVm = 0;

% possible X-Positions 2 and 55
minXPos = 1;
maxXPos = 55;

close all; % close figures

DURATION_TO_PLOT_BEFORE_FLASH = 0.5;
DURATION_TO_PLOT_PER_EPOCH = 1; % sec
DURATION_TO_BAR_PRESENTATION = 0.5; % sec

% total epoch is now 0.5 + 1 (0.5 nothing, 0.5 bar, 0.5 nothing)
BASELINE_START_TIME = DURATION_TO_PLOT_BEFORE_FLASH + DURATION_TO_BAR_PRESENTATION + 0.25; % sec
BASELINE_END_TIME = BASELINE_START_TIME + 0.25; % sec
DURATION_OF_BASELINE = 0.25; % sec

% timing for pre baseline period
PRE_BASELINE_START_TIME = DURATION_TO_PLOT_BEFORE_FLASH - DURATION_OF_BASELINE;
PRE_BASELINE_END_TIME = DURATION_TO_PLOT_BEFORE_FLASH; 

SEARCH_LENGTH = 0.3; % seconds
seqFramesToLookFor = SEARCH_LENGTH * settings.sampRate;

% intialize storage variables
aveSpikeCount = [];
semSpikeCount = [];
aveVoltage = [];
semVoltage = [];
allPosWithData = [];

aveBaselineSpikeCount = [];
semBaselineSpikeCount = [];

aveBaselineVoltage = [];
semBaselineVoltage = [];

avePreBaselineVoltage = [];
semPreBaselineVoltage = [];

aveVmChange =[];
semVmChange = [];

spikeCount = [];
preBaselineSpikeCount = [];
counter = 1;


for j = (minXPos : maxXPos)  % loop over all possible X-pos values and pull out traces for each that has data
    
    % find index where the pos value has changed
    change = find( diff (data.xPanelPos) ~= 0 ) + 1;
    
    % find which of those indexes the panels was at the currect position j
    indofJ = change (data.xPanelPos (change)  == j) ;
    
    % remove any indexs too close to the end of the trace to cause errors
    % below
    indofJ  = indofJ (indofJ < ( numel( data.xPanelPos )  - seqFramesToLookFor));
    
    epochStartInds = [];
    
    % only take those where there are not change in position values for 100 ms after the start
    for i = 1: length( indofJ)
        averagePosDecode = mean ( data.xPanelPos (indofJ(i) + 1 : indofJ(i) + seqFramesToLookFor ) );
        
        if(averagePosDecode == j) % if not other changes then take this index
            epochStartInds = [epochStartInds indofJ(i) ];
        end
        
    end
    
    % check that this epoch step was not empty for some reason
    if( ~ isempty(epochStartInds) )
        
        spikeCount = [];
        voltageAverage =[];
        baselineVoltageAverage = [];
        preBaselineVoltageAverage =[];
        currPositionVoltage = [];
        
        
        spikeFrameTimesOffset = {};
        allSpikeFrameTimesOffset = [];
        figure();
        set(gcf, 'Color', 'w');
        
        SECONDS_TO_IGNORE = 4;
        END_OF_TRACE_FOR_ANALYSIS = length( timeArray ) - ( SECONDS_TO_IGNORE * settings.sampRate ) ;
        
        for i = 1 : ( sum( epochStartInds < END_OF_TRACE_FOR_ANALYSIS ) ) % make sure does run over end of trace
            
            
            currEpochStart = epochStartInds(i) - (DURATION_TO_PLOT_BEFORE_FLASH * settings.sampRate);
            currEpochEnd = epochStartInds(i) + ( DURATION_TO_PLOT_PER_EPOCH * settings.sampRate);
            
            currBarDisplayStart = epochStartInds(i);
            currBarDisplayEnds = epochStartInds(i) + ( DURATION_TO_BAR_PRESENTATION * settings.sampRate);
            
            currTimeArray = timeArray( currEpochStart : currEpochEnd) ;
            currTimeArray = currTimeArray - currTimeArray(1); % offset to start at t = 0 for plot;
            
            currVoltage = dataTrace( currEpochStart : currEpochEnd);
            currVoltageOffset = currVoltage - (OFFSET_FOR_PLOT * (i - 1));% offset to see all the traces over eachother
            % Plot data trace
            plot(currTimeArray, currVoltageOffset); hold on;
            
            stimulusStep = data.xPanelPos( currEpochStart: currEpochEnd);
            
            % plot stimulus array to double check
            plot(currTimeArray, stimulusStep);
            
            ORDER_MEDFILTER = 800; % parameter than MM uses
            
            currVoltageFiltered = medfilt1(currVoltage, ORDER_MEDFILTER, 'truncate'); % Median filtering of the trace
            
            epochIndex = (DURATION_TO_PLOT_BEFORE_FLASH + DURATION_TO_BAR_PRESENTATION) > currTimeArray & currTimeArray > DURATION_TO_PLOT_BEFORE_FLASH ;
            baselineIndex = BASELINE_END_TIME > currTimeArray & currTimeArray > BASELINE_START_TIME ;
            preBaselineIndex = PRE_BASELINE_END_TIME > currTimeArray & currTimeArray > PRE_BASELINE_START_TIME ;
            
            currPositionVoltage(i,:) = currVoltageFiltered;
            
            voltageAverage(i) = mean( currVoltageFiltered (epochIndex) );
            baselineVoltageAverage(i) = mean( currVoltageFiltered (baselineIndex) );
            preBaselineVoltageAverage(i) = mean( currVoltageFiltered (preBaselineIndex) );
            
            % save spike count during bar:
            spikesInBarEpoch =  spikeIndex( currBarDisplayStart <= spikeIndex  & spikeIndex <= currBarDisplayEnds);
            spikeCount(i) = numel( spikesInBarEpoch );
            
            
            %save spike count before bar presentation from 250ms before -
            %to imidiately before           
            spikesInPreBaselinePeriod =  spikeIndex( (currEpochStart + PRE_BASELINE_START_TIME * settings.sampRate) <= spikeIndex  & spikeIndex <= ( currEpochStart + PRE_BASELINE_END_TIME * settings.sampRate) );
            preBaselineSpikeCount (i) = numel( spikesInPreBaselinePeriod );
            
            % save spike count for this trial to build PSTH
            spikeInWholeEpoch = spikeIndex( currEpochStart <= spikeIndex  & spikeIndex <= currEpochEnd);
            
            spikeFrameTimesOffset{ i } = spikeInWholeEpoch - currEpochStart;
            allSpikeFrameTimesOffset = [ allSpikeFrameTimesOffset; ( spikeInWholeEpoch - currEpochStart ) ];
        end
        
        % find average spike count for this stimulus
        aveSpikeCount(counter) = mean( spikeCount(:) ) * (1 / DURATION_TO_BAR_PRESENTATION) ; % average of all of these trials
        semSpikeCount(counter) = std( spikeCount(:) * (1 / DURATION_TO_BAR_PRESENTATION) ) / sqrt( numel( spikeCount(:)));
        allPosWithData (counter) = j;
        
        aveBaselineSpikeCount(counter) =  mean( preBaselineSpikeCount(:) ) * (1 / DURATION_OF_BASELINE) ; % average of all of these trials
        semBaselineSpikeCount(counter) = std( preBaselineSpikeCount(:) * (1 / DURATION_OF_BASELINE) ) / sqrt( numel( spikeCount(:)));
        
        % calc and store mean and sem of membrane voltage during bar presentation
        aveVoltage(counter) = mean( voltageAverage(:) ); % average of all of these trials;
        semVoltage(counter) = std( voltageAverage(:) ) / sqrt( numel( voltageAverage(:)));
        
        % calc and store mean and sem of membrane voltage in baseline
        % period
        aveBaselineVoltage(counter) = mean( baselineVoltageAverage(:) ) ; % average of all of these trials;
        semBaselineVoltage(counter) = std( baselineVoltageAverage(:)) / sqrt( numel( baselineVoltageAverage(:)));
        
        
        % calc and store mean and sem of membrane voltage in "pre-baseline"
        % period
        avePreBaselineVoltage(counter) = mean( preBaselineVoltageAverage(:) ) ; % average of all of these trials;
        semPreBaselineVoltage(counter) = std( preBaselineVoltageAverage(:)) / sqrt( numel( preBaselineVoltageAverage(:)));
        
        
        % calc and store mean and sem of delta Vm response - prebaseline
        VmChange =  voltageAverage(:) - preBaselineVoltageAverage(:);
        aveVmChange(counter) = mean( VmChange(:) ) ; % average of all of these trials;
        semVmChange(counter) = std( VmChange(:)) / sqrt( numel( VmChange(:)));
        
        
        % plot mean Vm filtered trace, offset upward for viewing
        aveCurrPosTrace = mean( currPositionVoltage(:,:) );
        plot( currTimeArray,  aveCurrPosTrace + OFFSET_FOR_aveVm, '-k', 'lineWidth', 1);
        
        %PLOT PSTH above the raw data using spike times
        spikeTimesSeconds = allSpikeFrameTimesOffset / (settings.sampRate);
        hold on;
        binWidthSec = .030; % sec
        binsEdges = 0 : binWidthSec : ( DURATION_TO_PLOT_BEFORE_FLASH + DURATION_TO_PLOT_PER_EPOCH );
        binsCenters = binsEdges(1 : end - 1) + ( binWidthSec / 2);
        %histogram ( spikeTimesSeconds , binsEdges, 'Normalization' , 'countdensity'  );
        n = histcounts ( spikeTimesSeconds , binsEdges  );
        
        % normalized by number of trials in this count and by bin width
        % (seconds)
        spikesPerSecond = n / ( numel(epochStartInds) * binWidthSec );
        bar(binsCenters,  spikesPerSecond, 'BarWidth', 1 , 'EdgeColor', 'none', 'FaceColor', [ 0.1 0.1 1]);
        
        % shade when bar was on the screen for visuallization
        startInd = DURATION_TO_PLOT_BEFORE_FLASH;
        endInd =  DURATION_TO_PLOT_BEFORE_FLASH + DURATION_TO_BAR_PRESENTATION;
        xcord = [startInd endInd endInd startInd];
        
        ytop = 0 ; %
        TYPICAL_Vm = -57;
        ybottom = ytop - (OFFSET_FOR_PLOT * (numel(epochStartInds)) - TYPICAL_Vm ) ;
        ycord = [ybottom ybottom ytop ytop];
        patch( xcord, ycord ,'b', 'FaceAlpha',.15, 'LineStyle', 'none');
        
        ylim( [ -50, -20 ] );
        barPositionDegrees = ( j - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;
                
        title ( [ 'bar location: ' num2str( barPositionDegrees , 3 ) ' deg, pos = ' num2str(j) ] );
        ylabel('mV  or  s/sec')
        xlabel('s')
        box off;
        
        counter  = counter + 1;
    end
end

%% Plot Receptive field summary
barPositionDegreesFromMidline = ( allPosWithData - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

figure;
set(gcf, 'Color', 'w');

errorbar(barPositionDegreesFromMidline, aveSpikeCount , semSpikeCount, 'LineWidth', 1.5, 'Color',[0.9,0.2,0]); hold on;
errorbar(barPositionDegreesFromMidline, aveBaselineSpikeCount , semBaselineSpikeCount, '-k', 'LineWidth', 1.5);

xlabel(' deg from midline ')
ylabel('spikes / sec')
xlim([-90 140]);
title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ' trial#: ' num2str(trialMeta.trialNum) ' stim: ' num2str(stimulus.name) ] )
box off;
hh = legend( ' response', 'baseline (250ms before)');
set(hh,'box','off')
niceaxes

% Plot of RF using ave Vm
barPositionDegreesFromMidline = ( allPosWithData - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

figure;
set(gcf, 'Color', 'w');
errorbar(barPositionDegreesFromMidline, aveVoltage , semVoltage, 'LineWidth', 1.5, 'Color',[0,0.7,0.9]); hold on;
errorbar(barPositionDegreesFromMidline, avePreBaselineVoltage , semPreBaselineVoltage, '-k', 'LineWidth', 1.5);

xlabel(' deg from midline ')
ylabel('mV')
xlim([-100 160]);
title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ' trial#: ' num2str(trialMeta.trialNum) ' stim: ' num2str(stimulus.name) ] )
box off;
hh = legend( ' response', 'baseline (250 before)');
set(hh,'box','off')
niceaxes


figure;
set(gcf, 'Color', 'w');
errorbar(barPositionDegreesFromMidline, aveVmChange, semVmChange , 'LineWidth', 1.5, 'Color',[0,0.3,0.9]); hold on;

xlabel(' deg from midline ')
ylabel('mV')
xlim([-100 160]);
title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ' trial#: ' num2str(trialMeta.trialNum) ' stim: ' num2str(stimulus.name) ] )
box off;
hh = legend( ' Delta Vm  (response - baseline)');
set(hh,'box','off')
niceaxes
