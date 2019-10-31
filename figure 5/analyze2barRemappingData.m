%% analyze2barRemappingData
% Extract data for 1bar, 2bar and 1 bar return data set for each fly and
% plot the data
%
%  Script to be used for making Figure 5 
%
% Yvette Fisher 11/13/18  
%% Critera for including data
%Open loop:
% 1) cells ave Vm for a trial must not be above -33mV (after junction
% potential correction)
% 2) include the 2 trials before CL switch, 2 trials after CL switch,
% unless there is a problem with one of those trials

% Closed loop
% 1) cells ave Vm for a trial must not be above -33mV (after junction
% potential correction)
% 2) include the 2 trials before CL switch when possible, and the last 2 trials of the 12
% min training period.   ***Flies must visit all heading locations at some
% point across both of these two trials-8min segments to make the data set

%% 1 BAR INITIAL OL Load all open loop data for initial 1 bar period
clear all
close all
ephysSettings;

% INPUT wanted trials here
possibleTrials_1barBefore = [4];


% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = false;
trialFilesList = extractTrialsWithCertainStimulusName( 'barRandLocON()', includeIfNameContainsString, possibleTrials_1barBefore);

%% Plot all selected trials to check if there are any issues and extract and calculate mean tuning
plotSelectedTrials( trialFilesList )

%% Trial section removal: input values for removal of a bad section of the trace is needed
traceRemovalOL(1).trial = 4;
traceRemovalOL(1).period =  [85,90];% seconds during the trial

[OL_oneBarTuningCurve , flyNum ] = analyzeOpenLoopTuning( trialFilesList , traceRemovalOL );

%% No trial section removal: extract trials data, and plot mean tuning 
[OL_oneBarTuningCurve , flyNum ] = analyzeOpenLoopTuning( trialFilesList );

%% 1 BAR INITIAL CLOSED LOOP for initial 1 bar period
possibleTrials_1barBefore = [6,7];

% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = true;
trialFilesListClosedLoop = extractTrialsWithCertainStimulusName( 'closedLoop_vertStripeON_270world', includeIfNameContainsString, possibleTrials_1barBefore);

%%
plotSelectedTrials( trialFilesListClosedLoop )
%% Trial section removal: input values for removal of a bad section of the trace is needed
traceRemovalCL(1).trial = 1;
traceRemovalCL(1).period =  [164, 174];% seconds during the trial

CL_before1BarTuningCurve  = analyzeClosedLoopTuning( trialFilesListClosedLoop , traceRemovalCL );

%% No trial section removal: extract trials data, and plot mean tuning 
CL_before1BarTuningCurve  = analyzeClosedLoopTuning( trialFilesListClosedLoop );


%% AFTER 2 BAR, OL Load all open loop data
% INPUT wanted trials here
possibleTrials_after2bar = [8];

% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = false;
trialFilesList = extractTrialsWithCertainStimulusName( 'barRandLocON()', includeIfNameContainsString, possibleTrials_after2bar);

%% Plot all selected trials to check if there are any issues and extract and calculate mean tuning
plotSelectedTrials( trialFilesList );

%% Trial section removal: input values for removal of a bad section of the trace
traceRemovalOL(1).trial = 12;
traceRemovalOL(1).period =  [7, 18];% seconds during the trial

[OL_after2BarTuningCurve , flyNum ] = analyzeOpenLoopTuning( trialFilesList , traceRemovalOL );
%% extract trials data, and plot mean tuning 
OL_after2BarTuningCurve = analyzeOpenLoopTuning( trialFilesList );

%% 2 BAR CLOSED LOOP period
possibleTrials_2bar = [6,7];

% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = true;
trialFilesListClosedLoop = extractTrialsWithCertainStimulusName( 'closedLoop_2bar', includeIfNameContainsString, possibleTrials_2bar);

%%
plotSelectedTrials( trialFilesListClosedLoop )

%% trial section removal:
traceRemovalCL(1).trial = 7;
traceRemovalCL(1).period =  [164, 175];% seconds during the trial


CL_2BarTuningCurve  = analyzeClosedLoopTuning( trialFilesListClosedLoop , traceRemovalCL );
%%
CL_2BarTuningCurve  = analyzeClosedLoopTuning( trialFilesListClosedLoop );

%% Plot the fly movement as a function of trial number if wanted to see this!!
% looks at fly movement for these trial
[ aveAngularSpeed,  aveForwardSpeed, trialNums ] = findAveFlyMovementByTrial( trialFilesListClosedLoop );

aveAcrossTrials_angularSpeed = mean( aveAngularSpeed );
aveAcrossTrials_forwardSpeed = mean( aveForwardSpeed );

%% Plot all 4 curves for this recording:
figure;
set(gcf, 'Color', 'w');

DEGREE_PER_LED_SLOT = 360 / 96;
POSSIBLE_BAR_LOCATIONS = 2:2:71;
MIDLINE_POSITION = 34;
barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

% CL 1 bar
subplot( 2, 2, 1)
plot( barPositionDegreesFromMidline, CL_before1BarTuningCurve , '-r' );
title(' 1 bar closed loop' ); 
niceaxes; box off
xlim([ -120 135])

% CL 2 bar
subplot(2, 2, 2)
plot( barPositionDegreesFromMidline, CL_2BarTuningCurve , '-r');
title(' 2 bar closed loop' );
niceaxes; box off
xlim([ -120 135])

% OL before 1 bar
subplot( 2, 2, 3)
plot( barPositionDegreesFromMidline, OL_oneBarTuningCurve );
title(' open loop (before)' ); 
niceaxes; box off
xlim([ -120 135])

% OL after 2 bar
subplot(2, 2, 4)
plot( barPositionDegreesFromMidline, OL_after2BarTuningCurve );
title(' open loop (after)' );
niceaxes; box off
xlim([ -120 135])

suptitle( ['flyNum' num2str( flyNum )] )

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
