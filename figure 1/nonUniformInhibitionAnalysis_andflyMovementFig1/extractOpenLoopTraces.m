% extractOpenLoopTraces
%
% Open openloop data set from specified date and trial number, 
% Extract traces and mean response Vm per trace
% Save extracted traces in one file
% save mean response Vm per bar presentation in another file
%
%
% Yvette Fisher 3/2019
%% 1) Load in open loop data set we want to work on for this fly
%clear all
close all
ephysSettings;

% INPUT wanted trials here
possibleTrials = [1];

% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = false;
trialFilesList = extractTrialsWithCertainStimulusName( 'barRandLocON()', includeIfNameContainsString, possibleTrials );

%

% 3) Extract Yaw and angular velocity variablies and voltage for each bar flash position
ephysSettings;

for fileNum = 1 : length ( trialFilesList )
    
    cd( trialFilesList( fileNum ).folder );
    % load current file for current trial
    load( trialFilesList(fileNum).name );
    
        
    LOWPASS_FILTER_CUTOFF= 25; % Hz
    THRESHOLD_ANGULAR_VELOCITY = 2500; % degrees / s  this is the max velocity that can be allowed into analysis

    % decode angular velocity and accumulated position
    [ angularVelocity , accumulatedPosition ] = ficTracSignalDecoding( data.ficTracAngularPosition, settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);

    % save TrialTime
    processedData(fileNum).trialStartTime = trialMeta.trialStartTime;
    % same Trial number
    processedData(fileNum).trialNum = trialMeta.trialNum; 
    processedData(fileNum).exptInfo = exptInfo;
    processedData(fileNum).stimulus = stimulus;
    
    minXPos = 1;
    maxXPos = 71;
    DURATION_TO_PLOT_BEFORE_FLASH = 0.5;
    DURATION_TO_PLOT_PER_EPOCH = 1; % sec
    
    % set how long the stimulus step (data.xPanelPos) needs to be constant
    % after change to be included in data set
    SEARCH_LENGTH = 0.3; % seconds
    seqFramesToLookFor = SEARCH_LENGTH * settings.sampRate;
    
    INTERVAL_TO_LOOK_BEHIND = 0.75; % seconds
    numFrameToCheckBehind = INTERVAL_TO_LOOK_BEHIND * settings.sampRate;
    
    FINAL_STIMULUS_PERIOD =  2; %seconds
    numOfFramesToIgnoreAtEnd = FINAL_STIMULUS_PERIOD * settings.sampRate;
    
    START_STIMULUS_PERIOD =  2; %seconds
    numOfFramesToIgnoreAtStart = START_STIMULUS_PERIOD * settings.sampRate;
    
    degreeBarMoved = [];
    voltageResponse = NaN( 35, 5);
    voltageDiff= NaN( 35, 5);
    
    counter = 1;
    for j = ( minXPos : maxXPos)  % loop over all possible X-pos values and pull out traces for each that has data
        % find index where the pos value has changed
        change = find( diff (data.xPanelPos) ~= 0 ) + 1;
        
        % find which of those indexes the panels was at the currect position j
        indofJ = change (data.xPanelPos (change)  == j) ;
        
        % remove any indexs too close to the end of the trace to cause errors
        indofJ  = indofJ (indofJ < ( numel( data.xPanelPos )  - numOfFramesToIgnoreAtEnd));
        
        % remove any indexes to close to start of the trace that might
        % cause errors
        indofJ  = indofJ (indofJ >  numOfFramesToIgnoreAtStart) ;
        
        
        epochStartInds = [];
        
        % only take those where there are not change in position values for 300 ms after the start
        for ii = 1: length( indofJ )
            averagePosDecode = mean ( data.xPanelPos (indofJ(ii) + 1 : indofJ(ii) + seqFramesToLookFor ) );
            
            if(averagePosDecode == j) % if not other changes then take this index
                epochStartInds = [epochStartInds indofJ(ii) ];
            end
        end
        

        
        % check that this epoch step was not empty for some reason
        if( ~ isempty(epochStartInds) )
            
            stimulusStep= [];
            currentVoltageByPosition = [];
            currentBallPositionByBarPosition = [];
            currentBallVelocityByBarPosition = [];
            currBarPosition = [];
            lastBarPosition = [];
            numOfBarPositionJumped = [];
            degreeBarJumped = [];
            lastBarPositionShown = []; % prespecified
            currBarPositionShown = []; % prespecified

            
            for i = 1 : ( numel(epochStartInds) )
                
                currEpochStart = epochStartInds(i) - (DURATION_TO_PLOT_BEFORE_FLASH * settings.sampRate);
                currEpochEnd = epochStartInds(i) + ( DURATION_TO_PLOT_PER_EPOCH * settings.sampRate);
                
                % store voltage and ball position and ball velocity trace
                currVoltage = data.scaledVoltage( currEpochStart : currEpochEnd);
                currentVoltageByPosition(i, :) = currVoltage; % save trace
                currTimeArray = (1 : length(  currVoltage ) ) / settings.sampRate;
                
                
                % Extract average from each epoch period
                PRE_BASELINE_START_TIME = 0.25;
                PRE_BASELINE_END_TIME = 0.5;
                
                RESP_START_TIME = 0.75;
                RESP_END_TIME = 1;
                
                % save average values of Vm for summary plot
                preBaselinesIndexes = PRE_BASELINE_START_TIME < currTimeArray & currTimeArray < PRE_BASELINE_END_TIME ;
                responsePeriodIndexes = RESP_START_TIME < currTimeArray & currTimeArray < RESP_END_TIME ;
                
                % save pre baseline mean
                preBaselineVoltage = mean ( currVoltage(preBaselinesIndexes, :) , 1 );
                responseVoltage = mean( currVoltage (responsePeriodIndexes, :), 1);
                
                % save delta Vm and raw Vm response
                voltageResponse(counter, i) = responseVoltage;
                voltageDiff(counter, i) = responseVoltage - preBaselineVoltage;
                
                
                currentBallPositionByBarPosition(i, :) = accumulatedPosition( currEpochStart : currEpochEnd);
                currentBallVelocityByBarPosition(i , :) = angularVelocity( currEpochStart : currEpochEnd);
                
                % same stimulus step
                stimulusStep(i, :)  = data.xPanelPos( currEpochStart: currEpochEnd );
                
                % solve for how far the bar jumped to get to this current location
                DEGREE_PER_LED_SLOT = 360 / 96;  %
                
                currBarPosition = data.xPanelPos( epochStartInds(i) );
                lastBarPosition = data.xPanelPos( epochStartInds(i) -  numFrameToCheckBehind );
                
                numOfBarPositionJumped(i) =  currBarPosition - lastBarPosition;% last - current
                
                % turn jump into degree
                degreeJumped = numOfBarPositionJumped(i) * DEGREE_PER_LED_SLOT;
                % if larger jump than 180 wrap other direction
                if( abs(degreeJumped) >180 && degreeJumped > 0 )
                    degreeJumped = degreeJumped - 360;
                    
                elseif( abs(degreeJumped) >180 && degreeJumped < 0 )
                    degreeJumped = degreeJumped + 360;
                end
                
                % save degree jumped value
                degreeBarJumped(i) = degreeJumped;
                lastBarPositionShown(i) = lastBarPosition;
                currBarPositionShown(i) = currBarPosition;
                
            end
             % save the data traces for this bar position:
            voltageByPosition{counter} = currentVoltageByPosition;
            
            ballPositionByBarPosition{counter} = currentBallPositionByBarPosition; 
            ballVelocityByBarPosition{counter} = currentBallVelocityByBarPosition;
            stimulusStepByPosition{counter} = stimulusStep;
            % 
            currentBarPosition{counter} = currBarPositionShown;
            if length( unique(currBarPositionShown) ) > 1
                error(' error bar flashes from different positions are grouped together')
            end
            currentBarPositionForBatch(counter) = unique(currBarPositionShown);
            
            previousBarPosition{counter} = lastBarPositionShown;
            barPositionJumped{counter} = numOfBarPositionJumped;
            degreeBarMoved{counter} = degreeBarJumped;
            
            counter = counter + 1;
        end
          
        
        
    end
    % same traces & stimulus steps for all bar positions
    processedData(fileNum).voltageByPosition =   voltageByPosition;
    
    % add in delta Vm and raw Vm response calculations
    processedData(fileNum).meanVoltageResp = voltageResponse;
    processedData(fileNum).meanVoltageDiff = voltageDiff;
    
    processedData(fileNum).ballPositionByBarPosition =   ballPositionByBarPosition; % Yaw
    processedData(fileNum).ballVelocityByBarPosition =   ballVelocityByBarPosition; % YAW
    processedData(fileNum).stimulusStep = stimulusStepByPosition;
    
    processedData(fileNum).previousBarPosition = previousBarPosition;
    processedData(fileNum).currentBarPositionByTrial = currentBarPosition;
    processedData(fileNum).currentBarPosition = currentBarPositionForBatch;
    processedData(fileNum).numBarPositionJumped =  barPositionJumped;
    processedData(fileNum).degreeBarJumped =  degreeBarMoved;
end

% save processed data into open loop folder
parentDir = '/Users/evettita/Dropbox (HMS)/EphysData/EP-G_recordings/openLoopMetaData/';
fileName = [processedData(1).exptInfo.dNum ...
    '_fly' num2str(processedData(1).exptInfo.flyNum)...
    '_cell' num2str( processedData(1).exptInfo.cellNum )...
    '_expt' num2str( processedData(1).exptInfo.cellExpNum ) '.mat' ];


fullFileName = [ parentDir fileName ];
save( fullFileName , 'processedData' , 'trialFilesList');

processedData;



