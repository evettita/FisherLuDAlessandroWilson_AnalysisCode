% analyzeVmResponseByPreviousBarLocation

% Yvette Fisher 3/5/19
% %% 1) Load in open loop data set we want to work on for this fly
% clear all
% close all
% ephysSettings;
% 
% % INPUT wanted trials here
% possibleTrials = [1:12];
% 
% % pull out file names for the trials where the wanted stimulus was shown:
% includeIfNameContainsString = false;
% trialFilesList = extractTrialsWithCertainStimulusName( 'barRandLocON()', includeIfNameContainsString, possibleTrials );
% % 
% % %% 2) Plot all selected trials to check how the raw traces look by eye
% % plotSelectedTrials( trialFilesList )
% % 
% % 
% % 3) Extract Yaw and angular velocity variablies/voltage for each bar flash position
% ephysSettings;
% 
% for fileNum = 1 : length ( trialFilesList )
%     
%     cd( trialFilesList( fileNum ).folder );
%     % load current file for current trial
%     load( trialFilesList(fileNum).name );
%     
%         
%     LOWPASS_FILTER_CUTOFF= 25; % Hz
%     THRESHOLD_ANGULAR_VELOCITY = 2500; % degrees / s  this is the max velocity that can be allowed into analysis
% 
%     % decode angular velocity and accumulated position
%     [ angularVelocity , accumulatedPosition ] = ficTracSignalDecoding( data.ficTracAngularPosition, settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);
% 
%     % save TrialTime
%     processedData(fileNum).trialStartTime = trialMeta.trialStartTime;
%     % same Trial number
%     processedData(fileNum).trialNum = trialMeta.trialNum; 
%     processedData(fileNum).exptInfo = exptInfo;
%     processedData(fileNum).stimulus = stimulus;
%     
%     minXPos = 1;
%     maxXPos = 71;
%     DURATION_TO_PLOT_BEFORE_FLASH = 0.5;
%     DURATION_TO_PLOT_PER_EPOCH = 2.5; % sec
%     
%     % set how long the stimulus step (data.xPanelPos) needs to be constant
%     % after change to be included in data set
%     SEARCH_LENGTH = 0.3; % seconds
%     seqFramesToLookFor = SEARCH_LENGTH * settings.sampRate;
%     
%     INTERVAL_TO_LOOK_BEHIND = 0.75; % seconds
%     numFrameToCheckBehind = INTERVAL_TO_LOOK_BEHIND * settings.sampRate;
%     
%     FINAL_STIMULUS_PERIOD =  4; %seconds
%     numOfFramesToIgnoreAtEnd = FINAL_STIMULUS_PERIOD * settings.sampRate;
%     
%     START_STIMULUS_PERIOD =  4; %seconds
%     numOfFramesToIgnoreAtStart = START_STIMULUS_PERIOD * settings.sampRate;
%     
%     degreeBarMoved = [];
%     
%     for j = ( minXPos : maxXPos)  % loop over all possible X-pos values and pull out traces for each that has data
%         % find index where the pos value has changed
%         change = find( diff (data.xPanelPos) ~= 0 ) + 1;
%         
%         % find which of those indexes the panels was at the currect position j
%         indofJ = change (data.xPanelPos (change)  == j) ;
%         
%         % remove any indexs too close to the end of the trace to cause errors
%         indofJ  = indofJ (indofJ < ( numel( data.xPanelPos )  - numOfFramesToIgnoreAtEnd));
%         
%         % remove any indexes to close to start of the trace that might
%         % cause errors
%         indofJ  = indofJ (indofJ >  numOfFramesToIgnoreAtStart) ;
%         
%         
%         epochStartInds = [];
%         
%         % only take those where there are not change in position values for 300 ms after the start
%         for ii = 1: length( indofJ )
%             averagePosDecode = mean ( data.xPanelPos (indofJ(ii) + 1 : indofJ(ii) + seqFramesToLookFor ) );
%             
%             if(averagePosDecode == j) % if not other changes then take this index
%                 epochStartInds = [epochStartInds indofJ(ii) ];
%             end
%         end
%         
% 
%         
%         % check that this epoch step was not empty for some reason
%         if( ~ isempty(epochStartInds) )
%             
%             stimulusStep= [];
%             currentVoltageByPosition = [];
%             currentBallPositionByBarPosition = [];
%             currentBallVelocityByBarPosition = [];
%             currBarPosition = [];
%             lastBarPosition = [];
%             numOfBarPositionJumped = [];
%             degreeBarJumped = [];
%             lastBarPositionShown =[];
%             
%             for i = 1 : ( numel(epochStartInds) )
%                 
%                 currEpochStart = epochStartInds(i) - (DURATION_TO_PLOT_BEFORE_FLASH * settings.sampRate);
%                 currEpochEnd = epochStartInds(i) + ( DURATION_TO_PLOT_PER_EPOCH * settings.sampRate);
%                 
%                 % store voltage and ball position and ball velocity trace
%                 currentVoltageByPosition(i, :) = data.scaledVoltage( currEpochStart : currEpochEnd);
%                 currentBallPositionByBarPosition(i, :) = accumulatedPosition( currEpochStart : currEpochEnd);
%                 currentBallVelocityByBarPosition(i , :) = angularVelocity( currEpochStart : currEpochEnd);
%                 
%                 % same stimulus step
%                 stimulusStep(i, :)  = data.xPanelPos( currEpochStart: currEpochEnd );
%                 
%                 % solve for how far the bar jumped to get to this current location
%                 DEGREE_PER_LED_SLOT = 360 / 96;  %
%                 
%                 currBarPosition = data.xPanelPos( epochStartInds(i) );
%                 lastBarPosition = data.xPanelPos( epochStartInds(i) -  numFrameToCheckBehind );
%                 
%                 numOfBarPositionJumped(i) =  currBarPosition - lastBarPosition;% last - current
%                 
%                 % turn jump into degree
%                 degreeJumped = numOfBarPositionJumped(i) * DEGREE_PER_LED_SLOT;
%                 % if larger jump than 180 wrap other direction
%                 if( abs(degreeJumped) >180 && degreeJumped > 0 )
%                     degreeJumped = degreeJumped - 360;
%                     
%                 elseif( abs(degreeJumped) >180 && degreeJumped < 0 )
%                     degreeJumped = degreeJumped + 360;
%                 end
%                 
%                 % save degree jumped value
%                 degreeBarJumped(i) = degreeJumped;
%                 lastBarPositionShown(i) = lastBarPosition;
%                 
%             end
%              % save the data traces for this bar position:
%             voltageByPosition{j} = currentVoltageByPosition;
%             ballPositionByBarPosition{j} = currentBallPositionByBarPosition; 
%             ballVelocityByBarPosition{j} = currentBallVelocityByBarPosition;
%             stimulusStepByPosition{j} = stimulusStep;
%             % 
%             previousBarPosition(j, 1:length( lastBarPositionShown ) ) = lastBarPositionShown;
%             barPositionJumped(j, 1:length(numOfBarPositionJumped) ) = numOfBarPositionJumped;
%             degreeBarMoved(j, 1:length(degreeBarJumped) ) = degreeBarJumped;
%         end
%              
%     end
%     % same traces & stimulus steps for all bar positions
%     processedData(fileNum).voltageByPosition =   voltageByPosition;
%     processedData(fileNum).ballPositionByBarPosition =   ballPositionByBarPosition; % Yaw
%     processedData(fileNum).ballVelocityByBarPosition =   ballVelocityByBarPosition; % YAW
%     processedData(fileNum).stimulusStep = stimulusStepByPosition;
%     processedData(fileNum).previousBarPosition = previousBarPosition;
%     processedData(fileNum).numBarPositionJumped =  barPositionJumped;
%     processedData(fileNum).degreeBarJumped =  degreeBarMoved;
% end

%% Plot raw voltage trace and raw bar movmements alinged for the same trials--
ephysSettings; 
close all;
% this only works for a single trial worth of data, choice which here
trialIndexToPlot = 2;

DURATION_TO_PLOT_BEFORE_FLASH = 0.5;
DURATION_TO_BAR_PRESENTATION = 0.5;
MIDLINE_POSITION = 34;

% possible X-Positions 2 and 55
minXPos = 1;
maxXPos = 70;

for j = (minXPos : maxXPos)
    
    % only plot if data exists for that bar position
    if( ~ isempty( processedData(trialIndexToPlot).voltageByPosition{j} ) )
        
        figure('Position',[50, 50, 400, 400]);
        set(gcf, 'Color', 'w');
        
        %subplot(3, 1, 1);
        % shade when bar was on the screen for visuallization
        startInd = DURATION_TO_PLOT_BEFORE_FLASH;
        endInd =  DURATION_TO_PLOT_BEFORE_FLASH + DURATION_TO_BAR_PRESENTATION;
        xcord = [startInd endInd endInd startInd];
        
        ytop = 0 ; %
        TYPICAL_Vm = -57;
        ybottom = ytop + TYPICAL_Vm ;
        ycord = [ybottom ybottom ytop ytop];
        patch( xcord, ycord ,'b', 'FaceAlpha',.15, 'LineStyle', 'none'); hold on;
        
        % shade second bar presentation location on the screen:
        startInd = 1.5;
        endInd = 2;
        xcord = [startInd endInd endInd startInd];
        
        ytop = 0 ; %
        TYPICAL_Vm = -57;
        ybottom = ytop + TYPICAL_Vm ;
        ycord = [ybottom ybottom ytop ytop];
        patch( xcord, ycord ,'b', 'FaceAlpha',.15, 'LineStyle', 'none'); hold on;
        
        
        barPositionDegrees = ( j - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;
        
        title ( [ 'bar location: ' num2str( barPositionDegrees , 3 ) ' deg, pos = ' num2str(j) ] );
        
        % membrane votlage
        timeArray = (1 : length(processedData(trialIndexToPlot).voltageByPosition{j}(1,:) ) ) / settings.sampRate; % seconds
        plot( timeArray,  processedData(trialIndexToPlot).voltageByPosition{j}' ); hold on
        ylabel(' mV ')
        
        ylim( [ -50, -25 ] );
        box off; 
        
        
%         subplot(3, 1, 2)
%         
%         % yaw position
%         plot( timeArray, processedData(trialIndexToPlot).ballPositionByBarPosition{j}' )
%         ylabel(' yaw position (deg) ');
%         box off; 
%         
%         subplot(3, 1, 3)
%         % yaw angular velocity
%         plot( timeArray,  processedData(trialIndexToPlot).ballVelocityByBarPosition{j}' )
%         ylabel(' yaw angular velocity (deg/s) ')
%         box off;
    end
end

