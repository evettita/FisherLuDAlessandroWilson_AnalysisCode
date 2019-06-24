%% QuickPlotting_flyOnTheBall
% plot fly on the ball data with visual stimuli/ or not....
% 
% Yvette Fisher 12/2017


% other analysis to do:  
% 1) Vm and/or spiking as a function of bar position
% 2) Analysis of ball angular velocity and how it correlates with Vm....
% 3)


ephysSettings

close all;
FigHand = figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');

timeArray = (1  :  length(data.current) ) / settings.sampRate; % seconds

% check if scaled voltage exists
% voltage = data.voltage   or voltage = data.scaleVoltage
if(isfield(data,'scaledVoltage'))
    voltage = data.scaledVoltage; % for current clamp traces
else
    voltage = data.voltage;
end
% plot voltage trace
ax(1) = subplot(4,1,1);
plot( timeArray, voltage, 'k'); hold on;
ylabel('mV');
xlim([0 timeArray(end) ])
box off

% plot panel data if it was aquired
if(isfield( data, 'xPanelPos') )
    % plot x and y panel pos
    ax(2) = subplot(4,1,2);

    % Ajust into degrees 230 degree panels 
    TOTAL_DEGREES_IN_PATTERN = 360;
    DEGREE_PER_LED_SLOT = 360 / 96;  % fixed from YEF math error on 1/3/18
    MIDLINE_POSITION = 34;
    MAXBARPOS = 72;%96;
    
    patternOffsetPosition = decodePatternOffset(  stimulus.panelParams.patternNum );
    
    xBarPostion =  data.xPanelPos + patternOffsetPosition;
    xBarPostion = mod( xBarPostion , MAXBARPOS );
    
    barPositionDegreesFromMidline = (  xBarPostion - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;
    
    EDGE_OF_SCREEN_POSITION = 72; % last LED slot in postion units
    EDGE_OF_SCREEN_DEG = (EDGE_OF_SCREEN_POSITION -  MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

    plot( timeArray , EDGE_OF_SCREEN_DEG * ones(1, length( timeArray )) , '-k' ,'LineStyle', '--', 'DisplayName' , 'right edge of screen '  ); hold on;
    
    plot( timeArray, barPositionDegreesFromMidline, 'DisplayName' , 'bar position (deg)' ); hold on;
    %plot( timeArray, data.xPanelPos , 'DisplayName', 'bar position (LED position)' ); hold on;
    ylabel('bar center (deg) ');
    xlim([0 timeArray(end) ])
    ylim([0 - ( MIDLINE_POSITION * DEGREE_PER_LED_SLOT),  TOTAL_DEGREES_IN_PATTERN - ( MIDLINE_POSITION * DEGREE_PER_LED_SLOT)]); %ranges from 0 -10 Volts
    
    %TICK_STEPS_DEG = 30;
    %yticks([ (0 - ( MIDLINE_POSITION * DEGREE_PER_LED_SLOT)) : TICK_STEPS_DEG : TOTAL_DEGREES_IN_PATTERN - ( MIDLINE_POSITION * DEGREE_PER_LED_SLOT) ])
    %legend('show')
    
    
    box off;
end


 % plot yaw postion of the bar over time so show how the fly was
 % walking/turning
    ficTracAngularPos = data.ficTracAngularPosition;
    ficTracIntx = data.ficTracIntx;
    ficTracInty = data.ficTracInty;
    
    % plot position data from fictrac (0 - 10 Volts)
    ax(3) = subplot(4,1,3);
    
    plot( timeArray, ficTracAngularPos , 'DisplayName' , 'heading' ); hold on;
    plot( timeArray, ficTracIntx , 'DisplayName' , 'X' ); hold on;
    plot( timeArray, ficTracInty , 'DisplayName' , 'Y' ); hold on;
    title('ficTrac position signal (10V = 1 ball revolution)');
    xlabel('time(s)')
    ylabel('V ');
    xlim([0 timeArray(end) ])
    ylim([0 10]); %ranges from 0 -10 Volts
    
    legend('show')
    box off
    
    % plot ball forward and angular velocity

    
    LOWPASS_FILTER_CUTOFF= 25; % Hz
    THRESHOLD_ANGULAR_VELOCITY = 2500; % degrees / s  this is the max velocity that can be allowed into analysis
    THRESHOLD_FORWARD_VELOCITY = 2500; % degrees / s  this is the max velocity that can be allowed into analysis
    [ angularVelocity , ~ ] = ficTracSignalDecoding( data.ficTracAngularPosition, settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);
        % decode forward velocity and accumulated X position
    [ forwardVelocity , ~ ] = ficTracSignalDecoding( data.ficTracIntx , settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_FORWARD_VELOCITY);
    
    ax(4) = subplot(4,1,4);
    plot( timeArray, angularVelocity, 'DisplayName' , 'angular velocity' ); hold on;
    plot( timeArray, forwardVelocity, 'DisplayName' , 'forward velocity' ); hold on;
    legend('show')
    box off

linkaxes(ax,'x');
suptitle( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ' trial#: ' num2str(trialMeta.trialNum) ' stim: ' num2str(stimulus.name) ] )

%%
    % save current plot
    dir =  'F:\Dropbox (HMS)\FisherLuWilson ms\figure 2\';
    fileName = [ dir 'ExampleCLTrace_flyNum_179.eps' ];
    print( fileName, '-dpdf');


%% Plot Vm as a function of closed loop ball position
binSizeForBarPosition = 20;

barPositions =  min(barPositionDegreesFromMidline): binSizeForBarPosition : max(barPositionDegreesFromMidline);

ORDER_MEDFILTER = 800; % parameter than MM uses    
medianFilteredVoltage = medfilt1( voltage, ORDER_MEDFILTER); % Median filtering of the trace

meanVoltageByPos = [];
stdVoltageByPos = [];
currIndex = [];

for i = 1: (length ( barPositions ) -1)
    
    currIndex = ( barPositions(i) < barPositionDegreesFromMidline) & (barPositionDegreesFromMidline < barPositions (i + 1)) ;
   
    
    meanVoltageByPos(i) = mean( medianFilteredVoltage( currIndex ) );
    stdVoltageByPos (i) = std( medianFilteredVoltage( currIndex ) );
     %
    
end

middleOfBarPositionValues = mean([barPositions(1:end-1);barPositions(2:end)]);

figure;
set(gcf, 'Color', 'w');
errorbar( middleOfBarPositionValues, meanVoltageByPos, stdVoltageByPos,'LineWidth', 1.5, 'Color',[0,0.7,0.9]); hold on;
ylabel( 'Vm')
xlabel( 'bar position (deg)');
box off
niceaxes;
xlim( [-80, 200]); 
TICK_STEPS_DEG = 20;
xticks([ (0 - ( MIDLINE_POSITION * DEGREE_PER_LED_SLOT)) : TICK_STEPS_DEG : TOTAL_DEGREES_IN_PATTERN - ( MIDLINE_POSITION * DEGREE_PER_LED_SLOT) ])

title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ' trial#: ' num2str(trialMeta.trialNum) ' stim: ' num2str(stimulus.name) ] )

