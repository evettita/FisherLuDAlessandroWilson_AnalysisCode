function [ aveAngularSpeed,  aveForwardSpeed, trialNums ] = findAveFlyMovementByTrial( trialFilesList )
%FINDAVEFLYMOVEMENTBYTRIAL extract the average angular speed and forward
%velocity for each trial and plot them for the user
%   Yvette Fisher 2/2018

ephysSettings;
LOWPASS_FILTER_CUTOFF= 25; % Hz
THRESHOLD_ANGULAR_VELOCITY = 2500; % degrees / s  this is the max velocity that can be allowed into analysis
THRESHOLD_FORWARD_VELOCITY = 2500; % degrees / s  this is the max velocity that can be allowed into analysis

%loop over all the files and extract the ficTrac bar movement information
for fileNum = 1 : length ( trialFilesList )
    
    cd( trialFilesList( fileNum ).folder );
    % load current file for current trial
    load( trialFilesList(fileNum).name , 'data' );
    
    % decode angular velocity and accumulated position
    [ angularVelocity , ~ ] = ficTracSignalDecoding( data.ficTracAngularPosition, settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);
    % save ave angular speed
    aveAngularSpeed( fileNum ) = mean ( abs( angularVelocity ) );
    
    
    % decode forward velocity and accumulated X position
    [ forwardVelocity , ~ ] = ficTracSignalDecoding( data.ficTracIntx , settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_FORWARD_VELOCITY);
    aveForwardSpeed( fileNum ) = mean ( abs( forwardVelocity ) );
    
    %load and store trial number
    load( trialFilesList(fileNum).name , 'trialMeta' );
    trialNums(fileNum) = trialMeta.trialNum;
end

end

