%% decodingVisualCuePosition
%   Bayesian analysis to look at decoding of visual cue position based on
%   the EPG population activity depending on bar position
%   -with lots of help from Jan Drugowitsch and Stephen Holtz
%
% definitions:
% N = number of neurons in populuation data set
% K = number of presenations of each bar position to each neuron
% M = number of Theta positions where the bar was flashed
%
%
% Yvette Fisher 3/2019

%% Extract responses from meta file to be used for baysian decoding analysis
clear all
NUM_OF_TRIALS_REQUIRED = 3; % if a recording doesn't have this many open loop trials, it will not be included
NUM_OF_BAR_PRESENTATIONS_PER_TRIALS = 4; % some bar locations have 5 samples these will be cut out

openLoopMetaDir = '/Users/evettita/Dropbox (HMS)/EphysData/EP-G_recordings/openLoopMetaData';
cd( openLoopMetaDir );

% check how many files in open loop folder
fileList = dir( fullfile(openLoopMetaDir, '*.mat' ));

counter = 1;
% loop over files
for i = 1: numel( fileList )
    % open the current
    load( fileList(i).name )
    
    cellVoltageResponse = [];
    if( length (trialFilesList) >= NUM_OF_TRIALS_REQUIRED ) % if the cell has enough trials,
        
        % loop over the number of trials to use
        for trialNum = 1: NUM_OF_TRIALS_REQUIRED
        
        %  take out the meanVoltageResp data and
        % only keep first number decided from each trial
        currResp = processedData(trialNum).meanVoltageResp(: , 1:NUM_OF_BAR_PRESENTATIONS_PER_TRIALS);
        
        % save values
        cellVoltageResponse = [ cellVoltageResponse currResp]; % save 35 X (NUM_OF_TRIALS_REQUIRED * NUM_OF_BAR_PRESENTATIONS_PER_TRIALS) matrix here
        
        end
        % responses to be used for this data set
        cellRespEnoughTrials{counter} = cellVoltageResponse; % all responses will enough trials
        counter = counter + 1;
    end
    cellAllResponses{i} = cellVoltageResponse;  % all response with blank if not enough trials
end


%% loop over cell array and build 'data' as N X M X K metrix
%Initial data set: matrix of voltage responses N X M X K in size (mV)
N = numel( cellRespEnoughTrials ); % number of neurons
M = 35; % theta positions
K = NUM_OF_TRIALS_REQUIRED * NUM_OF_BAR_PRESENTATIONS_PER_TRIALS; % responses measured per theta positions

data = NaN( [N, M, K] );

for n = 1:length( cellRespEnoughTrials )
    data( n, :, :) = cellRespEnoughTrials{n};    
end


%% Decoding cue position analysis

% loop over each k to be held out
for k = 1 : K
    
    % Pick 1 out of K trials to hold out for all N and M
    heldOutReponse =data( :, theta_m, k); % N x M
    % reduce the dimentionality of heldOutResponse
    heldOutReponse = squeeze( heldOutReponse );
    
    notHeldOutResponse = []; %%data( n , m, NOT k) ;  % N x M x (K - 1)
    
    
    % Use the rest of the (K-1) trials to build gaussian responses model
    % for each neuron (n) as each theta position (m)
    [ mean, variance ] = buildGaussianModel( notHeldOutResponse );
    
    % Choice a response (r_m) vector to calculate the posterior for
    for r_m = 1: M
        
        
        % loop over each theta position and calculate probabily of
        % p(reponse | theta) - likelihood
        % and p( theta | response) - posterior
        for theta_m = 1: M
            
            % use the held out reponses compared with the gaussian model for
            % each neuron to calculate the likelihood:
            
            probabilityOfEachResp = [];
            % calculate the probabily of each r(n)
            for n = 1: N
                
                % double check indexing!!
                probabilityOfEachResp(n) = probabilyResponseGivenGaussianModel( heldOutReponse( n, r_m ), mean(n,theta_m), variance(n,theta_m) );
                
            end
            
            % multiply together all the probabiliy for all N neurons
            sumOfProbabiliyOfResp = sum( probabilityOfEachResp );
            
            % store this probabily at a single value within the likelihood
            likelihood(theta_m) = exp ( -0.5 * sumOfProbabiliyOfResp );
            
        end
        
        % save the current likelyhood and posterior prob-function for this k hold out and this
        % particular r_m
        
        
        
    end
    
    
    % save all the posterior prob functions for this k hold out across all
    % r_m values
    
    
    
end


function [ mean, variance ] = buildGaussianModel ( input )
% input - N x M x (K-1) matrix 
% calculate mean and variances for each n,m 

mean = input;
variance  = input;
% both outputs should be  (N x M)

end

function [ prob ] = probabilyResponseGivenGaussianModel( response, modelMean, modelVariance )
% add description

prob = ( (response - modelMean)^2 ) / modelVariance + log( modelVariance ) ;

end
