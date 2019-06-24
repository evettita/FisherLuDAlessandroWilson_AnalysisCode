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
%% load in data set
clear all
% navigate the directory with data
% cd('C:\Users\Wilson\Documents\GitHub\EphyCode-Yvette\FisherLuWilson_ms_FigureMakingCode\figure 1\nonUniformInhibitionAnalysis_andflyMovementFig1\for Jan')
load('populationData.mat')

%% Decoding cue position analysis
NUM_OF_TRIALS_REQUIRED = 3; % if a recording doesn't have this many open loop trials, it will not be included
NUM_OF_BAR_PRESENTATIONS_PER_TRIALS = 4; % some bar locations have 5 samples these will be cut out

N = length( data(:,1,1)); % number of neurons
M = 35; % theta positions
K = NUM_OF_TRIALS_REQUIRED * NUM_OF_BAR_PRESENTATIONS_PER_TRIALS; % responses measured per theta positions


% loop over each k to be held out
for k = 1 : K
    
    % Pick 1 out of K trials to hold out for all N and M
    heldOutReponse =data( :, :, k); % N x M

    % reduce the dimentionality of heldOutResponse
    heldOutReponse = squeeze( heldOutReponse );
    
    possible_k = 1:K;
    not_k = possible_k( possible_k~= k);
    % all data except the help out epoch
    notHeldOutResponse = data( : , :, not_k) ;  % size is N x M x (K - 1)
    
    
    % Use the rest of the (K-1) trials to build gaussian responses model
    % for each neuron (n) as each theta position (m)
    [ meanNM, varianceNM ] = buildGaussianModel( notHeldOutResponse );
    
    % Choice a response (r_m) vector to calculate the posterior for
    for r_m = 1 : M      
        

        likelihood = [];
        % loop over each theta position and calculate probabily of
        % p(reponse | theta) aka likelihood
        % and p( theta | response) aka  posterior
        for theta_m = 1: M
            
            % use the held out reponses compared with the gaussian model for
            % each neuron to calculate the likelihood:
            comparisonToGaussian = [];
            for n = 1: N
                
                comparisonToGaussian(n) = compareResponseToGaussianModel( heldOutReponse( n, r_m ), meanNM(n , theta_m), varianceNM(n , theta_m) );
                
            end
            
            % sum together all the probabiliy across all N neurons
            sumOfProbabiliyOfResp = sum( comparisonToGaussian );
            
            % store this probabily (likelihood) for theta M
            likelihood(theta_m) = exp( -0.5 * sumOfProbabiliyOfResp );
            
        end
        
        % save the current likelihood and posterior prob-function for this k hold out and this
        % particular r_m
        likelihoodFunction{ k , r_m } = likelihood;
        
        % normalized this function to have all thetaM proabilities sum to 1
        posteriorFunction{ k , r_m } = likelihood / sum( likelihood ); 
        
    end
    
end

%% plot out each posterior function obtained from each held out population response
close all;
figure('Position',[50, 50, 1000, 800]); set(gcf, 'Color', 'w');

DEGREE_PER_LED_SLOT = 360 / 96; 
MIDLINE_POSITION = 34;%
POSSIBLE_BAR_LOCATIONS = 2:2:70;
barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

for m = 1: M
   subplot( 7, 5, m);
   
   currPosteriorFunc = [];
   for k = 1:K
   plot( barPositionDegreesFromMidline,  posteriorFunction{ k , m }, 'Color', [0.5 0.5 0.5]);
   hold on;
    
   currPosteriorFunc( k, :) = posteriorFunction{ k , m };
   end
   
   plot( barPositionDegreesFromMidline,  mean( currPosteriorFunc ), '-k' ); hold on; 
   xlim( [ min(barPositionDegreesFromMidline), max(barPositionDegreesFromMidline) ])
   box off;
   title( ['r: ' num2str( round( barPositionDegreesFromMidline(m),0)) char(0176) ] );
   
end


%% Sub functions
function [ meanNM, varianceNM ] = buildGaussianModel ( input )
% input - N x M x (K-1) matrix 
% calculate mean and variances for each n,m along dimentation K 
K_DIM = 3;

% both outputs should be  (N x M)
meanNM =  mean( input , K_DIM);
varianceNM  = var( input, [], K_DIM);

end

function [ out ] = compareResponseToGaussianModel( response, modelMean, modelVariance )

out = ( (response - modelMean)^2 ) / modelVariance + log( modelVariance ) ;

end
