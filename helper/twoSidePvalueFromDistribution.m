function [ pValue ] = twoSidePvalueFromDistribution( distribution , dataPoint )
% probabiliyOfPointWithinDisritubtion - plots probabily density funciton and cumlative
% distriution function for distributions, and then solved for 
% the pvalue and CI location of another data point within that distribution
%
% Yvette Fisher 12/2018

% make the probability density function curve by building a histgram of the values 
% figure; set(gcf, 'Color', 'w');
H = histogram( distribution, 100, 'Normalization','pdf'); hold on;
pdf = H.Values;

% loop over and build the cdf - cumlative distribution function
cumlativeSum = cumsum( pdf );

% normalize
cdf = cumlativeSum /  max(cumlativeSum);
binMiddles = H.BinEdges(1:end-1) + ( H.BinWidth / 2 );%

%  Step 1) calculate distance from mean(distribution) to dataPoint
meanToDataPointDist = abs( mean( distribution ) - dataPoint );

% Step 2) find bottom probability.  Take mean - meanToDataPoint, look up where
% this values falls on CDF - record that values as the bottom probablity
bottomCutOff = mean( distribution ) - meanToDataPointDist;
plot([bottomCutOff, bottomCutOff], [0,1]);

bottomPval = max( cdf( binMiddles<= bottomCutOff ));
if( isempty( bottomPval ) )
    bottomPval = 0;
end

% Step 3) find top probability.  Take mean + meanToDataPoint, look up where 
% this values falls on CDF - record probaility as 1-CDF location
topCutOff = mean( distribution ) + meanToDataPointDist;
plot([topCutOff, topCutOff], [0,1]);

topPval = 1 - max( cdf( binMiddles<= topCutOff ));
if( isempty( topPval ) )
    topPval = 0;
end

% 4) add top and bottom probabilites together to get two-sided pvalue
pValue = bottomPval + topPval;

end

