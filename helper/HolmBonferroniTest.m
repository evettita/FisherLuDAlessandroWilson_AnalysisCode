function [ significanceLogicalArray , numberOfSignValues] = HolmBonferroniTest( pvaluesFromData , alpha )
%HolmBonferroniTest run Bonferroni-Holm test on a collection of pvalues
%that were all tested together to account for multiple comparisons.
%
%       INPUTS:
%       pvaluesFromData - list of pvalues to be tested together with multiple
%                       comparisons correction
%
%       alpha - significance level to reject the Null hypothesis at,
%                       typically 0.05
%
%       OUTPUT:  significanceLogicalArray 0 = not significant, 1 =
%       significant @ alpha (reject the null)
%
%
% Sourse: https://www.statisticshowto.datasciencecentral.com/holm-bonferroni-method/
%
% Yvette Fisher 12/2018
nSamples = length ( pvaluesFromData );

% Step 1: Order the p-values from smallest to greatest:
[orderedpValues , sortedIndex] = sort ( pvaluesFromData );
% save a way to reverse the sorting for later
[~ ,  indexToUndoSort] = sort( sortedIndex );

orderedSignArray = zeros(1 , nSamples);

% Loop over the ordered pvalues and test them against alpha
for i = 1 : nSamples
    
    %Step 2: Work the Holm-Bonferroni formula for the i rank:
    HB = alpha / (nSamples - i + 1);
    
    % compared the ith pvalue to the i rank HB level
    if( orderedpValues(i) < HB)
        orderedSignArray(i) = 1;
    else
        break;%The testing stops when you reach the first non-rejected hypothesis.
        %All subsequent hypotheses are non-significant (i.e. not rejected).
    end
    
end

significanceLogicalArray = orderedSignArray( indexToUndoSort );

numberOfSignValues = sum( significanceLogicalArray );
end

