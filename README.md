# FisherLuDAlessandroWilson_AnalysisCode
-Analysis code from Fisher et al. 2019 manuscript. 

## Figure 1:
### analyzeOpenLoopTuningScript.m
- initial processing of visual tuning curves 
### plotOrderedReceptiveFields.m
- Plots heatmaps, histograms and population sum for E-PG visual responses data set
### plotOrderedReceptiveFields_comparingResponseToPost.m
- Compares heatmaps, histograms and population sum for E-PG visual responses when the visual stimulus was present vs 250ms after the visual stimulus was removed.
### plotVisualTuningVsDendriticLocation.m
- Scatter plot of max inhibition vs location of E-PG dendrite in the Ellipsoid body
- Circular correlation coefficient analysis of interactions between visual tuning and dendrite location
### analyzedDeltaYawbyBarPosition.m
- initial processing of fly movement (delta Yaw) relative to bar position
- obtain a pvalue by comparing individual fly movement responses vs a bootstrap distribution for each bar position response
- analyze significance of p values from all flies across all cue positions using Bonferroni-Holms analysis
### plotFlyPopulationDeltaYawByBarPosition.m
- compares fly movement response as a function of bar position.  Calculates 95% confidence interval using bootstrap distribution from randomly drawn position responses
### analyzedDeltaYawbyBarJumpDistance.m
- initial processing of fly movement (delta Yaw) relative to the distance the visual bar jumped
- obtain a pvalue by comparing individual fly movement responses vs a bootstrap distribution for each bar jump response
- compares grand fly movement response as a function of the distance the visual bar jumped. Calculates 95% confidence interval using bootstrap distribution

## Figure 2:
### analyzeClosedLoopTuningScript.m
- initial processing of closed loop heading tuning curves 
### plotOpenLoopVsClosednLoopTuningCurves.m  
- Plots open loop and closed loop tuning curves
- analyzes open loop vs closed loop correlation coef. and compares to shuffled data
- Plots true and shuffled correlation values

## Figure 3: 
### plotBarRandLoc_ringNeuron_180degScreens.m 
- Loads R neuron visual responses, finds spikes, plots tuning curve
### plotChrimsonResponseAmp.m
- Plots E-PG voltage response to chrimson stimuluation
- Scatter plot of mean amplitude vs controls

## Figure 4:
### PB_data_analysis.m
- Imaging and behavior data processing and analysis
### PB_ROI_analysis.m
- ROI analysis function

## Figure 5: 
### analyze2barRemappingData.m 
- script used for checking data and consolidated tuning curves for a full remapping data set. 
### plot2barRemappingDataSet.m
- plot tuning curves for remapping data set, analyzes relationships between receptive field shape changes, absolute changes and modulation of heading tuning during 2 bar training.
### plot1barControlDataSet.m
- plot tuning curves for control remapping data set

## Helpers/dependencies:
- analysisClosedLoopTuning.m
- analyzeOpenLoopTuning.m
- HolmBonferroniTest.m
- twoSidePvalueFromDistribution.m
- ephySettings.m   (Rig parameters)
- niceaxes.m
- findAveFlyMovementByTrial.m
- bluewhitered.m    (color map, MatLab File Exchanges, by Nathan Childress)
- circ_corrcc.m    (Circular Statistics toolbox for Matlab, by Philipp Berens)

