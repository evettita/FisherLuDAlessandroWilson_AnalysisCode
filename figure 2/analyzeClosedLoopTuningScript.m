%% analyzeClosedLoopTuningScript
% Analysis code for plotting data obtained using barRandLocON stimulus
% displays both all trials from an expeirment in order as well as
% averaged tuning curves for the full expeirment
% 
% Written to look at data for Figure 1 of the manuscript, modified form
% analyzeClosedLoopVsOpenLoopTuning script
%   
% Yvette Fisher 11/8/18
%% %% Load in closed loop data set 
clear all; close all;
ephysSettings;

% INPUT wanted trials here
possibleTrials_1barBefore = 1:7;

% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = true;
trialFilesListClosedLoop = extractTrialsWithCertainStimulusName( 'closedLoop_vertStripeON_270world', includeIfNameContainsString, possibleTrials_1barBefore);

%% Plot all selected trials 
plotSelectedTrials( trialFilesListClosedLoop )


% At this stage in the script exclude any trials from analysis that do not meet quality
% criterea:
%-ClosedLoop_vertStripe270 stimulus trials are only taken from those collected before any stimuli with inverted CL or 2-3 bars, especially 2 bar closed loop stimuli was shown to the fly
%-Ave Vm for each trial is within 15 mv of first trial 
%-Ave Vm for each trial is below -20mV
%-Spike amplitude within 50% of initial size =  not dramatic evidence of decrease in Ra
%-Remove trials with large hyper pol/ depot events that look like that crazy thing
% Fly visits all bar locations in 4min trial for the trial to be included
% Must have at least 2 ClosedLoop_vertStripe270 trials to be in the data
% set

%% loop over all the files and extract the medianFiltered version of barPosition data
MEDIAN_FILTER_WIDTH_SEC = 0.04; % sec
ORDER_MEDFILTER = MEDIAN_FILTER_WIDTH_SEC * settings.sampRate;

for fileNum = 1 : length ( trialFilesListClosedLoop )
    
    cd( trialFilesListClosedLoop( fileNum ).folder );
    % load current file for current trial
    load( trialFilesListClosedLoop(fileNum).name ); 
    xPanelPos{ fileNum } =  data.xPanelPos ;
    medFilteredVoltage { fileNum } = medfilt1( data.voltage, ORDER_MEDFILTER , 'truncate' ); % Median filtering of the trace

    %load/store trial number
    trialNums( fileNum ) = trialMeta.trialNum;
    
    trial_stimPatternNum ( fileNum ) = stimulus.panelParams.patternNum;
    
    %load/store stimulus name
    trialStimulusName{ fileNum } = stimulus.name;
end

% Ajust into degrees 270 degree panels
TOTAL_DEGREES_IN_PATTERN = 360;
DEGREE_PER_LED_SLOT = 360 / 96;  % 
MIDLINE_POSITION = 34;%
EDGE_OF_SCREEN_POSITION = 72; % last LED slot in postion units
EDGE_OF_SCREEN_DEG = (EDGE_OF_SCREEN_POSITION -  MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

POSSIBLE_BAR_LOCATIONS = 2:2:71;
barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

 MINBARPOS = 1;
 MAXBARPOS = 72;
% makes barPositions match the edges of the interval used for open loop
% presentation
binIntervalSize_degrees = mean(  barPositionDegreesFromMidline(2 : end) - barPositionDegreesFromMidline(1 : end-1) );
barPositions = [ barPositionDegreesFromMidline - (binIntervalSize_degrees / 2), barPositionDegreesFromMidline(end) + (binIntervalSize_degrees / 2) ];

% extract data for heat map

meanVoltageByPosAcrossTrials=[];
numberOfSamplesByPosAcrossTrials=[];

for fileNum = 1 : length ( trialFilesListClosedLoop )
    % 
    patternOffsetPosition = decodePatternOffset( trial_stimPatternNum( fileNum) );
    
    xBarPostion =  xPanelPos{fileNum} + patternOffsetPosition;
    xBarPostion = mod( xBarPostion , MAXBARPOS );
    
    barPositionDegreesVsMidline = (  xBarPostion - MIDLINE_POSITION) * DEGREE_PER_LED_SLOT;
    
    meanVoltageByPos = [];
    semVoltageByPos = [];
    currLogical = [];
    numberOfSamplesInduced = [];
    secondsSpendInPosition = [];
    
    for i = 1: (length ( barPositions ) -1 )
        currLogical = ( barPositions(i) < barPositionDegreesVsMidline) & (barPositionDegreesVsMidline < barPositions (i + 1)) ;
        
        % store mean voltage values, empty = NaN
        meanVoltageByPos(i) = mean( medFilteredVoltage{fileNum}( currLogical ) );
        semVoltageByPos (i) = std( medFilteredVoltage{fileNum}( currLogical ) ) / sqrt( length (medFilteredVoltage{fileNum})) ;
        
        % store number of samples that exist for current bin
        numberOfSamplesInduced(i) = sum( currLogical );
        secondsSpendInPosition(i) = numberOfSamplesInduced(i) / settings.sampRate ; % convert samples to seconds
    end
    
    % store voltage values
    meanVoltageByPosAcrossTrials( fileNum, : ) = meanVoltageByPos; 
    numberOfSamplesByPosAcrossTrials ( fileNum, : ) = numberOfSamplesInduced;
end
  
middleOfBarPositionValues = round( mean( [barPositions( 1:end - 1) ; barPositions( 2:end ) ]) );

% Plot Heatmap
figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');

imagesc( meanVoltageByPosAcrossTrials ,'AlphaData', ~isnan(meanVoltageByPosAcrossTrials)); hold on;
set(gca,'color',0*[1 1 1]); 

c = colorbar;
c.Label.String = ' Vm  (mV)';


set(gca, 'xtick', 1: length(middleOfBarPositionValues) );
set(gca, 'xticklabel', middleOfBarPositionValues )
xlabel(' Pattern position (deg)')

set(gca, 'ytick', 1:length(trialNums) )
% build tick labels
tickLabels = [];
for i = 1: length (trialNums)
    tickLabels{i} = [ 'trial ' num2str( trialNums(i)) ',' trialStimulusName{i} ];
end
set(gca, 'yticklabel', tickLabels )
niceaxes

title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ] )

% Figure out to get the average closed loop tuning curve without sampling different bar postions different amounts 
% perhaps pull the whole data set about before this analysis


% or used the information about samples in each bin to make sure the
% average is correctly weighted

% 1)take response from each trial the matrix multiple it by the number of samples in
% that bin
% 2) sum the Vm signals for each bar position,
% 3) divide the number of the total number of samples for the bar position
% in the data set
voltageWeightedBySampleNum = meanVoltageByPosAcrossTrials .* numberOfSamplesByPosAcrossTrials;
sumVoltageByBarPos =  nansum(voltageWeightedBySampleNum);
sumSampleByBarPos = nansum(numberOfSamplesByPosAcrossTrials);

voltageByBarPosition = sumVoltageByBarPos ./ sumSampleByBarPos;

figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');
imagesc( voltageByBarPosition );
c = colorbar;
c.Label.String = ' Vm  (mV)';


set(gca, 'xtick', 1: length(middleOfBarPositionValues) );
set(gca, 'xticklabel', middleOfBarPositionValues )
xlabel(' Pattern position (deg)')

set(gca, 'ytick', 1:length(trialNums) )
% build tick labels
tickLabels = [];
for i = 1: length (trialNums)
    tickLabels{i} = [ 'trial ' num2str( trialNums(i)) ',' trialStimulusName{i} ];
end
set(gca, 'yticklabel', tickLabels )
niceaxes;

title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum) ] )


%% Fig S2, Plot 2 EPG cells from the same fly:
% Ajust into degrees 270 degree panels
TOTAL_DEGREES_IN_PATTERN = 360;
DEGREE_PER_LED_SLOT = 360 / 96;  % 
MIDLINE_POSITION = 34;%
EDGE_OF_SCREEN_POSITION = 72; % last LED slot in postion units
EDGE_OF_SCREEN_DEG = (EDGE_OF_SCREEN_POSITION -  MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

POSSIBLE_BAR_LOCATIONS = 2:2:71;
barPositionDegreesFromMidline = ( POSSIBLE_BAR_LOCATIONS - MIDLINE_POSITION ) * DEGREE_PER_LED_SLOT;

% closed loop data 4/3/18 fly 113 2 fills
%twoRFs = [-36.0261827000000,-36.5012248900000,-36.4153547600000,-35.9661848300000,-37.5929094500000,-36.4996730700000,-39.2506758600000,-38.1202957100000,-38.5558632100000,-40.0056465600000,-40.6671602500000,-39.6547487900000,-41.2346608200000,-40.7858873700000,-41.4825391700000,-38.7292944800000,-35.3811916100000,-35.3661376600000,-33.5831814000000,-33.3607651700000,-33.2677294600000,-32.5916925200000,-32.8526500100000,-34.6333442000000,-36.0963436500000,-35.4182440700000,-35.0301238800000,-35.4748362600000,-33.7960856000000,-35.9654613100000,-37.0080588700000,-36.5090180000000,-39.8503516000000,-38.1774360900000,-36.8306262100000;-39.3625295100000,-39.5551826400000,-39.8609980200000,-40.7053161700000,-39.9523413000000,-40.8327924200000,-39.6940546800000,-40.5158011800000,-40.9846813100000,-40.4161667000000,-39.1655837200000,-38.7031752500000,-34.9429731300000,-41.4424115800000,-43.1507576300000,-41.0354601900000,-39.9665824400000,-36.4790737400000,-37.3108664500000,-34.2678216000000,-36.6918558900000,-35.0959459300000,-31.1766547700000,-31.0051825000000,-32.7382163400000,-33.6792377000000,-33.3999530600000,-33.6347358700000,-34.2088637300000,-33.0789394600000,-33.4910933100000,-33.1769654700000,-37.3927192000000,-37.0544063200000,-39.1644789200000];
% open loop fly 113
%twoRFs = [0.376362840000000,1.07270135300000,-0.623160552000000,-2.31035009500000,-2.06939229300000,-1.97138053400000,-0.845985682000000,-1.31160489200000,-0.519532785000000,-2.55616170300000,-2.60169263300000,-3.49347374400000,-1.96180256500000,-2.22315754800000,-2.33074791500000,-1.52437522400000,-1.10797164000000,0.559284946000000,2.40233722900000,2.33174510400000,1.83682456700000,2.71983070000000,2.46654431700000,2.40099126600000,-0.265080213000000,1.83364309000000,0.747280084000000,-0.0629772210000000,1.46412385700000,-0.637148314000000,0.152201025000000,-0.283898980000000,-1.56413479000000,-2.02194767500000,-2.55219002200000;-2.67730324900000,-2.26321198600000,-4.19545741600000,-4.82090191600000,-4.64421568500000,-4.21277224200000,-3.60416903400000,-3.05096199000000,-3.61545548800000,-2.91573140500000,-4.40956447500000,-4.59536392600000,-3.66596057500000,-1.76768763200000,-0.172912466000000,-1.71478934700000,-2.47906137700000,0.397644168000000,0.671504199000000,2.26425443900000,1.15906374200000,2.35714597500000,0.543872484000000,0.574566409000000,1.21770652200000,0.111790856000000,0.787858347000000,-0.124986525000000,0.819751325000000,-0.569442342000000,1.76171040900000,1.71558940000000,-0.948208647000000,-1.43893349800000,-1.16245587200000];

%2/7/19 fly 338 2 fills
%twoRFs = [-40.4375998900000,-41.1347307000000,-41.4080825400000,-41.7101146500000,-42.8750545800000,-42.4922226200000,-43.7252330700000,-42.7791902400000,-43.5967490000000,-43.0526073400000,-43.4039094200000,-43.5522891500000,-41.9896228500000,-41.7307299600000,-41.2237459600000,-40.0204012200000,-41.6520508200000,-39.5551465200000,-37.1746989200000,-35.8578015300000,-36.8091109200000,-36.7974689900000,-38.3656485600000,-35.7362432400000,-36.9770036000000,-36.7699834900000,-36.5547738500000,-36.7564302800000,-35.1912170600000,-36.3539991700000,-36.1489171900000,-36.9264266200000,-36.4948614700000,-37.1958689100000,-39.2886923400000;-40.9902799400000,-41.9201494300000,-43.2754989400000,-42.6920281300000,-44.5393966300000,-44.5042278500000,-43.8047090200000,-44.3988748200000,-42.4256673300000,-42.5891144500000,-41.4841968300000,-41.2610366700000,-41.5996774200000,-39.9206595300000,-39.5924644100000,-39.0522401700000,-38.4376230200000,-39.8081298300000,-38.9589659200000,-38.4520557700000,-38.2641887400000,-38.3616556300000,-38.4779518100000,-38.9860586500000,-40.1874334900000,-39.3310184800000,-39.1288955500000,-39.2901578100000,-40.7501745400000,-38.5849918900000,-39.1649918900000,-40.1589488500000,-39.4306901500000,-38.8502373200000,-40.7525797400000];
% open loop fly 338
%twoRFs = [-2.84466021300000,-1.03647838000000,-4.33353210800000,-2.91784324100000,-3.43794260000000,-4.92181037500000,-6.25341629000000,-3.02741737100000,-4.65095831900000,-3.40357628300000,-4.15300402900000,-5.21781887000000,-2.91887711300000,-2.79376490600000,-1.70586212100000,-0.693222535000000,-1.45929521600000,0.603711482000000,4.03875743300000,3.20971354600000,-0.109078099000000,0.629506530000000,1.20484363600000,0.698285773000000,0.720744686000000,3.12351132400000,3.79443353300000,2.43874916100000,2.89172064100000,2.78104744900000,3.79241692500000,4.26852744000000,1.21700174200000,-0.575347008000000,0.435750297000000;-2.14365042200000,-1.77862217400000,-2.90493611000000,-3.21618428600000,-4.03567675000000,-3.96421653500000,-4.05372870700000,-3.75499455000000,-3.55563605700000,-1.96025568100000,-1.98747098600000,-1.61513205200000,-0.364561054000000,-0.0694789460000000,0.0871327660000000,-0.771268168000000,0.337768123000000,-0.380036734000000,0.556456271000000,-0.850977497000000,-2.05227545500000,-1.73622920300000,-1.99174281900000,-2.09964527800000,-3.51435556500000,-3.61677379700000,-0.883314324000000,-2.18635033200000,-2.68519702600000,-2.23475666300000,-2.11194602900000,-2.27342096300000,-2.08897596400000,-1.40854330700000,-1.76954347000000];




% open loop fly 127
%twoRFs =[-0.0144623810000000,0.481359068000000,-0.355629614000000,1.25204459100000,0.994433288000000,2.83326213600000,2.85937263400000,2.39794429500000,2.58044341400000,2.51288775500000,2.84729798900000,2.73584631000000,2.05194535300000,1.47671963800000,0.301398429000000,1.74341249700000,2.02636216700000,1.71024893200000,1.46739416800000,0.261692395000000,-0.360441026000000,-0.938429611000000,-1.32558645500000,-1.56568774600000,-1.84585378500000,-1.13733809500000,-2.72906921300000,-1.18947866900000,-3.64267668000000,-2.41645668800000,-3.52825312400000,-3.69835156600000,-2.32389523600000,-2.59595468200000,-2.73833759400000;-0.829577220000000,-1.88917958000000,-1.89897960500000,-1.14610038400000,-1.61634806500000,-1.98791386900000,-2.13273285100000,-2.08547233600000,-1.49273117300000,-1.96590326100000,-0.926236053000000,-1.14593838900000,-0.340214322000000,-1.09006856200000,-2.23757017400000,-3.04111847400000,-2.05961470500000,-1.83240178900000,-1.89202446500000,-2.05867551700000,-2.24242428500000,-2.22678572700000,-1.91544160400000,-2.61023354200000,-2.34784862100000,-1.42161415000000,-1.77922097900000,-0.905887921000000,-0.763107041000000,-1.05546748600000,-1.47670612000000,-1.52815499500000,-1.46326435500000,-0.772325211000000,-1.04446651600000];
% open loop fly 171
%twoRFs =[0.626404132000000,-0.0705138140000000,-1.15684519800000,0.250383823000000,-1.02626554300000,0.851326784000000,1.01339056500000,1.36150600200000,1.77735820900000,1.52004526400000,0.300733744000000,0.981699315000000,-1.46476551800000,-0.523092365000000,-0.790013838000000,-1.43247572000000,-0.184499785000000,-0.858152506000000,-1.29677256300000,-2.21017433500000,-1.53571504000000,-1.45368703700000,-3.19143909500000,-3.87883550300000,-1.83702749000000,-3.11942807600000,-3.89865262200000,-2.54385376100000,-4.07594874100000,-1.98032919700000,-3.03441416400000,-1.87997664000000,-3.00650106800000,-1.44162927000000,-1.54435452900000;-0.981051208000000,-1.16251429300000,-0.759451616000000,-1.11568411200000,-1.05150907400000,-0.859705427000000,-1.16943035200000,-0.866548448000000,-0.514010391000000,-0.210724971000000,-0.517304420000000,-0.728854379000000,-0.843288667000000,-0.668752644000000,-0.840448489000000,-0.495271243000000,-0.310424891000000,-0.732768927000000,-0.810503791000000,-0.816253259000000,-0.458172564000000,-1.11289314000000,-1.49231605100000,-0.837652759000000,-1.03728364900000,-1.34063093600000,-1.52886799700000,-0.805152343000000,-1.37558923500000,-0.989756330000000,-0.821865410000000,-0.434078586000000,-1.17218964400000,-1.12150037200000,-1.55082386500000];

% open loop fly 179
% twoRFs =[-1.98073294400000,-2.70584195300000,-2.80326268100000,-2.56139742100000,-2.01188160500000,-2.41250340400000,-1.99371758600000,-2.19360831900000,-1.88753660100000,-1.59547754100000,-1.38990798300000,-1.74302759600000,-0.952492396000000,-1.30111286000000,-1.02908966800000,-0.872764710000000,-1.21731469700000,-0.460326259000000,0.0339315380000000,0.0304058910000000,-0.931247114000000,-1.01099124600000,-0.768654744000000,-0.525629025000000,-1.52099756600000,-0.679230673000000,-0.168597934000000,-0.637457225000000,-0.0656118490000000,-0.607688469000000,0.0214190420000000,0.0510006750000000,0.271502756000000,-1.00146730100000,0.0683936720000000;-2.31715423800000,-2.81385833300000,-2.97992202000000,-3.02658693900000,-3.36099012700000,-2.35980119600000,-2.80838706800000,-2.37325483800000,-2.88134749900000,-2.96866455000000,-2.71137004100000,-2.09168690000000,-2.30531878300000,-1.84410231100000,-1.37882943000000,-1.10083554200000,-0.986112196000000,-0.601205262000000,-0.157805322000000,-0.127356594000000,-0.904493672000000,-0.422123906000000,0.691011144000000,1.23081085600000,0.739640296000000,1.58073802300000,2.41051828900000,1.07345241300000,1.14697299400000,1.47616456800000,1.62455549600000,1.38146984800000,0.666975586000000,0.161404658000000,1.02599724100000]
% open loop fly 274
twoRFs = [-3.53733399900000,-3.47889562500000,-2.32725721300000,-2.64032577000000,-2.80650158400000,-2.71923976700000,-3.12297426300000,-1.33035233600000,-2.52255231700000,-1.73926981500000,-0.670432841000000,-0.539017042000000,-0.259020722000000,-0.232384370000000,1.36356333900000,1.30281398200000,-0.227066234000000,1.48214782500000,1.58788058800000,0.128826097000000,-1.06903555100000,-0.182369040000000,1.23775248900000,0.434236096000000,1.14831797900000,1.04702353700000,1.43534852900000,1.57776662100000,1.87378928700000,1.29155996200000,2.27529270500000,1.86424724200000,1.08412063500000,2.39957890000000,1.33444424400000;-1.10218677700000,-1.20357294000000,-1.02526056500000,-1.43557260500000,-1.04285668700000,-0.180808566000000,-0.0102932090000000,-0.437258309000000,0.710337094000000,-0.182872197000000,0.238579616000000,0.656919889000000,2.62837615900000,1.54395851300000,1.05495720200000,1.47720711300000,0.885076972000000,2.17177519600000,1.23751539300000,0.111021565000000,1.01216177900000,0.458112954000000,0.257614284000000,1.62502488700000,-0.301846954000000,1.62042807100000,-0.932525427000000,-0.0835083520000000,-0.664087766000000,-0.404883456000000,0.523168129000000,0.323125842000000,0.511484236000000,0.795770200000000,-0.0945815260000000];

% Smooth the data with median filter if user wants to:
% order for median filter smooth
ORDER_MEDFILTER =2; %    
filteredRFs = medfilt1( twoRFs , ORDER_MEDFILTER, 'truncate', 2); % Median filtering of the trace
 
figure('Position',[50, 50, 1000, 800]); set(gcf, 'Color', 'w');
plot( barPositionDegreesFromMidline, filteredRFs ,'Color', [ 0 , 0, 0] );        
niceaxes; box off
xlim([ -120 135])
%% % % % % %     % save current plot
         dir = 'F:\Dropbox (HMS)\FisherLuWilson ms\revision\Figures\';
         fileName = [ dir 'flyNum' num2str( 274 ) '_2EPGcell_visualTuning.eps' ];
         print( fileName, '-dpdf');

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