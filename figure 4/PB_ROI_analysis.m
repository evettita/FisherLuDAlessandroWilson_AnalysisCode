function roi_data = PB_ROI_analysis(parentDir, sid, tid_list)

%% ROI Analysis Function
% Objectives: From ROIs given, create an ROI data structure with following
% information:
%
% roi_analysis.baseline_f
% roi_analysis.dff
% roi_analysis.zscore
% roi_analysis.num
% roi_analysis.xi
% roi_analysis.yi
% roi_analysis.zplane
% roi_analysis.leftPB
% roi_analysis.rightPB
% roi_analysis.glomerulus
% roi_analysis.pixels
%

%% Look for session Files post registration. Find if there is one or two
% channels.
global slash;
if isunix() == 1
    slash = '/';
else
    slash = '\';
end
imagingDir = [parentDir slash '2p' slash];
num_tids = length(tid_list);

for i = 1:num_tids
    
    %% Get files
    tid = tid_list(i);
    imagingTrialDir = [imagingDir 'sid_' num2str(sid) '_tid_' num2str(tid) slash];
    if ~isdir(imagingTrialDir)
        disp('Dir does not exist');
        return;
    end
    expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*'];
    expression2 = ['rigid_*'];
    roiDir = [imagingDir slash 'ROI' slash];
    expression3 = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*'];
    imagingFile = dir(fullfile(imagingTrialDir, expression2));
    roiFile = dir(fullfile(roiDir, expression3));

    numChannels = size(imagingFile, 1);

    %% Load the GCaMP stack of each channel.
    summedStack = [];
    numSlices = [];
    numVols = [];

    if numChannels == 1
        tic;
        load(fullfile(imagingTrialDir, imagingFile.name));
        numSlices = size(regProduct, 3);
        numVols = size(regProduct, 4);
        summedStack = max(regProduct,[], 3);
        regProductGCamp = regProduct;
        toc;
    else % numChannels == 2
        load(fullfile(imagingTrialDir, imagingFile(1).name));
        summedStack = max(regProduct,[], 3);
        numSlices = size(regProduct, 3);
        numVols = size(regProduct, 4);

        regProductGCamp = regProduct;
        tic;
        load(fullfile(imagingTrialDir, imagingFile(2).name));
        regProductTdTom = regProduct;
        toc;
    end

    tic;
    %% Load the ROI file
    load(fullfile(roiDir, roiFile(1).name)); % loads variable roi

    toc;
    %% For each ROI
    roi_nums = [roi(:).num];
    num_rois = max(roi_nums);

    roi_data = [];

    tic;
    for i = 1:num_rois
       %% Look for the rois that match the roi num
       indices = find([roi(:).num]==i);
       summed = [];
       pixels = zeros(1, length(indices));
       
       %% Sum the values across all the rois
       for j = 1:length(indices)
           index = indices(j);
           roi_im = bsxfun(@times, squeeze(regProductGCamp(:,:,roi(index).z,:)), roi(index).BW);
           summed = [summed squeeze(sum(sum(roi_im,1)))];
           mask = roi(index).BW;
           pixels(j) = sum(mask(:));
       end
       %% Calculate the trace as the mean value in the rois (by pixel number)
       trace = summed*pixels'./sum(pixels);
       sorted = sort(trace);
       timepoints = size(trace,1);
       
       %% Detect if "frames lost"
       if sorted(end)/sorted(1) >= 30 % this is my current threshold
           T = clusterdata(trace, 2);
           if mode(T) == 1
               bad_cluster = 2;
           else
               bad_cluster = 1;
           end
           trace(T == bad_cluster) = nan;
           sorted = sort(trace(T == 2));
           timepoints = size(sorted, 1);
       end
       
       %% Baseline is the bottom fifth percentile
       fifthpercentile = floor(.05*timepoints);
       baseline_f = mean(sorted(1:fifthpercentile));
       
       %% Calculate dff and zscore
       dff = (trace-baseline_f)./baseline_f;
       z = zscore(trace);
       
       %% Create roi_analysis object for this ROI
       roi_analysis = [];
       roi_analysis.baseline_f = baseline_f;
       roi_analysis.dff = dff;
       roi_analysis.zscore = z;
       roi_analysis.num = i;
       roi_analysis.xi = {roi(indices).xi};
       roi_analysis.yi = {roi(indices).yi};
       roi_analysis.zplane = [roi(indices).z];
       roi_analysis.leftPB = roi(indices(1)).leftPB;
       roi_analysis.rightPB = roi(indices(1)).rightPB;
       roi_analysis.glomerulus = roi(indices(1)).glomerulus;
       mask = roi(i).BW;
       roi_analysis.pixels = sum(mask(:));

       %% Never updated this part b/c did not do two-color imaging
       if numChannels == 2
           summed_tdTom = [];
           pixels_tdTom = zeros(1, length(indices));
           for j = 1:length(indices)
               index = indices(j);
               roi_im = bsxfun(@times, squeeze(regProductTdTom(:,:,roi(index).z,:)), roi(index).BW);
               summed_tdTom = [summed_tdTom squeeze(sum(sum(roi_im,1)))];
               mask = roi(index).BW;
               pixels_tdTom(j) = sum(mask(:));
           end
           trace_tdTom = summed_tdTom*pixels_tdTom'./sum(pixels_tdTom);
           sorted_tdTom = sort(trace_tdTom);
           timepoints_tdTom = size(trace_tdTom,1);
           fifthpercentile_tdTom = floor(.05*timepoints_tdTom);
           baseline_f_tdTom = mean(sorted_tdTom(1:fifthpercentile_tdTom));
           dff_tdTom = (trace_tdTom-baseline_f_tdTom)./baseline_f_tdTom;
           z_tdTom = zscore(trace_tdTom);
           roi_analysis.baseline_f_tdTom = baseline_f_tdTom;
           roi_analysis.dff_tdTom = dff_tdTom;
           roi_analysis.zscore_tdTom = z_tdTom;
       end
       roi_data = [roi_data roi_analysis];
    end

    toc;
    
    tic;
    % Save data
    roi_analysis_path = [parentDir slash '2p' slash 'ROI_analysis'];
    if(~exist(roi_analysis_path, 'dir'))
    mkdir(roi_analysis_path);
    end
    filename = [roi_analysis_path slash ['ROI_data_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']];
    
    save(filename, 'roi_data');
    toc;
end
end