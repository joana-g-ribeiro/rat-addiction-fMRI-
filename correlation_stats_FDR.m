%% Extracting timeseries directly from scans .nii

% MIND THE CLEAR below
clear

%%%%%%%%%%%%%% CHOOSE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Animals:
animals = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};

% ROIs to test
ROI_labels = {''};

% slice indication for each ROI just to confirm the position in the
% highlighted figure, see below
z_ROI = [15 14 14 14 14 11 10 9 10 11];

main_directory = 'C:\';
masks_directory = 'C:\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose treatments to analyze
% Pre-CPP "_pre"; Post-CPP "_pos"; Morphine "_mor"; Saline "_sal"
% label treatments
% MIND the order below, do not skip treatmens as they maybe
% swapped during the loop
treatments = {'_pre', '_pos', '_mor', '_sal'};
treat_label = {'Pre-CPP', 'Post-CPP', 'Morphine', 'Saline'};

% pre-fill the final cell array with all_ROI_Bold in each group treatment
% it is a cell arrays so each cell can have different variables
% Prefill Results if needed (include the predicted number of variables to
% store
Results(1,:) = {treat_label{1},animals};
Results(2,:) = {treat_label{2},animals};
Results(3,:) = {treat_label{3},animals};
Results(4,:) = {treat_label{4},animals};
Results(5,:) = {'Animals', animals};

% outer loop through treatments

for t = 1:length(treatments)
    
    for a = 1:length(animals)
        
    % load swscan nifti file obtained in SPM, total scan
    % in the imagesc figure row = Y, column = X, 3d = Z, 4d = scan
    filename = [main_directory animals{a} treatments{t}];
    tscan = niftiread(filename);
    numslices = length(tscan(1,1,:,1));
    numvol = length(tscan(1,1,1,:));
    
    % prefill matrix
    all_ROI_Bold = zeros(length(ROI_labels), numvol);
    
for m = 1:length(ROI_labels)
    
    ROI_filename = [masks_directory ROI_labels{m}];
            
    % load mask for ROI .nii file
    roimask = niftiread(ROI_filename);
    
    % convert roimask uint8 format into int16 as tscan
    % and multiply tscan by mask to obtain only tscan values inside the mask
    roimask = cast(roimask,'like',tscan);
    tscan_masked = tscan.*roimask;
    
    % find the number of nonzero elements in volume 1
    NZ =nnz(tscan_masked(:,:,:,1));
    
    % Create timeseries for current mask
    % Each volume value is obtained by the sum of all slices / NZ
    sumSlice = zeros(1,numslices);
    Bold = zeros(1,numvol);
    
    for x = 1: numvol
        
        for k = 1:numslices
            sumSlice(k) = sum(sum(tscan_masked(:,:,k,x)));
        end
        totalVolume = sum(sumSlice);
        
        Bold(x) = totalVolume/NZ;
        
    end
    
    % fill the all ROI Bold matrix with curent BOLD
    all_ROI_Bold(m, :)= Bold;
    
    
    % Turn OFF this plot for now since too many figures will be created
    %scans = 1:numvol;
    %{
    % plot timeseries
    figure (a*10000 + t*m)
    plot(scans, Bold)
    title(['Treatment' treatments{t} ' ' 'Animal# ' animals{a} ' - ' ROI_labels{m}], 'Interpreter', 'none')
    xlabel('Scans')
    ylabel('BOLD')
    %}
end

% create tables with the ROI labels for entire and mean timeseries
labeled_all_ROI_timeseries = table(ROI_labels', all_ROI_Bold);

% create correlation and p value matrices
% maybe these p values are ignored now since a t-test will be performed
% later inter-treatment only with correlation coeficients Z-transformed
% note the need to transpose all_ROI_Bold so columns represent variables 
[R, P] = corrcoef(all_ROI_Bold(:,1:numvol)');

% Fisher Z-transform the correlation coeficients
R = 0.5*(log(1+R)-log(1-R));

% create table so ROIs are labeled
labeled_correlations = table(ROI_labels', R);
% MAYBE we ignore these p values for now but save anyway
labeled_p_values = table(ROI_labels', P);

Results{t,2}{a} = {labeled_all_ROI_timeseries, labeled_correlations, labeled_p_values};
    end
end

% SAVE .mat file with all relevant variables for each animal scan
% REMEMBER the order of treatments and variables:
% labeled_all_ROI_timeseries, labeled_correlations, labeled_p_values
% eg. Results{1,2}{1} is the labeled_correlations 1st treatment "Pre-CPP";
% eg. Results{2,2}{1} is the labeled_correlations for 2nd treatment "Post-CPP";

save('Results', 'Results');

%% Show ROI mask highlighted for last treatment last animal 

% choose volume/scan to view
    volview = 1;

for m = 1:length(ROI_labels)
 % show tscan with mask highlighted
    % Create figure for the last animal last treatment only, for quality
    % control
    
    filename = [main_directory animals{a} treatments{t}];
    tscan = niftiread(filename);
    
    ROI_filename = [masks_directory ROI_labels{m}];
    
    % load mask for ROI .nii file
    roimask = niftiread(ROI_filename);
    
    % convert roimask uint8 format into int16 as tscan
    % and multiply tscan by mask to obtain only tscan values inside the mask
    roimask = cast(roimask,'like',tscan);
    tscan_masked = tscan.*roimask;
    
    % choose slice and volume/scan to view ONLY
    sliceview = z_ROI(m);
     
    % make matrix with highlighted values in the mask region
    roimask_highlight = roimask + 1;
    tscan_highlight = tscan.*roimask_highlight;
    
    figure (a*1000 + t*m)
    imagesc(tscan_highlight(:,:,sliceview,volview));
    colorbar
    title(['Treatment' treatments{t} ' ' 'Animal# ' animals{a} ' ' ROI_labels{m}], 'Interpreter', 'none')
end

% Extract data like this:
%  = Results{treatment index,2}{animal index}{variable index}.variable name;
% Example: extract R matrix from table "labeled_correlations" (index 2), 4th treatment, 3rd animal:
%  = Results{4,2}{3}{2}.R;

%% Paired t-test
% Test the null hypothesis that the pairwise difference between data vectors has a mean equal to zero
% h = 1 indicates that ttest rejects the null hypothesis at the chosen significance (Alpha) level

% Preallocate matrices for loop
stats_matrix = ones(length(ROI_labels));
mean_matrix_Mor = ones(length(ROI_labels));
mean_matrix_Sal = ones(length(ROI_labels));
mean_matrix_Pre = ones(length(ROI_labels));
mean_matrix_Pos = ones(length(ROI_labels));
stats_matrix_FDR = zeros(length(ROI_labels));

t_Stats(1,:) = {'Morphine vs Saline',stats_matrix};
t_Stats(2,:) = {'Pre- vs Post-CPP',stats_matrix};
t_Stats(3,:) = {'Morphine vs Post-CPP',stats_matrix};
t_Stats(4,:) = {'Saline vs Post-CPP',stats_matrix};
t_Stats(5,:) = {'Morphine vs Saline FDR adjusted',stats_matrix};
t_Stats(6,:) = {'Pre- vs Post-CPP FDR adjusted',stats_matrix};
t_Stats(7,:) = {'Morphine vs Post-CPP FDR adjusted',stats_matrix};
t_Stats(8,:) = {'Saline vs Post-CPP FDR adjusted',stats_matrix};


% preallocale vectors for animals for t-test
vector1 = 1:length(animals);
vector2 = 1:length(animals);

% t-test for Morphine vs Saline

% loop through matrices' rows
for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % loop through animals for Morphine and Saline
        for a = 1:length(animals)
            % Morphine
            matrix1 = Results{3,2}{a}{2}.R;
            vector1(a) = matrix1(m,n);
            % Saline
            matrix2 = Results{4,2}{a}{2}.R;
            vector2(a) = matrix2(m,n);
        end
        
        % mean now for Morphine and Saline and store
        mean_matrix_Mor(m,n) = mean(vector1);
        mean_matrix_Sal(m,n) = mean(vector2);
        
        % t-test now between Morphine and Saline
        [~,p_ttest] = ttest(vector1,vector2,'Alpha',0.05);
        % store p_value in the correct position in stats_matrix       
        stats_matrix(m,n) = p_ttest;
               
    end
        
end


% Store current stats_matrix in t_Stats correct position cell array
% Morphine vs Saline on row 1
t_Stats{1,2} = stats_matrix;

% FDR adjusted p-values
% extract lower triangle from matrix excluding diagonal
low_tri = tril(stats_matrix, -1);
% reshape to a vector excluding zeros and sort in ascending order
rank_vector = sort(reshape(low_tri(low_tri~=0),[],1));
% calculate total number of tests
total_tests = (length(ROI_labels)*length(ROI_labels)-length(ROI_labels))/2;
% for FDR adjusted p-value = p-value*(total number of hypotheses tested)/(rank of the p-value)
% loop through matrices' rows
for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % if NaN in original matrix leave the same value
        % only assign new value if a value exists
        if ~isnan(stats_matrix(m,n))
            stats_matrix_FDR(m,n) = stats_matrix(m,n)*(total_tests/find(rank_vector == stats_matrix(m,n)));
        end
     end
end

% Store current stats_matrix in t_Stats correct position cell array
% Morphine vs Saline FDR on row 4
t_Stats{5,2} = stats_matrix_FDR;


% t-test for Pre-CPP vs Post-CPP

% loop through matrices' rows
for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % loop through animals for Pre-CPP and Post-CPP
        for a = 1:length(animals)
            % Pre-CPP
            matrix1 = Results{1,2}{a}{2}.R;
            vector1(a) = matrix1(m,n);
            % Post-CPP
            matrix2 = Results{2,2}{a}{2}.R;
            vector2(a) = matrix2(m,n);
        end
        % mean now for Pre-CPP and Post-CPP and store
        mean_matrix_Pre(m,n) = mean(vector1);
        mean_matrix_Pos(m,n) = mean(vector2);
        % t-test now between Morphine and Saline
        [~,p_ttest] = ttest(vector1,vector2,'Alpha',0.05);
        % store p_value in the correct position in stats_matrix       
        stats_matrix(m,n) = p_ttest;
        
                      
    end
        
end
    
% Store updated stats_matrix in t_Stats correct position cell array
% Pre vs Post on row 2
t_Stats{2,2} = stats_matrix;

% FDR adjusted p-values
% extract lower triangle from matrix excluding diagonal
low_tri = tril(stats_matrix, -1);
% reshape to a vector excluding zeros and sort in ascending order
rank_vector = sort(reshape(low_tri(low_tri~=0),[],1));
% calculate total number of tests
total_tests = (length(ROI_labels)*length(ROI_labels)-length(ROI_labels))/2;
% for FDR adjusted p-value = p-value*(total number of hypotheses tested)/(rank of the p-value)
% loop through matrices' rows
for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % if NaN in original matrix leave the same value
        % only assign new value if a value exists
        if ~isnan(stats_matrix(m,n))
            stats_matrix_FDR(m,n) = stats_matrix(m,n)*(total_tests/find(rank_vector == stats_matrix(m,n)));
        end
     end
end

% Store current stats_matrix in t_Stats correct position cell array
% Pre vs Post FDR on row 5
t_Stats{6,2} = stats_matrix_FDR;


% t-test for Morphine vs Post-CPP

for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % loop through animals for Morphine and Post-CPP
        for a = 1:length(animals)
            % Morphine
            matrix1 = Results{3,2}{a}{2}.R;
            vector1(a) = matrix1(m,n);
            % Post-CPP
            matrix2 = Results{2,2}{a}{2}.R;
            vector2(a) = matrix2(m,n);
        end
        
        % mean now for Morphine and Post-CPP and store
        mean_matrix_Mor(m,n) = mean(vector1);
        mean_matrix_Pre(m,n) = mean(vector2);
        % t-test now between Morphine and Post-CPP
        [~,p_ttest] = ttest(vector1,vector2,'Alpha',0.05);
        % store p_value in the correct position in stats_matrix       
        stats_matrix(m,n) = p_ttest;
               
    end
        
end

% Store current stats_matrix in t_Stats correct position cell array
% Morphine vs Post-CPP on row 3
t_Stats{3,2} = stats_matrix;

% FDR adjusted p-values
% extract lower triangle from matrix excluding diagonal
low_tri = tril(stats_matrix, -1);
% reshape to a vector excluding zeros and sort in ascending order
rank_vector = sort(reshape(low_tri(low_tri~=0),[],1));
% calculate total number of tests
total_tests = (length(ROI_labels)*length(ROI_labels)-length(ROI_labels))/2;
% for FDR adjusted p-value = p-value*(total number of hypotheses tested)/(rank of the p-value)
% loop through matrices' rows
for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % if NaN in original matrix leave the same value
        % only assign new value if a value exists
        if ~isnan(stats_matrix(m,n))
            stats_matrix_FDR(m,n) = stats_matrix(m,n)*(total_tests/find(rank_vector == stats_matrix(m,n)));
        end
     end
end

% Store current stats_matrix in t_Stats correct position cell array
% Morphine vs Post FDR on row 6
t_Stats{7,2} = stats_matrix_FDR;


% t-test for Saline vs Post-CPP

for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % loop through animals for Saline and Post-CPP
        for a = 1:length(animals)
            % Saline
            matrix1 = Results{4,2}{a}{2}.R;
            vector1(a) = matrix1(m,n);
            % Post-CPP
            matrix2 = Results{2,2}{a}{2}.R;
            vector2(a) = matrix2(m,n);
        end
        
        % mean now for Morphine and Post-CPP and store
        mean_matrix_Sal(m,n) = mean(vector1);
        mean_matrix_Pre(m,n) = mean(vector2);
        % t-test now between Morphine and Post-CPP
        [~,p_ttest] = ttest(vector1,vector2,'Alpha',0.05);
        % store p_value in the correct position in stats_matrix       
        stats_matrix(m,n) = p_ttest;
               
    end
        
end

% Store current stats_matrix in t_Stats correct position cell array
% Morphine vs Post-CPP on row 3
t_Stats{4,2} = stats_matrix;

% FDR adjusted p-values
% extract lower triangle from matrix excluding diagonal
low_tri = tril(stats_matrix, -1);
% reshape to a vector excluding zeros and sort in ascending order
rank_vector = sort(reshape(low_tri(low_tri~=0),[],1));
% calculate total number of tests
total_tests = (length(ROI_labels)*length(ROI_labels)-length(ROI_labels))/2;
% for FDR adjusted p-value = p-value*(total number of hypotheses tested)/(rank of the p-value)
% loop through matrices' rows
for m=1:length(ROI_labels)
    % loop through matrices' columns
    for n=1:length(ROI_labels)
        % if NaN in original matrix leave the same value
        % only assign new value if a value exists
        if ~isnan(stats_matrix(m,n))
            stats_matrix_FDR(m,n) = stats_matrix(m,n)*(total_tests/find(rank_vector == stats_matrix(m,n)));
        end
     end
end

% Store current stats_matrix in t_Stats correct position cell array
% Morphine vs Post FDR on row 6
t_Stats{8,2} = stats_matrix_FDR;


%%
% Create heatmaps for final correlation matrices' p values

% heatmap Morphine vs Saline
figure(1)
h1 = heatmap(t_Stats{1,2});
h1.XDisplayLabels = ROI_labels;
h1.YDisplayLabels = ROI_labels;
h1.GridVisible = 'off';
h1.Title = 'Morphine vs Saline (p-value)';

% heatmap Pre-CPP vs Post-CPP
figure(2)
h2 = heatmap(t_Stats{2,2});
h2.XDisplayLabels = ROI_labels;
h2.YDisplayLabels = ROI_labels;
h2.GridVisible = 'off';
h2.Title = 'Pre-CPP vs Post-CPP (p-value)';

% heatmap Morphine vs Post-CPP
figure(3)
h3 = heatmap(t_Stats{3,2});
h3.XDisplayLabels = ROI_labels;
h3.YDisplayLabels = ROI_labels;
h3.GridVisible = 'off';
h3.Title = 'Morphine vs Post-CPP (p-value)';

% heatmap Saline vs Post-CPP
figure(4)
h4 = heatmap(t_Stats{4,2});
h4.XDisplayLabels = ROI_labels;
h4.YDisplayLabels = ROI_labels;
h4.GridVisible = 'off';
h4.Title = 'Saline vs Post-CPP (p-value)';

% heatmap Morphine vs Saline FDR ADJUSTED p-values
figure(5)
h5 = heatmap(t_Stats{5,2});
h5.XDisplayLabels = ROI_labels;
h5.YDisplayLabels = ROI_labels;
h5.GridVisible = 'off';
h5.Title = 'Morphine vs Saline (FDR adjusted p-value)';

% heatmap Pre-CPP vs Post-CPP FDR ADJUSTED p-values
figure(6)
h6 = heatmap(t_Stats{6,2});
h6.XDisplayLabels = ROI_labels;
h6.YDisplayLabels = ROI_labels;
h6.GridVisible = 'off';
h6.Title = 'Pre-CPP vs Post-CPP (FDR adjusted p-value)';

% heatmap Morphine vs Post-CPP FDR ADJUSTED p-values
figure(7)
h7 = heatmap(t_Stats{7,2});
h7.XDisplayLabels = ROI_labels;
h7.YDisplayLabels = ROI_labels;
h7.GridVisible = 'off';
h7.Title = 'Morphine vs Post-CPP (FDR adjusted p-value)';

% heatmap Saline vs Post-CPP FDR ADJUSTED p-values
figure(8)
h8 = heatmap(t_Stats{8,2});
h8.XDisplayLabels = ROI_labels;
h8.YDisplayLabels = ROI_labels;
h8.GridVisible = 'off';
h8.Title = 'Saline vs Post-CPP (FDR adjusted p-value)';

% Create heatmaps for final Mean Correlations

% heatmap Morphine
figure(9)
h9 = heatmap(mean_matrix_Mor);
h9.XDisplayLabels = ROI_labels;
h9.YDisplayLabels = ROI_labels;
h9.GridVisible = 'off';
h9.Title = 'Morphine (mean)';

% heatmap Saline
figure(10)
h10 = heatmap(mean_matrix_Sal);
h10.XDisplayLabels = ROI_labels;
h10.YDisplayLabels = ROI_labels;
h10.GridVisible = 'off';
h10.Title = 'Saline (mean)';

% heatmap Pre-CPP
figure(11)
h11 = heatmap(mean_matrix_Pre);
h11.XDisplayLabels = ROI_labels;
h11.YDisplayLabels = ROI_labels;
h11.GridVisible = 'off';
h11.Title = 'Pre-CPP (mean)';

% heatmap Post-CPP
figure(12)
h12 = heatmap(mean_matrix_Pos);
h12.XDisplayLabels = ROI_labels;
h12.YDisplayLabels = ROI_labels;
h12.GridVisible = 'off';
h12.Title = 'Post-CPP (mean)';









