% Scatter for t-test and Scatter for Correlations with CPP_score

% MIND THE CLEAR below
clear

% Load results file in the current folder
% MIND to run correlation stats first to update Results file
load('Results.mat')
animals = Results{5,2};
treat_labels = Results(1:4,1);
ROI_labels = Results{1,2}{1}{1}.Var1;
numvol = length(Results{1,2}{1}{1}.all_ROI_Bold);
% Remember variables/tables to index are: labeled_all_ROI_timeseries, labeled_correlations, labeled_p_values

%%%%%%%%%%%%%%%%%%%%% CHOOSE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remember files p values order on Stats file is:
% [column factor (i.e. treatment), row factor (i.e. time), interaction]

% Select pairs of ROIs to scatter from ROI_labels
ROI_pairs = [8 3; 10 3; 4 3; 8 4];

% indicate CPP_score for each animal by the same order of animal list
%%CPP_score = [36.45404162, -6.873307148, 8.777353458, -3.563133323, 4.804733728, -27.1401917, -14.48187159, 6.21526452, 24.68508093, -9.997979938, 8.953079018, 89.72724838];

CPP_score = [66, -2, 62, 41, -29, -159, -7, -27, 113, -174, 30, 275];

% Scatter options t-test
% choose font size for axis labels
fontsize = 12;
% choose font weight 'normal' or 'bold'
fontweight = 'normal';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract data like this:
%  = Results{treatment index,2}{animal index}{variable index}.variable name;
% Example: extract all_ROI_Bold matrix from table "labeled_all_ROI_timeseries" (index 1), 4th treatment, 3rd animal:
%  = Results{4,2}{3}{1}.all_ROI_Bold;

% Create cell array with selected ROI pairs
ROI_pair_Labels = cell(1,size(ROI_pairs,1));
for u = 1: size(ROI_pairs,1)
    ROI_pair_Labels{u} = [ROI_labels{ROI_pairs(u,1)} ' - ' ROI_labels{ROI_pairs(u,2)}];
end

% preallocate cell array to store all ROI pair Stats
% loop all ROI pair Labels
ROI_pair_ScatterStats = cell(length(ROI_pair_Labels),2);
for r = 1:length(ROI_pair_Labels)
    ROI_pair_ScatterStats(r,:) = {ROI_pair_Labels(r),[1 1 1]};
end

% preallocate cell array to store comparisons (Morphine vs Saline) and
% (Pre- vs Post-CPP)
ScatterStats(1,:) = {'Morphine vs Saline',ROI_pair_ScatterStats};
ScatterStats(2,:) = {'Pre- vs Post-CPP',ROI_pair_ScatterStats};

% preallocale vectors for animals for scatter and t-test
vector1 = 1:length(animals);
vector2 = 1:length(animals);

% Scatter and t-test for Pre-CPP vs Post-CPP and Correlation with CPP_score

% loop through all ROI_pairs
for q=1:length(ROI_pair_Labels)
    % use ROI_pairs to index the position in correlation matrix
    % loop through all animals for Pre-CPP (treatment 1) and Post-CPP
    % (treatment 2)
    for a = 1:length(animals)
        % Pre-CPP
        matrix1 = Results{1,2}{a}{2}.R;
        vector1(a) = matrix1(ROI_pairs(q,1), ROI_pairs(q,2));
        % Post-CPP
        matrix2 = Results{2,2}{a}{2}.R;
        vector2(a) = matrix2(ROI_pairs(q,1), ROI_pairs(q,2));
    end
    % t-test now between Pre-CPP and Post-CPP
    [~,p_ttest] = ttest(vector1,vector2,'Alpha',0.05);
    % scatter now for Pre-CPP and Post-CPP
    x(:,1) = vector1';
    x(:,2) = vector2';
    
    %%nCats = 2;
    %%nDatas = 12;
    %%b = animals;
    %%f = cellstr(b);
    figure(1000+q);
    boxplot(x, 'BoxStyle', 'outline')
    set(gca,'XTickLabel',{'Pre-CPP'; 'Post-CPP'},'fontsize',fontsize,'FontWeight',fontweight')
    title(ROI_pair_Labels{q},'Interpreter', 'none')
    ylabel('Correlation Coeficient (z-score)')
    hold on
    scatter(linspace(0.9999,1.0001,size(x,1)), x(:,1))
    scatter(linspace(1.9999,2.0001,size(x,1)), x(:,2))
    %%line(repmat([(1:nCats).';NaN], [nDatas,1]), reshape(x(1:nDatas,[1:nCats, 1]).', [], 1), 'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10);
    %%uit = uitable('Data', x(:,1), 'Position', [440 80 100 280]);

    % coordinates for text box [x_begin y_begin length height]
    % By default, the units are normalized to the figure 0 to 1
    annotation('textbox',[0.43 0.5 0.3 0.3],'String',['  p = ', num2str(p_ttest)],'FitBoxToText','on');
    hold off
    % CPP_score vs ROI_pair correlation coeficients
    [R1,P1] = corrcoef(CPP_score, vector1);
    [R2,P2] = corrcoef(CPP_score, vector2);
    % Plot or scatter
    figure(2000+q);
    scatter(CPP_score, vector1)
    lsline
    hold on
    scatter(CPP_score, vector2)
    lsline
    %legend('Pre-CPP', '', 'Post-CPP', 'Location', 'northeastoutside')
    annotation('textbox',[0.7 0.6 0.3 0.3],'String',{'Pre-CPP', ['R = ', num2str(R1(1,2))], ['P = ', num2str(P1(1,2))]},'FitBoxToText','on', 'Color',  [0 0.4470 0.7410], 'EdgeColor', 'none');
    annotation('textbox',[0.7 0.45 0.3 0.3],'String',{'Post-CPP', ['R = ', num2str(R2(1,2))], ['P = ', num2str(P2(1,2))]},'FitBoxToText','on', 'Color',  [0.8500 0.3250 0.0980], 'EdgeColor', 'none');
    title(ROI_pair_Labels{q},'Interpreter', 'none')
    ylabel('Correlation Coeficient (CPP score vs z-score)')
    xlabel('CPP score')
    hold off
    % store current vectors of z-scores for all animals Pre-CPP (1st row) and Post-CPP (2nd row)
    %ROI_pair_ScatterStats{q,2}(1,:) = vector1;
    %ROI_pair_ScatterStats{q,2}(2,:) = vector2;
    
end

% store updated ROI_pair_ScatterStats into ScatterStats cell array Pre-CPP vs Post-CPP on row 2
ScatterStats{2,2} = ROI_pair_ScatterStats;

% Scatter and t-test for Morphine vs Saline and Correlation with CPP_score
%%
% loop through all ROI_pairs
for q=1:length(ROI_pair_Labels)
    % use ROI_pairs to index the position in correlation matrix
    % loop through all animals for Morphine (treatment 3) and Saline
    % (treatment 4)
    for a = 1:length(animals)
        % Morphine
        matrix1 = Results{3,2}{a}{2}.R;
        vector1(a) = matrix1(ROI_pairs(q,1), ROI_pairs(q,2));
        % Saline
        matrix2 = Results{4,2}{a}{2}.R;
        vector2(a) = matrix2(ROI_pairs(q,1), ROI_pairs(q,2));
    end
    % t-test now between Morphine and Saline
    [~,p_ttest] = ttest(vector1,vector2,'Alpha',0.05);
    % scatter now for Morphine and Saline
    x(:,1) = vector1';
    x(:,2) = vector2';
    figure(10000+q);
    boxplot(x, 'BoxStyle', 'outline')
    set(gca,'XTickLabel',{'Morphine'; 'Saline'},'fontsize',fontsize,'FontWeight',fontweight')
    title(ROI_pair_Labels{q},'Interpreter', 'none')
    ylabel('Correlation Coeficient (z-score)')
    hold on
    scatter(linspace(0.9999,1.0001,size(x,1)), x(:,1))
    scatter(linspace(1.9999,2.0001,size(x,1)), x(:,2))
    % coordinates for text box [x_begin y_begin length height]
    % By default, the units are normalized to the figure 0 to 1
    annotation('textbox',[0.43 0.5 0.3 0.3],'String',['  p = ', num2str(p_ttest)],'FitBoxToText','on');
    hold off
    % CPP_score vs ROI_pair correlation coeficients
    [R1,P1] = corrcoef(CPP_score, vector1);
    [R2,P2] = corrcoef(CPP_score, vector2);
    % Plot or scatter
    figure(20000+q);
    scatter(CPP_score, vector1)
    lsline
    hold on
    scatter(CPP_score, vector2)
    lsline
    %legend('Morphine', '', 'Saline', 'Location', 'northeastoutside')
    annotation('textbox',[0.7 0.6 0.3 0.3],'String',{'Morphine', ['R = ', num2str(R1(1,2))], ['P = ', num2str(P1(1,2))]},'FitBoxToText','on', 'Color',  [0 0.4470 0.7410], 'EdgeColor', 'none');
    annotation('textbox',[0.7 0.45 0.3 0.3],'String',{'Saline', ['R = ', num2str(R2(1,2))], ['P = ', num2str(P2(1,2))]},'FitBoxToText','on', 'Color',  [0.8500 0.3250 0.0980], 'EdgeColor', 'none');
    title(ROI_pair_Labels{q},'Interpreter', 'none')
    ylabel('Correlation Coeficient (CPP score vs z-score)')
    xlabel('CPP score')
    hold off
    % store current vectors of z-scores for all animals Morphine (1st row) and Saline (2nd row)
    %ROI_pair_ScatterStats{q,2}(1,:) = vector1;
    %ROI_pair_ScatterStats{q,2}(2,:) = vector2;
    
end

% store updated ROI_pair_ScatterStats into ScatterStats cell array Morphine
% vs Saline on row 1
ScatterStats{1,2} = ROI_pair_ScatterStats;

