% load nifti file, total scan
% in the imagesc figure row = Y, column = X, 3d = Z

%animals = {'zFC_001_pre' 'zFC_002_pre' 'zFC_003_pre' 'zFC_004_pre' 'zFC_005_pre' 'zFC_006_pre' 'zFC_007_pre' 'zFC_008_pre' 'zFC_009_pre' 'zFC_010_pre' 'zFC_011_pre' 'zFC_012_pre'};


animals = {'zFC_001_pos' 'zFC_002_pos' 'zFC_003_pos' 'zFC_004_pos' 'zFC_005_pos' 'zFC_006_pos' 'zFC_007_pos' 'zFC_008_pos' 'zFC_009_pos' 'zFC_010_pos' 'zFC_011_pos' 'zFC_012_pos'};

% ROIs to test
main_directory = 'C:\';

Results(1,:) = {animals};


for a = 1:length(animals)
        
    % load zFC nifti files obtained in RESTplus
    % in the imagesc figure row = Y, column = X, 3d = Z, 4d = scan
    filename = [main_directory animals{a}];
    tscan = niftiread(filename);
    numslices = length(tscan(1,1,:,1));
    
          
    % load mask for ROI .nii file
ROI_filename = 'C:\';
roimask = niftiread(ROI_filename);

% convert roimask uint8 format into int16 as tscan
% and multiply tscan by mask to obtain only tscan values inside the mask
roimask = cast(roimask,'like',tscan);
tscan_masked = tscan.*roimask;

% make matrix with highlighted values in the mask region
roimask_highlight = roimask + 1;
tscan_highlight = tscan.*roimask_highlight;

% choose slice to view ONLY
sliceview = 6;

figure (1)
imagesc(tscan(:,:,sliceview));
colorbar
title(['Slice ' num2str(sliceview)], 'Interpreter', 'none')

figure (2)
imagesc(tscan_highlight(:,:,sliceview));
colorbar
title('Mask Highlighted', 'Interpreter', 'none')



% CALCULATE AVERAGE zFC FOR ENTIRE BRAIN

% For whole brain no need to multiply by mask (is already masked)
% find the number of nonzero elements in all slices
NZ =nnz(tscan(:,:,:));
% The average is obtained by the sum of all slices / NZ
% prefill sumslice matrices
sumSlice = zeros(1,numslices);

% loop for all slices
 for k = 1:numslices
     sumSlice(k) = sum(sum(tscan(:,:,k)));
 end

    totalSlices = sum(sumSlice);
    meanBrain = totalSlices/NZ
      
    Results{a} = {meanBrain};
    
      
end






    
    
    
            