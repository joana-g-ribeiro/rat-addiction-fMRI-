% Extracting timeseries directly from smoothed scans .nii
% MIND THE CLEAR below
clear

% load swscan nifti file obtained in SPM, total scan
% in the imagesc figure row = Y, column = X, 3d = Z, 4d = scan
filename = 'C:\';
filenamefilt = 'C:\';
tscan = niftiread(filename);
tscanfilt = niftiread(filenamefilt);

% load mask for ROI .nii file and name the ROI
roimask = niftiread('C:\');
ROIname = ('');

% choose slice and volume/scan to view ONLY
sliceview = 5;
volview = 1;
%set window of timeseries to show
window = 1:;

numslices = length(tscan(1,1,:,1));
numvol = length(tscan(1,1,1,:));

figure (1)
imagesc(tscan(:,:,sliceview,volview));
colorbar
title(['Slice' ' ' num2str(sliceview) '/' num2str(numslices) '   ' 'Volume' ' ' num2str(volview) '/' num2str(numvol)])

% converts roimask uint8 format into int16 as tscan
% also create roimask for filtered into single as tscanfilt
roimask = cast(roimask,'like',tscan);
tscan_masked = tscan.*roimask;

roimaskfilt = cast(roimask, 'like', tscanfilt);
tscan_masked_filt = tscanfilt.*roimaskfilt;

% make matrix with highlighted values in the mask region
roimask_highlight = roimask + 1;
tscan_highlight = tscan.*roimask_highlight;

figure (2)
imagesc(tscan_highlight(:,:,sliceview,volview));
colorbar
title(ROIname,'Interpreter','none')


% find the number of nonzero elements in volume 1
NZ =nnz(tscan_masked(:,:,:,1));

% Create timeseries for current mask
% Each volume value is obtained by the sum of all slices / NZ
sumSlice = zeros(1,numslices);
Bold = zeros(1,numvol);
Bold_filt = zeros(1,numvol);

for x = 1: numvol
    
    for k = 1:numslices
    sumSlice(k) = sum(sum(tscan_masked(:,:,k,x)));
    end
    totalVolume = sum(sumSlice);
    Bold(x) = totalVolume/NZ;
   
end


for x = 1: numvol
    
    for k = 1:numslices
    sumSlice(k) = sum(sum(tscan_masked_filt(:,:,k,x)));
    end
    totalVolume = sum(sumSlice);
    Bold_filt(x) = totalVolume/NZ;
   
end



scans = 1:numvol;

figure (3)
title([filename '       ' ROIname],'Interpreter','none')
xlabel('Scans')
yyaxis left
plot(scans(window),Bold(window))
ylabel('BOLD')

figure (4)
title([filename '       ' ROIname],'Interpreter','none')
xlabel('Scans')
yyaxis left
plot(scans(window),Bold_filt(window))
ylabel('BOLD filtered')




