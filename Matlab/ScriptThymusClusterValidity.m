% ScriptThymusClusterValidity

% Read inputs (and undersample).

dataPath ='~/workspace/Preprocessing/ThymusPreprocessingShadingCorrection.04222010/';

% subjectID = '20091208_111335MRI-06014';
% NS = load_nii([dataPath subjectID '/20091208_111335WIPeTHRIVENOSUPPSENSEMRI-06014s501a1005_cropped_n3_corrected.nii']);
% WS = load_nii([dataPath subjectID '/20091208_111335WIPeTHRIVEWSSENSEMRI-06014s601a1006_cropped_n3_corrected_NS_aligned.nii']);
% FS = load_nii([dataPath subjectID '/20091208_111335WIPeTHRIVEFSSENSEMRI-06014s701a1007_cropped_n3_corrected_NS_aligned.nii']);

% subjectID = '20091207_094748MRI-089-01';
% NS = load_nii([dataPath subjectID '/20091207_094748WIPeTHRIVENOSUPPSENSEMRI-089-01s401a1004_cropped_n3_corrected.nii']);
% WS = load_nii([dataPath subjectID '/20091207_094748WIPeTHRIVEWSSENSEMRI-089-01s501a1005_cropped_n3_corrected_NS_aligned.nii']);
% FS = load_nii([dataPath subjectID '/20091207_094748WIPeTHRIVEFSSENSEMRI-089-01s601a1006_cropped_n3_corrected_NS_aligned.nii']);

subjectID = '20091214_143047MRI-09501';
NS = load_nii([dataPath subjectID '/20091214_143047WIPeTHRIVENOSUPPSENSEMRI-09501s901a1009_cropped_n3_corrected.nii']);
WS = load_nii([dataPath subjectID '/20091214_143047WIPeTHRIVEWSSENSEMRI-09501s701a1007_cropped_n3_corrected_NS_aligned.nii']);
FS = load_nii([dataPath subjectID '/20091214_143047WIPeTHRIVEFSSENSEMRI-09501s801a1008_cropped_n3_corrected_NS_aligned.nii']);


NS = NS.img;
NS = NS(:,:,25:30);
NS_vector = NS(:);

WS = WS.img;
WS = WS(:,:,25:30);
WS_vector = WS(:);

FS =  FS.img;
FS = FS(:,:,25:30);
FS_vector = FS(:);

DM = double([NS_vector, FS_vector, WS_vector]);

figure, imagesc( NS(:,:,3)'), colormap(gray), colorbar, axis image
saveas(gcf, ['NS_original.png']);

% Run clustering and validity metrics.

maxNClusters = 10;

% D = 3
infoString = [sprintf('3D case\n')];
infoString = [infoString  sprintf('# of clusters, PC, CE, SC, S, XB \n')];

for nClusters=2:maxNClusters
    
    result3D{nClusters} = ThymusClusterValidity(DM, nClusters);
    
    % Display original and the label images.
    [~, IDX] =  max(result3D{nClusters}.data.f');
    L = reshape( IDX, size(NS) );
    
    figure, scatter3(DM(:,1), DM(:, 2), DM(:, 3), 4, IDX), colormap(jet),
    xlabel('Non-suppressed'), ylabel('Fat-suppressed'), zlabel('Water-suppressed')
    saveas(gcf, ['scatterplot', '_', num2str(nClusters), '_', num2str(size(DM,2)), 'D', '.fig']);
    
    figure, imagesc( L(:,:,3)'), colormap(jet), colorbar, axis image
    saveas(gcf, ['labeled_image', '_', num2str(nClusters), '_', num2str(size(DM,2)), 'D', '.png']);
    
    infoString = [infoString  sprintf( '%d, %d, %d, %d, %d, %d\n', ...
        nClusters, ...
        result3D{nClusters}.validity.PC, result3D{nClusters}.validity.CE, ...
        result3D{nClusters}.validity.SC, result3D{nClusters}.validity.S, ...
        result3D{nClusters}.validity.XB)];
    
end


% D = 2

DM = DM(:, 2:3);
infoString = [infoString sprintf('2D case\n')];
infoString = [infoString  sprintf('# of clusters, PC, CE, SC, S, XB \n')];

for nClusters=2:maxNClusters
    
    result2D{nClusters} = ThymusClusterValidity(DM, nClusters);
    
    % Display original and the label images.
    [~, IDX] =  max(result2D{nClusters}.data.f');
    L = reshape( IDX, size(NS) );
    
    figure, scatter(DM(:,1), DM(:, 2), 4, IDX), colormap(jet),
    xlabel('Fat-suppressed'), ylabel('Water-suppressed')
    saveas(gcf, ['scatterplot', '_', num2str(nClusters), '_', num2str(size(DM,2)), 'D', '.fig']);
    
    figure, imagesc( L(:,:,3)'), colormap(jet), colorbar, axis image
    saveas(gcf, ['labeled_image', '_', num2str(nClusters), '_', num2str(size(DM,2)), 'D', '.png']);
    
    infoString = [infoString  sprintf( '%d, %d, %d, %d, %d, %d\n', ...
        nClusters, ...
        result2D{nClusters}.validity.PC, result2D{nClusters}.validity.CE, ...
        result2D{nClusters}.validity.SC, result2D{nClusters}.validity.S, ...
        result2D{nClusters}.validity.XB)];
    
end


% D = 1

DM = DM(:, 2);
infoString = [infoString sprintf('1D case\n')];
infoString = [infoString  sprintf('# of clusters, PC, CE, SC, S, XB \n')];

for nClusters=2:maxNClusters
    
    result1D{nClusters} = ThymusClusterValidity(DM, nClusters);
    
    % Display original and the label images.
    [~, IDX] =  max(result1D{nClusters}.data.f');
    L = reshape( IDX, size(NS) );
    
%     figure, scatter(DM(:,1), 4, IDX), colormap(jet),
%     xlabel('Water-suppressed')
%     saveas(gcf, ['scatterplot', '_', num2str(nClusters), '_', num2str(size(DM,2)), 'D', '.png']);
    
    figure, imagesc( L(:,:,3)'), colormap(jet), colorbar, axis image
    saveas(gcf, ['labeled_image', '_', num2str(nClusters), '_', num2str(size(DM,2)), 'D', '.png']);
    
    infoString = [infoString  sprintf( '%d, %d, %d, %d, %d, %d\n', ...
        nClusters, ...
        result1D{nClusters}.validity.PC, result1D{nClusters}.validity.CE, ...
        result1D{nClusters}.validity.SC, result1D{nClusters}.validity.S, ...
        result1D{nClusters}.validity.XB)];
    
end


% Store voxel counts in csv format.
imagefileprefix = [ subjectID, '_', 'ClusterValidation_'];

fid=fopen([imagefileprefix '_ThymusQuantification.csv'], 'wb', 'l');
fprintf(fid, infoString);
fclose(fid);
