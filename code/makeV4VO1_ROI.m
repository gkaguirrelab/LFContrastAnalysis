% initialize
subjId = 'LZ23';
analysisParams = getSubjectParams(subjId);
% load
mapNameV4  = ['/Users/mbarnett/labDropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/wangAtlas/dtseries/' subjId '_V4.dtseries.nii'];
mapNameVO1 = ['/Users/mbarnett/labDropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/wangAtlas/dtseries/' subjId '_VO1.dtseries.nii'];

mapValsV4 = loadCIFTI(mapNameV4);
mapValsVO1 = loadCIFTI(mapNameVO1);
% theshold due to linear interp step in mapping process
map_v4 = zeros(size(mapNameV4));
map_vo1 = zeros(size(mapNameVO1));

map_v4(mapValsV4 >= .5 ) = 1;
map_v4(mapValsV4 < .5  ) = 0;
map_vo1(mapValsVO1 >= .5 ) = 1;
map_vo1(mapValsVO1 < .5 ) = 0;

% save
dropBoxPath     = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.projectName);
templateFile   = fullfile(dropBoxPath,'surfaceMaps','templates','template.dscalar.nii');
outMapNameV4   = ['/Users/mbarnett/labDropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/' subjId '_wang_V4.dscalar.nii'];
outMapNameVO1  = ['/Users/mbarnett/labDropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/' subjId '_wang_VO1.dscalar.nii'];
makeWholeBrainMap(map_v4, [], templateFile, outMapNameV4)
makeWholeBrainMap(map_vo1, [], templateFile, outMapNameVO1)