function [analysisParams] = loadHRF(analysisParams)
% Loads a hemodynamic response function for a subject.
%
% Syntax:
%    [analysisParams] = loadHRF(analysisParams)
%
% Description:
%    Loads a hemodynamic response function created by the forwardModel gear
%    on flywheel for a specific subject. 
%
% Inputs:
%    analysisParams           - a struct with relevant meta data for
%                               analyzing the LFContrast data. 
%
% Outputs:
%    analysisParams           - the same struct with relevant meta data for
%                               analyzing the LFContrast data but updated 
%                               include a subfield for the HRF.
%                               - analysisParams.HRF.timebase 
%                               - analysisParams.HRF.values
% Optional key/value pairs:
%    - none for now

% MAB 12/22/19 created it

load(fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'),'LFContrastAnalysis','subjectHRFs',analysisParams.expSubjID,[analysisParams.expSubjID '_eventGain_results.mat']));
xBase = zeros(1,analysisParams.expLengthTR);
xBase(1:length(results.hrf')) = results.hrf';
analysisParams.HRF.values = xBase;
analysisParams.HRF.timebase =   analysisParams.timebase*1000;
scaleVal = trapz(analysisParams.HRF.timebase,analysisParams.HRF.values);
analysisParams.HRF.values = analysisParams.HRF.values./scaleVal;

end