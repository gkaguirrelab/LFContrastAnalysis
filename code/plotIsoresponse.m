function [] = plotIsoresponse(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit)
% Takes in a text file name and retuns a cell of the lines of the text file
%
% Syntax:
%   filesCell = textFile2cell(inFile)
%
% Description:
%    This function takes in a file name for a trext file and returns a cell
%    that is composed of the lines of the text file. Example of this would
%    be a text file of file names the output is a cell of files names.
%
% Inputs:
%    inFile            - File name of a text file. (string)
%
% Outputs:
%    fileCell          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    none

% MAB 09/09/18
%fix to match upsampling 


    
LminusMbetas = semIAMPBetas(1:5)- semIAMPBetas(end-1); 
LplusMbetas = semIAMPBetas(1:5)- semIAMPBetas(end-1); 
LIsoBetas = semIAMPBetas(1:5)- semIAMPBetas(end-1); 
MIsoBetas = semIAMPBetas(1:5)- semIAMPBetas(end-1); 

LminusMcontrast = contrastCoding.*maxContrastPerDir(1);
LplusMcontrast= contrastCoding.*maxContrastPerDir(2);
LIsocontrast= contrastCoding.*maxContrastPerDir(3);
MIsocontrast= contrastCoding.*maxContrastPerDir(4);

IAMPBetas = {LminusMbetas,LplusMbetas,LIsoBetas,MIsoBetas};
contrastLevels = {LminusMcontrast,LplusMcontrast,LIsocontrast,MIsocontrast};
directions = {[1,-1],[1,1],[1,0],[0,1]};
thresh = 0.25;
hdl = plotIsorespContour(paramsQCMFit,IAMPBetas,contrastLevels,directions,thresh,[],'r');
thresh = 0.5;
hdl = plotIsorespContour(paramsQCMFit,IAMPBetas,contrastLevels,directions,thresh,hdl,'g');
thresh = 0.75;
hdl = plotIsorespContour(paramsQCMFit,IAMPBetas,contrastLevels,directions,thresh,hdl,'b');
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);