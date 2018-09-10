function fileCell = textFile2cell(inFile)
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

% Open file
fid = fopen(inFile);

% Read lines of file
C = textscan(fid,'%s');

% extract cell of lines
fileCell = C{1};

% Close file
fclose(fid);
end