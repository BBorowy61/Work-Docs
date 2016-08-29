%% Import data from text file.
% Bogdan Borowy
%% Initialize variables.
filename = 'C:\Users\Bogdan.Borowy\Documents\McHenry\Performance Analysis\PerfTest.csv';
delimiter = ',';
startRow = 2;

%% Format string for each line of text:
%   column1: double (%f)
%	column2: datetimes (%{HH:mm:ss}D)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%{HH:mm:ss}D%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Interval = dataArray{:, 1};
TIME = dataArray{:, 2};
AGCMW = dataArray{:, 3};
RawECOBP = dataArray{:, 4};
EnergyRRatAGCMWMWmin = dataArray{:, 5};
RRMW10sec = dataArray{:, 6};
RampedEcoBP = dataArray{:, 7};
FleetTREG = dataArray{:, 8};
FleetRegulationControlSignal = dataArray{:, 9};

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TIME=datenum(TIME);


%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;