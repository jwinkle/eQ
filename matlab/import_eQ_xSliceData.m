function xSlice = import_eQ_xSliceData(filename, dataLines)
%IMPORTFILE Import data from a text file
%  XSLICE = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the numeric data.
%
%  XSLICE = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  xSlice = importfile("/home/winkle/Link to images/1578439556-17_25_55-Jan__7_2020/xSlice.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 07-Jan-2020 21:30:05

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["u", "Var2", "Var3", "Points0", "Var5", "Var6"];
opts.SelectedVariableNames = ["u", "Points0"];
opts.VariableTypes = ["double", "string", "string", "double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var2", "Var3", "Var5", "Var6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var3", "Var5", "Var6"], "EmptyFieldRule", "auto");

% Import the data
xSlice = readtable(filename, opts);

%% Convert to output type
xSlice = table2array(xSlice);
end