basedir = pwd;
addpath([basedir '/tools_mexcode/']);
addpath([basedir '/tools_splines/']);
addpath([basedir '/tools_misc/']);
addpath([basedir '/GLMcode/']);
addpath([basedir '/GLMcode/nlfuns/']);
addpath([basedir '/MIDCode/']);
addpath([basedir '/MIDCode/mexcode/']);

global RefreshRate;  % Stimulus frame rate (Hz)

initialize_mexcode;