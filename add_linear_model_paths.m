% add plant/ (model loader and step function used by the MUSCLE3 actor) and utils/
% (launchers and test tools) to the MATLAB path. The generator launchers in utils/ add
% generators/ themselves. Uses the location of this file as root, so the name of the
% cloned folder is irrelevant.

linear_model_root = fileparts(mfilename('fullpath'));
addpath(fullfile(linear_model_root,'plant'), fullfile(linear_model_root,'utils'));
clear linear_model_root
