function setup()
%SETUP Add the VRRK source code to the MATLAB path.

repositoryRoot=fileparts(mfilename('fullpath'));
sourceFolder=fullfile(repositoryRoot,'src');
resultsFolder=fullfile(repositoryRoot,'results');

addpath(repositoryRoot);
addpath(sourceFolder);

if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end

fprintf('VRRK-AS repository initialized.\n');
fprintf('Source directory: %s\n',sourceFolder);
fprintf('Run run_smoke_test to verify the installation.\n');
end
