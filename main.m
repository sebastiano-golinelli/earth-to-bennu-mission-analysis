%% Earth -> Bennu mission analysis
% Preliminary design of a mission from a geostationary transfer orbit to
% the near-Earth asteroid 101955 Bennu, in three steps:
%   1) geocentric transfer from the launch orbit to a parking orbit,
%   2) direct two-impulse heliocentric transfer Earth -> Bennu,
%   3) patched conics: Earth escape and capture at Bennu.
% Edit config/missionConfig.m to change data and options.

clear; close all; clc;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'src')));

cfg = missionConfig();
results = runMission(cfg);
