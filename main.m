clear all
close all

%% specify desired pattern
pattern = 'testPattern';

%% retrive model and optimisation settings
litoset = lito_settings( pattern );

%% pattern-specifics
% load the desired pattern
addpath(genpath('./pattern'))
load(pattern)

litoset.Zm = data;
litoset.xlim = xlim; % grid size (µm)
litoset.ylim = ylim;

% subdomain division
litoset.Nsdx = 2; 
litoset.Nsdy = 2;

% penalty parameter
litoset.gamma = 0.0001;

% initial value
litoset.wwInit = 0.3;

%% call the computational routine
tic;
result = lithography( litoset );
result.time = toc;

%% plot
plotter(litoset, result);

%% L2 norm
rms_error = getRMS( litoset, result );
fprintf('\nRMS-error: %5.10f\n',  rms_error)