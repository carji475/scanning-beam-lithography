function litoset = lito_settings( varargin )

% pattern-independent settings

%% Physical model parameters
litoset.tr = 0.5; % Threshold value
litoset.a  = 5; % Steepness of sigmoid used for threshold
% set the bandwidth according to 'Iterative deconvolution for inverse
% maskless lithography'
wx0 = 570e-3;
wy0 = 560e-3;
phi = 2.2*pi/180;
h11 = cos(phi)^2/wx0^2+sin(phi)^2/wy0^2;
h12 = -0.5*sin(2*phi)/wx0^2+0.5*sin(2*phi)/wy0^2;
h22 = sin(phi)^2/wx0^2+cos(phi)^2/wy0^2;
litoset.H = (4*[h11 h12; h12 h22])\eye(2); % bandwidth, B=K( H^(-1/2)*[x; y]) )

%% effective kernel support
litoset.tau = 10;

%% optimisation parameters
litoset.lambda = 0.0; % exposure time penalty term
litoset.gamma  = 0.1; % dosage penalty term
litoset.wwInit = 1; % initial exposure

%% reparametrisation
litoset.par.repar = 'none'; % reparametrisation (none / pos / possq / posb), see ww_reparametrised.m
litoset.par.wmax = 2; % upper bound when repar='posb'

%% lbfgsb parameters
litoset.memoryBound = 1; % memory limit
litoset.maxIter     = 2e2; % maximum nr of iterations
litoset.maxTotIter  = 5*litoset.maxIter; % including line search
litoset.printFreq   = 5; % printing every printFreq:th iteration
litoset.projGradTol = 1e-8; % infinity norm of proj gradient
litoset.factr       = 0; % rel red in f <= factr*machine_eps
