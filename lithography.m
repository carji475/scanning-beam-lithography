function result = lithography( litoset )
%{
    In this version we combine the LBFGSB algorithm (Stephen Becker's code)
    with FFT for calculating the dosage/gradient.
%}

%% add path
addpath(genpath('./lbfgsb/l-bfgs-b-c-master/Matlab')); % path to lbfgsb files

%% Simulation parameters
Nxm = size(litoset.Zm,1);
Nym = size(litoset.Zm,2);
x  = linspace(0,litoset.xlim,Nxm); % Spatial dimension
y  = linspace(0,litoset.ylim,Nym);
xres = x(2)-x(1); % x grid resolution
yres = y(2)-y(1); % y grid resolution

% H^(-1/2)
Hsqn = sqrtm(litoset.H)\eye(2);

%% set the subdomains
support_x = ceil(litoset.tau*sqrt(max(eig(litoset.H)))/xres);
support_y = ceil(litoset.tau*sqrt(max(eig(litoset.H)))/yres);
Nxstd = ceil(Nxm/litoset.Nsdx);
Nxext = min(support_x, Nxm-Nxstd);
Nystd = ceil(Nym/litoset.Nsdy);
Nyext = min(support_y, Nym-Nystd);

% sanity check
if Nxstd*(litoset.Nsdx-1)>Nxm
    error('Too many x-domains!')
elseif Nystd*(litoset.Nsdy-1)>Nym
    error('Too many y-domains!')
end

%% allocate all variables
wwTot = zeros(Nxm,Nym);

%% loop
for nnn = 1:litoset.Nsdx*litoset.Nsdy
    fprintf('\n========== Subdomain %d ==========\n', nnn)
    
    % subdomain divsion indices
    [Nx, ZxInd, wwxInd, wwTotxInd, Ny, ZyInd, wwyInd, wwTotyInd] =...
        subdomain(nnn, Nxm, Nxstd, Nxext, Nym, Nystd, Nyext, litoset);
    
    % indices specifying what to extract in this subdomain
    ZxmaxInd = min(ZxInd+Nx-1,size(litoset.Zm,1));
    ZxminInd = max(ZxInd,1);
    ZymaxInd = min(ZyInd+Ny-1,size(litoset.Zm,2));
    ZyminInd = max(ZyInd,1);
    
    % effective numbers
    Nxeff = ZxmaxInd - ZxminInd + 1;
    Nyeff = ZymaxInd - ZyminInd + 1;
    
    %% fft initialisation following the Gramacki brothers
    fft_exposure(Nxeff, Nyeff, support_x, support_y , xres, yres, Hsqn); % initilaise
    
    %% Define Feature
    Z = litoset.Zm(ZxminInd:ZxmaxInd, ZyminInd:ZymaxInd);
    
    %% reduce size for lbfgsb
    indices = Z(:)>1e-4; % only variables within desired pattern are optimised for
    nVar = sum(indices);
    ww = litoset.wwInit * ones(nVar,1); % initial guess is nonzero for any x,y's inside the desired pattern
    
    %% normalise the cost function
    litoset.scale = 1/nVar;
    
    %% set limits and defaults
    if strcmp(litoset.par.repar,'none')
        lb = 0; % lower bound
    else
        lb=-inf;
    end
    ub = inf; % upper bound
    l = lb*ones(nVar,1); u = ub*ones(nVar,1);
   
    if strcmp(litoset.par.repar,'posb') && litoset.wwInit>0.5*litoset.par.wmax
        % if bounded, we set the initial value to half of the max value
        ww = 0.5*litoset.par.wmax * ones(size(ww));
    end
    if isempty(ww) % go to next loop if there are no variables to optimise in this region
        continue
    end
    
    %% initialise peristenet variables
    exposureFunction_lbfgs_fft(litoset.a, litoset.tr, Nxeff, Nyeff, ...
        litoset.scale, litoset.lambda, litoset.gamma, litoset.par);
    
    %% LBFGSB
    opts = struct( 'x0', ww_reparametrised(ww, litoset.par, 1) ); % set start guess
    opts.printEvery = litoset.printFreq; % controls how often we print output from .m file wrapper
    opts.m  = litoset.memoryBound; % memory bound
    opts.maxIts = litoset.maxIter; % max number of iterations
    opts.maxTotalIts = litoset.maxTotIter; % including line search
    opts.pgtol = litoset.projGradTol; % infinity norm
    opts.factr = litoset.factr; % rel red in f <= factr*machine_eps
    
    start=tic; % measure time
    [ww,~,info] = lbfgsb( @(ww) exposureFunction_lbfgs_fft(ww,Z,indices), l, u, opts ); % optimise
    stop=toc(start);
    fprintf('Time: %5.5f\tTime/iteration: %5.5f\n' , [stop, stop/info.iterations])
    
    % save opimisation info
    result.info{nnn} = info;
    
    % reshape variables to rectangular grid
    ww = ww_reparametrised(ww, litoset.par);
    ww = squareUp(ww,indices,Nxeff,Nyeff);
    
    % place at correct place in total grid
    wwTot( max(wwTotxInd, 1):min(wwTotxInd+Nxstd-1, Nxm),...
        max(wwTotyInd, 1):min(wwTotyInd+Nystd-1, Nym) ) =...
        ww(max(wwxInd,1):min(wwxInd+Nxstd-1,Nxm), max(wwyInd,1):min(wwyInd+Nystd-1,Nym));
end

%% add total solution to result structure
result.wwTot=wwTot;

%% include some more variables in the result structure
result.Nxm = Nxm;
result.Nym = Nym;
result.xres = xres;
result.yres = yres;
result.x = x;
result.y = y;
result.support_x = support_x;
result.support_y = support_y;
result.Hsqn = Hsqn;

%% remove path
rmpath(genpath('./lbfgsb/l-bfgs-b-c-master/Matlab'));
end

% function that returns subdomain indices
function [Nx, ZxInd, wwxInd, wwTotxInd, Ny, ZyInd, wwyInd, wwTotyInd] =...
    subdomain(nnn, Nxm, Nxstd, Nxext, Nym, Nystd, Nyext, litoset)
if litoset.Nsdx==1
    Nx = Nxm;
    ZxInd     = 1;
    wwxInd    = 1;
    wwTotxInd = 1;
else
    if nnn<=litoset.Nsdy % top grid row
        Nx = Nxstd+Nxext;
        ZxInd     = 1;
        wwTotxInd = 1;
    elseif nnn>litoset.Nsdy*(litoset.Nsdx-1) % bottom grid row
        Nx = Nxstd+Nxext;
        ZxInd     = Nxm-Nx+1;
        wwTotxInd = Nxm-Nxstd+1;
    else
        theRem = floor((nnn-0.5)/litoset.Nsdy)+1;
        Nx = Nxstd+2*Nxext;
        ZxInd     = (theRem-1)*Nxstd-Nxext;
        wwTotxInd = (theRem-1)*Nxstd+1;
    end
    ZxInd = max(ZxInd,1);
    wwxInd = wwTotxInd-ZxInd+1;
end
if litoset.Nsdy==1
    Ny = Nym;
    ZyInd     = 1;
    wwyInd    = 1;
    wwTotyInd = 1;
else
    if rem(nnn-1,litoset.Nsdy)==0 % left grid col
        Ny = Nystd+Nyext;
        ZyInd     = 1;
        wwTotyInd = 1;
    elseif rem(nnn,litoset.Nsdy)==0 % right grid col
        Ny = Nystd+Nyext;
        ZyInd     = Nym-Ny+1;
        wwTotyInd = Nym-Nystd+1;
    else
        theRem=rem(nnn,litoset.Nsdy);
        Ny = Nystd+2*Nyext;
        ZyInd    = (theRem-1)*Nystd-Nyext;
        wwTotyInd = (theRem-1)*Nystd+1;
    end
    ZyInd = max(ZyInd,1);
    wwyInd = wwTotyInd-ZyInd+1;
end
end
