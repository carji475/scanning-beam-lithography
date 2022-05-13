clear all
close all
rng default

addpath(genpath('../'))

%% specify desired pattern
pattern = 'testPattern';

%% retrive model and optimisation settings
litoset = lito_settings( pattern );
litoset.projGradTol=0;
litoset.maxIter=50;

%% load pattern
load(pattern)

%% to loop over
nrloops = 3;
repars = {'none', 'none', 'pos'};
memBounds = [5 1 1];

%% loop
Ns = round(logspace(2.7,4,14));
nvars = zeros(length(Ns),nrloops); % tot. number of variables
times = zeros(length(Ns),nrloops); % optimisation time
errors = zeros(length(Ns),litoset.maxIter,nrloops); % cost function eval.
for q=1:length(Ns)
    for p=1:nrloops
        N=Ns(q);
        disp([p q N])  
        
        litoset.Zm = imresize(data, [N,N], 'nearest');
        litoset.xlim = xlim; % grid size (µm)
        litoset.ylim = ylim;
        
        % loop params
        litoset.memoryBound = memBounds(p);
        litoset.par.repar = repars{p};
        
        % subdomain division
        litoset.Nsdx = 1;
        litoset.Nsdy = 1;
        
        %% call the computational routine
        start=tic;
        result = lithography( litoset );
        times(q,p) = toc(start);
        
        nvars(q,p) = N*N;
        
        errors(q,:,p) = result.info{1}.err(:,1);
    end
end

%% plot
% time evaluation
figure;
fs=14;
markers = {'ro','c>','b^'};
hs={};
for p=1:nrloops
    hs{p}=loglog(nvars(:,p), times(:,p), markers{p}, 'linewidth',2); hold on
end
th=nvars(:,1).*log(nvars(:,1));
hs{end+1}=loglog(nvars(:,1),th/th(1)*times(1)*1.5,'-k','linewidth',1.5);
hs{end+1}=loglog(nvars(:,1),nvars(:,1)/nvars(1)*times(1)*1.5,'--k','linewidth',1.5);
grid on
xlabel('$N$', 'interpreter', 'latex','FontSize',14)
ylabel('Time [s]', 'interpreter', 'latex','FontSize',14)
set(hs{1},'DisplayName','L-BFGS-B ($m=5$)')
set(hs{2},'DisplayName','GD (bounded)')
set(hs{3},'DisplayName','GD (var. trans.)')
set(hs{end-1},'DisplayName','$\propto N\log N$')
set(hs{end},'DisplayName','$\propto N$')
leg=legend([hs{:}]);
set(leg,'interpreter','latex')
set(leg,'FontSize',fs)
set(leg,'Location','NorthWest')
axis tight
saveas(gcf,'opti_times.epsc')

% cost evaluation
figure;
styles = {'r-','c--','b-.'};
hs={};
for p=1:nrloops
    hs{p} = semilogy(1:litoset.maxIter, errors(end,:,p),styles{p}, 'linewidth',2); hold on
end
set(hs{1},'DisplayName','L-BFGS-B ($m=5$)')
set(hs{2},'DisplayName','GD (bounded)')
set(hs{3},'DisplayName','GD (var. trans.)')
leg=legend([hs{:}]);
set(leg,'interpreter','latex')
set(leg,'FontSize',fs)
set(leg,'Location','NorthEast')
axis tight
ylabel('Cost', 'interpreter', 'latex','FontSize',14)
xlabel('Iteration', 'interpreter', 'latex','FontSize',14)
grid on
saveas(gcf,'opti_cost.epsc')