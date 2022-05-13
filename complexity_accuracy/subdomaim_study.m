clear all
close all
rng default

addpath(genpath('../'))

%% specify desired pattern
pattern = 'testPattern';

%% retrive model and optimisation settings
litoset = lito_settings( pattern );
litoset.maxIter=5e1;

%% pattern-specifics
% load the desired pattern
load(pattern)
litoset.Zm = data;
litoset.xlim = xlim; % grid size (µm)
litoset.ylim = ylim;

%% subdomain division
Nsdms = [1, 1:5];
litoset.tau = 1e10; % exact FFT for the first computation
rms_values = zeros(length(Nsdms),1);
times = zeros(length(Nsdms),1);
for q=1:length(Nsdms)
    disp(q) 
    Nsdm = Nsdms(q);
    if q>1
        litoset.tau=10;
    end
    litoset.Nsdx = Nsdm; 
    litoset.Nsdy = Nsdm;

    %% call the computational routine
    tic;
    result = lithography( litoset );
    result.time = toc;

    %% L2 norm
    rms_values(q) = getRMS( litoset, result );
    times(q) = result.time;
end

%% plot
figure;
yyaxis left
plot(1:length(Nsdms), rms_values/rms_values(1), 'bo', 'linewidth', 2)
ylabel('RMS error','interpreter','latex')
grid on
yyaxis right
plot(1:length(Nsdms), times/times(1), '*', 'linewidth', 2)
ylabel('Time','interpreter','latex')
ax=gca;
ax.XTick = 1:length(Nsdms);
for q=1:length(Nsdms)
    ax.XTickLabel{q}=['$' num2str(Nsdms(q)) '\times ' num2str(Nsdms(q)) '$'];
end
ax.XAxis.TickLabelInterpreter = 'latex';
xlabel('Subdomain division','interpreter','latex')