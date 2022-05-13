clear all
close all
rng default

addpath(genpath('../'))

%% specify desired pattern
pattern = 'testPattern';

%% retrive model and optimisation settings
litoset = lito_settings( pattern );

%% load the desired pattern
load(pattern)

%% loop
usetaus=[0 1];
Ns = round(logspace(2.5,4,15));
nvars = zeros(length(Ns),1);
times = zeros(length(Ns),2);
Nrtimes=10;
for t=1:Nrtimes
    for q=1:length(Ns)
        N=Ns(q);
        litoset.Zm = imresize(data, [N,N], 'nearest');
        litoset.xlim = xlim; % grid size (µm)
        litoset.ylim = ylim;
        
        % subdomain division
        litoset.Nsdx = 1;
        litoset.Nsdy = 1;
        
        % Simulation parameters
        Nxm = size(litoset.Zm,2);
        Nym = size(litoset.Zm,1);
        x  = linspace(0,litoset.xlim,Nxm); % Spatial dimension
        y  = linspace(0,litoset.ylim,Nym);
        xres = x(2)-x(1); % x grid resolution
        yres = y(2)-y(1); % y grid resolution
        
        % H^(-1/2)
        Hsqn = sqrtm(litoset.H)\eye(2);
        
        % random exposure times
        ww = ones(size(litoset.Zm));
        for usetau=usetaus
            % reduce convolutional cost
            litoset.tau=10;
            if usetau
                support_x = ceil(litoset.tau*sqrt(max(eig(litoset.H)))/xres);
                support_y = ceil(litoset.tau*sqrt(max(eig(litoset.H)))/yres);
            else
                support_x = Nxm;
                support_y = Nym;
            end
            
            fft_exposure(Nxm, Nym, support_x, support_y , xres, yres, Hsqn); % initilaise
            
            par.repar='none';
            exposureFunction_lbfgs_fft(litoset.a, litoset.tr, Nxm, Nym, ...
                litoset.scale, litoset.lambda, litoset.gamma, par);
            
            start=tic;
            [c,g]=exposureFunction_lbfgs_fft(ww(:), litoset.Zm, ww(:)>=0);
            done=toc(start);
            disp(c)
            times(q,usetau+1) = times(q,usetau+1) + done/Nrtimes;
        end
        nvars(q) = numel(ww);
        
        disp([usetau q nvars(q)])
    end
end

%% plot
figure;
colours = {'b','r'};
markers = {'*','o'};
lines = {'--', '-'};
fs=14;
est_exact = polyfit(log(nvars),log(times(:,1)),1);
est_appr = polyfit(log(nvars),log(times(:,2)),1);
ests={est_exact, est_appr};
hs={};
for q=1:2
    hs{q}=loglog(nvars, times(:,q), [markers{q} colours{q}], 'linewidth',2); hold on
end
th=nvars.*log(nvars);
hs{end+1}=loglog(nvars,th/th(1)*times(1,2)*1.5,'-k','linewidth',1.5);
hs{end+1}=loglog(nvars,nvars/nvars(1)*times(1,2)*1.5,'--k','linewidth',1.5);
grid on
xlabel('$N$', 'interpreter', 'latex','FontSize',14)
ylabel('Time [s]', 'interpreter', 'latex','FontSize',14)
set(hs{1},'DisplayName','$\tau=\infty$')
set(hs{2},'DisplayName','$\tau=10$')
set(hs{end-1},'DisplayName','$\propto N\log N$')
set(hs{end},'DisplayName','$\propto N$')
leg=legend([hs{:}]);
set(leg,'interpreter','latex')
set(leg,'FontSize',fs)
set(leg,'Location','NorthWest')
axis tight
saveas(gcf,'cost_scaling.epsc')