function plotter(litoset, result)

x=result.x;
y=result.y;
Nx = length(x);
Ny = length(y);

% initialise
fft_exposure(result.Nxm, result.Nym, result.support_x, result.support_y, ...
    result.xres, result.yres, result.Hsqn);

% compute initial quantities
ww_i = zeros(Nx, Ny);
ww_i( litoset.Zm>1e-4 ) = litoset.wwInit;
XX_i = fft_exposure(ww_i); % Dosage
fx_i = 1./( 1 + exp(-litoset.a*(XX_i-litoset.tr)) ); % feature

% compute final quantities
XX_f = fft_exposure(result.wwTot); % Dosage
fx_f = 1./( 1 + exp(-litoset.a*(XX_f-litoset.tr)) ); % feature

% colormap
cmap = 'jet';

% plotting margins
xinsetl = 0.05;
xinsetr = 0.05;
yinsetb = 0.05;
yinsett = 0.05;
xgap = 0.000;
ygap = 0.10;
width = (1- xinsetl -xinsetr - 2*xgap)/3;
height = (1 - yinsetb -yinsett - 2*ygap)/2;%/3;
l1 = xinsetl;
l2 = l1+width+xgap;
l3 = l2+width+xgap;
% b3 = yinsetb;
% b2 = b3 + height+ygap;
b2 = ygap;
b1 = b2+ygap+height;

%% plot exposure/dosage/feature
FigHandle = figure;
set(FigHandle, 'Position', [100 100 1500 1000]);

% Plot Exposure
ax1=subplot('position',[l1,b1,width,height]);
imagesc(x,y,ww_i); % exposure plot
axis square
title('Initial exposure')
ylabel('$x$ [$\mu$m]','interpreter','latex','FontSize',14)
xlabel('$y$ [$\mu$m]','interpreter','latex','FontSize',14)
colormap(ax1,cmap);
caxis([0 max(ww_i(:))])
colorbar

% Plot Dosage
ax2=subplot('position',[l2,b1,width,height]);
imagesc(x,y,XX_i);
axis square
title('Initial dosage')
ylabel('$x$ [$\mu$m]','interpreter','latex','FontSize',14)
xlabel('$y$ [$\mu$m]','interpreter','latex','FontSize',14)
colormap(ax2,cmap);
colorbar

% Plot 3D Feature
ax3=subplot('position',[l3,b1,width,height]);
imagesc(x,y,fx_i)
axis square
colormap(ax3,cmap)
caxis([0 1])
colorbar
ylabel('$x$ [$\mu$m]','interpreter','latex','FontSize',14)
xlabel('$y$ [$\mu$m]','interpreter','latex','FontSize',14)
title('Initial feature')

% Final Exposure
ax4=subplot('position',[l1,b2,width,height]);
imagesc(x,y,result.wwTot); % exposure plot
axis square
title('Final exposure')
ylabel('$x$ [$\mu$m]','interpreter','latex','FontSize',14)
xlabel('$y$ [$\mu$m]','interpreter','latex','FontSize',14)
colormap(ax4,cmap);
caxis([0 max(result.wwTot(:))])
colorbar

% Final Dosage
ax5=subplot('position',[l2,b2,width,height]);
imagesc(x,y,XX_f);
axis square
title('Final dosage')
ylabel('$x$ [$\mu$m]','interpreter','latex','FontSize',14)
xlabel('$y$ [$\mu$m]','interpreter','latex','FontSize',14)
colormap(ax5,cmap);
caxis([0 max(XX_f(:))])
colorbar

% Final Feature
ax6=subplot('position',[l3,b2,width,height]);
imagesc(x,y,fx_f)
axis square
colormap(ax6,cmap)
caxis([0 1])
colorbar
ylabel('$x$ [$\mu$m]','interpreter','latex','FontSize',14)
xlabel('$y$ [$\mu$m]','interpreter','latex','FontSize',14)
title('Final feature')

drawnow

%% plot optimisation cost(s)
figure;
for qq = 1:litoset.Nsdx*litoset.Nsdy
    subplot( litoset.Nsdy, litoset.Nsdx, qq);
    try
        semilogy(result.info{qq}.err(:,1),'LineWidth',1.5)
        grid on
        axis square
        title 'Cost decrease'
        xlabel 'Iteration'
        ylabel('$f(\mathbf{w})$','interpreter','latex','FontSize',14)
    catch
        title 'No variables'
    end
end
end