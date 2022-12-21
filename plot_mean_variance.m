%% Marco Iglesias, Universtity of Nottingham, 2022
% plots the truth and the prior/posterior mean/variance for the
% storage/loss modulus, as well as the prior/posterior probability for the
% cyst

clear 
close all
load 'Results_first'
load Truth
Grid=Model.Grid;
Storage_truth = reshape(Truth.Storage,Truth.Grid.Nx,Truth.Grid.Ny);
Loss_truth = reshape(Truth.Loss,Truth.Grid.Nx,Truth.Grid.Ny);



min_p_k=min(Storage_truth(:)/1e3);
max_p_k=max(Storage_truth(:)/1e3);
Lx=Model.Lx;
Ly=Model.Ly;
X=linspace(0,Lx,Model.Grid.Nx);
Y=linspace(0,Ly,Model.Grid.Nx);

figure('Position', [10 10 1250 750]);
ha = tight_subplot(2,3,[.02 .05],[.06 .06],[.06 .06]);
axes(ha(1))
imagesc(X,Y,Storage_truth/1e3);shading flat; colormap jet
clim([min_p_k,max_p_k]);
axis square
c=colorbar;
c.FontSize=15;
c.Label.String = '[kPa]';
c.Label.Rotation=0;
c.Label.Position=[0.6 ,  4.6   ,      0];
set(gca,'xtick',[])
set(gca,'ytick',[])
title('$${\bf True ~storage}$$','interpreter','latex','fontsize',25);


Storage_mean=reshape(meanSto,Model.Grid.Nx,Grid.Ny);

axes(ha(2))
imagesc(X,Y,Storage_mean/1e3);shading flat; colormap jet
clim([min_p_k,max_p_k]);
axis square
set(gca,'xtick',[])
set(gca,'ytick',[])
title('$${\bf Prior~mean~storage}$$','interpreter','latex','fontsize',25);
c=colorbar;
c.FontSize=15;
c.Label.String = '[kPa]';
c.Label.Rotation=0;
c.Label.Position=[0.6 ,  4.6   ,      0];



min_p_loss=min(Loss_truth(:)/1e3);
max_p_loss=1.002*max(Loss_truth(:)/1e3);
axes(ha(4))
imagesc(X,Y,Loss_truth/1e3);shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])
title('$${\bf True~ loss}$$','interpreter','latex','fontsize',25);
clim([min_p_loss,max_p_loss]);
c=colorbar;
c.FontSize=15;
c.Label.String = '[kPa]';
c.Label.Rotation=0;
c.Label.Position=[0.6 ,  3.6   ,      0];
axis square



Loss_mean=reshape(meanLoss,Model.Grid.Nx,Grid.Ny);
axes(ha(5))
imagesc(X,Y,Loss_mean/1e3);shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])
title('$${\bf Prior~mean~loss}$$','interpreter','latex','fontsize',25);
clim([min_p_loss,max_p_loss]);
c=colorbar;
c.FontSize=15;
c.Label.String = '[kPa]';
c.Label.Rotation=0;
c.Label.Position=[0.6 ,  3.6   ,      0];
axis square


load Results_final.mat
Storage_mean=reshape(meanSto,Model.Grid.Nx,Grid.Ny);
axes(ha(3))
imagesc(X,Y,Storage_mean/1e3);shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])
title('$${\bf Posterior~mean~storage}$$','interpreter','latex','fontsize',25);
c=colorbar;
c.FontSize=15;
c.Label.String = '[kPa]';
c.Label.Rotation=0;
c.Label.Position=[1.5 ,  4.6   ,      0];
clim([min_p_k,max_p_k]);
axis square


Loss_mean=reshape(meanLoss,Model.Grid.Nx,Grid.Ny);
axes(ha(6))
imagesc(X,Y,Loss_mean/1e3);shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])

title('$${\bf Posterior~mean~loss}$$','interpreter','latex','fontsize',25);
clim([min_p_loss,max_p_loss]);
c=colorbar;
c.FontSize=15;
c.Label.String = '[kPa]';
c.Label.Rotation=0;
c.Label.Position=[0.6 ,  3.6   ,      0];
axis square
drawnow;

pause(0.01);




figure('Position', [10 10 1250 750]);
ha = tight_subplot(2,3,[.02 .05],[.06 .06],[.06 .06]);

axes(ha(4))
imagesc(X,Y,reshape(log(VarSto/1e6),Model.Grid.Nx,Model.Grid.Ny));shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])
title('$${\bf (log)~Posterior~storage~variance}$$','interpreter','latex','fontsize',19);
clim([-9.5,0.0]);
axis square
c=colorbar;
c.FontSize=15;
c.Label.String = '[log(kPa^2)]';
c.Label.Rotation=0;
c.Label.Position=[2.4,1.0,0];

axes(ha(5))
imagesc(X,Y,reshape(log(VarLoss/1e6),Model.Grid.Nx,Model.Grid.Ny));shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])

title('$${\bf (log)~Posterior~loss~variance}$$','interpreter','latex','fontsize',19);
clim([-9.0,-1.2]);
c=colorbar;
c.FontSize=15;
c.Label.String = '[log(kPa^2)]';
c.Label.Rotation=0;
c.Label.Position=[2.4375,-0.5,0];
axis square
drawnow;



load Results_first.mat
axes(ha(1))
imagesc(X,Y,reshape(log(VarSto/1e6),Model.Grid.Nx,Model.Grid.Ny));shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])

title('$${\bf (log)~Prior~storage~variance}$$','interpreter','latex','fontsize',19);
clim([-9.5,0.0]);
c=colorbar;
c.FontSize=15;
c.Label.String = '[log(kPa^2)]';
c.Label.Rotation=0;
c.Label.Position=[2.4,1.0,0];
axis square



axes(ha(2))
imagesc(X,Y,reshape(log(VarLoss/1e6),Model.Grid.Nx,Model.Grid.Ny));shading flat; colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])
title('$${\bf (log)~Prior~loss~variance}$$','interpreter','latex','fontsize',19);
clim([-9.0,-1.2]);
c=colorbar;
c.FontSize=15;
c.Label.String = '[log(kPa^2)]';
c.Label.Rotation=0;
c.Label.Position=[2.4375,-0.5,0];
axis square
drawnow;

axes(ha(3))
mask=Model.mask;
U_temp=(( mask==1)|(mask==2)).*reshape(prob,Model.Grid.Nx,Model.Grid.Ny);
U_phys=10.*( mask==0)+U_temp.*( mask~=0);
U_phys(U_phys==10)=NaN;
imagescwithnan(X,Y,U_phys,cool,[1 1 1]);
title('$${\bf Prior~cyst~probability}$$','interpreter','latex','fontsize',19);
set(gca,'xtick',[])
set(gca,'ytick',[])
clim([0,1]);
c=colorbar;
c.FontSize=15;
axis square;



load Results_final.mat
axes(ha(6))
mask=Model.mask;
U_temp=(( mask==1)|(mask==2)).*reshape(prob,Model.Grid.Nx,Model.Grid.Ny);
U_phys=10.*( mask==0)+U_temp.*( mask~=0);
U_phys(U_phys==10)=NaN;
imagescwithnan(X,Y,U_phys,cool,[1 1 1]);
title('$${\bf Posterior~cyst~probability}$$','interpreter','latex','fontsize',19);
axis square;
set(gca,'xtick',[])
set(gca,'ytick',[])
c=colorbar;
c.FontSize=15;
colormap(ha(1),jet);
colormap(ha(2),jet);
colormap(ha(4),jet);
colormap(ha(5),jet);




