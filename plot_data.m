%% Marco Iglesias, Universtity of Nottingham, 2022
function plot_data(data_noisy, Truth)
%%plot the noisy (displacements) data

l=length(data_noisy)/4;
data_res=reshape(data_noisy,l,4);
fig=figure('Position', [10 10 700 600]);
ha = tight_subplot(2,2,[.05 .1],[.05 .1],[.1 .1]);
M=length(Truth.Grid.X_meas); 
axes(ha(1))
imagesc(Truth.Grid.X_meas(:), Truth.Grid.Y_meas(:),  reshape(data_res(:,1)/1e-6,M,M));colormap jet;colorbar;shading flat
xlabel('$$X$$ [m]','interpreter','latex','fontsize',15);
ylabel('$$Y$$ [m]','interpreter','latex','fontsize',15);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',14)
c=colorbar;
c.Label.String = '$${\bf d_{re}^{x}}[\mu m]$$';
c.Label.Interpreter = 'latex';
c.FontSize=20;
c.Label.Rotation=0;
c.Label.Position=[2.25, 2.5,0];
axes(ha(2))
imagesc(Truth.Grid.X_meas(:), Truth.Grid.Y_meas(:), reshape(data_res(:,2)/1e-6,M,M));colormap jet;colorbar;shading flat
axis square
xlabel('$$X$$ [m]','interpreter','latex','fontsize',15);
ylabel('$$Y$$ [m]','interpreter','latex','fontsize',15);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',14)
c=colorbar;
c.Label.String = '$${\bf d_{re}^{y}}[\mu m]$$';
c.Label.Interpreter = 'latex';
c.FontSize=20;
c.Label.Rotation=0;
c.Label.Position=[2.0,1.225,0];
axes(ha(3))
imagesc(Truth.Grid.X_meas(:), Truth.Grid.Y_meas(:),  reshape(data_res(:,3)/1e-6,M,M));colormap jet;colorbar;shading flat
axis square
xlabel('$$X$$ [m]','interpreter','latex','fontsize',15);
ylabel('$$Y$$ [m]','interpreter','latex','fontsize',15);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',14)
c=colorbar;
c.Label.String = '$${\bf d_{im}^{x}}[\mu m]$$';
c.Label.Interpreter = 'latex';
c.FontSize=20;
c.Label.Rotation=0;
c.Label.Position=[2.300   1.0         0];
axes(ha(4))
imagesc(Truth.Grid.X_meas(:), Truth.Grid.Y_meas(:),  reshape(data_res(:,4)/1e-6,M,M));colormap jet;colorbar;shading flat
xlabel('$$X$$ [m]','interpreter','latex','fontsize',15);
ylabel('$$Y$$ [m]','interpreter','latex','fontsize',15);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',14)
axis square
c=colorbar;
c.Label.String = '$${\bf d_{im}^{y}}[\mu m]$$';
c.Label.Interpreter = 'latex';
c.FontSize=18;
c.Label.Rotation=0;
c.Label.Position=[2.0000 ,  0.3   ,      0];