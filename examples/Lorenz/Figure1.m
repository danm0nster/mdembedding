fontSize = 18;
close all

data = load('lorenz_3d_timeseries.txt');

figure1a = figure
set(gcf,'color','white')
axes1a = axes('Parent',figure1a);
hold(axes1a,'on');
plot3(data(:,1),data(:,2),data(:,3),'k')
xlabel('x')
ylabel('y')
zlabel('z')
view(axes1a,[67 25]);
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
print('Figure1a','-dpng')

figure
set(gcf,'color','white')
subplot(3,1,1), plot(data(:,1),'k')
ylabel('x')
axis([0 2000 min(data(:,1)) max(data(:,1))])
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
subplot(3,1,2), plot(data(41:end,1),'k')
ylabel('x')
axis([0 2000 min(data(:,1)) max(data(:,1))])
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
subplot(3,1,3), plot(data(81:end,1),'k')
ylabel('x')
xlabel('time')
axis([0 2000 min(data(:,1)) max(data(:,1))])
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')

figure1e = figure
set(gcf,'color','white')
axes1e = axes('Parent',figure1e);
hold(axes1e,'on');
plot3(data(1:end-20,1),data(11:end-10,1),data(21:end,1),'k')
xlabel('x')
ylabel('y')
zlabel('z')
view(axes1e,[67 25]);
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
