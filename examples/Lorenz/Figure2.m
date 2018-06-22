%
% This file contains examples of how to use the functions mdDelay and mdFnn
% on example data from the Lorenz equations. This script was used to
% produce Figure 2 in the article.
%

%
% Set the font size for plots
%
fontSize = 18;

%
% Time delay and plot using all variables
%
figure()
data = load('lorenz_3d_timeseries.txt');
tau = mdDelay(data, 'maxLag', 25, 'plottype', 'all');
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
disp('xyz: tau = ' + string(tau))
print('Figure2a','-dpng')

%
% Estimate the embedding dimension
%
figure()
[fnnPercent, embeddingDimension] = mdFnn(data, round(tau));
print('Figure2b','-dpng')

%
% Time delay and plot only the x-variable
%
figure()
tau = mdDelay(data(:,1), 'maxLag', 25, 'plottype', 'all');
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
disp('x: tau = ' + string(tau))

%
% Time delay and plot only the y-variable
%
figure()
tau = mdDelay(data(:,2), 'maxLag', 25, 'plottype', 'all');
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
disp('y: tau = ' + string(tau))

%
% Time delay and plot only the z-variable
%
figure()
tau = mdDelay(data(:,3), 'maxLag', 25, 'plottype', 'all');
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
disp('z: tau = ' + string(tau))

%
% Time delay and plot x and y variables
%
figure()
tau = mdDelay(data(:,1:2), 'maxLag', 25, 'plottype', 'all');
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
disp('xy: tau = ' + string(tau))

%
% Time delay and plot x and z variables
%
figure()
tau = mdDelay(data(:,[1,3]), 'maxLag', 25, 'plottype', 'all');
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
disp('xz: tau = ' + string(tau))

%
% Time delay and plot y and z variables
%
figure()
tau = mdDelay(data(:,2:3), 'maxLag', 25, 'plottype', 'all');
set(gca,'FontSize',fontSize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontSize,'fontWeight','normal')
disp('yz: tau = ' + string(tau))