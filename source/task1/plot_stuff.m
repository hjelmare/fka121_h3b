%% potential plot
clear all
clc
clf
textStorlek = 14;
legendStorlek = 11;

data = dlmread('grid.data','\t');
exactData = dlmread('../phi_exact.data', '\t');

data = data(:,1:end-1);

middleIndex = (length(data)-1)/2 + 1;

xData = linspace(0,1,length(data));
xExact = linspace(0,1, length(exactData));

surf(xData,xData,data)
xlabel('Y', 'FontSize', textStorlek);
ylabel('X', 'FontSize', textStorlek);
zlabel('\Phi(x,y)', 'FontSize', textStorlek);
shading flat
%view(45,0)
figure
hold on
plot(xData, data(:,middleIndex))
plot(xExact, exactData, 'r');
xlabel('X', 'FontSize', textStorlek);
ylabel('\Phi(x,1/2)', 'FontSize', textStorlek);
text=legend('Simulated result', 'Exact result');
set(text, 'FontSize', legendStorlek);

saveas(gcf,'task1.png','png')

axis([0.58 0.62 0.3 1.42])

saveas(gcf,'task1_zoom.png','png')

