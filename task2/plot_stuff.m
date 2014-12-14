%% potential plot

clear all
clc
clf

data = dlmread('grid.data','\t');
data = data(:,1:end-1);
data = data(:,fix(end/2)+1);

xData = linspace(0,1,length(data));

exactData = dlmread('../phi_exact.data','\t');
xExactData = linspace(0,1,length(exactData));

plot(xExactData,exactData,'r')
hold on
plot(xData,data,'b')
hold off


%%

nPoints = length(data);

surf([0:nPoints-1],[0:nPoints-1],data)
shading flat
view(45,0)

%% depth

clear all
clc
clf

data = dlmread('log.data','\t');

semilogy(data)


%% times

times = [0.05,0.45,3.01,17.1];
times2 = [0.05,0.44,2.99,16.9];

avgTimes = times + times2;
avgTimes = avgTimes/2;

plot(avgTimes,'r')
hold on
plot(0.05*[1 2 4 8].^2.8)
hold off

%% Felsökning
clear all
clc
clf
close all

data = dlmread('grid.data','\t');
dataStart = dlmread('start.data','\t');

data = data(:,1:end-1);
dataStart = dataStart(:,1:end-1);

nPoints = length(data);
nPointsStart = length(dataStart);

surf([0:nPoints-1],[0:nPoints-1],data)
view(45,0)

figure
surf([0:nPointsStart-1],[0:nPointsStart-1],dataStart)
view(45,0)
