%% potential plot

clear all
clc
clf

cc = [ 0 1 1 ; 1 0 1 ; 0 1 0 ; 0 0 1 ; 1 0 0];
textStorlek = 14;
legendStorlek = 11;

exactData = dlmread('../phi_exact.data','\t');
xExactData = linspace(0,1,length(exactData));
plot(xExactData,exactData,'k')

hold on

for i = 1:5
    filename = ['grid' num2str(2^(i+2)) num2str(1) '.data'];

    data = dlmread(filename,'\t');
    data = data(:,1:end-1);
    data = data(:,fix(end/2)+1);

    xData = linspace(0,1,length(data));

    plot(xData,data,'Color',cc(i,:))

end

xlabel('x','FontSize',textStorlek)
ylabel('\Phi(x,1/2)','FontSize',textStorlek)

h = legend('Exact solution','81','161','321','641','1281')
set(h,'FontSize',legendStorlek);
hold off

saveas(gcf,'task3.png','png')

axis([0.58 0.62 0.3 1.42])

saveas(gcf,'task3_zoom.png','png')


%% depth

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data = dlmread('log161.data','\t');

semilogy(data,'x-')
set(gca, 'YTick', [0 11 21 41 81 161 321 641 1281]);

xlabel('Iterations','FontSize',textStorlek)
ylabel('Number of points','FontSize',textStorlek)

h = legend('Number of grid points');
set(h,'FontSize',legendStorlek);

saveas(gcf,'task3_depth.png','png')

%% times - not used right now

times = [0.05,0.45,3.01,17.1];
times2 = [0.05,0.44,2.99,16.9];

avgTimes = times + times2;
avgTimes = avgTimes/2;

plot(avgTimes,'r')
hold on
plot(0.05*[1 2 4 8].^2.8)
hold off

