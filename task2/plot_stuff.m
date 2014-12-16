%% potential plot - v-cycle

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
    filename = ['grid' num2str(2^(i+2)) num2str(1) 'v.data'];

    data = dlmread(filename,'\t');
    data = data(:,1:end-1);
    data = data(:,fix(end/2)+1);

    xData = linspace(0,1,length(data));

    plot(xData,data,'Color',cc(i,:))
end

xlabel('x','FontSize',textStorlek)
ylabel('\Phi(x,1/2)','FontSize',textStorlek)

h = legend('Exact solution','81 V-cycle','161 V-cycle','321 V-cycle','641 V-cycle','1281 V-cycle')
set(h,'FontSize',legendStorlek);
hold off

saveas(gcf,'task2v.png','png')

axis([0.58 0.62 0.3 1.42])

saveas(gcf,'task2v_zoom.png','png')

%% potential plot - w-cycle

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
    filename = ['grid' num2str(2^(i+2)) num2str(1) 'w.data'];

    dataV = dlmread(filename,'\t');
    dataV = dataV(:,1:end-1);
    dataV = dataV(:,fix(end/2)+1);

    xData = linspace(0,1,length(dataV));

    plot(xData,dataV,'Color',cc(i,:))
end

xlabel('x','FontSize',textStorlek)
ylabel('\Phi(x,1/2)','FontSize',textStorlek)

h = legend('Exact solution','81 W-cycle','161 W-cycle','321 W-cycle','641 W-cycle','1281 W-cycle')
set(h,'FontSize',legendStorlek);
hold off

saveas(gcf,'task2w.png','png')

axis([0.58 0.62 0.3 1.42])

saveas(gcf,'task2w_zoom.png','png')


%% depth

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data = dlmread('log81w.data','\t',1,0);
savename = 'task2_w_depth.png';
data = data(:,1);
data = data( data ~= 0 );

semilogy(data,'x-');
set(gca, 'YTick', [0 11 21 41 81 161 321 641 1281]);

xlabel('Step nr','FontSize',textStorlek)
ylabel('Number of points','FontSize',textStorlek)

h = legend('Number of grid points');
set(h,'FontSize',legendStorlek);


saveas(gcf,savename,'png')

%% nIterations

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;
for i = 1:5
    filename = ['log' num2str(2^(i+2)) num2str(1) 'v.data'];
    
    data = dlmread(filename,'\t',1,0);
    data = data(:,2);
    data = data( data ~= 0 );
    yV(i) = sum(data);
    
    filename = ['log' num2str(2^(i+2)) num2str(1) 'w.data'];
    
    data = dlmread(filename,'\t',1,0);
    data = data(:,2);
    data = data( data ~= 0 );
    yW(i) = sum(data);   
    
    filename = ['../task3/log' num2str(2^(i+2)) num2str(1) '.data'];

    data = dlmread(filename,'\t',1,0);
    data = data(:,2);
    yFMG(i) = sum(data);
    
end

xData = [81 161 321 641 1281];
yBehavior = xData.^2;
yBehavior = yBehavior * yV(1)/yBehavior(1);
yBehavior = yBehavior ;


plot(xData, yV,'r')
hold on
plot(xData, yW,'b')
plot(xData, yFMG,'g')
plot(xData, yBehavior, 'k')
hold off

xlabel('Grid size','FontSize',textStorlek)
ylabel('Number of GS iterations','FontSize',textStorlek)

h = legend('V-cycle', 'W-cycle', 'Full multigrid', '(Grid size)^2');
set(h,'FontSize',legendStorlek);

saveas(gcf,'task2_its.png','png')
