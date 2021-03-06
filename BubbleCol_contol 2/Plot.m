clear
results = load('results_sample6.mat');
%results = struct();
%results.results = res;

sx = size(results.results);

%for sample5
hrs = 960;
samp = 961;
plottime = [0:.25:samp/4-.25];
plottime2 = [0:.25:hrs/4-.25];


%for sample4
% hrs = 401;
% samp = 401;
% plottime = [1:1:hrs];

Tsamp = 1;
slt = .84;

Ce = zeros(sx(1), hrs);
Ca = zeros(sx(1), hrs);
Cz = zeros(sx(1), hrs);

D = zeros(sx(1), samp);

T11 = zeros(sx(1), samp);
T12 = zeros(sx(1), samp);
T13 = zeros(sx(1), samp);
T21 = zeros(sx(1), samp);
T22 = zeros(sx(1), samp);
T23 = zeros(sx(1), samp);
T31 = zeros(sx(1), samp);
T32 = zeros(sx(1), samp);
T33 = zeros(sx(1), samp);

st = zeros(sx(1),1);

for i =1:sx(1)
    Ce(i,:) = results.results(i).Ce;
    Ca(i,:) = results.results(i).Ca;
    Cx(i,:) = results.results(i).Cx;

    D(i,:) = results.results(i).D;
    ug(i,:) = results.results(i).u_g;

    T11(i,:) = results.results(i).T11;
    T21(i,:) = results.results(i).T21;
    T31(i,:) = results.results(i).T31;
    T12(i,:) = results.results(i).T12;
    T22(i,:) = results.results(i).T22;
    T32(i,:) = results.results(i).T32;
    T13(i,:) = results.results(i).T13;
    T23(i,:) = results.results(i).T23;
    T33(i,:) = results.results(i).T33;
    
    st(i) = results.results(i).solve_time;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

leg = {'line 1'};
for ii = 2:sx(1)
   leg = [leg, "line "+num2str(ii)] ;
end


figure(1)
hold on
for i = 1:sx(1)
    plot(plottime2,Ce(i,:))
end
xlabel('Time [hr]')
ylabel('Ethanol [g/L]')
%xlim([0 hrs/4])
clickableLegend(leg)

figure(2)
hold on
for i = 1:sx(1)
    plot(plottime2,Ca(i,:))
end
xlabel('Time [hr]')
ylabel('acetate [g/L]')
%xlim([0 hrs/4])
clickableLegend(leg)


figure(3)
hold on
for i = 1:sx(1)
    plot(plottime2,Cx(i,:))
end
xlabel('Time [hr]')
ylabel('Biomass [g/L]')
%xlim([0 hrs/4])
clickableLegend(leg)


figure(4)
hold on
for i = 1:sx(1)
    plot(plottime,D(i,:))
end
xlabel('Time [hr]')
ylabel('Dilution ')
%xlim([0 hrs/4])
ylim([.0099 .1001])
clickableLegend(leg)


figure(5)
hold on
for i = 1:sx(1)
    plot(plottime,ug(i,:))
end
xlabel('Time [hr]')
ylabel('Gas velocity ')
xlim([0 hrs/4])
%ylim([.0099 .1001])
clickableLegend(leg)


figure(6)
hold on
for i = 1:sx(1)
    plot(plottime2,Ce(i,:)./Ca(i,:))
end
plot([0:hrs-1],ones(hrs,1)*slt)
xlabel('Time [hr]')
ylabel('Selectivity [g/L]')
xlim([0 hrs/4])
ylim([0 3.5])
clickableLegend(leg)


figure(7)
sgtitle('rate term, i = \theta_{i,1}D + \theta_{i2}u_g + \theta_{i3}')
subplot(3,3,1)
hold on 
    for i = 1:sx(1)
        plot(plottime,T11(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{1,1} - v_e ')
%    xlim([0 hrs/4])

subplot(3,3,4)
hold on
    for i = 1:sx(1)
        plot(plottime,T21(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{2,1} - v_a ')
%    xlim([0 hrs/4])

subplot(3,3,7)
hold on
    for i = 1:sx(1)
        plot(plottime,T31(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{3,1} - \mu ')
    xlim([0 hrs/4])

subplot(3,3,2)
hold on
    for i = 1:sx(1)
        plot(plottime,T12(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{1,2} - v_e ')
%    xlim([0 hrs/4])

subplot(3,3,5)
hold on
    for i = 1:sx(1)
        plot(plottime,T22(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{2,2} - v_a ')
%    xlim([0 hrs/4])

subplot(3,3,8)
hold on
    for i = 1:sx(1)
        plot(plottime,T32(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{3,2} - \mu ')
%    xlim([0 hrs/4])
    
    
 subplot(3,3,3)
hold on
    for i = 1:sx(1)
        plot(plottime,T13(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta{1,3} - v_e ')
%    xlim([0 hrs/4])

subplot(3,3,6)
hold on
    for i = 1:sx(1)
        plot(plottime,T23(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{2,3} - v_a ')
%    xlim([0 hrs/4])

subplot(3,3,9)
hold on
    for i = 1:sx(1)
        plot(plottime,T33(i,:))
    end
    xlabel('Time [hr]')
    ylabel('\theta_{3,3} - \mu ')
%    xlim([0 hrs/4])

    


%%

for i = 1:sx(1)
    [tot,viol] =  score(Ce(i,:),Ca(i,:));
    "run " + num2str(i) + "= score " + num2str(tot) + "; violations " + num2str(viol) 
end

avg_solve_time =mean(st)









