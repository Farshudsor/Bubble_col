clear
results = load('results_sample3.mat');
%results = struct();
%results.results = res;

sx = size(results.results);
hrs = 401;
samp = 401;
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold on
for i = 1:sx(1)
    plot([1:hrs],Ce(i,:))
end
xlabel('Time [hr]')
ylabel('Ethanol [g/L]')
xlim([0 hrs])
legend('1','2','3','4')

figure(2)
hold on
for i = 1:sx(1)
    plot([0:hrs-1],Ca(i,:))
end
xlabel('Time [hr]')
ylabel('acetate [g/L]')
xlim([0 hrs])

figure(3)
hold on
for i = 1:sx(1)
    plot([0:hrs-1],Cx(i,:))
end
xlabel('Time [hr]')
ylabel('Biomass [g/L]')
xlim([0 hrs])

figure(4)
hold on
for i = 1:sx(1)
    plot([0:hrs-1],D(i,:))
end
xlabel('Time [hr]')
ylabel('Dilution ')
xlim([0 hrs])
ylim([.0099 .1001])

figure(5)
hold on
for i = 1:sx(1)
    plot([0:hrs-1],ug(i,:))
end
xlabel('Time [hr]')
ylabel('Gas velocity ')
xlim([0 hrs])
%ylim([.0099 .1001])


figure(6)
hold on
for i = 1:sx(1)
    plot([0:hrs-1],Ce(i,:)./Ca(i,:))
end
plot([0:hrs-1],ones(hrs,1)*slt)

xlabel('Time [hr]')
ylabel('Selectivity [g/L]')
xlim([0 hrs])
ylim([0 3.5])

figure(7)
subplot(3,3,1)
hold on 
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T11(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 1,1 - v_e ')
    xlim([0 hrs])

subplot(3,3,4)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T21(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 2,1 - v_a ')
    xlim([0 hrs])

subplot(3,3,7)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T31(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 3,1 - \mu ')
    xlim([0 hrs])

subplot(3,3,2)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T12(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 1,2 - v_e ')
    xlim([0 hrs])

subplot(3,3,5)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T22(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 2,2 - v_a ')
    xlim([0 hrs])

subplot(3,3,8)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T32(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 3,2 - \mu ')
    xlim([0 hrs])
    
    
 subplot(3,3,3)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T13(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 1,3 - v_e ')
    xlim([0 hrs])

subplot(3,3,6)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T23(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 2,3 - v_a ')
    xlim([0 hrs])

subplot(3,3,9)
hold on
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T33(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 3,3 - \mu ')
    xlim([0 hrs])


%%

score_1 = score(Ce(1,:),Ca(1,:))
score_2 = score(Ce(2,:),Ca(2,:))
score_3 = score(Ce(3,:),Ca(3,:))
score_4 = score(Ce(4,:),Ca(4,:))
score_5 = score(Ce(5,:),Ca(5,:))
score_6 = score(Ce(6,:),Ca(6,:))









