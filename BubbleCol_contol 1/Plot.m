clear
res = load('results.mat');
results = struct();
results.results = res.results(3);

sx = size(results.results);
hrs = 501;
samp = 501;
Tsamp = 1;
slt = .84;

Ce = zeros(sx(2), hrs);
Ca = zeros(sx(2), hrs);
Cz = zeros(sx(2), hrs);

D = zeros(sx(2), samp);

T11 = zeros(sx(2), samp);
T12 = zeros(sx(2), samp);
T21 = zeros(sx(2), samp);
T22 = zeros(sx(2), samp);
T31 = zeros(sx(2), samp);
T32 = zeros(sx(2), samp);

for i =1:sx
    Ce(i,:) = results.results(i).Ce;
    Ca(i,:) = results.results(i).Ca;
    Cx(i,:) = results.results(i).Cx;

    D(i,:) = results.results(i).D;

    T11(i,:) = results.results(i).T11;
    T21(i,:) = results.results(i).T21;
    T31(i,:) = results.results(i).T31;
    T12(i,:) = results.results(i).T12;
    T22(i,:) = results.results(i).T22;
    T32(i,:) = results.results(i).T32;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold on
for i = 1:sx(1)
    plot([0:hrs-1],Ce(i,:))
end
xlabel('Time [hr]')
ylabel('Ethanol [g/L]')
xlim([0 hrs])

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
    plot([0:hrs-1],Ce(i,:)./Ca(i,:))
end
plot([0:hrs-1],ones(hrs,1)*slt)

xlabel('Time [hr]')
ylabel('Selectivity [g/L]')
xlim([0 hrs])
ylim([0 3.5])

figure(6)
hold on 
subplot(3,2,1)
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T11(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 1,1 - v_e ')
    xlim([0 hrs])

subplot(3,2,3)
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T21(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 2,1 - v_a ')
    xlim([0 hrs])

subplot(3,2,5)
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T31(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 3,1 - \mu ')
    xlim([0 hrs])

subplot(3,2,2)
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T12(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 1,2 - v_e ')
    xlim([0 hrs])

subplot(3,2,4)
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T22(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 2,2 - v_a ')
    xlim([0 hrs])

subplot(3,2,6)
    for i = 1:sx(1)
        plot([1:Tsamp:hrs],T32(i,:))
    end
    xlabel('Time [hr]')
    ylabel('Theta 3,2 - \mu ')
    xlim([0 hrs])














