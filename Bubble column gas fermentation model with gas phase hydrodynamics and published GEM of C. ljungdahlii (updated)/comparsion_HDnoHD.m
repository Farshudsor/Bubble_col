function [] = comparsion_HDnoHD( )

clear all
close all
clc

load 20190319_HD_p.mat;
biomass1 = yfin1;
cl1 = yfin2;
cg1 = yfin3_2;
hl1 = yfin4;
hg1 = yfin5_2;
c2l1 = yfin6_2;
c2g1 = yfin7_2;
acetate1 = yfin8;
ethanol1 = yfin9;
p1 = yfin10;
jg1 = yfin11;
db1 = yfin12;
ub1 = yfin13;
eg1 = yfin14;
klac1 = klac;
zt_1 = zt;
mu_1 = mu_spatial;
vc_1 = vc_spatial;
vc2_1 = vc2_spatial;
va_1 = va_spatial;
ve_1 = ve_spatial;
biomass1 = biomass1*Qmedia;
acetate1 = acetate1*Qmedia;
ethanol1 = ethanol1*Qmedia;
cg1 = cg1*28.01*3600/1000;
hg1 = hg1*2*3600/1000;

load 20190319_noHD_p.mat;
biomass2 = yfin1;
cl2 = yfin2;
cg2 = yfin3;
hl2 = yfin4;
hg2 = yfin5;
c2l2 = yfin6;
c2g2 = yfin7;
acetate2 = yfin8;
ethanol2 = yfin9;
mu_2 = mu_spatial;
vc_2 = vc_spatial;
vc2_2 = vc2_spatial;
va_2 = va_spatial;
ve_2 = ve_spatial;

zt_2 = zt_1;

ug = 81.7902;
for i = 1:N
    eg2(i) = 0.1246;
    jg2(i) = ug/3600;
    db2(i) = 1.2249;
    p2(i) = pL + zs*(N-i)*9.81*993*(1-eg2(i));
    klac2(i) = 6*eg2(i)/(1-eg2(i))/(db2(i)/1000)*1e-4*3600;
end

A = 3;

biomass2 = biomass2*Qmedia;
acetate2 = acetate2*Qmedia;
ethanol2 = ethanol2*Qmedia;
cg2 = cg2*ug*A*28.01/1000;
hg2 = hg2*ug*A*2/1000;
(biomass2(end) - biomass1(end))/biomass2(end)

figure(1);
subplot(2,3,1)
hold on
plot(zt_1,biomass1,'-b','linewidth',2)
plot(zt_2,biomass2,'-r','linewidth',2)
legend1 = legend('Recycle (Hydrodynamics)','Recycle (No Hydrodynamics)','location','northwest');
x1 = xlabel('Location [m]','fontsize',14);
y1 = ylabel('Biomass [kg/h]','fontsize',14);
get(x1)
get(y1)
set(legend1,'fontsize',12)
set(gca,'fontsize',16)
title('a')
xlim([0 L])
ylim([4 6])

subplot(2,3,2)
hold on
plot(zt_1,cl1,'-b','linewidth',2)
plot(zt_2,cl2,'-r','linewidth',2)
x2 = xlabel('Location [m]','fontsize',14);
y2 = ylabel('dissolved CO concentration [mmol/L]','fontsize',14);
get(x2)
get(y2)
set(gca,'fontsize',16)
title('b')
xlim([0 L])

subplot(2,3,3)
hold on
plot(zt_1,cg1,'-b','linewidth',2)
plot(zt_2,cg2,'-r','linewidth',2)
x3 = xlabel('Location [m]','fontsize',14);
y3 = ylabel('gas phase CO flowrate [kg/h]','fontsize',14);
get(x3)
get(y3)
set(legend1,'fontsize',12)
set(gca,'fontsize',16)
title('c')
xlim([0 L])

subplot(2,3,4)
hold on
plot(zt_1,hg1,'-b','linewidth',2)
plot(zt_2,hg2,'-r','linewidth',2)
x6 = xlabel('Location [m]','fontsize',14);
y6 = ylabel('H2 in gas [kg/h]','fontsize',14);
get(x6)
get(y6)
set(gca,'fontsize',16)
title('d')
xlim([0 L])
% ylim([0 1.6])

subplot(2,3,5)
hold on
plot(zt_1,acetate1,'-b','linewidth',2)
plot(zt_2,acetate2,'-r','linewidth',2)
x4 = xlabel('Location [m]','fontsize',14);
y4 = ylabel('Acetate [kg/h]','fontsize',14);
get(x4)
get(y4)
set(gca,'fontsize',16)
title('e')
xlim([0 L])
% ylim([10 14])

subplot(2,3,6)
hold on
plot(zt_1,ethanol1,'-b','linewidth',2)
plot(zt_2,ethanol2,'-r','linewidth',2)
x5 = xlabel('Location [m]','fontsize',14);
y5 = ylabel('Ethanol [kg/h]','fontsize',14);
get(x5)
get(y5)
set(gca,'fontsize',16)
xlim([0 L])
title('f')
% ylim([7 11])

figure(2)

subplot(2,3,1)
hold on
plot(zt_1,jg1*3600,'-b','linewidth',2)
plot(zt_2,jg2*3600,'-r','linewidth',2)
legend2 = legend('Recycle (Hydrodynamics)','Recycle (No Hydrodynamics)','location','northwest');
x7 = xlabel('Location [m]','fontsize',14);
y7 = ylabel('Superficial gas velocity [m/h]','fontsize',14);
get(x7)
get(y7)
set(legend2,'fontsize',12)
set(gca,'fontsize',16)
title('a')
xlim([0 L])

subplot(2,3,2)
hold on
plot(zt_1,db1,'-b','linewidth',2)
plot(zt_2,db2,'-r','linewidth',2)
x8 = xlabel('Location [m]','fontsize',14);
y8 = ylabel('Bubble diameter [mm]','fontsize',14);
get(x8)
get(y8)
set(gca,'fontsize',16)
title('b')
xlim([0 L])

subplot(2,3,3)
hold on
plot(zt_1,eg1,'-b','linewidth',2)
plot(zt_2,eg2,'-r','linewidth',2)
x9 = xlabel('Location [m]','fontsize',14);
y9 = ylabel('Gas holdup [-]','fontsize',14);
get(x9)
get(y9)
set(gca,'fontsize',16)
title('c')
xlim([0 L])
% ylim([0.18 0.26])

subplot(2,3,4)
hold on
plot(zt_1,klac1,'-b','linewidth',2)
plot(zt_2,klac2,'-r','linewidth',2)
x10 = xlabel('Location [m]','fontsize',14);
y10 = ylabel('Volumetric mass transfer coefficient [1/h]','fontsize',14);
get(x10)
get(y10)
set(gca,'fontsize',16)
title('d')
xlim([0 L])
% ylim([350 520])

subplot(2,3,5)
hold on
plot(zt_1,mu_1,'-b','linewidth',2)
plot(zt_2,mu_2,'-r','linewidth',2)
x11 = xlabel('Location [m]','fontsize',14);
y11 = ylabel('Growth rate [1/h]','fontsize',14);
get(x11)
get(y11)
set(gca,'fontsize',15)
title('e')
xlim([0 L])
ylim([0 0.12])

subplot(2,3,6)
hold on
plot(zt_1,vc_1,'-b','linewidth',2)
plot(zt_2,vc_2,'-r','linewidth',2)
x12 = xlabel('Location [m]','fontsize',14);
y12 = ylabel('CO uptake rate [mmol/gDW/h]','fontsize',14);
get(x12)
get(y12)
set(gca,'fontsize',16)
title('f')
xlim([0 L])

end

