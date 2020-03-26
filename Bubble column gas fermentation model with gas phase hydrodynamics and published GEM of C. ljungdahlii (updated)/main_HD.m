%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Jose A. Gomez and Kai H�ffner                                %
%                                                                         % 
% This code can only be used for academic purposes. When using this code  %
% please cite:                                                            %
%                                                                         %
% Gomez, J.A., H�ffner, K. and Barton, P. I.                              %
% DFBAlab: A fast and reliable MATLAB code for Dynamic Flux Balance       %
% Analysis. Submitted.                                                    %
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This example is based on the High-rate algal pond (HRAP) proposed in
% Buhr, H.O., Miller, S.B. A dynamic model of the high-rate algal-bacterial
% wastewater treatment pond. Water Res. 17, 29-37 (1983), and in Yang, A.
% Modeling and evaluation of CO2 supply and utilization in algal ponds.
% Industrial and Engineering Chemistry Research 50, 11181-11192 (2011). In
% this case, a monoculture was simulated (instead of a coculture). 
% The algae model used is iRC1080 from Chang, R.L. et al. Metabolic network
% reconstruction of Chlamydnomnas offers insight into light-driven algal
% metabolism. Molecular Systems Biology 7(518) (2011).


% 20150119 use modified solveModel.m

clear all
close all
clc
path(pathdef)

% Add paths to LP solver and DFBAlab

% addpath(genpath('C:\Users\lixia_000\Documents\MATLAB\DFBAlab'))
% addpath(genpath('C:\LXA\Gurobi\win64'))

% Set number of models and simulation time
tic    
nmodel = 51;
tspan = [0,500];
N = nmodel;                     % number of reactor discretization points
ns = 11;                        % number of state variables

% Load models. These should be .mat files generated by the COBRA toolbox. 
% When generating these files using the COBRA toolbox, a big number is used
% as infinity. This number should be fed to the DB vector (Default bound).
% INFO.DB = DB;

for i = 1:nmodel
    load model.mat 
    model{i} = CL_Auto; 
    DB(i) = 1000;
end

%% exID array
% You can either search the reaction names by name or provide them directly
% in the exID array.
% RxnNames = {'EX_glc(e)', 'EX_ac(e)', 'biomass'};
% for i = 1:length(RxnNames)
%    [a,exID(i)] = ismember(RxnNames(i),model.rxns);
% end
% INFO.exID = exID;
% Lower bounds and upper bounds for these reactions should be provided in
% the RHS code. 

for i = 1:nmodel
    exID{i}=[20, 21, 43, 1, 26, 4]; % This order determines the order in RHS
end

% This codes solves the LPs in standard form. Bounds on exchange fluxes in 
% the exID array can be modified directly on the first 2*n rows where n is 
% the number of exchange fluxes. Order will be lower bound, upper bound, 
% lower bound, upper bound in the same order as exID. 
%
% NOTE: All bounds on fluxes in the exID arrays are relaxed to -Inf and 
% + Inf. These bounds need to be updated if needed in the RHS file.

%% Cost vectors
% Usually the first cost vector will be biomass maximization, but it can
% be any other objective. The CPLEX objects will minimize by default. 
% Report only nonzero elements. 
% The structure should be:
% C{model} = struct
% Element C{k}(i) of C, is the cost structure i for model k. 
% C{k}(i).sense = +1 for minimize, or -1 for maximize.
% C{k}(i).rxns = array containing the reactions in this objective. 
% C{k}(i).wts = array containing coefficients for reactions reported in 
% rxns. Both arrays should have the same length. 
% Example, if:
% C{k}(i).rxns = [144, 832, 931];
% C{k}(i).wts = [3, 1, -1];
% Then the cost vector for this LP will be:
% Cost{k}(i) = 3*v_144 + v_832 - v_931 (fluxes for model k). 
% This cost vector will be either maximized or minimized depending on the
% value of C{k}(i).sense.

% In SBML files, usually production fluxes are positive and uptake fluxes
% are negative. Keep in mind that maximizing a negative flux implies 
% minimizing its absolute value.
% Different models can have different number of objectives. 
% INFO.C = C; 

minim = 1;
maxim = -1;

for i=1:nmodel
    
    % Maximize growth
    C{i}(1).sense = maxim;
    C{i}(1).rxns = [1];
    C{i}(1).wts = [1];
    % Maximize CO uptake
    C{i}(2).sense = minim;
    C{i}(2).rxns = [20];
    C{i}(2).wts = [1];
    % Maximize H2 uptake
    C{i}(3).sense = minim;
    C{i}(3).rxns = [43];
    C{i}(3).wts = [1];
    % Minimize CO2 production
    C{i}(4).sense = minim;
    C{i}(4).rxns = [21];
    C{i}(4).wts = [1];
    % Maximize acetate uptake
    C{i}(5).sense = minim;
    C{i}(5).rxns = [4];
    C{i}(5).wts = [1];
    % Maximize ethanol uptake
    C{i}(6).sense = minim;
    C{i}(6).rxns = [26];
    C{i}(6).wts = [1];

end

% import steady-state solutions as initial conditions
load Initial_Conditions.mat
biomass_initial = yfin1i;
clsi = yfin2i;
cgi = yfin3i;
hlsi = yfin4i;
hgi = yfin5i;
c2lsi = yfin6i;
c2gi = yfin7i;
acetate_initial = yfin8i;
ethanol_initial = yfin9i;
p = yfin10i;
jgi = yfin11i;
dbi = yfin12i;
ubi = yfin13i;
egi = yfin14i;

%%
% Pass parameters
INFO.nmodel = nmodel;
INFO.DB = DB;
INFO.exID = exID;
INFO.C = C; 

% Set operating conditions
L = 10;                         % reactor length, m 
Area = 3;                       % reactor cross-sectional area, m^2
dR = sqrt(4*Area/pi);           % reactor diameter, m  (0.6308 m)
zs = L/(N-1);                   % spatial step size, m
dl = 0.25;                      % liquid dispersion coefficient, m^2/h
tr = 310.15;                    % temperature, K
pL = 1.013e5;                   % pressure at top of column, Pa
pc = 0.5;                       % partial pressure of CO %(60+2*(mm-1))/100;
Hc = 8.0e-4;                    % Henry's constant for CO in water, mol/(L*atm)
ph = 0.2;                       % partial pressure of H2
Hh = 6.6e-4;                    % Henry's constant for H2 in water, mol/(L*atm) (NIST)
pn = 0.3;                       % partial pressure of N2
pc2 = 0;                        % partial pressure of CO2
Hc2 = 2.5e-2;                   % Henry's constant for CO2 in water, mol/(L*atm)
g = 9.81;                       % gravitation, m/s^2
densityL = 993.34;              % density of liquid phase (roughly water at 310K), kg/m^3
viscosityL = 0.0009242;         % viscosity of liquid phase, Pa*s
D = 0.06;                       % dilution rate
Qmedia = D*L*Area;              % volumetric flow rate
p_atm = 1;                      % 1 atmospheric pressure

% Feed gas composition (inlet condition)
eg0 = 0.1;                      % gas holdup
ul = -50;                       % superficial liquid velocity, m/h 
ug0 = 150/3600;                 % superficial gas velocity, m/s
db0 = 1.5;                      % bubble diameter, mm
el0 = 1-eg0;                    % liquid phase holdup

% Set mass transfer coefficients
klc = (1e-4)*3600;              % CO gas-liquid mass transfer coefficient, m/h  
klh = klc;                      % H2 gas-liquid mass transfer coefficient, m/h
klc2 = klc;                     % CO2 gas-liquid mass transfer coefficient, m/h
klo= 0.54;                      % O2 gas-liquid mass transfer coeff m/h

% Form condition vector
condit = [klc,Hc,klh,Hh,klc2,Hc2,tr,zs,N,ns,ul,0,dl,pc,ph,pc2,...
    dR,pL,eg0,el0,Qmedia,D,Area,g,db0,ug0,densityL,viscosityL,pn,p_atm];

% Set uptake parameters
vcm = 35;
vhm = 70;
Kmc = 0.02;
Kmh = 0.02;
Kic = 0.601;
Amax = 20;
Emax = 60;
param = [
     vcm                 % maximum CO uptake rate, mmol/g/h; 34.36 in Mohammodi 2014
     Kmc                 % CO uptake saturation constant, mmol/L; 0.021atm in Mohamodi 2014 
     vhm                 % maximum H2 uptake rate, mmol/g/h 
     Kmh                 % H2 uptake saturation constant, mmol/L;
     vcm                 % maximum CO2 uptake rate, mmol/g/h 
     Kmc                 % CO2 uptake saturation constant, mmol/L;
     Kic                 % CO inhibition constant; atm
     Amax                % maximum acetate concentration. g/L;
     Emax                % maximum ethanol concentration, g/L;
    ];

% Pass parameters
INFO.param = param;
INFO.ns = ns;
INFO.N = N;
INFO.condit = condit;

yo = [];
for i=1:N
%     p(i) = pL+1000*9.81*zs*(N-i)*el0;
%     cgi(i) = pc*p(i)/8.314/tr;
%     hgi(i) = ph*p(i)/8.314/tr;
%     c2gi(i) = pc2*p(i)/8.314/tr;
%     clsi(i) = cgi(i)*8.314*tr*Hc*1000/1.013e5;
%     hlsi(i) = hgi(i)*8.314*tr*Hh*1000/1.013e5;
%     c2lsi(i) = c2gi(i)*8.314*tr*Hc2*1000/1.013e5;
%     jgi(i) = ug0*(pL/p(i));
%     r(i) = (cgi(i)+hgi(i)).*jgi(i)/(cgi(1)+hgi(1))/jgi(1);
%     dbi(i) = nthroot(pL/p(i)*r(i),3)*db0;
%     ubi(i) = 1/18*dbi(i)*dbi(i)*1e-6*densityL*g/viscosityL;
%     egi(i) = 0.3;
    yz = [clsi(i) cgi(i) hlsi(i) hgi(i) c2lsi(i) c2gi(i) p(i) jgi(i) dbi(i) ubi(i) egi(i)]; 
    yo = [yo yz];
end

yo = [yo, 0.1, 0, 0, 0];
Y0 = yo;

% Mass Matrix
M=[];
for i=1:N
    Mi = zeros(ns,ns); %11x11
    Mi(1,1) = 1;   % CO liquid concentration
    Mi(2,2) = 1;   % CO gas concentration
    Mi(3,3) = 1;   % H2 liquid concentration
    Mi(4,4) = 1;   % H2 gas concentration
    Mi(5,5) = 1;   % CO2 liquid concentration
    Mi(6,6) = 1;   % CO2 gas concentration
    M = blkdiag(M,Mi); % adding all the node points
end

M = blkdiag(M,1);
M = blkdiag(M,1);
M = blkdiag(M,1);
K = eye(1,1); % for the penalty term
M = blkdiag(M,K);

% CPLEX Objects construction parameters

INFO.LPsolver = 1; % CPLEX = 0, Gurobi = 1.
                   % CPLEX works equally fine with both methods.
                   % Gurobi seems to work better with Method = 1, and 
                   % Mosek with Method = 0.
INFO.tol = 1E-7; % Feasibility, optimality and convergence tolerance for Cplex (tol>=1E-9). 
                 % It is recommended it is at least 2 orders of magnitude
                 % tighter than the integrator tolerance. 
                 % If problems with infeasibility messages, tighten this
                 % tolerance.
INFO.tolPh1 = INFO.tol; % Tolerance to determine if a solution to phaseI equals zero.
                   % It is recommended to be the same as INFO.tol. 
INFO.tolevt = 10*INFO.tol; % Tolerance for event detection. Has to be greater 
                   % than INFO.tol.

% You can modify the integration tolerances here.
% If some of the flows become negative after running the simulation once
% you can add the 'Nonnegative' option.

options = odeset('AbsTol',1E-5,'RelTol',1E-5,'Events',@evts,'Mass',M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[model,INFO] = ModelSetupM(model,Y0,INFO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if INFO.LPsolver == 0
    [INFO] = LexicographicOpt(model,INFO);
elseif INFO.LPsolver == 1
    [INFO] = LexicographicOptG(model,INFO);
else
    display('Solver not currently supported.');
end

tint = 0;
TF = [];
YF = [];

while tint<tspan(2)
    
% Look at MATLAB documentation if you want to change solver.
% ode15s is more or less accurate for stiff problems. 
%     [T,Y] = ode15s(@DRHS,tspan,Y0,options,INFO);

    [T,Y] = ode15s(@DRHS_dispersion_new_HD,tspan,Y0,options,INFO);
    TF = [TF;T];
    YF = [YF;Y];
    tint = T(end);
    tspan = [tint,tspan(2)];
    Y0 = Y(end,:);
    if tint == tspan(2)
        break;
    end
    
% Determine model with basis change

    value = evts(tint,Y0,INFO);
    [jjj,j] = min(value);
    ct = 0;
    k = 0;
    while j>ct
        k = k + 1;
        ct = ct + size(model{k}.A,1);
    end
    INFO.flagbasis = k;
    fprintf('Basis change at time %d. \n',tint);
    
% Update b vector 

[INFO] = bupdate(tint,Y0,INFO);

% Perform lexicographic optimization

if INFO.LPsolver == 0
    [INFO] = LexicographicOpt(model,INFO);
elseif INFO.LPsolver == 1
    [INFO] = LexicographicOptG(model,INFO);
else
    display('Solver not currently supported.');
end

end

elapsedTime = toc;
disp(['End Time = ',num2str(elapsedTime)]);

T = TF;
Y = YF;
%%
for k=1:length(T)
    [~, coeff(k,:),mu(k,:),vc(k,:),vh(k,:),vc2(k,:),va(k,:),ve(k,:)]=DRHS_dispersion_new_HD(T(k),Y(k,:),INFO);
end

%----------------------The first node point--------------------------------
figure(1)
% 
subplot(4,4,1)
hold on
plot(T,Y(:,N*ns+1))
xlabel('Time [h]')
ylabel('Biomass [g/L]')

subplot(4,4,2)
hold on
plot(T,Y(:,1))
xlabel('Time [h]')
ylabel('CO in liquid [mmol/L]')

subplot(4,4,3)
hold on
plot(T,Y(:,2))
xlabel('Time [h]')
ylabel('CO in gas [mmol/L]')

subplot(4,4,4)
hold on
plot(T,Y(:,3))
xlabel('Time [h]')
ylabel('H2 in liquid [mmol/L]')

subplot(4,4,5)
hold on
plot(T,Y(:,4))
xlabel('Time [h]')
ylabel('H2 in gas [mmol/L]')

subplot(4,4,6)
hold on
plot(T,Y(:,5))
xlabel('Time [h]')
ylabel('CO2 in liquid [mmol/L]')

subplot(4,4,7)
hold on
plot(T,Y(:,6))
xlabel('Time [h]')
ylabel('CO2 in gas [mmol/L]')

subplot(4,4,8)
hold on
plot(T,Y(:,N*ns+2))
xlabel('Time [h]')
ylabel('Acetate [g/L]')

subplot(4,4,9)
hold on
plot(T,Y(:,N*ns+3))
xlabel('Time [h]')
ylabel('Ethanol [g/L]')

subplot(4,4,10)
hold on
plot(T,Y(:,7))
xlabel('Time [h]')
ylabel('Pressure [Pa]')

subplot(4,4,11)
hold on
plot(T,Y(:,8))
xlabel('Time [h]')
ylabel('Superficial gas velocity [m/s]')

subplot(4,4,12)
hold on
plot(T,Y(:,9))
xlabel('Time [h]')
ylabel('Bubble diameter [mm]')

subplot(4,4,13)
hold on
plot(T,Y(:,10))
xlabel('Time [h]')
ylabel('Bubble rising velocity [m/s]')

subplot(4,4,14)
hold on
plot(T,Y(:,11))
xlabel('Time [h]')
ylabel('Gas holdup [-]')


%--------------The middle node point---------------------------------------
mz = round(N/2);

figure(2)

subplot(4,4,1)
hold on
plot(T,Y(:,N*ns+1))
xlabel('Time [h]')
ylabel('Biomass [g/L]')

subplot(4,4,2)
hold on
plot(T,Y(:,ns*mz+1))
xlabel('Time [h]')
ylabel('CO in liquid [mmol/L]')

subplot(4,4,3)
hold on
plot(T,Y(:,ns*mz+2))
xlabel('Time [h]')
ylabel('CO in gas [mmol/L]')

subplot(4,4,4)
hold on
plot(T,Y(:,ns*mz+3))
xlabel('Time [h]')
ylabel('H2 in liquid [mmol/L]')

subplot(4,4,5)
hold on
plot(T,Y(:,ns*mz+4))
xlabel('Time [h]')
ylabel('H2 in gas [mmol/L]')

subplot(4,4,6)
hold on
plot(T,Y(:,ns*mz+5))
xlabel('Time [h]')
ylabel('CO2 in liquid [mmol/L]')

subplot(4,4,7)
hold on
plot(T,Y(:,ns*mz+6))
xlabel('Time [h]')
ylabel('CO2 in gas [mmol/L]')

subplot(4,4,8)
hold on
plot(T,Y(:,N*ns+2))
xlabel('Time [h]')
ylabel('Acetate [g/L]')

subplot(4,4,9)
hold on
plot(T,Y(:,N*ns+3))
xlabel('Time [h]')
ylabel('Ethanol [g/L]')

subplot(4,4,10)
hold on
plot(T,Y(:,ns*mz+7))
xlabel('Time [h]')
ylabel('Pressure [Pa]')

subplot(4,4,11)
hold on
plot(T,Y(:,ns*mz+8))
xlabel('Time [h]')
ylabel('Superficial gas velocity [m/s]')

subplot(4,4,12)
hold on
plot(T,Y(:,ns*mz+9))
xlabel('Time [h]')
ylabel('Bubble diameter [mm]')

subplot(4,4,13)
hold on
plot(T,Y(:,ns*mz+10))
xlabel('Time [h]')
ylabel('Bubble rising velocity [m/s]')

subplot(4,4,14)
hold on
plot(T,Y(:,ns*mz+11))
xlabel('Time [h]')
ylabel('Gas holdup [-]')


%-------------------Outlet-------------------------------------------------
figure(3)

subplot(4,4,1)
hold on
plot(T,Y(:,N*ns+1))
xlabel('Time [h]')
ylabel('Biomass [g/L]')

subplot(4,4,2)
hold on
plot(T,Y(:,ns*(N-1)+1))
xlabel('Time [h]')
ylabel('CO in liquid [mmol/L]')

subplot(4,4,3)
hold on
plot(T,Y(:,ns*(N-1)+2))
xlabel('Time [h]')
ylabel('CO in gas [mmol/L]')

subplot(4,4,4)
hold on
plot(T,Y(:,ns*(N-1)+3))
xlabel('Time [h]')
ylabel('H2 in liquid [mmol/L]')

subplot(4,4,5)
hold on
plot(T,Y(:,ns*(N-1)+4))
xlabel('Time [h]')
ylabel('H2 in gas [mmol/L]')

subplot(4,4,6)
hold on
plot(T,Y(:,ns*(N-1)+5))
xlabel('Time [h]')
ylabel('CO2 in liquid [mmol/L]')

subplot(4,4,7)
hold on
plot(T,Y(:,ns*(N-1)+6))
xlabel('Time [h]')
ylabel('CO2 in gas [mmol/L]')

subplot(4,4,8)
hold on
plot(T,Y(:,N*ns+2))
xlabel('Time [h]')
ylabel('Acetate [g/L]')

subplot(4,4,9)
hold on
plot(T,Y(:,N*ns+3))
xlabel('Time [h]')
ylabel('Ethanol [g/L]')

subplot(4,4,10)
hold on
plot(T,Y(:,ns*(N-1)+7))
xlabel('Time [h]')
ylabel('Pressure [Pa]')

subplot(4,4,11)
hold on
plot(T,Y(:,ns*(N-1)+8))
xlabel('Time [h]')
ylabel('Superficial gas velocity [m/s]')

subplot(4,4,12)
hold on
plot(T,Y(:,ns*(N-1)+9))
xlabel('Time [h]')
ylabel('Bubble diameter [mm]')

subplot(4,4,13)
hold on
plot(T,Y(:,ns*(N-1)+10))
xlabel('Time [h]')
ylabel('Bubble rising velocity [m/s]')

subplot(4,4,14)
hold on
plot(T,Y(:,ns*(N-1)+11))
xlabel('Time [h]')
ylabel('Gas holdup [-]')

%---------------spatial profiles versus time-------------------------------
for i=1:(nmodel-1)
    zt(i)=i*zs;
end

for i=1:nmodel
    yfin2(i)=Y(end,ns*(i-1)+1);
    yfin3(i)=Y(end,ns*(i-1)+2);
    yfin4(i)=Y(end,ns*(i-1)+3);
    yfin5(i)=Y(end,ns*(i-1)+4);
    yfin6(i)=Y(end,ns*(i-1)+5);
    yfin7(i)=Y(end,ns*(i-1)+6);
    yfin10(i) = Y(end,ns*(i-1)+7);
    yfin11(i) = Y(end,ns*(i-1)+8);
    yfin12(i) = Y(end,ns*(i-1)+9);
    yfin13(i) = Y(end,ns*(i-1)+10);
    yfin14(i) = Y(end,ns*(i-1)+11);
    yfin1(i)=Y(end,ns*N+1);
    yfin8(i)=Y(end,ns*N+2);
    yfin9(i)=Y(end,ns*N+3);
    klac(i) = 6*Y(end,ns*(i-1)+11)./(1-Y(end,ns*(i-1)+11))./(Y(end,ns*(i-1)+9)/1000)*klc;
end

zt = [0 zt];

figure(4)
subplot(4,4,1)
hold on
plot(zt,yfin1)
xlabel('Location [m]')
ylabel('Biomass [g/L]')
% xlim([0 25])
xlim([0 L])

subplot(4,4,2)
hold on
plot(zt,yfin2)
xlabel('Location [m]')
ylabel('CO in liquid [mmol/L]')
xlim([0 L])

subplot(4,4,3)
hold on
plot(zt,yfin3)
xlabel('Location [m]')
ylabel('CO in gas [mmol/L]')
xlim([0 L])

subplot(4,4,4)
hold on
plot(zt,yfin4)
xlabel('Location [m]')
ylabel('H2 in liquid [mmol/L]')
xlim([0 L])

subplot(4,4,5)
hold on
plot(zt,yfin5)
xlabel('Location [m]')
ylabel('H2 in gas [mmol/L]')
xlim([0 L])

subplot(4,4,6)
hold on
plot(zt,yfin6)
xlabel('Location [m]')
ylabel('CO2 in liquid [mmol/L]')
xlim([0 L])

subplot(4,4,7)
hold on
plot(zt,yfin7)
xlabel('Location [m]')
ylabel('CO2 in gas [mmol/L]')
xlim([0 L])

subplot(4,4,8)
hold on
plot(zt,yfin8)
xlabel('Location [m]')
ylabel('Acetate [g/L]')
xlim([0 L])

subplot(4,4,9)
hold on
plot(zt,yfin9)
xlabel('Location [m]')
ylabel('Ethanol [g/L]')
xlim([0 L])

subplot(4,4,10)
hold on
plot(zt,yfin10)
xlabel('Location [m]')
ylabel('Pressure [Pa]')

subplot(4,4,11)
hold on
plot(zt,yfin11)
xlabel('Location [m]')
ylabel('Superficial gas velocity [m/s]')

subplot(4,4,12)
hold on
plot(zt,yfin12)
xlabel('Location [m]')
ylabel('Bubble diameter [mm]')

subplot(4,4,13)
hold on
plot(zt,yfin13)
xlabel('Location [m]')
ylabel('Bubble rising velocity [m/s]')

subplot(4,4,14)
hold on
plot(zt,yfin14)
xlabel('Location [m]')
ylabel('Gas holdup [-]')

subplot(4,4,15)
hold on
plot(zt,klac)
xlabel('Location [m]')
ylabel('Volumetric mass transfer rate [1/h]')

%--------------------------flow rate plot----------------------------------

for i=1:nmodel
     Y(:,ns*(i-1)+2) = Y(:,ns*(i-1)+2).*Area.*Y(:,ns*(i-1)+8);   %mol/s
     Y(:,ns*(i-1)+4) = Y(:,ns*(i-1)+4).*Area.*Y(:,ns*(i-1)+8);
     Y(:,ns*(i-1)+6) = Y(:,ns*(i-1)+6).*Area.*Y(:,ns*(i-1)+8);
end
% liquid phase
vl_flowrate = Qmedia;
for i=1:nmodel
     Y(:,ns*(i-1)+1) = Y(:,ns*(i-1)+1)*vl_flowrate;
     Y(:,ns*(i-1)+3) = Y(:,ns*(i-1)+3)*vl_flowrate;
     Y(:,ns*(i-1)+5) = Y(:,ns*(i-1)+5)*vl_flowrate;
end


for i=1:nmodel
    yfin2_2(i)=Y(end,ns*(i-1)+1);
    yfin3_2(i)=Y(end,ns*(i-1)+2);
    yfin4_2(i)=Y(end,ns*(i-1)+3);
    yfin5_2(i)=Y(end,ns*(i-1)+4);
    yfin6_2(i)=Y(end,ns*(i-1)+5);
    yfin7_2(i)=Y(end,ns*(i-1)+6);
    yfin10_2(i) = Y(end,ns*(i-1)+7);
    yfin11_2(i) = Y(end,ns*(i-1)+8);
    yfin12_2(i) = Y(end,ns*(i-1)+9);
    yfin13_2(i) = Y(end,ns*(i-1)+10);
    yfin14_2(i) = Y(end,ns*(i-1)+11);
    yfin1_2(i)=Y(end,ns*N+1);
    yfin8_2(i)=Y(end,ns*N+2);
    yfin9_2(i)=Y(end,ns*N+3);
end

figure(5)
subplot(4,4,1)
hold on
plot(zt,yfin1_2)
xlabel('Location [m]')
ylabel('Biomass [g/L]')
% xlim([0 25])
xlim([0 L])

subplot(4,4,2)
hold on
plot(zt,yfin2_2)
xlabel('Location [m]')
ylabel('CO in liquid [mol/s]')
xlim([0 L])

subplot(4,4,3)
hold on
plot(zt,yfin3_2)
xlabel('Location [m]')
ylabel('CO in gas [mol/s]')
xlim([0 L])

subplot(4,4,4)
hold on
plot(zt,yfin4_2)
xlabel('Location [m]')
ylabel('H2 in liquid [mol/s]')
xlim([0 L])

subplot(4,4,5)
hold on
plot(zt,yfin5_2)
xlabel('Location [m]')
ylabel('H2 in gas [mol/s]')
xlim([0 L])

subplot(4,4,6)
hold on
plot(zt,yfin6_2)
xlabel('Location [m]')
ylabel('CO2 in liquid [mol/s]')
xlim([0 L])

subplot(4,4,7)
hold on
plot(zt,yfin7_2)
xlabel('Location [m]')
ylabel('CO2 in gas [mol/s]')
xlim([0 L])

subplot(4,4,8)
hold on
plot(zt,yfin8_2)
xlabel('Location [m]')
ylabel('Acetate [g/L]')
xlim([0 L])

subplot(4,4,9)
hold on
plot(zt,yfin9_2)
xlabel('Location [m]')
ylabel('Ethanol [g/L]')
xlim([0 L])

subplot(4,4,10)
hold on
plot(zt,yfin10_2)
xlabel('Location [m]')
ylabel('Pressure [Pa]')

subplot(4,4,11)
hold on
plot(zt,yfin11_2)
xlabel('Location [m]')
ylabel('Superficial gas velocity [m/s]')

subplot(4,4,12)
hold on
plot(zt,yfin12_2)
xlabel('Location [m]')
ylabel('Bubble diameter [mm]')

subplot(4,4,13)
hold on
plot(zt,yfin13_2)
xlabel('Location [m]')
ylabel('Bubble rising velocity [m/s]')

subplot(4,4,14)
hold on
plot(zt,yfin14_2)
xlabel('Location [m]')
ylabel('Gas holdup [-]')

subplot(4,4,15)
hold on
plot(zt,klac)
xlabel('Location [m]')
ylabel('Volumetric mass transfer rate [1/h]')

figure(6)

for i=1:nmodel
    mu_spatial(i)=mu(end,i);
    vc_spatial(i)=vc(end,i);
    vh_spatial(i)=vh(end,i);
    vc2_spatial(i)=vc2(end,i);
    va_spatial(i)=va(end,i);
    ve_spatial(i)=ve(end,i);
end

subplot(2,3,1)
hold on
plot(zt,mu_spatial)
xlabel('Location [m]')
ylabel('Growth rate [1/h]')
xlim([0 L])

subplot(2,3,2)
hold on
plot(zt,vc_spatial)
xlabel('Location [m]')
ylabel('CO uptake rate [mmol/g/h]')
xlim([0 L])

subplot(2,3,3)
hold on
plot(zt,vh_spatial)
xlabel('Location [m]')
ylabel('H2 uptake rate [mmol/g/h]')
xlim([0 L])

subplot(2,3,4)
hold on
plot(zt,vc2_spatial)
xlabel('Location [m]')
ylabel('CO2 uptake rate [mmol/g/h]')
xlim([0 L])

subplot(2,3,5)
hold on
plot(zt,va_spatial)
xlabel('Location [m]')
ylabel('acetate synthesis rate [mmol/g/h]')
xlim([0 L])

subplot(2,3,6)
hold on
plot(zt,ve_spatial)
xlabel('Location [m]')
ylabel('ethanol synthesis rate [mmol/g/h]')
xlim([0 L])
