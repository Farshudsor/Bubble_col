%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Jose A. Gomez and Kai Höffner                                %
%                                                                         % 
% This code can only be used for academic purposes. When using this code  %
% please cite:                                                            %
%                                                                         %
% Gomez, J.A., Höffner, K. and Barton, P. I.                              %
% DFBAlab: A fast and reliable MATLAB code for Dynamic Flux Balance       %
% Analysis. Submitted.                                                    %
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dy, coeff, mu, vc, vh, vc2, va, ve] = DRHS_dispersion_new(t, y, INFO)

% Print time

t

% Assign values

nmodel = INFO.nmodel;
N = INFO.N;
ns = INFO.ns;
condit = INFO.condit;

% Set conditions
klc = condit(1);
Hc = condit(2);
klh = condit(3);
Hh = condit(4);
klc2 = condit(5);
Hc2 = condit(6);
tr = condit(7);
zs = condit(8);
N = condit(9);
ns = condit(10);
ul = condit(11);
ug = condit(12);
dl = condit(13);
pc = condit(14);
ph = condit(15);
pc2 = condit(16);
dR = condit(17);
pL = condit(18);
eg = condit(19);
el = condit(20);
Qmedia = condit(21);
D = condit(22);
cgi = condit(23);
hgi = condit(24);
c2gi = condit(25);
db = condit(26);
Area = condit(27);
pn = condit(28);
p_atm = condit(29);

% Set other constants
Ma = 60.0/1000;       % Ma molecular weight g/mmol
Me = 46.1/1000;       % CO molecular weight g/mmol
klac = 6*eg/(1-eg)/(db/1000)*klc;
klah = klac;
klac2 = klac;

% preallocating
cl = zeros(1,N);
cg = zeros(1,N);
hl = zeros(1,N);
hg = zeros(1,N);
c2l = zeros(1,N);
c2g = zeros(1,N);
mu = zeros(1,nmodel);
vc = zeros(1,nmodel);
vh = zeros(1,nmodel);
vc2 = zeros(1,nmodel);
va = zeros(1,nmodel);
ve = zeros(1,nmodel);

% Define extracellular state variables
% same value, that of Y0, at all node points
for i=1:N
    cl(i) = y(1+(i-1)*ns);   % liquid CO concentration, mmol/L
    cg(i) = y(2+(i-1)*ns);   % gas CO concentration, mmol/L
    hl(i) = y(3+(i-1)*ns);   % liquid H2 concentration, mmol/L
    hg(i) = y(4+(i-1)*ns);   % gas H2 concentration, mmol/L
    c2l(i) = y(5+(i-1)*ns);  % liquid CO2 concentration, mmol/L
    c2g(i) = y(6+(i-1)*ns);  % gas CO2 concentration, mmol/L
end

Z = y(N*ns+1);
A = y(N*ns+2);
E = y(N*ns+3);

% The elements of the flux matrix have the sign given to them by the
% coefficients in the Cost vector in main. 
% Example, if:
% C{k}(i).rxns = [144, 832, 931];
% C{k}(i).wts = [3, 1, -1];
% Then the cost vector for this LP will be:
% flux(k,i) = 3*v_144 + v_832 - v_931 
% The penalty is an array containing the feasibility LP problem optimal
% objective function value for each model. 

%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);
%%

for i = 1:nmodel
    mu(i) = flux(i,1);      % growth rate, h-1
    vc(i) = flux(i,2);      % CO uptake rate, mmol/g/h
    vh(i) = flux(i,3);      % H2 uptake rate, mmol/g/h
    vc2(i) = flux(i,4);     % CO2 uptake rate, mmol/g/h  
    va(i) = flux(i,5);      % acetate synthesis rate, mmol/g/h 
    ve(i) = flux(i,6);      % ethanol synthesis rate, mmol/g/h
end

amu = mean(mu);
ava = mean(va);
ave = mean(ve);

%% Dynamics

dy = [];

i=1:N;
    
% Saturation gas concentrations
    p(i) = pL + 1000*9.81*zs*(N-i)*el;
    cls(i) = cg(i)*8.314*tr*Hc/1.013e5*1000.*(p(i)/1.013e5);
    hls(i) = hg(i)*8.314*tr*Hh/1.013e5*1000.*(p(i)/1.013e5);
    c2ls(i) = c2g(i)*8.314*tr*Hc2/1.013e5*1000.*(p(i)/1.013e5);
   
% Extracellular balances
        
i=1;
        cgd(i) = (cg(i) - cgi)/zs;
        hgd(i) = (hg(i) - hgi)/zs;
        c2gd(i) = (c2g(i)-c2gi)/zs;
        dcg(i) = -klac*(cls(i)-cl(i))/eg - ug*cgd(i)/eg;
        dhg(i) = -klah*(hls(i)-hl(i))/eg - ug*hgd(i)/eg;
        dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg - ug*c2gd(i)/eg;

        cld(i) = (cl(i+1)-cl(i))/zs;
        hld(i) = (hl(i+1)-hl(i))/zs;
        c2ld(i) = (c2l(i+1)-c2l(i))/zs;
        cld2(i) = (cl(i+1)-cl(i))/zs^2;
        hld2(i) = (hl(i+1)-hl(i))/zs^2;
        c2ld2(i) = (c2l(i+1)-c2l(i))/zs^2;
        dcl(i) = vc(i)*Z + klac*(cls(i)-cl(i))/el - ul*cld(i)/el + dl*cld2(i);
        dhl(i) = vh(i)*Z + klah*(hls(i)-hl(i))/el - ul*hld(i)/el + dl*hld2(i);
        dc2l(i) = vc2(i)*Z + klac2*(c2ls(i)-c2l(i))/el - ul*c2ld(i)/el + dl*c2ld2(i);
       
i=2;
        cgd(i) = (cg(i) - cg(i-1))/zs;
        hgd(i) = (hg(i) - hg(i-1))/zs;
        c2gd(i) = (c2g(i)-c2g(i-1))/zs;
        dcg(i) = -klac*(cls(i)-cl(i))/eg - ug*cgd(i)/eg;
        dhg(i) = -klah*(hls(i)-hl(i))/eg - ug*hgd(i)/eg;
        dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg - ug*c2gd(i)/eg;
       
        cld(i) = (cl(i+1)-cl(i))/zs;
        hld(i) = (hl(i+1)-hl(i))/zs;
        c2ld(i) = (c2l(i+1)-c2l(i))/zs;
        cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
        hld2(i) = (hl(i+1)-2*hl(i)+hl(i-1))/zs^2;
        c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
        dcl(i) = vc(i)*Z + klac*(cls(i)-cl(i))/el - ul*cld(i)/el + dl*cld2(i);
        dhl(i) = vh(i)*Z + klah*(hls(i)-hl(i))/el - ul*hld(i)/el + dl*hld2(i);
        dc2l(i) = vc2(i)*Z + klac2*(c2ls(i)-c2l(i))/el - ul*c2ld(i)/el + dl*c2ld2(i);
       
i=3;
        cgd(i) = (cg(i) - cg(i-1))/zs;
        hgd(i) = (hg(i) - hg(i-1))/zs;
        c2gd(i) = (c2g(i)-c2g(i-1))/zs;
        
        dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;   
        dhg(i) = -klah*(hls(i)-hl(i))/eg-ug*hgd(i)/eg; 
        dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg-ug*c2gd(i)/eg;

        cld(i) = (cl(i+1)-cl(i))/zs;
        hld(i) = (hl(i+1)-hl(i))/zs;
        c2ld(i) = (c2l(i+1)-c2l(i))/zs;
        cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
        hld2(i) = (hl(i+1)-2*hl(i)+hl(i-1))/zs^2;
        c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
        dcl(i) = vc(i)*Z + klac*(cls(i)-cl(i))/el - ul*cld(i)/el + dl*cld2(i);
        dhl(i) = vh(i)*Z + klah*(hls(i)-hl(i))/el - ul*hld(i)/el + dl*hld2(i);
        dc2l(i) = vc2(i)*Z + klac2*(c2ls(i)-c2l(i))/el - ul*c2ld(i)/el + dl*c2ld2(i);

i=4:N-2;
        cgd(i) = (cg(i) - cg(i-1))/zs;
        hgd(i) = (hg(i) - hg(i-1))/zs;
        c2gd(i) = (c2g(i)-c2g(i-1))/zs;
        dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;     
        dhg(i) = -klah*(hls(i)-hl(i))/eg-ug*hgd(i)/eg;
        dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg-ug*c2gd(i)/eg;

        cld(i) = (cl(i+1)-cl(i))/zs;
        hld(i) = (hl(i+1)-hl(i))/zs;
        c2ld(i) = (c2l(i+1)-c2l(i))/zs;
        cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
        hld2(i) = (hl(i+1)-2*hl(i)+hl(i-1))/zs^2;
        c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
        dcl(i) = vc(i)*Z + klac*(cls(i)-cl(i))/el - ul*cld(i)/el + dl*cld2(i);
        dhl(i) = vh(i)*Z + klah*(hls(i)-hl(i))/el - ul*hld(i)/el + dl*hld2(i);
        dc2l(i) = vc2(i)*Z + klac2*(c2ls(i)-c2l(i))/el - ul*c2ld(i)/el + dl*c2ld2(i);

i=N-1;
        cgd(i) = (cg(i) - cg(i-1))/zs;
        hgd(i) = (hg(i) - hg(i-1))/zs;
        c2gd(i) = (c2g(i)-c2g(i-1))/zs;
        dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;     
        dhg(i) = -klah*(hls(i)-hl(i))/eg-ug*hgd(i)/eg;
        dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg-ug*c2gd(i)/eg;

        cld(i) = (cl(i+1)-cl(i))/zs;
        hld(i) = (hl(i+1)-hl(i))/zs;
        c2ld(i) = (c2l(i+1)-c2l(i))/zs;
        cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
        hld2(i) = (hl(i+1)-2*hl(i)+hl(i-1))/zs^2;
        c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
        dcl(i) = vc(i)*Z + klac*(cls(i)-cl(i))/el-ul*cld(i)/el+dl*cld2(i);
        dhl(i) = vh(i)*Z + klah*(hls(i)-hl(i))/el-ul*hld(i)/el+dl*hld2(i);
        dc2l(i) = vc2(i)*Z + klac2*(c2ls(i)-c2l(i))/el-ul*c2ld(i)/el+dl*c2ld2(i);

i=N;

        cgd(i) = (cg(i)-cg(i-1))/zs;      
        hgd(i) = (hg(i)-hg(i-1))/zs; 
        c2gd(i) = (c2g(i)-c2g(i-1))/zs;
             
        clr = cl(1)*(abs(ul)*Area-Qmedia)/(abs(ul)*Area);  % cl(N+1) calculate from CO reci
        hlr = hl(1)*(abs(ul)*Area-Qmedia)/(abs(ul)*Area);
        c2lr = c2l(1)*(abs(ul)*Area-Qmedia)/(abs(ul)*Area);  % c2l(N+1)

        cln = (ul*zs*clr/el/dl-cl(i))/(ul*zs/el/dl-1);    % cl(N)
        hln = (ul*zs*hlr/el/dl-hl(i))/(ul*zs/el/dl-1);
        c2ln = (ul*zs*c2lr/el/dl-c2l(i))/(ul*zs/el/dl-1); % c2l(N)
        coeff=[cln, hln, c2ln];
       
        cld(i) = (cln-cl(i))/zs;
        hld(i) = (hln-hl(i))/zs;
        c2ld(i) = (c2ln-c2l(i))/zs;
        cld2(i) = (cln-2*cl(i)+cl(i-1))/zs^2;
        hld2(i) = (hln-2*hl(i)+hl(i-1))/zs^2;
        c2ld2(i) = (c2ln-2*c2l(i)+c2l(i-1))/zs^2;     
        
        dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;    
        dhg(i) = -klah*(hls(i)-hl(i))/eg-ug*hgd(i)/eg;
        dc2g(i) = -klac2*(c2ls(i)-c2ln)/eg-ug*c2gd(i);

        dcl(i) = vc(i)*Z + klac*(cls(i)-cl(i))/el-ul*cld(i)/el+dl*cld2(i);
        dhl(i) = vh(i)*Z + klah*(hls(i)-hl(i))/el-ul*hld(i)/el+dl*hld2(i);
        dc2l(i) = vc2(i)*Z + klac2*(c2ls(i)-c2l(i))/el-ul*c2ld(i)/el+dl*c2ld2(i);

        
i=1:N;
dyTmp = [dcl(i); dcg(i); dhl(i); dhg(i); dc2l(i); dc2g(i)];
dyTmp1 = reshape(dyTmp,[],1);
dy=[dy dyTmp1'];

dZ =  amu*Z - D*Z;
dA =  Ma*ava*Z - D*A;               
dE =  Me*ave*Z - D*E; 

dy = [dy dZ dA dE];

dylast=0; % for penalty value
for i=1:nmodel
    dylast = dylast + penalty(i);
end
dy = [dy dylast];
dy = dy';
end