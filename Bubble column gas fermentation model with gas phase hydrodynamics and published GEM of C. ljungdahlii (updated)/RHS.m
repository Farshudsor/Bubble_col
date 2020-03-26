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

function [lb,ub] = RHS( t,y,INFO )

% This subroutine updates the upper and lower bounds for the fluxes in the
% exID arrays in main. The output should be two matrices, lb and ub. The lb matrix
% contains the lower bounds for exID{i} in the ith row in the same order as
% exID. The same is true for the upper bounds in the ub matrix.
% Infinity can be used for unconstrained variables, however, it should be 
% fixed for all time. 

% Assign value

nmodel = INFO.nmodel;
ns = INFO.ns;
param = INFO.param;

v_cm = param(1);
Kc   = param(2);
v_hm = param(3);
Kh   = param(4);
v_c2m = param(5);
Kc2   = param(6);
Kic = param(7);
Amax = param(8);
Emax = param(9);

% Conpute uptake rate bounds

for j = 1:nmodel
    
    % CO
    lb(j,1) = -v_cm*max([y(1+ns*(j-1)) 0])/(Kc + y(1+ns*(j-1)) + y(1+ns*(j-1))*y(1+ns*(j-1))/Kic)*max([(1 - y(2+ns*nmodel)/Amax) 0])*max([(1 - y(3+ns*nmodel)/Emax) 0]);
    ub(j,1) = 0;
    
    % CO2 
    lb(j,2) = -v_c2m*max([y(5+ns*(j-1)) 0])/(Kc2 + y(5+ns*(j-1)))*max([(1 - y(2+ns*nmodel)/Amax) 0])*max([(1 - y(3+ns*nmodel)/Emax) 0]);
    ub(j,2) = Inf;

    % H2
    lb(j,3) = -v_hm*max([y(3+ns*(j-1)) 0])/(Kh + y(3+ns*(j-1)))*max([(1 - y(2+ns*nmodel)/Amax) 0])*max([(1 - y(3+ns*nmodel)/Emax) 0]);
    ub(j,3) = 0;
    
    % Biomass
    lb(j,4) = 0;
    ub(j,4) = Inf;
    
    % Ethanol
    lb(j,5) = 0;
    ub(j,5) = Inf;
    
    % Acetate
    lb(j,6) = 0;
    ub(j,6) = Inf;
 
end

end

