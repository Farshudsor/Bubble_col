
import casadi.*

%bubble colomn opperates for [run_time*Tsamp] hrs
run_time = 100; % Time horizon
n_pred = 100; % prediction horizon
n_ctrl = 100; % number of control intervals
Tsamp = 10;   % timestemps between control actions

hrs = run_time*Tsamp;

n_st = 3;   n_ip = 2;   n_par = 3;

res_xk = zeros(run_time+1,n_st);
res_theta = zeros(hrs+1,n_par);
res_uk = zeros(run_time+1,n_ip);

% Declare model variables
x1 = SX.sym('x1');  % Ethanol
x2 = SX.sym('x2');  % Acetate
x3 = SX.sym('x3');  % Biomass

u1 = SX.sym('u1');  % gas flow rate
u2 = SX.sym('u2');  % Dilution rate

theta1 = SX.sym('theta1');
theta2 = SX.sym('theta2');
theta3 = SX.sym('theta3');

xk = [x1; x2; x3];  uk = [u1; u2];  theta = [theta1; theta2; theta3];

vk = SX.sym('vk',n_st);    % measurment noise

% Parameters
Ma = 60.0/1000;       % Ma molecular weight g/mmol
Me = 46/1000;       % CO molecular weight g/mmol

theta0 = [1, 1, 1];     x0 = [0.0; 0.0; 0.1];   u0 = [12.8; 0.06];

% Constraints
xlb=[0,0,0];     xub=[inf,inf,inf];
ulb=[11,0];     uub=[14,1];

% Measurment noise
vkp = rand([run_time,1])*.01;

% Model equations
dx1 = Me*(theta1*u1)*x3 - u2*x1;    dtheta1 = 0*theta1;
dx2 = Ma*(theta2*u1)*x3 - u2*x2;    dtheta2 = 0*theta2;
dx3 =    (theta3*u1)*x3 - u2*x3;    dtheta3 = 0*theta3;

sys_ode = [dx1; dx2; dx3];
mdl_ode = [dx1; dx2; dx3; dtheta1; dtheta2; dtheta3];

%mdl_ode2 =vertcat(de1+wk2[0],de2+wk2[1],de3+wk2[2],de4+wk2[3],de5+wk2[4],de6+wk2[5],de7+wk2[6],dthetah1,dthetah2)
%L2   = Function('L2' ,[xk2,uk2,thetah2,wk2],[jacobian((vertcat(xk2,thetah2)+mdl_ode2),wk2)])

% Jacobians for predictions
Jx = Function('Jx',{xk,uk,theta},{jacobian((sys_ode),xk)});
Jz = Function('Jz',{xk,uk,theta},{jacobian((vertcat(xk,theta)+mdl_ode),vertcat(xk,theta))});
Ju = Function('Ju',{xk,uk,theta},{jacobian((vertcat(xk,theta)+mdl_ode),uk)});

% Integrators
ode = struct('x',xk, 'p',vertcat(theta,uk), 'ode',sys_ode);
F_ode = integrator('F_ode', 'cvodes', ode);

m_ode = struct('x',vertcat(xk,theta), 'p',uk, 'ode',mdl_ode);
M_ode = integrator('M_ode', 'cvodes', m_ode);

% Output equation (Observations)
C = eye(n_st);
fy = C*xk + vk;
h = Function('h',{xk,vk}, {fy});
Ch = Function('Ch',{vertcat(xk,theta)},{jacobian(fy,vertcat(xk,theta))});
Chx = Function('Chx',{xk},{jacobian(fy,xk)});

% Initial values
uk_opt = [12.3, .06];   xkp = [0, 0, .1];
theta_par = [.0001, .0001, 0.001];

%EKF covar matricies
Qz = eye(6);    R = eye(3);     Sigmak_p = Qz;

xkh0 = xkp;
zkh0 =vertcat(xkh0',theta_par');

% Storing initial conditions
res_uk(1,:) = uk_opt;
res_xk(1,:) = xkp;
res_theta(1,:) = theta_par;

Yend = []; %full (

for k = 1:run_time

    time(k) = k;
    
    if false    %for batch sys
        time_remain = run_time - k;
        if time_remain < n_pred
            n_pred = time_remain;
        end
    else        %for cont. sys
        time_remain = run_time; 
    end
        
    %Compute the measurement
    ykp = h(xkp,vkp(k,:)');
    
    %EKF update step/measurment
    Czh = Ch(zkh0);
    Kkh = mtimes(Sigmak_p,mtimes(Czh',inv(mtimes(Czh,mtimes(Sigmak_p,Czh')) + R)));
    zkh0 = zkh0 + mtimes(Kkh,(ykp - h(zkh0(1:n_st), zeros(n_st,1))));
    xkh0 = zkh0(1:n_st);
    theta_par = zkh0(n_st+1:end);
    Sigmak = mtimes((eye(n_st+n_par) - mtimes(Kkh,Czh)),Sigmak_p);

    %Generates predictions 
    [Jce, qu_ce, lbq, ubq, g, lbg, ubg, qu_init] = prediction(F_ode,n_pred,n_ctrl,n_st,n_par,n_ip,ulb,uub,xlb,xub,xk,theta_par,uk,Tsamp,xkh0,uk_opt);
    %Formulate into nlp,
    prob = struct('x',qu_ce, 'f',Jce, 'g',g);
    solver_mpc = nlpsol('solver_mpc', 'ipopt', prob);
    %Solve
    res_mpc = solver_mpc('x0',qu_init, 'lbx',lbq, 'ubx',ubq, 'lbg',lbg, 'ubg',ubg);
    
    %Take first optimal control action
    uk_ce = res_mpc.x;
    uk_opt = uk_ce(n_st+1:n_st+n_ip)';

    %Simulate the system
    [xki, Yend] = main(uk_opt, xkp, Yend);
    xkp = xki(end, :);   
   
    %EKF prediction
    Az = Jz(xkh0,uk_opt,theta_par);
    %Jw = L2(xkh0,uk_opt,theta_par,wkp(k,:));
    z_end = M_ode('x0',zkh0, 'p',uk_opt);
    zkh0 = z_end.xf.full();
    %Sigmak_p = mtimes(Az,mtimes(Sigmak,Az')) + Jw*Qz*Jw';
    Sigmak_p = mtimes(Az,mtimes(Sigmak,Az')) + Qz;

    % Save results
    res_uk(k+1,:) =full(uk_opt) ;
    res_xk((k-1)*Tsamp+1:k*Tsamp,:) =xki ;
    res_theta(k+1,:) =full(theta_par);

end
    

%solve_time[1] =  Time.time() - start











