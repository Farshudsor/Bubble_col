
import casadi.*

%bubble colomn opperates for [run_time*Tsamp] hrs
run_time = 500; % Time horizon
n_pred = 10; % prediction horizon
n_ctrl = 10; % number of control intervals
Tsamp = 1;   % timestemps between control actions

hrs = run_time*Tsamp;

n_st = 3;   n_ip = 1;   n_par = 6;

res_xk = zeros(hrs+1,n_st);
res_theta = zeros(run_time+1,n_par);
res_uk = zeros(run_time+1,n_ip);

% Declare model variables
x1 = SX.sym('x1');  % Ethanol
x2 = SX.sym('x2');  % Acetate
x3 = SX.sym('x3');  % Biomass

uk = SX.sym('uk');  % Dilution rate

theta11 = SX.sym('theta1');      theta12 = SX.sym('theta2');  
theta21 = SX.sym('theta3');      theta22 = SX.sym('theta1');  
theta31 = SX.sym('theta2');      theta32 = SX.sym('theta3');  

xk = [x1; x2; x3];  
theta = [theta11; theta21; theta31; theta12; theta22; theta32];

vk = SX.sym('vk',n_st);    % measurment noise

% Parameters
slt = .84;          %selectivity
slt_p = 10;         %power on selectivity soft constraint
Ma = 60.0/1000;       % Ma molecular weight g/mmol
Me = 46/1000;       % CO molecular weight g/mmol

% Constraints
xlb=[0,0,0];     xub=[inf,inf,inf];
ulb=[0.01];     uub=[.1];

% Measurment noise
vkp = rand([run_time,1])*.00;

% Model equations
dx1 = Me*(theta11*uk+theta12)*x3 - uk*x1;    
dx2 = Ma*(theta21*uk+theta22)*x3 - uk*x2;    
dx3 =    (theta31*uk+theta32)*x3 - uk*x3;    
dtheta11 = 0*theta11;    dtheta21 = 0*theta21;    dtheta31 = 0*theta31;
dtheta12 = 0*theta12;    dtheta22 = 0*theta22;    dtheta32 = 0*theta32;

sys_ode = [dx1; dx2; dx3];
mdl_ode = [dx1; dx2; dx3; dtheta11; dtheta21; dtheta31; dtheta12; ...
            + dtheta22; dtheta32];

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
uk_opt = [.06];   xkp = [0, 0, .1];
theta_par = [1.316, -0.597, -0.001, 1.8753, .83, .064];

%EKF covar matricies
Qz = eye(n_st+n_par);    R = eye(n_st)*0;     Sigmak_p = Qz;

xkh0 = xkp;
zkh0 =vertcat(xkh0',theta_par');

% Storing initial conditions
res_uk(1,:) = uk_opt;
res_xk(1,:) = xkp;
res_theta(1,:) = theta_par;

Yend = []; %all HF model states for next integration step

% Set IPOPT conditions
opts = struct;
opts.ipopt.max_iter = 1000;
opts.ipopt.print_level = 5;
opts.ipopt.output_file = 'Main_out.txt';

tic()

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
    [Jce, qu_ce, lbq, ubq, g, lbg, ubg, qu_init] = prediction(F_ode,...
                            + n_pred,n_ctrl,n_st,n_par,n_ip,ulb,uub,...
                            + xlb,xub,xk,theta_par,slt,Tsamp,xkh0,uk_opt,slt_p);
    %Formulate nlp & solver
    prob = struct('x',qu_ce, 'f',Jce, 'g',g);
    solver = nlpsol('solver_mpc', 'ipopt', prob, opts);
    %Solve
    res_mpc = solver('x0',qu_init,'lbx',lbq,'ubx',ubq,'lbg',lbg,'ubg',ubg);
    
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
    res_xk((k-1)*Tsamp+1:k*Tsamp,:) =xki(:,:) ;
    res_theta(k+1,:) =full(theta_par);
    
    '************'
    (k*Tsamp)
    '************'
    xkp
    '************'
    
    if k == 100
        'break'
    end
    

end

solve_time = toc();

%% Save results
if 1
    results_loaded = load('results.mat');
    var_names = {'Ce','Ca','Cx', 'D', 'T11', 'T21', 'T31','T12', 'T22', 'T32', 'v', 'solve_time'};

    if 0
        results = table(res_xk(:,1)',res_xk(:,2)',res_xk(:,3)',res_uk',...
                        +res_theta(:,1)', res_theta(:,2)', res_theta(:,3)',...
                        +res_theta(:,4)', res_theta(:,5)', res_theta(:,6)',...
                        +vkp', solve_time,'VariableNames',var_names);
        results_loaded.results = table2struct(results);

    else
        add = struct('Ce',res_xk(:,1)','Ca',res_xk(:,2)','Cx',...
            + res_xk(:,3)','D',res_uk','T11',res_theta(:,1)','T21', ...
            + res_theta(:,2)','T31', res_theta(:,3)','T12',...
            + res_theta(:,4)','T22', res_theta(:,5)','T32', ...
            + res_theta(:,6)', 'v', vkp','solve_time', solve_time);
        results_loaded.results = [results_loaded.results;add];
    end

    save('results.mat', '-struct', 'results_loaded');
end

%% Plotting

if 1
    figure(1)
    plot([0:hrs],res_xk(:,3))
    xlabel('Time [hr]')
    ylabel('Biomass [g/L]')
    xlim([0 hrs])

    figure(2)
    plot([0:hrs],res_xk(:,2))
    xlabel('Time [hr]')
    ylabel('acetate [g/L]')
    xlim([0 hrs])

    figure(3)
    plot([0:hrs],res_xk(:,1))
    xlabel('Time [hr]')
    ylabel('Ethanol [g/L]')
    xlim([0 hrs])

    figure(4)
    plot([1:Tsamp:hrs],res_uk(2:end))
    xlabel('Time [hr]')
    ylabel('Dilution ')
    xlim([0 hrs])

    figure(5)
    hold on 
    subplot(3,2,1)
        plot([1:Tsamp:hrs],res_theta(2:end,1))
        xlabel('Time [hr]')
        ylabel('Theta 1,1 - v_e ')
        xlim([0 hrs])

    subplot(3,2,3)
        plot([1:Tsamp:hrs],res_theta(2:end,2))
        xlabel('Time [hr]')
        ylabel('Theta 2,1 - v_a ')
        xlim([0 hrs])

    subplot(3,2,5)
        plot([1:Tsamp:hrs],res_theta(2:end,3))
        xlabel('Time [hr]')
        ylabel('Theta 3,1 - \mu ')
        xlim([0 hrs])
        
    subplot(3,2,2)
        plot([1:Tsamp:hrs],res_theta(2:end,4))
        xlabel('Time [hr]')
        ylabel('Theta 1,2 - v_e ')
        xlim([0 hrs])

    subplot(3,2,4)
        plot([1:Tsamp:hrs],res_theta(2:end,5))
        xlabel('Time [hr]')
        ylabel('Theta 2,2 - v_a ')
        xlim([0 hrs])

    subplot(3,2,6)
        plot([1:Tsamp:hrs],res_theta(2:end,6))
        xlabel('Time [hr]')
        ylabel('Theta 3,2 - \mu ')
        xlim([0 hrs])
        
    figure(6)
    hold on
    plot([0:hrs],res_xk(:,1)./res_xk(:,2))
    plot([0:hrs], ones(hrs+1,1)*slt)
    xlabel('Time [hr]')
    ylabel('selectivity [-]')
    xlim([0 hrs])
end








