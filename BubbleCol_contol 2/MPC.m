import casadi.*

T = 599; % Time horizon
N = 6; % number of control intervals

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x = [x1; x2; x3];

u1 = SX.sym('u1');
u2 = SX.sym('u2');
u = [u1; u2];

theta1 = SX.sym('theta1');
theta2 = SX.sym('theta2');
theta3 = SX.sym('theta3');
theta = [theta1; theta2; theta3];

% Parameters
Me = 10;
Ma = 10;
theta0 = [10, 10, 10];
x0 = [0.2; 0.3; .2];
u0 = [12.8; 0.4];

% Model equations
ve = theta1*u1;
va = theta2*u2;
mu = theta3*u3;

dx1 = Me*ve*x3 - u2*x1;
dx2 = Ma*va*x3 - u2*x2;
dx3 = mu*x3 - u2*x3;
xdot = [dx1; dx2; dx3];

% Objective term
L = -x1^2 ;% + x2^2 + u^2;

% Continuous time dynamics
f = Function('f', {x, theta, u}, {xdot, L});

% Formulate discrete time dynamics
% CVODES from the SUNDIALS suite
dae = struct('x',vertcat(x,theta),'p',u,'ode',xdot,'quad',L);
opts = struct('tf',T/N);
F = integrator('F', 'cvodes', dae, opts);

% Evaluate at a test point
Fk = F('x0',[0.2; 0.3, .2],'p',[12.8; 0.4]);
disp(Fk.xf)
disp(Fk.qf)

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
Xk = MX.sym('X0', 3);
w = {w{:}, Xk};
lbw = [lbw; xlb];
ubw = [ubw; xub];
w0 = [w0; x0];

% Formulate the NLP
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],2);
    w = {w{:}, Uk};
    lbw = [lbw; ulb];
    ubw = [ubw; uub];
    w0 = [w0;  u0];

    % Integrate till the end of the interval
    Fk = F('x0', Xk, 'p', Uk);
    Xk_end = Fk.xf;
    J=J+Fk.qf;

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], 3);
    w = [w, {Xk}];
    lbw = [lbw; xlb];
    ubw = [ubw; xub];
    w0 = [w0; x0];

    % Add equality constraint
    g = [g, {Xk_end-Xk}];
    lbg = [lbg; 0; 0; 0];
    ubg = [ubg; 0; 0; 0];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

%% Plot the solution
x1_opt = w_opt(1:3:end);
x2_opt = w_opt(2:3:end);
u_opt = w_opt(3:3:end);
tgrid = linspace(0, T, N+1);
clf;

hold on
plot(tgrid, x1_opt, '--')
plot(tgrid, x2_opt, '-')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x1','x2','u')
