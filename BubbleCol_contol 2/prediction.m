function [Jce, qu_ce, lbq, ubq, g, lbg, ubg, qu_init] = prediction(F_ode,n_pred,n_ctrl,n_st,n_par,n_ip,ulb,uub,xlb,xub,xk,theta_par,uk,Tsamp,xkh0,uk_opt)

import casadi.*

Jce = 0.0;
qu_ce = [];

g = [];
lbg = [];
ubg = [];
%lbq = [];
%ubq=[];

X0 = MX.sym('X0', n_st);
qu_ce = X0;
lbq = xkh0;
ubq = xkh0;
qu_init = xkh0;
Xk = X0;

for i = 1:n_pred
    if i ~= 1
        Uk_ = Uk;
    else
        Uk_ = uk_opt;
    end
    
    if mod(i-1,Tsamp) == 0
        Uk = MX.sym(['U_', num2str(i)],n_ip);
        qu_ce = vertcat(qu_ce, Uk);
        lbq = vertcat(lbq, ulb');
        ubq = vertcat(ubq, uub');
        qu_init = vertcat(qu_init, uk_opt');
    end


    x_end = F_ode('x0',Xk, 'p',vertcat(theta_par,Uk));
    xk_end = x_end.xf;
    Jce = Jce + 0.5*xk_end(1)^2;

    Xk = MX.sym(['X_' num2str(i+1)],n_st);
    qu_ce = vertcat(qu_ce, Xk);
    lbq = vertcat(lbq, xlb');
    ubq = vertcat(ubq, xub');
    qu_init = vertcat(qu_init, xkh0);

    g = vertcat(g, (xk_end-Xk));
    lbg = [lbg, xlb];
    ubg = [ubg, xub];
end

