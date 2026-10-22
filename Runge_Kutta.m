function u_n1 = Runge_Kutta(u_n,t,dt,g,restr)
% This function let the solution propagate of a step dt
% for the equation du/dt = g(t,u)

% u_n function of all positions at time t_n = t
% u_n1 function of all positions at time t_{n+1} = t+dt

if exist('restr','var') && ~isempty(restr)
    
    u_n1 = u_n;
    u_tmp = u_n;
    
    k1 = g(t,u_n);

    u_tmp(restr,:,:) =  u_n(restr,:,:)+0.5*dt*k1;
    k2 = g(t+0.5*dt, u_tmp);

    u_tmp(restr,:,:) =  u_n(restr,:,:)+0.5*dt*k2;
    k3 = g(t+0.5*dt, u_tmp);

    u_tmp(restr,:,:) =  u_n(restr,:,:)+dt*k3;
    k4 = g(t, u_tmp);
    u_n1(restr,:,:) = u_n(restr,:,:) + (dt/6)*(k1+2*k2+2*k3+k4);

else    
    k1 = g(t,u_n);
    k2 = g(t+0.5*dt, u_n+0.5*dt*k1);
    k3 = g(t+0.5*dt, u_n+0.5*dt*k2);
    k4 = g(t, u_n+k3*dt);
    % u_n1 = u_n;
    % u_n1(:,2:end-1) = u_n(:,2:end-1) + (dt/6)*(k1(:,2:end-1)+2*k2(:,2:end-1)+2*k3(:,2:end-1)+k4(:,2:end-1));
    u_n1 = u_n + (dt/6)*(k1+2*k2+2*k3+k4);
end

