
% function solving the ODE equations
% y(1) -> TC, y(2) -> FC

function [Sim_t, Sim_y] = fnSim_NCCR(input_NCCR)

T_ini = input_NCCR.T;
P_ini = input_NCCR.P;
F_ini = input_NCCR.F;
kf = input_NCCR.kf;
kr = input_NCCR.kr;
tspan = input_NCCR.t;
y_ini = input_NCCR.y0;

[Sim_t, Sim_y] = ode45(@(t, y) odeEqn_NCCR(t, y, T_ini, P_ini, F_ini, kf, kr), tspan, y_ini );

end


