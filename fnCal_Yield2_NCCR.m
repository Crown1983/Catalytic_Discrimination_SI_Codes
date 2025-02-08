
% Help
% 
% Description  -> Calculate the yield of a single Toehold-exchange reaction
% 
%         Rxn-1       T       +      CP      <->     TC      +      P           ;    K
%                   t-y(1)        1-y(1)-y(2)       y(1)        p+y(1)+y(2)

%         Rxn-2       F       +      TC      <->     FC      +      T           ;   1/K
%                   f-y(2)          y(1)            y(2)          t-y(1)
% 
%                   K = exp(-(dG+ddG)/(R*T) )
% 
% Fn_Input     -> Input(1) = [T0]:  <Array> Conc. of T for NCCR rxn
%              -> Input(2) = ddG:  <Array> SNV free energy differences
% 
% Fn_Parameter -> Param(1) = dG:  <Scalar> Free energy of TER rxn
%              -> Param(2) = p:  <Array> Dimensionless [free P]
%              -> Param(3) = Tmp:  <Scalar> Temperature in Celcius unit
% 
% Fn Output    -> Output(1) = y(1):  <Matrix> Yield of `TC`
%              -> Output(2) = y(2):  <Matrix> Yield of `FC`
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [y2_sim_TC, y2_sim_FC, y2_sim] = fnCal_Yield2_NCCR(Input2, Param2)

dG_TER = Param2.dG;
P_ini = Param2.p0;
F_ini = Param2.f0;
Tmp = Param2.Tmp;

R = 1.987*10^-3; % kcal*K^-1*mol^-1
T_kel = 273.15 + Tmp; % temperatue K

dG_TER_arr = reshape(dG_TER + Input2.ddG, [length(Input2.ddG) ,1 ] ); % Column of ddG values
K_TER_arr = exp(-dG_TER_arr./(R*T_kel) );

[T_ini_xx, K_TER_yy] = meshgrid(Input2.t0, K_TER_arr);

y2_sim_TC = zeros(size(T_ini_xx) );
y2_sim_FC = zeros(size(T_ini_xx) );
y2_sim = zeros(size(T_ini_xx) );

syms x l
for i = 1: numel(T_ini_xx)

    Eqn_CalYield_1 = x*(P_ini+x+l)/((T_ini_xx(i)-x)*(1-x-l) ) == K_TER_yy(i);
    Eqn_CalYield_2 = l*(T_ini_xx(i)-x)/(x*(F_ini-l) ) == 1/K_TER_yy(i);

    sol = vpasolve([Eqn_CalYield_1, Eqn_CalYield_2], [x, l]);
    y2_sim_TC(i) = sol.x(sol.x > 0 & sol.x < min(T_ini_xx(i), 1) & (sol.x+sol.l) < 1 );
    y2_sim_FC(i) = sol.l(sol.l > 0 & sol.l < min(F_ini, 1) & (sol.x+sol.l) < 1 );
    
    y2_sim(i) = y2_sim_TC(i) + y2_sim_FC(i);

end; clc

end


