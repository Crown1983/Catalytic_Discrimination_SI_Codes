
% ODEs of Toehold-Exchange Reaction Kinetic Model

% Rxn-1 ->         T    +   CP    <>  CT   +   P                   ; k_1, k_-1
%              (T_ini-x)  (20-x)  <>  (x)  (P_ini+x)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function y_sim = fnFit_TER(k_assume, inTER)

% ---------------------------------- ODE of Toehold-Exchange ---------------------------------------

T_ini = inTER.T;
P_ini = inTER.P;
TC_ini = inTER.TC;
tspan = inTER.t;

kf = k_assume(1); % assumed forward reaction rate
kr = k_assume(2); % kf/K, reverse reaction rate

f = @(t, y) kf*(T_ini-y)*(1-y) - kr*y*(P_ini+y);
[~, y_sim] = ode45(f, tspan, TC_ini);

end

