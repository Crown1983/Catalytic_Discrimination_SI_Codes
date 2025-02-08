
% ODEs of Noncovalent Catalysis Kinetic Model

% Rxn-1 ->            T       +      CP      <>     TC      +      P           ;    kf, kr
%                 T_ini-y(1)     1-y(1)-y(2)        y(1)     P_ini+y(1)+y(2)

% Rxn-2 ->            F       +      TC      <>     FC      +      T           ;    kr, kf
%                 F_ini-y(2)         y(1)           y(2)       T_ini-y(1)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function dydt = odeEqn_NCCR(t, y, T_ini, P_ini, F_ini, kf, kr)

% ---------------------------- ODEs of Noncovalent Catalysis Network -------------------------------

% T  = T_ini-y(1);
% CP = 1-y(1)-y(2);
% TC = y(1);
% P  = P_ini+y(1)+y(2);
% F  = F_ini-y(2);
% FC = y(2);
% 
% dydt_TC =  kf*T*CP - kr*TC*P - kr*F*TC + kf*FC*T;
% dydt_FC =  kr*F*TC - kf*FC*T;
% dydt = [dydt_TC; dydt_FC];

T  = T_ini-y(1);
CP = 1-y(1)-y(2);
TC = y(1);
P  = P_ini+y(1)+y(2);
F  = F_ini-y(2);
FC = y(2);

dydt_TC =  kf*T*CP - kr*TC*P - kr*F*TC + kf*FC*T;
dydt_FC =  kr*F*TC - kf*FC*T;
dydt = [dydt_TC; dydt_FC];

end

