
% Help
% 
% Description  -> Calculate the yield of a single Toehold-exchange reaction
% 
%                T    +    CP      <->    TC    +     P      ; K = exp(-(dG+ddG)/(R*T) ) 
%              (t-x)      (1-y)    <->    (y)       (p+x)
% 
% Fn_Input     -> Input(1) = dG: <Array> Free energies of TER rxn
%              -> Input(2) = ddG: <Array> SNV free energy differences
% 
% Fn_Parameter -> Param(1) = t: <Array> Dimensionless [T]
%              -> Param(2) = p: <Array> Dimensionless [free P]
%              -> Param(3) = Tmp: <Scalar> Temperature in Celcius unit
% 
% Fn_Output    -> Output(1) = y: <Matrix> Yield of `TC`
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function y_sim = fnCal_Yield_TER(Input, Param)

T_ini = Param.t0;
P_ini = Param.p0;
Tmp = Param.Tmp;

dG_rxn = reshape(Input.dG, [1, length(Input.dG)] );
ddG_snv = reshape(Input.ddG, [length(Input.ddG), 1] );

dG_mat = ones(length(Input.ddG), 1)*dG_rxn + ddG_snv*ones(1, length(Input.dG));

R = 1.987*10^-3; % kcal*K^-1*mol^-1
T_kel = 273.15 + Tmp; % temperatue K

K_mat = exp(-dG_mat./(R*T_kel) );

y_sim = zeros(size(dG_mat) );

% wb = waitbar(0,'Please wait...');
syms x
for i = 1: numel(dG_mat)
    Eqn_CalYield = x*(P_ini+x)/((T_ini-x)*(1-x) ) == K_mat(i);
    
    sol = vpasolve(Eqn_CalYield, x);
    y_sim(i) = sol(sol > 0 & sol < min(T_ini, 1) );
    
%     waitbar(i/numel(dG_mat), wb, sprintf('Calculating SimY_TER: %d out of %d',i,numel(dG_mat)) )
end; clc
% close(wb)

end


