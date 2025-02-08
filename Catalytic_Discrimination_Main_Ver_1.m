
% Date    -> 2024-01-02
% Project -> Two-dimensional Characterization of SNV discrimination
% Version -> V20231228
% 
% Description -> Chemical Reaction Model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Rxn-1 ->        T       +      CP      <>     TC      +      P           ;    kf, kr
%             T_ini-y(1)     1-y(1)-y(2)        y(1)     P_ini+y(1)+y(2)
% 
% Rxn-2 ->        F       +      TC      <>     FC      +      T           ;    kr, kf
%             F_ini-y(2)         y(1)           y(2)       T_ini-y(1)
% ==================================================================================================


%% Part - I. Thermodynamics Only

% -|> Descroption:
%   a. Given contstant Conc_P, changing Conc_F
%   b. To predict the performance of discriminating SNVs, changing ddG_SNV
% -|> Outputs:
%   1. Regular TER, SimY_TER and DF_TER
%   2. NCCR, SimY_NCCR and DF_NCCR

% ==================================================================================================
clear; clc
Temperature = 37; % Celcius

T_ini = 15.6/20; % dimensionlessed [target]
P_ini = 1-1; % dimensionlessed [free P]
F_ini = linspace(1,20,20)'; % dimensionlessed [F]
leng_F = length(F_ini);

range_dG = linspace(-5, 5, 101)';
range_ddG = linspace(0, 5, 51)'; % Assume the ddG of SNV
leng_dG = length(range_dG);
leng_ddG = length(range_ddG);

% Settings used for Calculating Fig (7 ~ )
arr_dG = linspace(-5, 0, 6)';
arr_T0 = 2.^linspace(-1.0342, 8.9658, 101)'/20;
arr_F0 = 19; 

% .......................................................................................
Input1_NCCR = struct('dG', range_dG, 'ddG', range_ddG);
Param_TER = struct('t0', T_ini, 'p0', P_ini, 'Tmp', Temperature );

SimY_TER = fnCal_Yield_TER(Input1_NCCR, Param_TER);

% .......................................................................................
% Calculate the Yields of NCCR, @ different [F]_0
SimY_NCCR_TC = cell(size(F_ini) );
SimY_NCCR_FC = cell(size(F_ini) );
SimY_NCCR = cell(size(F_ini) );

parfor j = 1: length(F_ini)
    Param1_NCCR = struct('t0', T_ini, 'p0', P_ini, 'f0', F_ini(j), 'Tmp', Temperature );
    [temp_SimY_TC, temp_SimY_FC, temp_SimYfam] = fnCal_Yield_NCCR(Input1_NCCR, Param1_NCCR);

    SimY_NCCR_TC{j} = temp_SimY_TC;
    SimY_NCCR_FC{j} = temp_SimY_FC;
    SimY_NCCR{j} = temp_SimYfam;

    fprintf('Progress > %d th of 20,  \n', j );
end; clc
% fprintf('Progress of Calculation - Part I is Done !  \n');
% save("Part_I_Thermodynamics.mat");


% ==================================================================================================
% Visualization
Save_1_SimY_TER = SimY_TER; % !!! why `SimY_TER*2' ???
Save_1_SimY_TER(Save_1_SimY_TER < 0.01) = 0.01;
Save_2_SimDF_TER = Save_1_SimY_TER(1,:)./Save_1_SimY_TER;

[axis_X1, axis_Y1] = meshgrid(range_dG, range_ddG);

% Figure (1). <TER> : Yield ~ (range_dG, range_ddG)
figure
Fig1_m_TER = pcolor(axis_X1, axis_Y1, Save_1_SimY_TER);
Fig1_c_TER = colorbar;
Fig1_c_TER.Label.String = 'Sim Yield';
set(Fig1_m_TER, 'LineStyle', 'none');
clim([0 1]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');

% Figure (2). <TER> : DF ~ (range_dG, range_ddG)
figure
Fig2_m_DF = pcolor(axis_X1, axis_Y1, Save_2_SimDF_TER);
Fig2_c_DF = colorbar;
Fig2_c_DF.Label.String = 'Sim DF';
set(Fig2_m_DF, 'LineStyle', 'none');
clim([1 20]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');


% --------------------------------------------------------------------------------------------------
Save_3_SimYfam_NCCR = zeros(leng_ddG, leng_dG, leng_F);
Save_3_SimYrox_NCCR = zeros(leng_ddG, leng_dG, leng_F);
Save_4_SimDF_NCCR = zeros(leng_ddG, leng_dG, leng_F);

for l = 1: leng_F
    % temp_SimYfam = SimY_NCCR_TC{l}*2;
    % temp_SimYrox = SimY_NCCR{l}*2;
    temp_SimYfam = SimY_NCCR_TC{l};
    temp_SimYrox = SimY_NCCR{l};

    temp_SimYfam(temp_SimYfam < 0.01) = 0.01;
    temp_SimYrox(temp_SimYrox < 0.01) = 0.01;
    Save_3_SimYfam_NCCR(:, :, l) = temp_SimYfam;
    Save_3_SimYrox_NCCR(:, :, l) = temp_SimYrox;

    temp_SimDF = temp_SimYfam(1, :)./temp_SimYfam;
    Save_4_SimDF_NCCR(:, :, l) = temp_SimDF;
end

% .......................................................................................
id_F0 = 10;

% Figure (3). <NCCR> : Yield-ROX ~ (range_dG, range_ddG) @ F = 10
figure
Fig3_m_NCCR = pcolor(axis_X1, axis_Y1, Save_3_SimYrox_NCCR(:,:,id_F0));
Fig3_c_NCCR = colorbar;
Fig3_c_NCCR.Label.String = 'Sim Yield';
set(Fig3_m_NCCR, 'LineStyle', 'none');
clim([0 1]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');

% Figure (4). <NCCR> : DF-FAM ~ (range_dG, range_ddG) @ F = 10
figure
Fig4_m_NCCR = pcolor(axis_X1, axis_Y1, Save_4_SimDF_NCCR(:,:,id_F0));
Fig4_c_NCCR = colorbar;
Fig4_c_NCCR.Label.String = 'Sim DF';
set(Fig4_m_NCCR, 'LineStyle', 'none');
clim([1 20]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');
clc

% .......................................................................................
% Figure (5). <NCCR> : Yield Improvement ~ (range_dG, range_ddG) @ F = 10
figure
Fig5_m_NCCR = pcolor(axis_X1, axis_Y1, ...
    Save_3_SimYrox_NCCR(:,:,id_F0)./Save_1_SimY_TER );
Fig5_c_NCCR = colorbar;
Fig5_c_NCCR.Label.String = 'Yield Multiplicity';
set(Fig5_m_NCCR, 'LineStyle', 'none');
clim([1 20]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');
clc

% Figure (6). <NCCR> : DF Improvement ~ (range_dG, range_ddG) @ F = 10
figure
Fig6_m_NCCR = pcolor(axis_X1, axis_Y1, ...
    Save_4_SimDF_NCCR(:,:,id_F0)./Save_2_SimDF_TER );
Fig6_c_NCCR = colorbar;
Fig6_c_NCCR.Label.String = 'DF Multiplicity';
set(Fig6_m_NCCR, 'LineStyle', 'none');
clim([1 20]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');
clc


% --------------------------------------------------------------------------------------------------
% .......................................................................................
% Calculate the Yields of NCCR, @ [F_0] = 19
SimY2_NCCR_TC = cell(size(arr_dG) );
SimY2_NCCR_FC = cell(size(arr_dG) );
SimY2_NCCR = cell(size(arr_dG) );

parfor l = 1: length(arr_dG)

    Input2_NCCR = struct('t0', arr_T0, 'ddG', range_ddG);
    Param2_NCCR = struct('dG', arr_dG(l), 'p0', P_ini, 'f0', arr_F0, 'Tmp', Temperature );

    [temp_SimY_TC, temp_SimY_FC, temp_SimYfam] = fnCal_Yield2_NCCR(Input2_NCCR, Param2_NCCR);

    SimY2_NCCR_TC{l} = temp_SimY_TC;
    SimY2_NCCR_FC{l} = temp_SimY_FC;
    SimY2_NCCR{l} = temp_SimYfam;

    fprintf('Progress > %d th of %i,  \n', l, length(arr_dG) );
end; clc

fprintf('Progress of Calculation - Part I is Done !  \n');

% .......................................................................................
[axis_X2, axis_Y2] = meshgrid(arr_T0, range_ddG );

Save_5_SimYfam_NCCR = zeros(length(range_ddG), length(arr_T0), length(arr_dG));
Save_5_SimYrox_NCCR = zeros(length(range_ddG), length(arr_T0), length(arr_dG));
Save_6_SimDF_NCCR = zeros(length(range_ddG), length(arr_T0), length(arr_dG));

for l = 1: length(arr_dG)

    temp_SimYfam = SimY2_NCCR_TC{l};
    temp_SimYrox = SimY2_NCCR{l};

    temp_SimYfam(temp_SimYfam < 0.01) = 0.01;
    temp_SimYrox(temp_SimYrox < 0.01) = 0.01;
    Save_5_SimYfam_NCCR(:, :, l) = temp_SimYfam;
    Save_5_SimYrox_NCCR(:, :, l) = temp_SimYrox;

    temp_SimDF = temp_SimYfam(1, :)./temp_SimYfam;
    Save_6_SimDF_NCCR(:, :, l) = temp_SimDF;
end

% .......................................................................................
id_dG_TER = 1;

% Figure (7). <NCCR> : Yield-ROX ~ ([T0], range_ddG) @ constant dG_TER
figure
Fig7_m_NCCR = pcolor(axis_X2*20, axis_Y2, Save_5_SimYrox_NCCR(:,:,id_dG_TER));
Fig7_c_NCCR = colorbar;
Fig7_c_NCCR.Label.String = 'Sim Yield';
set(Fig7_m_NCCR, 'LineStyle', 'none');
clim([0 1]);
colormap('jet');
xlabel('[T_0] / nM')
ylabel('ddG / kcal/mol' )
% set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');

% Figure (8). <NCCR> : Yield-FAM ~ ([T0], range_ddG) @ constant dG_TER
figure
Fig8_m_NCCR = pcolor(axis_X2*20, axis_Y2, Save_5_SimYfam_NCCR(:,:,id_dG_TER));
Fig8_c_NCCR = colorbar;
Fig8_c_NCCR.Label.String = 'Sim Yield';
set(Fig8_m_NCCR, 'LineStyle', 'none');
clim([0 1]);
colormap('jet');
xlabel('[T_0] / nM')
ylabel('ddG / kcal/mol' )
% set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');

% Figure (9). <NCCR> : DF-FAM ~ ([T0], range_ddG) @ constant dG_TER
figure
Fig9_m_NCCR = pcolor(axis_X2*20, axis_Y2, Save_6_SimDF_NCCR(:,:,id_dG_TER));
Fig9_c_NCCR = colorbar;
Fig9_c_NCCR.Label.String = 'Sim DF';
set(Fig9_m_NCCR, 'LineStyle', 'none');
clim([1 20]);
colormap('jet');
xlabel('[T_0] / nM')
ylabel('ddG / kcal/mol' )
% set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');
clc


% --------------------------------------------------------------------------------------------------
% Temporat
Save_TempData_NCCR1_FAM = Save_3_SimYfam_NCCR(:,:, 20 );
Save_TempData_NCCR1_ROX = Save_3_SimYrox_NCCR(:,:, 20 );

Save_TempData_NCCR2_FAM = Save_5_SimYfam_NCCR(:,:, 6 );
Save_TempData_NCCR2_ROX = Save_5_SimYrox_NCCR(:,:, 6 );


%% Part - II. Kinetics and Thermodynamics
% ==================================================================================================
% -|> Descroption:
%   a. Given contstant Conc_P, changing Conc_F
%   b. To predict the performance of discriminating SNVs, changing ddG_SNV
% -|> Outputs:
%   1. Regular TER, SimY_TER and DF_TER
%   2. NCCR, SimY_NCCR and DF_NCCR


% Input ---------------------------------------------------------------------------------

Conc_ref = 20; % nM
Exp_T0 = 15.6/20; % initial concentration of <Target> strand
freeP0 = 1-1; % initial concentration of <Free Protect> strand
TC0 = 0;
FC0 = 0;
Fuel0 = [1, 19, 1, 19]; % initial concentration of <Fuel> strand

FileName = 'DataForSim_SX26_v20231228.xlsx';
FAM_sheet = 1; 
ROX_sheet = 2;

rangetime = 'D7: D67';
rangeY = 'E7: J67';

Exp_time = readmatrix(FileName, 'Sheet', FAM_sheet, 'Range', rangetime);
ExpY_FAM = readmatrix(FileName, 'Sheet', FAM_sheet, 'Range', rangeY);
ExpY_ROX = readmatrix(FileName, 'Sheet', ROX_sheet, 'Range', rangeY);


%% S1. Fit the Rate Constants -------------------------------------------------------------------

%  S1.1 Calculate the equilibrium constant K_rxn /////////////////////////////

fn_Y2K = @(yield) (yield ).*(yield+freeP0 )./((Exp_T0-yield ).*(1-yield ) );

K_Exp_FAM = fn_Y2K(mean(ExpY_FAM(end-4:end, :) ) ); % Equilibrium constant
K_Exp_ROX = fn_Y2K(mean(ExpY_ROX(end-4:end, :) ) ); % Equilibrium constant


%  S1.2 Fit the kf and kr of TER, FAM Channel /////////////////////////////
%       Using data <Without Fuel> strand to fit out the <kf> and <kr> values of standard 
%       Toehold-Exchange reaction
%       Parameters need to be fitted ->   kf_wt, kr_wt, kf_snv, kr_snv

Data_TER_FAM = [ExpY_FAM(:, 1), ExpY_FAM(:, 4)]; % Col 1: F = 1; Col 2: F = 19
K_TER_FAM = [K_Exp_FAM(1), K_Exp_FAM(4)];

FitY_TER_FAM = zeros(size(Data_TER_FAM) );
Fitk_TER_FAM = zeros(2, length(K_TER_FAM) ); % row 1: the kf; row 2: the kr
FitK_TER_FAM = zeros(size(K_TER_FAM) );

kf_assume_FAM = [0.0064    2.2331]; % Col 1: F = 1; Col 2: F = 19
kr_assume_FAM = [0.1884    0.0309];

% input for the Fitting Model -> fnFit_TER
input_TER_FAM = struct('T', Exp_T0, 'P', freeP0, 'TC', 0, 't', Exp_time);

for i = 1: length(K_TER_FAM)
    k_assume = [kf_assume_FAM(i); kr_assume_FAM(i)];
    
    [pb_TER_FAM, rn_TER_FAM, res_TER_FAM, ef_TER_FAM, out_TER_FAM] = ...
        lsqcurvefit(@fnFit_TER, k_assume, input_TER_FAM, Data_TER_FAM(:, i) );
    
    Fitk_TER_FAM(:, i) = pb_TER_FAM;
    FitK_TER_FAM(i) = pb_TER_FAM(1)/pb_TER_FAM(2);
    FitY_TER_FAM(:, i) = fnFit_TER(pb_TER_FAM, input_TER_FAM);
    
    figure(1)
    subplot(2, 2, i)
    plot(Exp_time, Data_TER_FAM(:, i), 'o', Exp_time, FitY_TER_FAM(:, i), '-', 'Linewidth', 1.5)
    ylim([0 1] )
    title('Kinetic Fitting of Exp-TER', 'FAM Channel' )
end

%  S1.3 Fit the kf and kr of TER, ROX Channel /////////////////////////////

Data_TER_ROX = [ExpY_ROX(:, 1), ExpY_ROX(:, 4)];
K_TER_ROX = [K_Exp_ROX(1), K_Exp_ROX(4)];

FitY_TER_ROX = zeros(size(Data_TER_ROX) );
Fitk_TER_ROX = zeros(2, length(K_TER_ROX) ); % row 1: the fowward reaction rate; row 2: reverse reaction
FitK_TER_ROX = zeros(size(K_TER_ROX) );

kf_assume_ROX = [0.0055    0.7000];
kr_assume_ROX = [0.1699    0.0675];

% input for the Fitting Model -> fnFit_TER
inTER_ROX = struct('T', Exp_T0, 'P', freeP0, 'TC', 0, 't', Exp_time);

for i = 1: length(K_TER_ROX)
    k_assume = [kf_assume_ROX(i); kr_assume_ROX(i)];
    
    [pb_TER_ROX, rn_TER_ROX, res_TER_ROX, ef_TER_ROX, out_TER_ROX] = ...
        lsqcurvefit(@fnFit_TER, k_assume, inTER_ROX, Data_TER_ROX(:, i) );
    
    Fitk_TER_ROX(:, i) = pb_TER_ROX;
    FitK_TER_ROX(i) = pb_TER_ROX(1)/pb_TER_ROX(2);
    FitY_TER_ROX(:, i) = fnFit_TER(pb_TER_ROX, inTER_ROX);
    
    figure(1)
    subplot(2, 2, 2+i)
    plot(Exp_time, Data_TER_ROX(:, i), 'o', Exp_time, FitY_TER_ROX(:, i), '-', 'Linewidth', 1.5)
    ylim([0 1] )
    title('Kinetic Fitting of Exp-TER', 'ROX Channel' )
end

Header_TER = ["Experimental K"; "Simulational K"; "Fitted k_f"; "Fitted k_r"];
FitResults_TER_FAM = [K_TER_FAM; FitK_TER_FAM; Fitk_TER_FAM];
array2table([Header_TER, FitResults_TER_FAM], 'VariableNames', {'FAM Channel','WT','SNV'} )

FitResults_TER_ROX = [K_TER_ROX; FitK_TER_ROX; Fitk_TER_ROX];
array2table([Header_TER, FitResults_TER_ROX], 'VariableNames', {'Rox Channel','WT','SNV'} )

% exportgraphics(gcf, ...
%     'Fig 1. Fitting the kinetic parameters from TER data.jpg', ...
%     'Resolution', 600)
% close(gcf)


%% S2. Sim the Kinetic Curves of NCCR -----------------------------------------------------------

%  S2.1 Fit the Kinetics of NCCR, FAM Channel /////////////////////////////

Data_NCCR_FAM = [ExpY_FAM(:, 2:3), ExpY_FAM(:, 5:6)]; % Col 1: F = 1; Col 2: F = 19
K_NCCR_FAM = [K_Exp_FAM(2:3), K_Exp_FAM(5:6)];

kf_WT_FAM = Fitk_TER_FAM(1, 1);
kr_WT_FAM = Fitk_TER_FAM(2, 1);
kf_SNV_FAM = Fitk_TER_FAM(1, 2);
kr_SNV_FAM = Fitk_TER_FAM(2, 2);

kf_NCCR_FAM = [kf_WT_FAM, kf_WT_FAM, kf_SNV_FAM, kf_SNV_FAM];
kr_NCCR_FAM = [kr_WT_FAM, kr_WT_FAM, kr_SNV_FAM, kr_SNV_FAM];

FitY_TC_FAM = zeros(size(Data_NCCR_FAM) );
FitY_FC_FAM = zeros(size(Data_NCCR_FAM) );
FitY_NCCR_FAM = zeros(size(Data_NCCR_FAM) );

plot_title = ["WT", "WT", "SNV", "SNV"];
for i = 1: length(K_NCCR_FAM)

    input_NCCR_FAM = struct('T', Exp_T0, 'P', freeP0, 'F', Fuel0(i), ...
        'kf', kf_NCCR_FAM(i), 'kr', kr_NCCR_FAM(i), 't', Exp_time, 'y0', [0, 0]);
    
    [t, y] = fnSim_NCCR(input_NCCR_FAM);

    FitY_TC_FAM(:, i) = y(:, 1);
    FitY_FC_FAM(:, i) = y(:, 2);
    FitY_NCCR_FAM(:, i) = y(:, 1);
    
    figure(2)
    subplot(2, 4, i)
    plot(Exp_time, Data_NCCR_FAM(:, i), 'o', Exp_time, FitY_NCCR_FAM(:, i), '-', 'Linewidth', 1)
    ylim([0 1] )
    title('FAM Channel. ' + plot_title(i), ['F = ', num2str(Fuel0(i)*Conc_ref), ' nM'])
    
end


%  S2.2 Fit the Kinetics of NCCR, ROX Channel /////////////////////////////

Data_NCCR_ROX = [ExpY_ROX(:, 2:3), ExpY_ROX(:, 5:6)];
K_NCCR_ROX = [K_Exp_ROX(2:3), K_Exp_ROX(5:6)];

kf_WT_ROX = Fitk_TER_ROX(1, 1);
kr_WT_ROX = Fitk_TER_ROX(2, 1);
kf_SNV_ROX = Fitk_TER_ROX(1, 2);
kr_SNV_ROX = Fitk_TER_ROX(2, 2);

kf_NCCR_ROX = [kf_WT_ROX, kf_WT_ROX, kf_SNV_ROX, kf_SNV_ROX];
kr_NCCR_ROX = [kr_WT_ROX, kr_WT_ROX, kr_SNV_ROX, kr_SNV_ROX];

FitY_TC_ROX = zeros(size(Data_NCCR_ROX) );
FitY_FC_ROX = zeros(size(Data_NCCR_ROX) );
FitY_NCCR_ROX = zeros(size(Data_NCCR_ROX) );

for i = 1: length(K_NCCR_ROX)

    input_NCCR_ROX = struct('T', Exp_T0, 'P', freeP0, 'F', Fuel0(i), ...
        'kf', kf_NCCR_ROX(i), 'kr', kr_NCCR_ROX(i), 't', Exp_time, 'y0', [0, 0]);
    
    [t, y] = fnSim_NCCR(input_NCCR_ROX);
    
    FitY_TC_ROX(:, i) = y(:, 1);
    FitY_FC_ROX(:, i) = y(:, 2);
    FitY_NCCR_ROX(:, i) = y(:, 1) + y(:, 2);
    
    figure(2)
    subplot(2, 4, 4+i)
    plot(Exp_time, Data_NCCR_ROX(:, i), 'o', Exp_time, FitY_NCCR_ROX(:, i), '-', 'Linewidth', 1)
    ylim([0 1] )
    title('ROX Channel. ' + plot_title(i), ['F = ', num2str(Fuel0(i)*Conc_ref), ' nM'])
    
end

% exportgraphics(gcf, ...
%     'Fig 2. Fitting the kinetic parameters from Fuelled-TER data.jpg', ...
%     'Resolution', 600)
% close(gcf)


%% S3. Sim the Kinetic Curves of NCCR -----------------------------------------------------------
tic

Sim_Time = [linspace(0, 120, 121)'; 720; 1440];
TimeSlice = [1 2 6 11 31 61 121 122 123]'; % min

leng_T0 = 197;
leng_F0 = 197;

SimStd_T0 = 2.^linspace(-1.0342, 8.9658, leng_T0)'/20; % when = 99, T0 = 15.6 nM
SimStd_F0 = linspace(20, 1000, leng_F0)'/20; % (1) = 1; (73) = 19

SimStd_T0_nM = 20*SimStd_T0;
SimStd_F0_nM = 20*SimStd_F0;

SimStdY_TER_WT_FAM = zeros(length(Sim_Time), leng_T0 );
SimStdY_TER_SNV_FAM = zeros(length(Sim_Time), leng_T0 );


%  S3.1 Sim the Kinetics of TER, Only FAM /////////////////////////////
for i = 1: leng_T0
% parfor i = 1: leng_T0

    % WT Sim data
    input_TER_WT = struct('T', SimStd_T0(i), 'P', freeP0, 'F', 0, ...
        'kf', kf_WT_FAM, 'kr', kr_WT_FAM, 't', Sim_Time, 'y0', [0, 0]);

    [~, y_WT] = fnSim_NCCR(input_TER_WT);

    SimStdY_TER_WT_FAM(:, i) = y_WT(:, 1);
    
    % SNV Sim data
    input_TER_SNV = struct('T', SimStd_T0(i), 'P', freeP0, 'F', 0, ...
        'kf', kf_SNV_FAM, 'kr', kr_SNV_FAM, 't', Sim_Time, 'y0', [0, 0]);

    [~, y_SNV] = fnSim_NCCR(input_TER_SNV);

    SimStdY_TER_SNV_FAM(:, i) = y_SNV(:, 1);
    
end
fprintf('^^^^^^^^^^^^^^^^^^^ Progress-S3.1 : Done ^^^^^^^^^^^^^^^^^^^^^\n')


% ----------------------------------------------------------------------------------
%  S3.2.1 Sim the Kinetics of NCCR, FAM Channel /////////////////////////////
[tt0, ff0] = ndgrid(SimStd_T0, SimStd_F0);
allpair = [tt0(:), ff0(:)];

temp_T0 = allpair(:, 1);
temp_F0 = allpair(:, 2);

tempY_NCCR_WT_FAM = cell(1, size(allpair,1) );
tempY_NCCR_SNV_FAM = cell(1, size(allpair,1) );

for ii = 1: size(allpair, 1)
% parfor ii = 1: size(allpair, 1)
    
    % WT Sim data
    input_NCCR_FAM_WT = struct('T', temp_T0(ii), 'P', freeP0, 'F', temp_F0(ii), ...
        'kf', kf_WT_FAM, 'kr', kr_WT_FAM, 't', Sim_Time, 'y0', [0, 0]);

    [~, y_FAM_WT] = fnSim_NCCR(input_NCCR_FAM_WT);

    tempY_NCCR_WT_FAM{ii} = y_FAM_WT(:, 1);

    % SNV Sim data
    input_NCCR_FAM_SNV = struct('T', temp_T0(ii), 'P', freeP0, 'F', temp_F0(ii), ...
        'kf', kf_SNV_FAM, 'kr', kr_SNV_FAM, 't', Sim_Time, 'y0', [0, 0]);

    [~, y_FAM_SNV] = fnSim_NCCR(input_NCCR_FAM_SNV);

    tempY_NCCR_SNV_FAM{ii} = y_FAM_SNV(:, 1);

end
fprintf('^^^^^^^^^^^^^^^^^^^ Progress-S3.2.1 : Done ^^^^^^^^^^^^^^^^^^^\n')

SimStdY_NCCR_WT_FAM = reshape([tempY_NCCR_WT_FAM{:}], [length(Sim_Time), leng_T0, leng_F0] );
SimStdY_NCCR_SNV_FAM = reshape([tempY_NCCR_SNV_FAM{:}], [length(Sim_Time), leng_T0, leng_F0] );


% ............................................................................................
%  S3.2.2 Sim the Kinetics of NCCR, ROX Channel /////////////////////////////
tempY_NCCR_WT_ROX = cell(1, size(allpair,1) );
tempY_NCCR_SNV_ROX = cell(1, size(allpair,1) );

for ii = 1: size(allpair, 1)
% parfor ii = 1: size(allpair, 1)

    % WT Sim data
    inNCCR_ROX_WT = struct('T', temp_T0(ii), 'P', freeP0, 'F', temp_F0(ii), ...
        'kf', kf_WT_ROX, 'kr', kr_WT_ROX, 't', Sim_Time, 'y0', [0, 0]);

    [~, y_ROX_WT] = fnSim_NCCR(inNCCR_ROX_WT);

    tempY_NCCR_WT_ROX{ii} = y_ROX_WT(:, 1) + y_ROX_WT(:, 2);

    % SNV Sim data
    inNCCR_ROX_SNV = struct('T', temp_T0(ii), 'P', freeP0, 'F', temp_F0(ii), ...
        'kf', kf_SNV_ROX, 'kr', kr_SNV_ROX, 't', Sim_Time, 'y0', [0, 0]);

    [~, y_ROX_SNV] = fnSim_NCCR(inNCCR_ROX_SNV);

    tempY_NCCR_SNV_ROX{ii} = y_ROX_SNV(:, 1) + y_ROX_SNV(:, 2);

end
fprintf('^^^^^^^^^^^^^^^^^^^ Progress-S3.2.2 : Done ^^^^^^^^^^^^^^^^^^^\n')

SimStdY_NCCR_WT_ROX = reshape([tempY_NCCR_WT_ROX{:}], [length(Sim_Time), leng_T0, leng_F0] );
SimStdY_NCCR_SNV_ROX = reshape([tempY_NCCR_SNV_ROX{:}], [length(Sim_Time), leng_T0, leng_F0] );

clc; toc


%% S4. Visualization in Heatmap format ----------------------------------------------------------

% the Dimension of Matrices is [Time - T0 - F0]
figY_TER_WT_FAM = SimStdY_TER_WT_FAM;
figY_TER_SNV_FAM = SimStdY_TER_SNV_FAM;
figY_NCCR_WT_FAM = SimStdY_NCCR_WT_FAM;
figY_NCCR_SNV_FAM = SimStdY_NCCR_SNV_FAM;
figY_NCCR_WT_ROX = SimStdY_NCCR_WT_ROX;
figY_NCCR_SNV_ROX = SimStdY_NCCR_SNV_ROX;

figY_TER_WT_FAM(figY_TER_WT_FAM < 0.01) = 0.01;
figY_TER_SNV_FAM(figY_TER_SNV_FAM < 0.01) = 0.01;
figY_NCCR_WT_FAM(figY_NCCR_WT_FAM < 0.01) = 0.01;
figY_NCCR_SNV_FAM(figY_NCCR_SNV_FAM < 0.01) = 0.01;
figY_NCCR_WT_ROX(figY_NCCR_WT_ROX < 0.01) = 0.01;
figY_NCCR_SNV_ROX(figY_NCCR_SNV_ROX < 0.01) = 0.01;

Conc_T0 = Conc_ref*SimStd_T0;
Conc_F0 = Conc_ref*SimStd_F0;

[X1, Y1] = meshgrid(Conc_T0, Conc_F0 );


% ............................................................................................
figure(3)
% X axis -> [Target]/nM, Y axis -> Yield of TER
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    hold on
    PL_WT = plot(Conc_T0, figY_TER_WT_FAM(TimeSlice(i), :), '-', 'Linewidth', 3);
    PL_WT.Color = [0 0.4470 0.7410 (0.9*i/length(TimeSlice)+0.1)];

    PL_SNV = plot(Conc_T0, figY_TER_SNV_FAM(TimeSlice(i), :), '-', 'Linewidth', 3);
    PL_SNV.Color = [0.8500 0.3250 0.0980 (0.9*i/length(TimeSlice)+0.1)];

    fill_x = [Conc_T0; flipud(Conc_T0)];
    inBetween = [figY_TER_WT_FAM(TimeSlice(i), :), fliplr(figY_TER_SNV_FAM(TimeSlice(i), :))];
    fill(fill_x, inBetween, [224,255,255]/255, 'LineStyle', 'none');
    hold off

    set(gca, 'XScale', 'log')
    ylim([0 1.1] )
    legend('WT', 'SNV', 'Location', 'NorthWest')
    title('Std-TER', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end


% ............................................................................................
figure(4)
% Heatmap, Snapshots, FAM of WT for NCCR ~ ([Target], [Fuel])
% X axis -> [Target]/nM, Y axis -> [Fuel]/nM, Color -> Yield
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    
    reSimY_NCCR_WT_FAM = reshape( figY_NCCR_WT_FAM(TimeSlice(i), :, :), ...
        [leng_T0 leng_F0] )'; % IMPORTANT !!! the reshaping !!!
    
    m_NCCR_WT_FAM = pcolor(X1, Y1, reSimY_NCCR_WT_FAM);
    c_NCCR_WT_FAM = colorbar;
    c_NCCR_WT_FAM.Label.String = 'Yield-FAM';
    set(m_NCCR_WT_FAM, 'LineStyle', 'none');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colormap('jet')
    clim([0, 1])
    xlabel('[Target] (nM)')
    ylabel('[Fuel] (nM)')
    title('Std-NCCR, WT, FAM', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end

figure(5)
% Heatmap, Snapshots, FAM of SNV for NCCR ~ ([Target], [Fuel])
% X axis -> [Target]/nM, Y axis -> [Fuel]/nM, Color -> Yield
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    
    reSimY_NCCR_SNV_FAM = reshape( figY_NCCR_SNV_FAM(TimeSlice(i), :, :), ...
        [leng_T0 leng_F0] )';
    
    m_NCCR_SNV_FAM = pcolor(X1, Y1, reSimY_NCCR_SNV_FAM);
    c_NCCR_SNV_FAM = colorbar;
    c_NCCR_SNV_FAM.Label.String = 'Yield-FAM';
    set(m_NCCR_SNV_FAM, 'LineStyle', 'none');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colormap('jet')
    clim([0, 1])
    xlabel('[Target] (nM)')
    ylabel('[Fuel] (nM)')
    title('Std-NCCR, SNV, FAM', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end

figure(6)
% Heatmap, Snapshots, ROX of WT for NCCR ~ ([Target], [Fuel])
% X axis -> [Target]/nM, Y axis -> [Fuel]/nM, Color -> Yield
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    
    reSimY_NCCR_WT_ROX = reshape( figY_NCCR_WT_ROX(TimeSlice(i), :, :), ...
        [leng_T0 leng_F0] )';
    
    m_NCCR_WT_ROX = pcolor(X1, Y1, reSimY_NCCR_WT_ROX);
    c_NCCR_WT_ROX = colorbar;
    c_NCCR_WT_ROX.Label.String = 'Yield-FAM';
    set(m_NCCR_WT_ROX, 'LineStyle', 'none');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colormap('jet')
    clim([0, 1])
    xlabel('[Target] (nM)')
    ylabel('[Fuel] (nM)')
    title('Std-NCCR, WT, ROX', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end

figure(7)
% Heatmap, Snapshots, ROX of SNV for NCCR ~ ([Target], [Fuel])
% X axis -> [Target]/nM, Y axis -> [Fuel]/nM, Color -> Yield
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    
    reSimY_NCCR_SNV_ROX = reshape( figY_NCCR_SNV_ROX(TimeSlice(i), :, :), ...
        [leng_T0 leng_F0] )';
    
    m_NCCR_SNV_ROX = pcolor(X1, Y1, reSimY_NCCR_SNV_ROX);
    c_NCCR_SNV_ROX = colorbar;
    c_NCCR_SNV_ROX.Label.String = 'Yield-FAM';
    set(m_NCCR_SNV_ROX, 'LineStyle', 'none');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colormap('jet')
    clim([0, 1])
    xlabel('[Target] (nM)')
    ylabel('[Fuel] (nM)')
    title('Std-NCCR, SNV, ROX', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end


%% S5. Analytical Evaluation of Simulation Results ----------------------------------------------

% step 1: permute the Yield matrices of NCCR by [F0 - T0 - Time]
figYpm_NCCR_WT_FAM = permute(figY_NCCR_WT_FAM, [3 2 1] );
figYpm_NCCR_WT_ROX = permute(figY_NCCR_WT_ROX, [3 2 1] );

figYpm_NCCR_SNV_FAM = permute(figY_NCCR_SNV_FAM, [3 2 1] );
figYpm_NCCR_SNV_ROX = permute(figY_NCCR_SNV_ROX, [3 2 1] );

% step 2: calculating the Distance between Origin and Data points
Norm_NCCR_WT = (figYpm_NCCR_WT_FAM.^2 + figYpm_NCCR_WT_ROX.^2).^0.5;
Norm_NCCR_SNV = (figYpm_NCCR_SNV_FAM.^2 + figYpm_NCCR_SNV_ROX.^2).^0.5;

% step 3: calculating the Distance between SNV and WT points
dVec_NCCR = ( (figYpm_NCCR_SNV_FAM - figYpm_NCCR_WT_FAM).^2 + ...
            (figYpm_NCCR_SNV_ROX - figYpm_NCCR_WT_ROX).^2 ).^0.5;

% ............................................................................................
%  Visualization the Distances /////////////////////////////

figure(8)
% Heatmap, Snapshots, Original Distance of WT for NCCR ~ ([Target], [Fuel])
% X axis -> [Target]/nM, Y axis -> [Fuel]/nM, Color -> Yield
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    
    m_oL_NCCR_WT = pcolor(X1, Y1, Norm_NCCR_WT(:,:,TimeSlice(i) ) );
    c_oL_NCCR_WT = colorbar;
    c_oL_NCCR_WT.Label.String = 'Norm Length';
    set(m_oL_NCCR_WT, 'LineStyle', 'none');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colormap('jet')
    clim([0, 1.5])
    xlabel('[Target] (nM)')
    ylabel('[Fuel] (nM)')
    title('Std-NCCR, WT, Origin', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end


figure(9)
% Heatmap, Snapshots, ROX of SNV for NCCR ~ ([Target], [Fuel])
% X axis -> [Target]/nM, Y axis -> [Fuel]/nM, Color -> Yield
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    
    m_oL_NCCR_SNV = pcolor(X1, Y1, Norm_NCCR_SNV(:,:,TimeSlice(i) ) );
    c_oL_NCCR_SNV = colorbar;
    c_oL_NCCR_SNV.Label.String = 'Norm Length';
    set(m_oL_NCCR_SNV, 'LineStyle', 'none');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colormap('jet')
    clim([0, 1.5])
    xlabel('[Target] (nM)')
    ylabel('[Fuel] (nM)')
    title('Std-NCCR, SNV, Origin', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end


figure(10)
% Heatmap, Snapshots, ROX of SNV for NCCR ~ ([Target], [Fuel])
% X axis -> [Target]/nM, Y axis -> [Fuel]/nM, Color -> Yield
for i = 1: length(TimeSlice)
    subplot(3, 3, i)
    
    m_dL_NCCR = pcolor(X1, Y1, dVec_NCCR(:,:,TimeSlice(i) ) );
    c_dL_NCCR = colorbar;
    c_dL_NCCR.Label.String = 'Vector Distance';
    set(m_dL_NCCR, 'LineStyle', 'none');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colormap('jet')
    clim([0, 1.5])
    xlabel('[Target] (nM)')
    ylabel('[Fuel] (nM)')
    title('Std-NCCR, Distance', ['Time = ', num2str(Sim_Time(TimeSlice(i) ) ), ' min'])
end


%% S6. Creating 3D Gif plots --------------------------------------------------------------------

gifAxis_Y = struct('X', Conc_T0, 'Y', Conc_F0, 'scale', 'log', 'limit', [0 1], ...
    'frame', Sim_Time, 'face', 'flat');
gifAxis_dL = struct('X', Conc_T0, 'Y', Conc_F0, 'scale', 'log', 'limit', [0 1.5], ...
    'frame', Sim_Time, 'face', 'flat');


fnPlot_Gif(gifAxis_Y, figYpm_NCCR_WT_FAM, 'Yield_NCCR_Spurious_FAM');
close(gcf)

fnPlot_Gif(gifAxis_Y, figYpm_NCCR_WT_ROX, 'Yield_NCCR_Spurious_ROX');
close(gcf)

fnPlot_Gif(gifAxis_Y, figYpm_NCCR_SNV_FAM, 'Yield_NCCR_Correct_FAM');
close(gcf)

fnPlot_Gif(gifAxis_Y, figYpm_NCCR_SNV_ROX, 'Yield_NCCR_Correct_ROX');
close(gcf)

% ............................................................................................
fnPlot_Gif(gifAxis_dL, Norm_NCCR_WT, 'VecNorm_NCCR_Spurious');
close(gcf)

fnPlot_Gif(gifAxis_dL, Norm_NCCR_SNV, 'VecNorm_NCCR_Correct');
close(gcf)

fnPlot_Gif(gifAxis_dL, dVec_NCCR, 'VecDistance_NCCR');
close(gcf)


%% Part - III. Detection Window

R = 1.987*10^-3; % kcal*K^-1*mol^-1
T_kel = 273.15 + Temperature; % temperatue K

SimStd_P0 = linspace(0, 25, 26)';

dG_TER_SX26 = -4.01;
ddG_SNV_XS26 = 2.91;
K_TER_SX26_c = exp(-dG_TER_SX26./(R*T_kel) );
K_TER_SX26_s = exp(-(dG_TER_SX26+ddG_SNV_XS26)./(R*T_kel) );

SimY_TER_SX26_c = zeros(length(SimStd_T0), length(SimStd_P0) );
SimY_TER_SX26_s = zeros(length(SimStd_T0), length(SimStd_P0) );

syms x
for j = 1: length(SimStd_P0)
    for i = 1: length(SimStd_T0)

        Eqn_SX26_c = x*(SimStd_P0(j)+x)/((SimStd_T0(i)-x)*(1-x) ) == K_TER_SX26_c;
        sol_c = vpasolve(Eqn_SX26_c, x);
        SimY_TER_SX26_c(i,j) = sol_c(sol_c > 0 & sol_c < min(SimStd_T0(i), 1) );

        Eqn_SX26_s = x*(SimStd_P0(j)+x)/((SimStd_T0(i)-x)*(1-x) ) == K_TER_SX26_s;
        sol_s = vpasolve(Eqn_SX26_s, x);
        SimY_TER_SX26_s(i,j) = sol_s(sol_s > 0 & sol_s < min(SimStd_T0(i), 1) );

    end
    sprintf('Calculating SimY_TER: %d out of %d', j, length(SimStd_P0) )
end; clc

% .......................................................................................

figure(11)
for j = 1: length(SimStd_P0)

    subplot(5, 6, j)
    hold on
    plot(Conc_T0, SimY_TER_SX26_c(:,j), '-', 'Linewidth', 1)
    plot(Conc_T0, SimY_TER_SX26_s(:,j), '-', 'Linewidth', 1)
    set(gca, 'XScale', 'log')
    ylim([0 1] )
    title(['[P]_0 = ', num2str(SimStd_P0(j)*Conc_ref), ' nM'] )
    hold off

end


%% .......................................................................................

select_F0 = 37; %  1: 20nM, 17: 100nM, 37: 200nM,

figure(12) % https://ww2.mathworks.cn/help/matlab/ref/view.html
for i = 1: length(TimeSlice)
    subplot(3,3,i)

    hold on
    x1 = Conc_T0;
    y1 = figYpm_NCCR_SNV_FAM(select_F0, :, TimeSlice(i) );
    z1 = figYpm_NCCR_SNV_ROX(select_F0, :, TimeSlice(i) );

    x2 = Conc_T0;
    y2 = figYpm_NCCR_WT_FAM(select_F0, :, TimeSlice(i) );
    z2 = figYpm_NCCR_WT_ROX(select_F0, :, TimeSlice(i) );

    c = [0.3010 0.7450 0.9330];

    PL3_SNV = plot3(x1, y1, z1, '-', 'Linewidth', 1.5);
    PL3_WT = plot3(x2, y2, z2, '-', 'Linewidth', 1.5);
    p = fill3([x1;flip(x2)]', [y1,flip(y2)], [z1,flip(z2)], c, ...
        'FaceAlpha',.25,'LineStyle','none');

    hold off
    grid on
    box on
    view(-35, 20 )
    set(gca, 'XScale', 'log')
    set(gca, 'BoxStyle', 'full')
    ylim([0 1] )
    zlim([0 1] )
    yticks([0 0.5 1])
    zticks([0 0.5 1])
    title(['[Fuel] = ', num2str(Conc_F0(select_F0) ), ' nM; ', ...
        'Time = ', num2str(Sim_Time(TimeSlice(i) ) ) ])
end


%% .......................................................................................

% Correct 3D line
x1 = Conc_T0;
y1_2h = figYpm_NCCR_SNV_FAM(select_F0, :, TimeSlice(7) );
z1_2h = figYpm_NCCR_SNV_ROX(select_F0, :, TimeSlice(7) );
% Spurious 3D line
x2 = Conc_T0;
y2_2h = figYpm_NCCR_WT_FAM(select_F0, :, TimeSlice(7) );
z2_2h = figYpm_NCCR_WT_ROX(select_F0, :, TimeSlice(7) );

c = [0.3010 0.7450 0.9330];
c_ROX = [0.8500 0.3250 0.0980];
c_FAM = [0.4660 0.6740 0.1880];


figure(13)
hold on
PL3_SNV_2h = plot3(x1, y1_2h, z1_2h, '-', 'Linewidth', 1.5);
PL3_WT_2h = plot3(x2, y2_2h, z2_2h, '-', 'Linewidth', 1.5);

p = fill3([x1;flip(x2)]', [y1_2h,flip(y2_2h)], [z1_2h,flip(z2_2h)], c, ...
    'FaceAlpha',.25,'LineStyle','none');
p_ROX = fill3([x1;flip(x2)]', ones(size([y1_2h,flip(y2_2h)])), [z1_2h,flip(z2_2h)], c_ROX, ...
    'FaceAlpha',.25,'LineStyle','none');
p_FAM = fill3([x1;flip(x2)]', [y1_2h,flip(y2_2h)], zeros(size([z1_2h,flip(z2_2h)])), c_FAM, ...
    'FaceAlpha',.25,'LineStyle','none');

hold off
grid off
box on
view(-35, 20 )
set(gca, 'XScale', 'log')
set(gca, 'BoxStyle', 'full')
ylim([0 1] )
zlim([0 1] )
yticks([0 0.5 1])
zticks([0 0.5 1])
% title(['[Fuel] = ', num2str(Conc_F0(select_F0) ), ' nM; ', ...
%     'Time = ', num2str(Sim_Time(TimeSlice(7) ) ) ])




%% Figure Export Setup

% temp_fig = gcf;
% set(gcf, 'position' ,[1, 1, 1028, 775])
% set(gca, "FontSize", 11)
% 
% title('');
% xlabel('[T] / nM');
% ylabel('[F] / nM');
% xticks([10^0 10^1 10^2]);
% yticks([10^1 10^2 10^3]);
% 
% temp_c = colorbar;
% temp_c.Label.String = '';













