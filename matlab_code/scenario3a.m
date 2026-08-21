%% Scenario 04 - Disponibilidade limitada de HAPS e compensação por LEO
% Objetivo:
% Como a limitação do número de HAPS disponíveis afeta a seleção
% gulosa de transmissores
% Quantos transmissores das demais arquiteturas
% são necessários para atingir BCRB <= epsilon.
%
% Parte A: 
%   C/N0 = 50 dB-Hz, epsilon = 1 m
%   HAPS_max = 4, 2, 0
%   LEO = 8, MEO = 8, GEO = 6 disponíveis
%   Adaptado do scenario 03
%
% Parte B:
%   C/N0 = 50 dB-Hz, epsilon = 1 m
%   HAPS_max = 0
%   LEO_disponiveis = 8, 12, 16, 20
%   MEO = 8, GEO = 6 disponíveis
%   Compensação por LEO
%
% Latência calculada napenas o atraso de propagação de ida,
% aproximado por tau = d/c, 
% onde d é a distância geométrica entre transmissor e receptor.

close all;
clear;
clc;
format long;
rng(2);  % reprodutibilidade

%% ===================== Saída ==============================================

outDir = 'results_scenario_04';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ===================== Cenário do usuário =================================

lat = -22.8596582;      % [graus]
lon = -43.2303236;      % [graus]
h   = 10;               % [m]

UserPosition = llh2ecef(lat, lon, h);   % [x y z] ECEF [m]

%% ===================== Parâmetros globais =================================

c     = 3e8;        % [m/s]
beta  = 1.023e6;    % [Hz] largura de banda do código GPS C/A
T_coh = 0.02;       % [s] tempo de integração coerente

CN0_ref_dBHz = 50;  % C/N0 de referência [dB-Hz]
epsilon      = 1.0; % precisão-alvo [m]

scenarioNames = {'HAPS','LEO','MEO','GEO'};
numScenarios = numel(scenarioNames);

%% ===================== Tabela de candidatos base ==========================

Arquitetura = ["HAPS"; "LEO"; "MEO"; "GEO"];

% Número de dispositivos candidatos do caso base
NumDispositivosBase = [4; 8; 8; 6];

% Número máximo de LEO usado na análise de compensação qdo a meta
% não é atingida
N_LEO_max_compensacao = 20;
N_LEO_extra = N_LEO_max_compensacao - NumDispositivosBase(2);

if N_LEO_extra < 0
    error('N_LEO_max_compensacao deve ser >= ao número base de LEO.');
end

% Máscara mínima de elevação [graus]
MascaraElevacao_min_deg = [15; 5; 5; 5];

% Faixa de distância de propagação transmissor-receptor [km]
DistanciaPropagacao_min_km = [20; 800; 20000; 35786];
DistanciaPropagacao_max_km = [50; 1500; 27000; 41100];

% Faixa-alvo de PDOP para geração geométrica controlada
PDOP_min = [2.0; 2.0; 2.0; 4.0];
PDOP_max = [3.0; 3.0; 3.0; 6.0];

% Deslocamento relativo de C/N0 [dB]
Delta_CN0_rel_dB = [8; 4; 0; -4];

Tparam = table( ...
    Arquitetura, ...
    NumDispositivosBase, ...
    MascaraElevacao_min_deg, ...
    DistanciaPropagacao_min_km, ...
    DistanciaPropagacao_max_km, ...
    PDOP_min, ...
    PDOP_max, ...
    Delta_CN0_rel_dB);

disp('===== Tabela - Parâmetros dos candidatos base =====');
disp(Tparam);

writetable(Tparam, fullfile(outDir, ...
    'tabela_parametros_candidatos_base.csv'));

%% ===================== Parâmetros do Cenário 4 ============================

Tscenario4 = table( ...
    CN0_ref_dBHz, ...
    epsilon, ...
    sum(NumDispositivosBase), ...
    NumDispositivosBase(1), ...
    NumDispositivosBase(2), ...
    N_LEO_max_compensacao, ...
    NumDispositivosBase(3), ...
    NumDispositivosBase(4), ...
    'VariableNames', {'CN0_ref_dBHz','Precisao_alvo_m', ...
    'Total_candidatos_base','N_HAPS_cand','N_LEO_base', ...
    'N_LEO_max_compensacao','N_MEO_cand','N_GEO_cand'} ...
);

disp('===== Parâmetros do Cenário 4 =====');
disp(Tscenario4);

writetable(Tscenario4, fullfile(outDir, ...
    'tabela_parametros_cenario_04.csv'));

%% ===================== Geração do conjunto candidato base =================

SatByArch   = cell(numScenarios,1);
RangeByArch = cell(numScenarios,1);
PdopByArch  = zeros(numScenarios,1);

SatPool_base   = [];
archLabel_base = strings(0,1);
rangePool_base = [];

for sc = 1:numScenarios

    scenario = scenarioNames{sc};

    numSV       = NumDispositivosBase(sc);
    elevMaskDeg = MascaraElevacao_min_deg(sc);
    pdopTarget  = [PDOP_min(sc), PDOP_max(sc)];
    rangeLimits = [DistanciaPropagacao_min_km(sc), ...
                   DistanciaPropagacao_max_km(sc)] * 1e3; % [m]

    geom = generateGeometryPDOP(UserPosition, numSV, ...
        'ElevMaskDeg', elevMaskDeg, ...
        'PdopTarget', pdopTarget, ...
        'RangeLimits', rangeLimits, ...
        'MaxTries', 10000);

    SatByArch{sc}   = geom.SatPosition;
    RangeByArch{sc} = geom.range(:);
    PdopByArch(sc)  = geom.PDOP;

    SatPool_base   = [SatPool_base; geom.SatPosition];
    archLabel_base = [archLabel_base; repmat(string(scenario), numSV, 1)];
    rangePool_base = [rangePool_base; geom.range(:)];

    fprintf('%s -> PDOP gerado com %d candidatos: %.4f\n', ...
        scenario, numSV, geom.PDOP);
end

ID_base = (1:size(SatPool_base,1)).';
Tcand_base = table( ...
    ID_base, ...
    archLabel_base, ...
    rangePool_base/1e3, ...
    'VariableNames', {'ID','Arquitetura','Distancia_km'} ...
);

disp('===== Lista de candidatos base gerados =====');
disp(Tcand_base);

writetable(Tcand_base, fullfile(outDir, ...
    'tabela_candidatos_base_gerados.csv'));

Tpdop_arch = table( ...
    string(scenarioNames(:)), ...
    NumDispositivosBase, ...
    PdopByArch, ...
    'VariableNames', {'Arquitetura','Num_candidatos_base','PDOP_gerado'} ...
);

disp('===== PDOP dos conjuntos candidatos base por arquitetura =====');
disp(Tpdop_arch);

writetable(Tpdop_arch, fullfile(outDir, ...
    'tabela_pdop_conjuntos_candidatos_base.csv'));

%% ===================== Geração de LEO adicionais ==========================

SatLEO_extra   = [];
RangeLEO_extra = [];
PDOP_LEO_extra = NaN;

if N_LEO_extra > 0

    fprintf('\n===== Gerando %d candidatos LEO adicionais =====\n', N_LEO_extra);

    geomExtraLEO = generateGeometryPDOP(UserPosition, N_LEO_extra, ...
        'ElevMaskDeg', MascaraElevacao_min_deg(2), ...
        'PdopTarget', [PDOP_min(2), PDOP_max(2)], ...
        'RangeLimits', [DistanciaPropagacao_min_km(2), ...
                        DistanciaPropagacao_max_km(2)] * 1e3, ...
        'MaxTries', 10000);

    SatLEO_extra   = geomExtraLEO.SatPosition;
    RangeLEO_extra = geomExtraLEO.range(:);
    PDOP_LEO_extra = geomExtraLEO.PDOP;

    ID_extra = (1:N_LEO_extra).';
    Tcand_extra_leo = table( ...
        ID_extra, ...
        repmat("LEO_extra", N_LEO_extra, 1), ...
        RangeLEO_extra/1e3, ...
        'VariableNames', {'ID_extra','Arquitetura','Distancia_km'} ...
    );

    disp('===== Lista de candidatos LEO adicionais =====');
    disp(Tcand_extra_leo);

    fprintf('LEO adicional -> PDOP gerado com %d candidatos: %.4f\n', ...
        N_LEO_extra, PDOP_LEO_extra);

    writetable(Tcand_extra_leo, fullfile(outDir, ...
        'tabela_candidatos_LEO_adicionais.csv'));
end

%% ===================== Parte A: Limitação de HAPS ==========================

HAPS_max_vec_A = [4; 2; 0];
N_LEO_avail_A  = NumDispositivosBase(2);

NcasesA = numel(HAPS_max_vec_A);

N_HAPS_A   = zeros(NcasesA,1);
N_LEO_A    = zeros(NcasesA,1);
N_MEO_A    = zeros(NcasesA,1);
N_GEO_A    = zeros(NcasesA,1);
N_TOTAL_A  = zeros(NcasesA,1);

PDOP_final_A = zeros(NcasesA,1);
BCRB_final_A = zeros(NcasesA,1);
Reached_A    = false(NcasesA,1);

Dist_media_A_km = zeros(NcasesA,1);
Dist_max_A_km   = zeros(NcasesA,1);
Tau_media_A_ms  = zeros(NcasesA,1);
Tau_max_A_ms    = zeros(NcasesA,1);

IDs_A   = strings(NcasesA,1);
Archs_A = strings(NcasesA,1);

allHistories_A = cell(NcasesA,1);
allSelected_A  = cell(NcasesA,1);

fprintf('\n===== Parte A - Limitação do número de HAPS =====\n');

for kk = 1:NcasesA

    HAPS_max = HAPS_max_vec_A(kk);

    [SatPool_case, archLabel_case, rangePool_case] = buildCasePool( ...
        SatByArch, RangeByArch, SatLEO_extra, RangeLEO_extra, N_LEO_avail_A);

    sigmaPool_case = buildSigmaPoolFromLabels( ...
        CN0_ref_dBHz, archLabel_case, ...
        scenarioNames, Delta_CN0_rel_dB, c, beta, T_coh);

    out = runGreedySelectionHAPSLimit( ...
        UserPosition, ...
        SatPool_case, ...
        sigmaPool_case, ...
        archLabel_case, ...
        epsilon, ...
        HAPS_max);

    selectedRanges = rangePool_case(out.Selected); % [m]

    N_HAPS_A(kk)  = out.N_HAPS;
    N_LEO_A(kk)   = out.N_LEO;
    N_MEO_A(kk)   = out.N_MEO;
    N_GEO_A(kk)   = out.N_GEO;
    N_TOTAL_A(kk) = out.Total;

    PDOP_final_A(kk) = out.PDOP_final;
    BCRB_final_A(kk) = out.BCRB_final;
    Reached_A(kk)    = out.Reached;

    Dist_media_A_km(kk) = mean(selectedRanges)/1e3;
    Dist_max_A_km(kk)   = max(selectedRanges)/1e3;
    Tau_media_A_ms(kk)  = mean(selectedRanges ./ c)*1e3;
    Tau_max_A_ms(kk)    = max(selectedRanges ./ c)*1e3;

    IDs_A(kk)   = formatSelectedIDs(out.Selected);
    Archs_A(kk) = strjoin(archLabel_case(out.Selected).', "+");

    allHistories_A{kk} = out.History;
    allSelected_A{kk}  = out.Selected;

    fprintf('\nHAPS max = %d -> HAPS = %d, LEO = %d, MEO = %d, GEO = %d, Total = %d, BCRB = %.4f m, Meta = %d\n', ...
        HAPS_max, out.N_HAPS, out.N_LEO, out.N_MEO, out.N_GEO, ...
        out.Total, out.BCRB_final, out.Reached);
    fprintf('Distância média = %.2f km, Distância máxima = %.2f km\n', ...
        Dist_media_A_km(kk), Dist_max_A_km(kk));
    fprintf('Tau médio = %.4f ms, Tau máximo = %.4f ms\n', ...
        Tau_media_A_ms(kk), Tau_max_A_ms(kk));

    writetable(out.History, fullfile(outDir, ...
        sprintf('historico_parteA_HAPSmax_%d.csv', HAPS_max)));
end

T_HAPS_lim = table( ...
    HAPS_max_vec_A, ...
    N_HAPS_A, ...
    N_LEO_A, ...
    N_MEO_A, ...
    N_GEO_A, ...
    N_TOTAL_A, ...
    PDOP_final_A, ...
    BCRB_final_A, ...
    Reached_A, ...
    Dist_media_A_km, ...
    Dist_max_A_km, ...
    Tau_media_A_ms, ...
    Tau_max_A_ms, ...
    IDs_A, ...
    Archs_A, ...
    'VariableNames', {'HAPS_max','HAPS','LEO','MEO','GEO', ...
    'Total_dispositivos','PDOP_final','BCRB_final_m','Meta_atingida', ...
    'Dist_media_km','Dist_max_km','Tau_media_ms','Tau_max_ms', ...
    'IDs_selecionados','Arquiteturas_selecionadas'} ...
);

disp('===== Parte A - Sensibilidade à disponibilidade de HAPS =====');
disp(T_HAPS_lim);

writetable(T_HAPS_lim, fullfile(outDir, ...
    'tabela_cenario4_parteA_limitacao_HAPS.csv'));

%% ===================== Parte B: Compensação por LEO ========================

LEO_avail_vec_B = [8; 12; 16; 20];
HAPS_max_B      = 0;

NcasesB = numel(LEO_avail_vec_B);

N_HAPS_B   = zeros(NcasesB,1);
N_LEO_B    = zeros(NcasesB,1);
N_MEO_B    = zeros(NcasesB,1);
N_GEO_B    = zeros(NcasesB,1);
N_TOTAL_B  = zeros(NcasesB,1);

TotalCand_B = zeros(NcasesB,1);

PDOP_final_B = zeros(NcasesB,1);
BCRB_final_B = zeros(NcasesB,1);
Reached_B    = false(NcasesB,1);

Dist_media_B_km = zeros(NcasesB,1);
Dist_max_B_km   = zeros(NcasesB,1);
Tau_media_B_ms  = zeros(NcasesB,1);
Tau_max_B_ms    = zeros(NcasesB,1);

IDs_B   = strings(NcasesB,1);
Archs_B = strings(NcasesB,1);

allHistories_B = cell(NcasesB,1);
allSelected_B  = cell(NcasesB,1);

fprintf('\n===== Parte B - Compensação por aumento de LEO sem HAPS =====\n');

for kk = 1:NcasesB

    N_LEO_avail = LEO_avail_vec_B(kk);

    [SatPool_case, archLabel_case, rangePool_case] = buildCasePool( ...
        SatByArch, RangeByArch, SatLEO_extra, RangeLEO_extra, N_LEO_avail);

    sigmaPool_case = buildSigmaPoolFromLabels( ...
        CN0_ref_dBHz, archLabel_case, ...
        scenarioNames, Delta_CN0_rel_dB, c, beta, T_coh);

    out = runGreedySelectionHAPSLimit( ...
        UserPosition, ...
        SatPool_case, ...
        sigmaPool_case, ...
        archLabel_case, ...
        epsilon, ...
        HAPS_max_B);

    selectedRanges = rangePool_case(out.Selected); % [m]

    N_HAPS_B(kk)  = out.N_HAPS;
    N_LEO_B(kk)   = out.N_LEO;
    N_MEO_B(kk)   = out.N_MEO;
    N_GEO_B(kk)   = out.N_GEO;
    N_TOTAL_B(kk) = out.Total;

    TotalCand_B(kk) = size(SatPool_case,1);

    PDOP_final_B(kk) = out.PDOP_final;
    BCRB_final_B(kk) = out.BCRB_final;
    Reached_B(kk)    = out.Reached;

    Dist_media_B_km(kk) = mean(selectedRanges)/1e3;
    Dist_max_B_km(kk)   = max(selectedRanges)/1e3;
    Tau_media_B_ms(kk)  = mean(selectedRanges ./ c)*1e3;
    Tau_max_B_ms(kk)    = max(selectedRanges ./ c)*1e3;

    IDs_B(kk)   = formatSelectedIDs(out.Selected);
    Archs_B(kk) = strjoin(archLabel_case(out.Selected).', "+");

    allHistories_B{kk} = out.History;
    allSelected_B{kk}  = out.Selected;

    fprintf('\nHAPS max = 0, LEO disponíveis = %d -> HAPS = %d, LEO = %d, MEO = %d, GEO = %d, Total = %d, BCRB = %.4f m, Meta = %d\n', ...
        N_LEO_avail, out.N_HAPS, out.N_LEO, out.N_MEO, out.N_GEO, ...
        out.Total, out.BCRB_final, out.Reached);
    fprintf('Distância média = %.2f km, Distância máxima = %.2f km\n', ...
        Dist_media_B_km(kk), Dist_max_B_km(kk));
    fprintf('Tau médio = %.4f ms, Tau máximo = %.4f ms\n', ...
        Tau_media_B_ms(kk), Tau_max_B_ms(kk));

    writetable(out.History, fullfile(outDir, ...
        sprintf('historico_parteB_HAPSmax_0_LEOdisp_%d.csv', N_LEO_avail)));
end

T_LEO_comp = table( ...
    repmat(HAPS_max_B, NcasesB, 1), ...
    LEO_avail_vec_B, ...
    TotalCand_B, ...
    N_HAPS_B, ...
    N_LEO_B, ...
    N_MEO_B, ...
    N_GEO_B, ...
    N_TOTAL_B, ...
    PDOP_final_B, ...
    BCRB_final_B, ...
    Reached_B, ...
    Dist_media_B_km, ...
    Dist_max_B_km, ...
    Tau_media_B_ms, ...
    Tau_max_B_ms, ...
    IDs_B, ...
    Archs_B, ...
    'VariableNames', {'HAPS_max','LEO_disponiveis','Total_candidatos', ...
    'HAPS','LEO','MEO','GEO','Total_dispositivos','PDOP_final', ...
    'BCRB_final_m','Meta_atingida','Dist_media_km','Dist_max_km', ...
    'Tau_media_ms','Tau_max_ms','IDs_selecionados', ...
    'Arquiteturas_selecionadas'} ...
);

disp('===== Parte B - Compensação por LEO sem HAPS =====');
disp(T_LEO_comp);

writetable(T_LEO_comp, fullfile(outDir, ...
    'tabela_cenario4_parteB_compensacao_LEO.csv'));

%% ===================== Figuras - Parte A ==================================

figure;
hold on;
grid on;

bar(categorical(string(HAPS_max_vec_A)), ...
    [N_HAPS_A, N_LEO_A, N_MEO_A, N_GEO_A], ...
    'stacked');

xlabel('Número máximo de HAPS disponíveis', 'Interpreter', 'none');
ylabel('Número de transmissores selecionados', 'Interpreter', 'none');
title('Composição selecionada versus limite de HAPS', 'Interpreter', 'none');
legend('HAPS','LEO','MEO','GEO','Location','best');

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_8_composicao_limitacao_haps.png'), ...
    'Resolution', 300);

figure;
hold on;
grid on;

plot(HAPS_max_vec_A, BCRB_final_A, '-o', 'LineWidth', 1.5);
yline(epsilon, '--', sprintf('Meta %.1f m', epsilon), 'LineWidth', 1.2);

xlabel('Número máximo de HAPS disponíveis', 'Interpreter', 'none');
ylabel('BCRB final [m]', 'Interpreter', 'none');
title('BCRB final versus limite de HAPS', 'Interpreter', 'none');
legend('BCRB final','Precisão-alvo','Location','best');
xticks(sort(HAPS_max_vec_A));
xlim([min(HAPS_max_vec_A)-0.2, max(HAPS_max_vec_A)+0.2]);

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_9_bcrb_limitacao_haps.png'), ...
    'Resolution', 300);

figure;
hold on;
grid on;

plot(HAPS_max_vec_A, Tau_media_A_ms, '-o', 'LineWidth', 1.5);
plot(HAPS_max_vec_A, Tau_max_A_ms, '-s', 'LineWidth', 1.5);

xlabel('Número máximo de HAPS disponíveis', 'Interpreter', 'none');
ylabel('Atraso de propagação de ida [ms]', 'Interpreter', 'none');
title('Atraso de propagação versus limite de HAPS', 'Interpreter', 'none');
legend('Atraso médio','Atraso máximo','Location','best');
xticks(sort(HAPS_max_vec_A));
xlim([min(HAPS_max_vec_A)-0.2, max(HAPS_max_vec_A)+0.2]);

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_10_latencia_limitacao_haps.png'), ...
    'Resolution', 300);

%% ===================== Figuras - Parte B ==================================

figure;
hold on;
grid on;

bar(LEO_avail_vec_B, ...
    [N_HAPS_B, N_LEO_B, N_MEO_B, N_GEO_B], ...
    'stacked');

xlabel('Número de LEO disponíveis', 'Interpreter', 'none');
ylabel('Número de transmissores selecionados', 'Interpreter', 'none');
title('Composição selecionada sem HAPS versus disponibilidade de LEO', ...
    'Interpreter', 'none');
legend('HAPS','LEO','MEO','GEO','Location','best');
xticks(LEO_avail_vec_B);

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_11_composicao_compensacao_leo.png'), ...
    'Resolution', 300);

figure;
hold on;
grid on;

plot(LEO_avail_vec_B, BCRB_final_B, '-o', 'LineWidth', 1.5);
yline(epsilon, '--', sprintf('Meta %.1f m', epsilon), 'LineWidth', 1.2);

xlabel('Número de LEO disponíveis', 'Interpreter', 'none');
ylabel('BCRB final [m]', 'Interpreter', 'none');
title('BCRB final sem HAPS versus disponibilidade de LEO', ...
    'Interpreter', 'none');
legend('BCRB final','Precisão-alvo','Location','best');
xticks(LEO_avail_vec_B);

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_12_bcrb_compensacao_leo.png'), ...
    'Resolution', 300);

figure;
hold on;
grid on;

plot(LEO_avail_vec_B, Tau_media_B_ms, '-o', 'LineWidth', 1.5);
plot(LEO_avail_vec_B, Tau_max_B_ms, '-s', 'LineWidth', 1.5);

xlabel('Número de LEO disponíveis', 'Interpreter', 'none');
ylabel('Atraso de propagação de ida [ms]', 'Interpreter', 'none');
title('Atraso de propagação sem HAPS versus disponibilidade de LEO', ...
    'Interpreter', 'none');
legend('Atraso médio','Atraso máximo','Location','best');
xticks(LEO_avail_vec_B);

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_13_latencia_compensacao_leo.png'), ...
    'Resolution', 300);

%% ===================== Salvar variáveis ===================================

save(fullfile(outDir, 'scenario_04_results.mat'), ...
    'CN0_ref_dBHz', ...
    'epsilon', ...
    'scenarioNames', ...
    'NumDispositivosBase', ...
    'N_LEO_max_compensacao', ...
    'Delta_CN0_rel_dB', ...
    'Tparam', ...
    'Tscenario4', ...
    'Tcand_base', ...
    'Tpdop_arch', ...
    'SatByArch', ...
    'RangeByArch', ...
    'SatLEO_extra', ...
    'RangeLEO_extra', ...
    'T_HAPS_lim', ...
    'T_LEO_comp', ...
    'allHistories_A', ...
    'allSelected_A', ...
    'allHistories_B', ...
    'allSelected_B');

%% ===================== FUNÇÕES LOCAIS =====================================

function ecef = llh2ecef(lat_deg, lon_deg, h_m)
% LLH2ECEF  Converte coordenadas geodésicas WGS-84 para ECEF [x y z].

    a  = 6378137.0;
    f  = 1/298.257223563;
    e2 = f*(2-f);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    sinlat = sin(lat);  coslat = cos(lat);
    sinlon = sin(lon);  coslon = cos(lon);

    N = a / sqrt(1 - e2 * sinlat.^2);

    x = (N + h_m) * coslat * coslon;
    y = (N + h_m) * coslat * sinlon;
    z = (N * (1 - e2) + h_m) * sinlat;

    ecef = [x, y, z];
end

function [SatPool_case, archLabel_case, rangePool_case] = buildCasePool( ...
    SatByArch, RangeByArch, SatLEO_extra, RangeLEO_extra, N_LEO_avail)
% BUILDCASEPOOL --> Monta o conjunto candidato de um caso.
%
% HAPS, MEO e GEO são mantidos com a disponibilidade base.
% O número de LEO disponíveis pode variar, usando primeiro os LEO base
% e, em seguida, os LEO adicionais caso a meta não for atingida.

    SatHAPS = SatByArch{1};
    SatLEO_base = SatByArch{2};
    SatMEO = SatByArch{3};
    SatGEO = SatByArch{4};

    RangeHAPS = RangeByArch{1};
    RangeLEO_base = RangeByArch{2};
    RangeMEO = RangeByArch{3};
    RangeGEO = RangeByArch{4};

    SatLEO_all = [SatLEO_base; SatLEO_extra];
    RangeLEO_all = [RangeLEO_base; RangeLEO_extra];

    if N_LEO_avail > size(SatLEO_all,1)
        error('Número de LEO disponíveis maior que o número de candidatos LEO gerados.');
    end

    SatLEO = SatLEO_all(1:N_LEO_avail,:);
    RangeLEO = RangeLEO_all(1:N_LEO_avail);

    SatPool_case = [SatHAPS; SatLEO; SatMEO; SatGEO];

    archLabel_case = [
        repmat("HAPS", size(SatHAPS,1), 1);
        repmat("LEO",  size(SatLEO,1),  1);
        repmat("MEO",  size(SatMEO,1),  1);
        repmat("GEO",  size(SatGEO,1),  1)
    ];

    rangePool_case = [RangeHAPS; RangeLEO; RangeMEO; RangeGEO];
end

function sigmaPool = buildSigmaPoolFromLabels(CN0_ref_dBHz, archLabel, ...
    scenarioNames, Delta_CN0_rel_dB, c, beta, T_coh)
% BUILDSIGMAPOOL Fixa um sigma de pseudodistância por arquitetura.

    sigmaPool = zeros(numel(archLabel),1);

    for ii = 1:numel(archLabel)

        arch = char(archLabel(ii));
        sc = find(strcmp(scenarioNames, arch), 1);

        if isempty(sc)
            error('Arquitetura não reconhecida: %s', arch);
        end

        CN0_eff_dBHz = CN0_ref_dBHz + Delta_CN0_rel_dB(sc);
        sigmaPool(ii) = sigmaFromCN0(CN0_eff_dBHz, c, beta, T_coh);
    end
end

function sigma_rho = sigmaFromCN0(CN0_dBHz, c, beta, T_coh)

    CN0_linear = 10^(CN0_dBHz/10);
    sigma_rho = c / (2*pi*beta*sqrt(CN0_linear * T_coh));
end

function out = generateGeometryPDOP(UserPosition, M, varargin)
% GENERATEGEOMETRYPDOP Gera geometria sintética visível com PDOP controlado.

    p = inputParser;
    p.addParameter('ElevMaskDeg', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 90);
    p.addParameter('PdopTarget', [1.5 3.5], @(x) isnumeric(x) && numel(x) == 2 && x(1) <= x(2));
    p.addParameter('MaxTries', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('RangeLimits', [2.0e7 2.7e7], @(x) isnumeric(x) && numel(x) == 2 && x(1) > 0 && x(1) < x(2));
    p.addParameter('NumPRNPool', 32, @(x) isnumeric(x) && isscalar(x) && x >= M);
    p.parse(varargin{:});

    elevMaskDeg = p.Results.ElevMaskDeg;
    pdopTarget  = p.Results.PdopTarget;
    maxTries    = p.Results.MaxTries;
    rangeLimits = p.Results.RangeLimits;
    numPRNPool  = p.Results.NumPRNPool;

    UserPosition = UserPosition(:).';
    if numel(UserPosition) ~= 3
        error('UserPosition deve ter 3 elementos [x y z].');
    end
    if M < 4
        error('M deve ser >= 4 para cálculo de PDOP.');
    end

    x = UserPosition(1);
    y = UserPosition(2);
    z = UserPosition(3);

    lon = atan2(y, x);
    hyp = hypot(x, y);
    lat = atan2(z, hyp);

    e_hat = [-sin(lon),              cos(lon),             0          ];
    n_hat = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat)];
    u_hat = [ cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat)];
    Benu2ecef = [e_hat(:), n_hat(:), u_hat(:)];

    elMask = deg2rad(elevMaskDeg);

    best = struct('SatPosition', [], 'azel_deg', [], 'H', [], ...
                  'PDOP', [], 'u', [], 'range', []);
    bestDistanceToTarget = inf;

    for it = 1:maxTries
        az = 2*pi*rand(1, M);

        sinEl = sin(elMask) + (1 - sin(elMask)) * rand(1, M);
        el = asin(sinEl);

        ce = cos(el);
        sa = sin(az);
        ca = cos(az);
        v_enu = [ce .* sa; ce .* ca; sin(el)];

        v_ecef = Benu2ecef * v_enu;
        u_los = v_ecef.';

        R = rangeLimits(1) + (rangeLimits(2) - rangeLimits(1)) * rand(M,1);
        SatPosition_try = UserPosition + u_los .* R;

        [H_try, PDOP_try, u_try, range_try] = ...
            localGeometryAndPDOP(UserPosition, SatPosition_try);

        if isnan(PDOP_try)
            continue;
        end

        if PDOP_try >= pdopTarget(1) && PDOP_try <= pdopTarget(2)
            out = struct();
            out.SatPosition = SatPosition_try;
            out.SatPRN      = randperm(numPRNPool, M).';
            out.azel_deg    = [rad2deg(az(:)), rad2deg(el(:))];
            out.H           = H_try;
            out.PDOP        = PDOP_try;
            out.u           = u_try;
            out.range       = range_try;
            out.FoundTarget = true;
            return;
        end

        targetMid = mean(pdopTarget);
        distToTarget = abs(PDOP_try - targetMid);

        if distToTarget < bestDistanceToTarget
            bestDistanceToTarget = distToTarget;
            best.SatPosition = SatPosition_try;
            best.azel_deg    = [rad2deg(az(:)), rad2deg(el(:))];
            best.H           = H_try;
            best.PDOP        = PDOP_try;
            best.u           = u_try;
            best.range       = range_try;
        end
    end

    if isempty(best.PDOP)
        error('Não foi possível gerar uma geometria válida.');
    end

    out = struct();
    out.SatPosition = best.SatPosition;
    out.SatPRN      = randperm(numPRNPool, M).';
    out.azel_deg    = best.azel_deg;
    out.H           = best.H;
    out.PDOP        = best.PDOP;
    out.u           = best.u;
    out.range       = best.range;
    out.FoundTarget = false;
end

function [H, PDOP, u, d] = localGeometryAndPDOP(UserPosition, SatPosition)

    M = size(SatPosition,1);

    if M < 4
        H = [];
        PDOP = NaN;
        u = [];
        d = [];
        return;
    end

    r = SatPosition - UserPosition;
    d = sqrt(sum(r.^2, 2));

    if any(d <= eps)
        H = [];
        PDOP = NaN;
        u = [];
        d = [];
        return;
    end

    u = r ./ d;
    H = [-u, ones(M,1)];

    A = H.' * H;
    if rcond(A) < 1e-12
        H = [];
        PDOP = NaN;
        u = [];
        d = [];
        return;
    end

    C = A \ eye(4);
    PDOP = sqrt(trace(C(1:3,1:3)));
end

function [BCRB, PDOP, valid] = computeBCRB_fromSubset( ...
    UserPosition, SatPool, sigmaPool, idx)

    valid = false;
    BCRB = NaN;
    PDOP = NaN;

    idx = idx(:);
    M = numel(idx);

    if M < 4
        return;
    end

    SatPosition = SatPool(idx,:);
    sigma       = sigmaPool(idx);

    [H, PDOP_tmp] = localGeometryAndPDOP(UserPosition, SatPosition);

    if isempty(H) || isnan(PDOP_tmp)
        return;
    end

    R = diag(sigma.^2);
    J = H.' * (R \ H);

    if rcond(J) < 1e-12
        return;
    end

    Pcrb = J \ eye(4);

    BCRB = sqrt(trace(Pcrb(1:3,1:3)));
    PDOP = PDOP_tmp;
    valid = true;
end

function allowed = isHAPSAllowed(idx, archLabel, HAPS_max)

    selectedLabels = archLabel(idx);
    allowed = sum(selectedLabels == "HAPS") <= HAPS_max;
end

function out = runGreedySelectionHAPSLimit(UserPosition, SatPool, sigmaPool, ...
    archLabel, epsilon, HAPS_max)

    Mtotal = size(SatPool,1);
    minTx = 4;

    comb4 = nchoosek(1:Mtotal, minTx);

    bestBCRB = inf;
    bestComb = [];
    bestPDOP = NaN;

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

        if ~isHAPSAllowed(idx, archLabel, HAPS_max)
            continue;
        end

        [BCRB_k, PDOP_k, valid_k] = computeBCRB_fromSubset( ...
            UserPosition, ...
            SatPool, ...
            sigmaPool, ...
            idx);

        if ~valid_k
            continue;
        end

        if BCRB_k < bestBCRB
            bestBCRB = BCRB_k;
            bestComb = idx;
            bestPDOP = PDOP_k;
        end
    end

    if isempty(bestComb)
        error('Não foi possível encontrar subconjunto inicial válido com 4 dispositivos.');
    end

    selected   = bestComb(:).';
    candidates = setdiff(1:Mtotal, selected);

    currentBCRB = bestBCRB;
    currentPDOP = bestPDOP;

    histIter  = [];
    histNumTx = [];
    histAdded = strings(0,1);
    histArch  = strings(0,1);
    histBCRB  = [];
    histPDOP  = [];

    histIter(end+1,1)  = 0;
    histNumTx(end+1,1) = numel(selected);
    histAdded(end+1,1) = "Inicial";
    histArch(end+1,1)  = strjoin(archLabel(selected).', "+");
    histBCRB(end+1,1)  = currentBCRB;
    histPDOP(end+1,1)  = currentPDOP;

    iter = 1;
    stopReason = "Meta atingida";

    while currentBCRB > epsilon && ~isempty(candidates)

        bestAdd     = NaN;
        bestAddBCRB = inf;
        bestAddPDOP = NaN;

        for kk = 1:numel(candidates)

            cand = candidates(kk);
            trialSet = [selected, cand];

            if ~isHAPSAllowed(trialSet, archLabel, HAPS_max)
                continue;
            end

            [BCRB_trial, PDOP_trial, valid_trial] = ...
                computeBCRB_fromSubset( ...
                UserPosition, ...
                SatPool, ...
                sigmaPool, ...
                trialSet);

            if ~valid_trial
                continue;
            end

            if BCRB_trial < bestAddBCRB
                bestAdd     = cand;
                bestAddBCRB = BCRB_trial;
                bestAddPDOP = PDOP_trial;
            end
        end

        if isnan(bestAdd)
            stopReason = "Sem candidato adicional permitido ou válido";
            fprintf('Nenhum candidato adicional permitido/válido encontrado para HAPS_max = %d.\n', HAPS_max);
            break;
        end

        selected   = [selected, bestAdd];
        candidates = setdiff(candidates, bestAdd);

        currentBCRB = bestAddBCRB;
        currentPDOP = bestAddPDOP;

        histIter(end+1,1)  = iter;
        histNumTx(end+1,1) = numel(selected);
        histAdded(end+1,1) = "T" + string(bestAdd);
        histArch(end+1,1)  = archLabel(bestAdd);
        histBCRB(end+1,1)  = currentBCRB;
        histPDOP(end+1,1)  = currentPDOP;

        iter = iter + 1;
    end

    if currentBCRB > epsilon && isempty(candidates)
        stopReason = "Todos os candidatos disponíveis foram selecionados";
    end

    selectedLabels = archLabel(selected);

    out.Selected = selected;
    out.N_HAPS = sum(selectedLabels == "HAPS");
    out.N_LEO  = sum(selectedLabels == "LEO");
    out.N_MEO  = sum(selectedLabels == "MEO");
    out.N_GEO  = sum(selectedLabels == "GEO");
    out.Total  = numel(selected);
    out.BCRB_final = currentBCRB;
    out.PDOP_final = currentPDOP;
    out.Reached = currentBCRB <= epsilon;
    out.StopReason = stopReason;

    out.History = table( ...
        histIter, ...
        histNumTx, ...
        histAdded, ...
        histArch, ...
        histBCRB, ...
        histPDOP, ...
        'VariableNames', {'Iteracao','Num_dispositivos','Adicionado', ...
        'Arquitetura','BCRB_m','PDOP'} ...
    );
end

function s = formatSelectedIDs(selected)

    tmp = "T" + string(selected(:).');
    s = strjoin(tmp, ",");
end
