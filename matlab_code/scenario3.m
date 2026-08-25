%% Scenario 03 - Seleção gulosa variando C/N0
% Objetivo:
% Avaliar como a configuração selecionada muda em função do C/N0 de referência.
%
% Critério:
% Para cada valor de C/N0, selecionar dispositivos HAPS, LEO, MEO e GEO
% até atender BCRB <= epsilon, quando possível.

close all;
clear;
clc;
format long;
rng(2);  % reprodutibilidade do Scenario 03

%% ##################### Saída ##############################################

outDir = 'results_scenario_03';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ##################### Cenário do usuário #################################

lat = -22.8596582;      % [graus]
lon = -43.2303236;      % [graus]
h   = 10;               % [m]

UserPosition = llh2ecef(lat, lon, h);   % [x y z] ECEF [m]

%% ##################### Parâmetros globais #################################

c     = 3e8;        % [m/s]
beta  = 1.023e6;    % [Hz] largura de banda do código GPS C/A
T_coh = 0.02;       % [s] tempo de integração coerente

CN0_vec = 30:1:50;  % C/N0 de referência [dB-Hz]
epsilon = 1.0;      % precisão-alvo [m]

scenarioNames = {'HAPS','LEO','MEO','GEO'};
numScenarios = numel(scenarioNames);

%% ##################### Tabela de candidatos disponíveis ##################=

Arquitetura = ["HAPS"; "LEO"; "MEO"; "GEO"];

% Número de dispositivos candidatos por arquitetura
% Mantido igual ao Scenario 02 com disponibilidade limitada de HAPS.
NumDispositivos = [4; 8; 8; 6];

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
    NumDispositivos, ...
    MascaraElevacao_min_deg, ...
    DistanciaPropagacao_min_km, ...
    DistanciaPropagacao_max_km, ...
    PDOP_min, ...
    PDOP_max, ...
    Delta_CN0_rel_dB);

disp('===== Tabela - Parâmetros dos candidatos disponíveis =====');
disp(Tparam);

writetable(Tparam, fullfile(outDir, ...
    'tabela_parametros_candidatos_disponiveis.csv'));

%% ##################### Parâmetros do problema ###########################

Tsel_param = table( ...
    min(CN0_vec), ...
    max(CN0_vec), ...
    epsilon, ...
    sum(NumDispositivos), ...
    NumDispositivos(1), ...
    NumDispositivos(2), ...
    NumDispositivos(3), ...
    NumDispositivos(4), ...
    'VariableNames', {'CN0_min_dBHz','CN0_max_dBHz','Precisao_alvo_m', ...
    'Total_candidatos','N_HAPS','N_LEO','N_MEO','N_GEO'} ...
);

disp('===== Parâmetros da varredura em C/N0 =====');
disp(Tsel_param);

writetable(Tsel_param, fullfile(outDir, ...
    'tabela_parametros_varredura_cn0.csv'));

%% ##################### Geração do conjunto candidato #####################

SatPool   = [];
archLabel = strings(0,1);
rangePool = [];
pdopArch  = zeros(numScenarios,1);

for sc = 1:numScenarios

    scenario = scenarioNames{sc};

    numSV       = NumDispositivos(sc);
    elevMaskDeg = MascaraElevacao_min_deg(sc);
    pdopTarget  = [PDOP_min(sc), PDOP_max(sc)];
    rangeLimits = [DistanciaPropagacao_min_km(sc), ...
                   DistanciaPropagacao_max_km(sc)] * 1e3; % [m]

    geom = generateGeometryPDOP(UserPosition, numSV, ...
        'ElevMaskDeg', elevMaskDeg, ...
        'PdopTarget', pdopTarget, ...
        'RangeLimits', rangeLimits, ...
        'MaxTries', 10000);

    SatPosition_sc = geom.SatPosition;
    pdopArch(sc)   = geom.PDOP;

    SatPool   = [SatPool; SatPosition_sc];
    archLabel = [archLabel; repmat(string(scenario), numSV, 1)];
    rangePool = [rangePool; geom.range(:)];

    fprintf('%s -> PDOP gerado com %d candidatos: %.4f\n', ...
        scenario, numSV, geom.PDOP);
end

Mtotal = size(SatPool,1);

ID = (1:Mtotal).';
Tcand = table( ...
    ID, ...
    archLabel, ...
    rangePool/1e3, ...
    'VariableNames', {'ID','Arquitetura','Distancia_km'} ...
);

disp('===== Lista de candidatos gerados =====');
disp(Tcand);

writetable(Tcand, fullfile(outDir, ...
    'tabela_candidatos_gerados.csv'));

Tpdop_arch = table( ...
    string(scenarioNames(:)), ...
    NumDispositivos, ...
    pdopArch, ...
    'VariableNames', {'Arquitetura','Num_candidatos','PDOP_gerado'} ...
);

disp('===== PDOP dos conjuntos candidatos por arquitetura =====');
disp(Tpdop_arch);

writetable(Tpdop_arch, fullfile(outDir, ...
    'tabela_pdop_conjuntos_candidatos.csv'));

%% ##################### Varredura em C/N0 #############################

Npoints = numel(CN0_vec);

N_HAPS_vec   = zeros(Npoints,1);
N_LEO_vec    = zeros(Npoints,1);
N_MEO_vec    = zeros(Npoints,1);
N_GEO_vec    = zeros(Npoints,1);
N_TOTAL_vec  = zeros(Npoints,1);

BCRB_final_vec = zeros(Npoints,1);
PDOP_final_vec = zeros(Npoints,1);
Reached_vec    = false(Npoints,1);

allHistories = cell(Npoints,1);
allSelected  = cell(Npoints,1);

for kk = 1:Npoints

    CN0_ref_dBHz = CN0_vec(kk);

    sigmaPool = buildSigmaPool( ...
        CN0_ref_dBHz, ...
        NumDispositivos, ...
        Delta_CN0_rel_dB, ...
        c, beta, T_coh);

    out = runGreedySelection( ...
        UserPosition, ...
        SatPool, ...
        sigmaPool, ...
        archLabel, ...
        epsilon);

    N_HAPS_vec(kk)  = out.N_HAPS;
    N_LEO_vec(kk)   = out.N_LEO;
    N_MEO_vec(kk)   = out.N_MEO;
    N_GEO_vec(kk)   = out.N_GEO;
    N_TOTAL_vec(kk) = out.Total;

    BCRB_final_vec(kk) = out.BCRB_final;
    PDOP_final_vec(kk) = out.PDOP_final;
    Reached_vec(kk)    = out.Reached;

    allHistories{kk} = out.History;
    allSelected{kk}  = out.Selected;

    fprintf('\nC/N0 = %.1f dB-Hz -> Total = %d, BCRB = %.4f m, Meta = %d\n', ...
        CN0_ref_dBHz, out.Total, out.BCRB_final, out.Reached);
end

%% ##################### Tabela 7.6 - Resultado da varredura #############

T_CN0 = table( ...
    CN0_vec(:), ...
    N_HAPS_vec, ...
    N_LEO_vec, ...
    N_MEO_vec, ...
    N_GEO_vec, ...
    N_TOTAL_vec, ...
    PDOP_final_vec, ...
    BCRB_final_vec, ...
    Reached_vec, ...
    'VariableNames', {'CN0_ref_dBHz','HAPS','LEO','MEO','GEO', ...
    'Total_dispositivos','PDOP_final','BCRB_final_m','Meta_atingida'} ...
);

disp('===== Tabela 7.6 - Seleção gulosa versus C/N0 =====');
disp(T_CN0);

writetable(T_CN0, fullfile(outDir, ...
    'tabela_selecao_gulosa_vs_cn0.csv'));

T_CN0_viavel = T_CN0(Reached_vec,:);
disp('===== Configurações viáveis =====');
disp(T_CN0_viavel);

writetable(T_CN0_viavel, fullfile(outDir, ...
    'tabela_configuracoes_viaveis_vs_cn0.csv'));

%% ##################### Figura 7.5 - Total versus C/N0 #####################

N_TOTAL_feasible = N_TOTAL_vec;
N_TOTAL_feasible(~Reached_vec) = NaN;

N_TOTAL_infeasible = N_TOTAL_vec;
N_TOTAL_infeasible(Reached_vec) = NaN;

figure;
hold on;
grid on;

plot(CN0_vec, N_TOTAL_feasible, '-o', 'LineWidth', 1.5);
plot(CN0_vec, N_TOTAL_infeasible, 'x', ...
    'LineWidth', 1.5, 'MarkerSize', 8);

xlabel('C/N0 de referência [dB-Hz]', 'Interpreter', 'none');
ylabel('Número de dispositivos selecionados', 'Interpreter', 'none');
title('Número de dispositivos selecionados versus C/N0', ...
    'Interpreter', 'none');

legend('Meta atingida', 'Meta não atingida', 'Location', 'best');

ylim([0 Mtotal+2]);
yticks(0:2:Mtotal+2);

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_5_total_dispositivos_vs_cn0.png'), ...
    'Resolution', 300);

%% ##################### Figura 7.6 - Composição viável versus C/N0 ##########

if any(Reached_vec)

    figure;
    hold on;
    grid on;

    bar(CN0_vec(Reached_vec), ...
        [N_HAPS_vec(Reached_vec), ...
         N_LEO_vec(Reached_vec), ...
         N_MEO_vec(Reached_vec), ...
         N_GEO_vec(Reached_vec)], ...
         'stacked');

    xlabel('C/N0 de referência [dB-Hz]', 'Interpreter', 'none');
    ylabel('Número de dispositivos selecionados', 'Interpreter', 'none');
    title('Composição das configurações viáveis selecionadas', ...
        'Interpreter', 'none');

    legend('HAPS','LEO','MEO','GEO','Location','best');

    ylim([0 max(N_TOTAL_vec(Reached_vec))+2]);
    yticks(0:1:max(N_TOTAL_vec(Reached_vec))+2);

    exportgraphics(gcf, fullfile(outDir, ...
        'fig_7_6_composicao_viavel_vs_cn0.png'), ...
        'Resolution', 300);
else
    warning('Nenhuma configuração atingiu a meta. Figura de composição viável não gerada.');
end

%% ##################### Figura 7.7 - BCRB final versus C/N0 ################

BCRB_feasible = BCRB_final_vec;
BCRB_feasible(~Reached_vec) = NaN;

BCRB_infeasible = BCRB_final_vec;
BCRB_infeasible(Reached_vec) = NaN;

figure;
hold on;
grid on;

plot(CN0_vec, BCRB_feasible, '-o', 'LineWidth', 1.5);
plot(CN0_vec, BCRB_infeasible, 'x', ...
    'LineWidth', 1.5, 'MarkerSize', 8);

yline(epsilon, '--', sprintf('Meta %.1f m', epsilon), ...
    'LineWidth', 1.2);

xlabel('C/N0 de referência [dB-Hz]', 'Interpreter', 'none');
ylabel('BCRB final [m]', 'Interpreter', 'none');
title('BCRB final da configuração selecionada versus C/N0', ...
    'Interpreter', 'none');

legend('Meta atingida', 'Meta não atingida', 'Precisão-alvo', ...
    'Location', 'best');

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_7_bcrb_final_vs_cn0.png'), ...
    'Resolution', 300);

%% ##################### Figura 7.8 - Primeiro C/N0 viável ##################=

idxFirstFeasible = find(Reached_vec, 1, 'first');

if ~isempty(idxFirstFeasible)

    CN0_first = CN0_vec(idxFirstFeasible);

    fprintf('\n===== Primeiro C/N0 viável =====\n');
    fprintf('Primeiro C/N0 que atinge a meta: %.1f dB-Hz\n', CN0_first);
    fprintf('Total de dispositivos: %d\n', N_TOTAL_vec(idxFirstFeasible));
    fprintf('BCRB final: %.6f m\n', BCRB_final_vec(idxFirstFeasible));

    T_first = T_CN0(idxFirstFeasible,:);

    writetable(T_first, fullfile(outDir, ...
        'tabela_primeiro_cn0_viavel.csv'));
else
    fprintf('\nNenhum valor de C/N0 atingiu a meta.\n');
end

%% ##################### Salvar variáveis ###################################

save(fullfile(outDir, 'scenario_03_results.mat'), ...
    'CN0_vec', ...
    'epsilon', ...
    'scenarioNames', ...
    'SatPool', ...
    'archLabel', ...
    'rangePool', ...
    'Tparam', ...
    'Tsel_param', ...
    'Tcand', ...
    'Tpdop_arch', ...
    'T_CN0', ...
    'T_CN0_viavel', ...
    'allHistories', ...
    'allSelected', ...
    'N_HAPS_vec', ...
    'N_LEO_vec', ...
    'N_MEO_vec', ...
    'N_GEO_vec', ...
    'N_TOTAL_vec', ...
    'BCRB_final_vec', ...
    'PDOP_final_vec', ...
    'Reached_vec');

%% ##################### FUNÇÕES LOCAIS #####################################

function ecef = llh2ecef(lat_deg, lon_deg, h_m)
%  Converte coordenadas geodésicas WGS-84 para ECEF [x y z].

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

function sigmaPool = buildSigmaPool(CN0_ref_dBHz, NumDispositivos, ...
    Delta_CN0_rel_dB, c, beta, T_coh)

    sigmaPool = [];

    for sc = 1:numel(NumDispositivos)

        CN0_eff_dBHz = CN0_ref_dBHz + Delta_CN0_rel_dB(sc);

        sigma_rho = sigmaFromCN0(CN0_eff_dBHz, c, beta, T_coh);

        sigmaPool = [sigmaPool; sigma_rho * ones(NumDispositivos(sc),1)];
    end
end

function sigma_rho = sigmaFromCN0(CN0_dBHz, c, beta, T_coh)

    CN0_linear = 10^(CN0_dBHz/10);

    sigma_rho = c / (2*pi*beta*sqrt(CN0_linear * T_coh));
end

function out = generateGeometryPDOP(UserPosition, M, varargin)
%  Gera geometria sintética visível com PDOP controlado.

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

function out = runGreedySelection(UserPosition, SatPool, sigmaPool, ...
    archLabel, epsilon)

    Mtotal = size(SatPool,1);
    minTx = 4;

    comb4 = nchoosek(1:Mtotal, minTx);

    bestBCRB = inf;
    bestComb = [];
    bestPDOP = NaN;

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

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

    while currentBCRB > epsilon && ~isempty(candidates)

        bestAdd     = NaN;
        bestAddBCRB = inf;
        bestAddPDOP = NaN;

        for kk = 1:numel(candidates)

            cand = candidates(kk);
            trialSet = [selected, cand];

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
            warning('Nenhum candidato adicional válido encontrado.');
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