%% Scenario 02 - Seleção gulosa baseada em BCRB
% Seleção de transmissores HAPS, LEO, MEO e GEO
% Critério: menor subconjunto que atende BCRB <= epsilon

close all;
clear;
clc;
format long;
rng(2);  % reprodutibilidade do Scenario 02

%% ##################### Saída ##############################################

outDir = 'results_scenario_02';
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

CN0_ref_dBHz = 50;  % C/N0 de referência [dB-Hz]
epsilon      = 1.; % precisão-alvo [m]

scenarioNames = {'HAPS','LEO','MEO','GEO'};
numScenarios = numel(scenarioNames);

%% ##################### Tabela de candidatos disponíveis ##################

Arquitetura = ["HAPS"; "LEO"; "MEO"; "GEO"];

% Número de dispositivos candidatos por arquitetura
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

disp('===== Tabela 7.4 - Parâmetros dos candidatos disponíveis =====');
disp(Tparam);

writetable(Tparam, fullfile(outDir, ...
    'tabela_parametros_candidatos_disponiveis.csv'));

%% ##################### Parâmetros do problema de seleção ##############

Tsel_param = table( ...
    CN0_ref_dBHz, ...
    epsilon, ...
    sum(NumDispositivos), ...
    NumDispositivos(1), ...
    NumDispositivos(2), ...
    NumDispositivos(3), ...
    NumDispositivos(4), ...
    'VariableNames', {'CN0_ref_dBHz','Precisao_alvo_m','Total_candidatos', ...
    'N_HAPS','N_LEO','N_MEO','N_GEO'} ...
);

disp('===== Parâmetros do problema de seleção =====');
disp(Tsel_param);

writetable(Tsel_param, fullfile(outDir, ...
    'tabela_parametros_problema_selecao.csv'));

%% ##################### Geração do conjunto candidato #####################

SatPool   = [];
archLabel = strings(0,1);
sigmaPool = [];
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

    CN0_eff_dBHz = CN0_ref_dBHz + Delta_CN0_rel_dB(sc);
    CN0_linear   = 10^(CN0_eff_dBHz/10);

    sigma_rho = c / (2*pi*beta*sqrt(CN0_linear*T_coh));

    SatPool   = [SatPool; SatPosition_sc];
    archLabel = [archLabel; repmat(string(scenario), numSV, 1)];
    sigmaPool = [sigmaPool; sigma_rho * ones(numSV,1)];
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
    sigmaPool, ...
    'VariableNames', {'ID','Arquitetura','Distancia_km','Sigma_rho_m'} ...
);

disp('===== Lista de candidatos gerados =====');
disp(Tcand);

writetable(Tcand, fullfile(outDir, 'tabela_candidatos_gerados.csv'));

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

%% ##################### Melhor subconjunto inicial com 4 dispositivos #############

fprintf('\n===== Busca do melhor subconjunto inicial com 4 dispositivos =====\n');

comb4 = nchoosek(1:Mtotal, 4);

bestBCRB = inf;
bestPDOP = inf;
bestComb = [];

for k = 1:size(comb4,1)

    idx = comb4(k,:);

    [BCRB_k, PDOP_k, valid_k] = computeBCRB_fromSubset( ...
        UserPosition, SatPool, sigmaPool, idx);

    if ~valid_k
        continue;
    end

    if BCRB_k < bestBCRB
        bestBCRB = BCRB_k;
        bestPDOP = PDOP_k;
        bestComb = idx;
    end
end

if isempty(bestComb)
    error('Não foi possível encontrar subconjunto inicial válido com 4 dispositivos.');
end

selected   = bestComb(:).';
candidates = setdiff(1:Mtotal, selected);

currentBCRB = bestBCRB;
currentPDOP = bestPDOP;

fprintf('Melhor subconjunto inicial:\n');
fprintf('IDs selecionados: ');
fprintf('%d ', selected);
fprintf('\n');

fprintf('Arquiteturas: ');
fprintf('%s ', archLabel(selected));
fprintf('\n');

fprintf('BCRB inicial: %.6f m\n', currentBCRB);
fprintf('PDOP inicial: %.4f\n', currentPDOP);

%% ##################### Histórico do algoritmo #####################

Iteracao      = 0;
DispositivoAdicionado = 0;
ArquiteturaAdicionada = "Inicial";
TotalSelecionado = numel(selected);
BCRB_hist     = currentBCRB;
PDOP_hist     = currentPDOP;
N_HAPS_hist   = sum(archLabel(selected) == "HAPS");
N_LEO_hist    = sum(archLabel(selected) == "LEO");
N_MEO_hist    = sum(archLabel(selected) == "MEO");
N_GEO_hist    = sum(archLabel(selected) == "GEO");

%% ##################### Algoritmo guloso ###################################

iter = 0;

while currentBCRB > epsilon && ~isempty(candidates)

    iter = iter + 1;

    bestAdd     = NaN;
    bestAddBCRB = inf;
    bestAddPDOP = inf;

    for kk = 1:numel(candidates)

        cand = candidates(kk);
        testSet = [selected, cand];

        [BCRB_k, PDOP_k, valid_k] = computeBCRB_fromSubset( ...
            UserPosition, SatPool, sigmaPool, testSet);

        if ~valid_k
            continue;
        end

        if BCRB_k < bestAddBCRB
            bestAdd     = cand;
            bestAddBCRB = BCRB_k;
            bestAddPDOP = PDOP_k;
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

    fprintf('\nIteração %d:\n', iter);
    fprintf('Dispositivo adicionado: T%d (%s)\n', ...
        bestAdd, archLabel(bestAdd));
    fprintf('Total selecionado: %d\n', numel(selected));
    fprintf('BCRB atual: %.6f m\n', currentBCRB);
    fprintf('PDOP atual: %.4f\n', currentPDOP);

    Iteracao(end+1,1) = iter;
    DispositivoAdicionado(end+1,1) = bestAdd;
    ArquiteturaAdicionada(end+1,1) = archLabel(bestAdd);
    TotalSelecionado(end+1,1) = numel(selected);
    BCRB_hist(end+1,1) = currentBCRB;
    PDOP_hist(end+1,1) = currentPDOP;
    N_HAPS_hist(end+1,1) = sum(archLabel(selected) == "HAPS");
    N_LEO_hist(end+1,1)  = sum(archLabel(selected) == "LEO");
    N_MEO_hist(end+1,1)  = sum(archLabel(selected) == "MEO");
    N_GEO_hist(end+1,1)  = sum(archLabel(selected) == "GEO");
end

atingiu_meta = currentBCRB <= epsilon;

%% ##################### Resultados finais #####################

N_HAPS_final = sum(archLabel(selected) == "HAPS");
N_LEO_final  = sum(archLabel(selected) == "LEO");
N_MEO_final  = sum(archLabel(selected) == "MEO");
N_GEO_final  = sum(archLabel(selected) == "GEO");

fprintf('\n===== Configuração final =====\n');
fprintf('HAPS: %d\n', N_HAPS_final);
fprintf('LEO : %d\n', N_LEO_final);
fprintf('MEO : %d\n', N_MEO_final);
fprintf('GEO : %d\n', N_GEO_final);
fprintf('Total: %d\n', numel(selected));
fprintf('BCRB final: %.6f m\n', currentBCRB);
fprintf('Precisão-alvo: %.2f m\n', epsilon);

if atingiu_meta
    fprintf('Resultado: meta atingida.\n');
else
    fprintf('Resultado: meta não atingida.\n');
end

Tfinal = table( ...
    CN0_ref_dBHz, ...
    epsilon, ...
    N_HAPS_final, ...
    N_LEO_final, ...
    N_MEO_final, ...
    N_GEO_final, ...
    numel(selected), ...
    currentPDOP, ...
    currentBCRB, ...
    atingiu_meta, ...
    'VariableNames', {'CN0_ref_dBHz','Precisao_alvo_m', ...
    'HAPS','LEO','MEO','GEO','Total_dispositivos', ...
    'PDOP_final','BCRB_final_m','Meta_atingida'} ...
);

disp('===== Tabela 7.5 - Configuração selecionada =====');
disp(Tfinal);

writetable(Tfinal, fullfile(outDir, ...
    'tabela_configuracao_selecionada.csv'));

Tselected = table( ...
    selected(:), ...
    archLabel(selected), ...
    rangePool(selected)/1e3, ...
    sigmaPool(selected), ...
    'VariableNames', {'ID','Arquitetura','Distancia_km','Sigma_rho_m'} ...
);

disp('===== Dispositivos selecionados =====');
disp(Tselected);

writetable(Tselected, fullfile(outDir, ...
    'tabela_dispositivos_selecionados.csv'));

Titer = table( ...
    Iteracao, ...
    DispositivoAdicionado, ...
    ArquiteturaAdicionada, ...
    TotalSelecionado, ...
    BCRB_hist, ...
    PDOP_hist, ...
    N_HAPS_hist, ...
    N_LEO_hist, ...
    N_MEO_hist, ...
    N_GEO_hist, ...
    'VariableNames', {'Iteracao','Dispositivo_adicionado','Arquitetura_adicionada', ...
    'Total_selecionado','BCRB_m','PDOP','HAPS','LEO','MEO','GEO'} ...
);

disp('===== Histórico do algoritmo guloso =====');
disp(Titer);

writetable(Titer, fullfile(outDir, ...
    'tabela_historico_algoritmo_guloso.csv'));

%% ##################### Figura 7.3 - Evolução do BCRB #####################

figure;
hold on;
grid on;

plot(TotalSelecionado, BCRB_hist, '-o', 'LineWidth', 1.5);

yline(epsilon, '--', sprintf('Meta %.1f m', epsilon), ...
    'LineWidth', 1.2);

xlabel('Número de dispositivos selecionados', 'Interpreter', 'none');
ylabel('$\mathcal{B}_{\rm CRB}$ posicional [m]', 'Interpreter', 'latex');

title('$\mathcal{B}_{\rm CRB}$ no algoritmo guloso', ...
    'Interpreter', 'latex');

xticks(TotalSelecionado);
xlim([min(TotalSelecionado)-0.5, max(TotalSelecionado)+0.5]);

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_3_evolucao_bcrb_algoritmo_guloso.png'), ...
    'Resolution', 300);

%% ##################### Figura 7.4 - Composição selecionada ##############

figure;
hold on;
grid on;

counts = [N_HAPS_final, N_LEO_final, N_MEO_final, N_GEO_final];

cats = categorical({'HAPS','LEO','MEO','GEO'}, ...
    {'HAPS','LEO','MEO','GEO'}, ...
    'Ordinal', true);

bar(cats, counts, 0.55);

xlabel('Arquitetura', 'Interpreter', 'none');
ylabel('Número de dispositivos selecionados', 'Interpreter', 'none');
title('Composição da configuração selecionada', 'Interpreter', 'none');

ylim([0, max(counts)+1]);
yticks(0:1:max(counts)+1);

for ii = 1:numel(counts)
    text(ii, counts(ii), sprintf('%d', counts(ii)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10);
end

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_4_composicao_configuracao_selecionada.png'), ...
    'Resolution', 300);

%% ##################### Salvar variáveis ###################################

save(fullfile(outDir, 'scenario_02_results.mat'), ...
    'CN0_ref_dBHz', ...
    'epsilon', ...
    'scenarioNames', ...
    'SatPool', ...
    'archLabel', ...
    'sigmaPool', ...
    'rangePool', ...
    'selected', ...
    'Tparam', ...
    'Tsel_param', ...
    'Tcand', ...
    'Tpdop_arch', ...
    'Tfinal', ...
    'Tselected', ...
    'Titer', ...
    'BCRB_hist', ...
    'PDOP_hist', ...
    'TotalSelecionado');

%% ##################### FUNÇÕES LOCAIS ###################################==

function ecef = llh2ecef(lat_deg, lon_deg, h_m)
%   Converte coordenadas geodésicas WGS-84 para ECEF [x y z].

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

        [H_try, PDOP_try, u_try, range_try] = localGeometryAndPDOP(UserPosition, SatPosition_try);

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

function [BCRB, PDOP, valid] = computeBCRB_fromSubset(UserPosition, SatPool, sigmaPool, idx)

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

    C = J \ eye(4);

    BCRB = sqrt(trace(C(1:3,1:3)));
    PDOP = PDOP_tmp;
    valid = true;
end