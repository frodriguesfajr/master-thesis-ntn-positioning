%% Scenario 01 - Arquiteturas isoladas
% Comparação HAPS, LEO, MEO e GEO
% Métricas: PDOP, BCRB e RMSE por Monte Carlo

close all;
clear;
clc;
format long;
%rng(1);  % reprodutibilidade

%% ##################### Saída ##########################################

outDir = 'results_scenario_01';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ##################### Cenário do usuário ###############################

lat = -22.8596582;      % [graus]
lon = -43.2303236;      % [graus]
h   = 10;               % [m]

UserPosition = llh2ecef(lat, lon, h);   % [x y z] ECEF [m]

%% ##################### Parâmetros globais ###############################

CNosim = 30:1:50;       % C/N0 de referência [dB-Hz]

% Para teste --> 500 ou 1000.
% Para resultado final --> 50000 se possível
Nexpe  = 100000;

c     = 3e8;        % [m/s]
beta  = 1.023e6;    % [Hz] largura de banda do código GPS C/A
T_coh = 0.02;       % [s] tempo de integração coerente

scenarioNames = {'HAPS','LEO','MEO','GEO'};
numScenarios = numel(scenarioNames);

RMSE_WLS_all  = zeros(numScenarios, numel(CNosim));
BCRB_all      = zeros(numScenarios, numel(CNosim));
PDOP_all      = zeros(numScenarios, 1);
numSV_all     = zeros(numScenarios, 1);
offset_all    = zeros(numScenarios, 1);
PosErr_MC_all = zeros(numScenarios, numel(CNosim), Nexpe);

%% ##################### Tabela de parâmetros dos cenários ##########

Arquitetura = ["HAPS"; "LEO"; "MEO"; "GEO"];

% Número de transmissores disponíveis em cada cenário isolado
NumTransmissores = [8; 8; 8; 8];

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
    NumTransmissores, ...
    MascaraElevacao_min_deg, ...
    DistanciaPropagacao_min_km, ...
    DistanciaPropagacao_max_km, ...
    PDOP_min, ...
    PDOP_max, ...
    Delta_CN0_rel_dB);

disp('===== Tabela 7.1 - Parâmetros dos cenários isolados =====');
disp(Tparam);

writetable(Tparam, fullfile(outDir, 'tabela_parametros_cenarios_isolados.csv'));

%% ##################### Loop dos cenários ###############################

for sc = 1:numScenarios

    scenario = scenarioNames{sc};

    %% ########### Configuração do cenário ###########

    % A tabela Tparam é a fonte única dos parâmetros do cenário.
    CN0_offset  = Delta_CN0_rel_dB(sc);
    elevMaskDeg = MascaraElevacao_min_deg(sc);
    pdopTarget  = [PDOP_min(sc), PDOP_max(sc)];
    rangeLimits = [DistanciaPropagacao_min_km(sc), ...
                   DistanciaPropagacao_max_km(sc)] * 1e3; % [m]
    numSV       = NumTransmissores(sc);

    numSV_all(sc)  = numSV;
    offset_all(sc) = CN0_offset;

    %% ########### Geração da geometria ###########

    geom = generateGeometryPDOP(UserPosition, numSV, ...
        'ElevMaskDeg', elevMaskDeg, ...
        'PdopTarget', pdopTarget, ...
        'RangeLimits', rangeLimits, ...
        'MaxTries', 10000);

    SatPosition = geom.SatPosition;
    H           = geom.H;
    PDOP        = geom.PDOP;
    PDOP_all(sc)= PDOP;

    fprintf('%s -> PDOP (%d transmissores): %.4f\n', scenario, numSV, PDOP);

    dXYZ = SatPosition - UserPosition;
    distances = sqrt(sum(dXYZ.^2, 2)).';

    %% ########### Estrutura do estimador WLS ###########

    wls_input.UserPosition  = UserPosition;
    wls_input.SatPosition   = SatPosition;
    wls_input.numSV         = numSV;
    wls_input.numIterations = 10;

    PosErrWLS = zeros(numel(CNosim), Nexpe);

    %% ########### Loop de C/N0 ###########

    for CNo_idx = 1:numel(CNosim)

        CN0_dBHz = CNosim(CNo_idx);

        CN0_dBHz_eff = CN0_dBHz + CN0_offset;
        CN0_linear = 10^(CN0_dBHz_eff/10);

        sigma_base = c / (2*pi*beta*sqrt(CN0_linear * T_coh));

        for exp_idx = 1:Nexpe

            x0 = UserPosition(:) + randn(3,1) * 400;
            clockBias = 30 + 5 * randn;

            sigma_rho_vec = sigma_base * ones(1, numSV);

            % Ruído nominal gaussiano
            noise_nom = sigma_rho_vec .* randn(1, numSV);

            % Erro adicional de média zero dependente de C/N0
            p_out = min(0.20, 0.20 * exp(-(CN0_dBHz - 30)/4));
            sigma_extra = 4 * sigma_base * exp(-(CN0_dBHz - 30)/6);

            noise_extra = (rand(1, numSV) < p_out) .* ...
                          (sigma_extra .* randn(1, numSV));

            rho_meas = distances + clockBias + noise_nom + noise_extra;

            wls_input.x0 = x0;
            wls_input.rho_meas = rho_meas;
            wls_input.sigma_rho_vec = sigma_rho_vec;

            pos_est_wls = conv2stepsPVT_WLS(wls_input);

            PosErrWLS(CNo_idx, exp_idx) = norm(pos_est_wls - UserPosition);
        end

        RMSE_WLS_all(sc, CNo_idx) = sqrt(mean(PosErrWLS(CNo_idx,:).^2));

        % BCRB ideal gaussiano
        sigma_rho_vec_crb = sigma_base * ones(numSV,1);
        R = diag(sigma_rho_vec_crb.^2);
        J = H.' * (R \ H);
        C = J \ eye(4);

        BCRB_all(sc, CNo_idx) = sqrt(trace(C(1:3,1:3)));
    end

    PosErr_MC_all(sc,:,:) = PosErrWLS;
end

%% ##################### Tabela de PDOP obtido #####################

Tpdop = table( ...
    string(scenarioNames(:)), ...
    numSV_all, ...
    PDOP_all, ...
    offset_all, ...
    'VariableNames', {'Arquitetura','Num_transmissores','PDOP_obtido','Delta_CN0_rel_dB'} ...
);

disp('===== Tabela 7.2 - PDOP obtido por arquitetura =====');
disp(Tpdop);

writetable(Tpdop, fullfile(outDir, 'tabela_pdop_obtido.csv'));

%% ##################### Tabelas completas por arquitetura ##########

ratio_all = RMSE_WLS_all ./ BCRB_all;

for sc = 1:numScenarios

    fprintf('\n===== %s =====\n', scenarioNames{sc});

    T = table(CNosim(:), RMSE_WLS_all(sc,:).', BCRB_all(sc,:).', ...
        ratio_all(sc,:).', ...
        'VariableNames', {'CN0_ref_dBHz','RMSE_WLS_m','BCRB_m','RMSE_sobre_BCRB'} ...
    );

    disp(T);

    writetable(T, fullfile(outDir, ...
        sprintf('tabela_%s_rmse_bcrb.csv', lower(scenarioNames{sc}))));
end

%% ##################### Tabela resumida em C/N0 selecionados ##########

CN0_sel = [30 35 40 45 50];
idx_sel = arrayfun(@(x) find(CNosim == x, 1), CN0_sel);

nRows = numScenarios * numel(CN0_sel);

Resumo_Arquitetura = strings(nRows,1);
Resumo_NumTx       = zeros(nRows,1);
Resumo_PDOP        = zeros(nRows,1);
Resumo_CN0         = zeros(nRows,1);
Resumo_BCRB        = zeros(nRows,1);
Resumo_RMSE        = zeros(nRows,1);
Resumo_Ratio       = zeros(nRows,1);

row = 1;
for sc = 1:numScenarios
    for kk = 1:numel(CN0_sel)

        idx = idx_sel(kk);

        Resumo_Arquitetura(row) = string(scenarioNames{sc});
        Resumo_NumTx(row)       = numSV_all(sc);
        Resumo_PDOP(row)        = PDOP_all(sc);
        Resumo_CN0(row)         = CN0_sel(kk);
        Resumo_BCRB(row)        = BCRB_all(sc,idx);
        Resumo_RMSE(row)        = RMSE_WLS_all(sc,idx);
        Resumo_Ratio(row)       = ratio_all(sc,idx);

        row = row + 1;
    end
end

Tsummary = table( ...
    Resumo_Arquitetura, ...
    Resumo_NumTx, ...
    Resumo_PDOP, ...
    Resumo_CN0, ...
    Resumo_BCRB, ...
    Resumo_RMSE, ...
    Resumo_Ratio, ...
    'VariableNames', {'Arquitetura','Num_transmissores','PDOP','CN0_ref_dBHz', ...
                      'BCRB_m','RMSE_WLS_m','RMSE_sobre_BCRB'} ...
);

disp('===== Tabela 7.3 - Resumo em C/N0 selecionados =====');
disp(Tsummary);

writetable(Tsummary, fullfile(outDir, 'tabela_resumo_cn0_selecionados.csv'));

%% ##################### C/N0 mínimo de eficiência ###############################

effThreshold = 1.10;   % RMSE até 10% acima de BCRB

CN0_eff_min = NaN(numScenarios,1);
idx_eff_min = NaN(numScenarios,1);
ratio_eff   = NaN(numScenarios,1);

for sc = 1:numScenarios

    ratio_sc = ratio_all(sc,:);

    % Primeiro C/N0 a partir do qual todos os pontos seguintes
    % permanecem abaixo do limiar de eficiência.
    for k = 1:numel(CNosim)
        if all(ratio_sc(k:end) <= effThreshold)
            idx_eff_min(sc) = k;
            CN0_eff_min(sc) = CNosim(k);
            ratio_eff(sc)   = ratio_sc(k);
            break;
        end
    end
end

T_eff = table( ...
    string(scenarioNames(:)), ...
    PDOP_all, ...
    CN0_eff_min, ...
    ratio_eff, ...
    'VariableNames', {'Arquitetura','PDOP','CN0_min_eficiente_dBHz','RMSE_sobre_BCRB'} ...
);

disp('===== C/N0 mínimo para eficiência do estimador WLS =====');
disp(T_eff);

writetable(T_eff, fullfile(outDir, 'tabela_cn0_min_eficiencia_wls.csv'));

%% ##################### Figura 7.1 - RMSE vs BCRB por arquitetura ##########

figure;

for sc = 1:numScenarios

    subplot(2,2,sc);
    hold on;
    grid on;

    h1 = semilogy(CNosim, RMSE_WLS_all(sc,:), '--o', 'LineWidth', 1.2);
    h2 = semilogy(CNosim, BCRB_all(sc,:), '-s', 'LineWidth', 1.2);

    if ~isnan(CN0_eff_min(sc))
        xline(CN0_eff_min(sc), ':', ...
            sprintf('%.0f dB-Hz', CN0_eff_min(sc)), ...
            'LineWidth', 1.0);
    end

    xlabel('C/N0 de referencia [dB-Hz]');
    ylabel('Erro de posicao [m]');
    title(sprintf('%s', scenarioNames{sc}));

    legend([h1 h2], {'RMSE WLS', '$\mathcal{B}_{\rm CRB}$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast');
end

sgtitle('RMSE WLS vs $\mathcal{B}_{\rm CRB}$', ...
    'Interpreter', 'latex');

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_1_rmse_wls_vs_bcrb_2x2.png'), ...
    'Resolution', 300);

%% ##################### Figura 7.2 - CDF 2x2 em C/N0 de referência ##########

CN0_ref_CDF = 50;
idx_CDF = find(CNosim == CN0_ref_CDF, 1);

figure;

for sc = 1:numScenarios

    subplot(2,2,sc);
    hold on;
    grid on;

    erros = squeeze(PosErr_MC_all(sc, idx_CDF, :));
    erros = sort(erros(:));

    F = (1:numel(erros))' / numel(erros);

    plot(erros, F, 'LineWidth', 1.5);

    xlabel('Erro de posição [m]', 'Interpreter', 'none');
    ylabel('Probabilidade acumulada', 'Interpreter', 'none');

    title(sprintf('%s: C/N0 = %d dB-Hz', ...
        scenarioNames{sc}, CNosim(idx_CDF)), ...
        'Interpreter', 'none');
end

sgtitle('CDF do erro de posição para C/N0 de referência fixo', ...
    'Interpreter', 'none');

exportgraphics(gcf, fullfile(outDir, ...
    'fig_7_2_cdf_2x2_cn0_50_dbhz.png'), ...
    'Resolution', 300);

%% ##################### Tabela - Percentis em C/N0 de referência ##########

CN0_ref_CDF = 50;
idx_CDF = find(CNosim == CN0_ref_CDF, 1);

Erro_50_m = NaN(numScenarios,1);
Erro_90_m = NaN(numScenarios,1);
Erro_95_m = NaN(numScenarios,1);

for sc = 1:numScenarios

    erros = squeeze(PosErr_MC_all(sc, idx_CDF, :));

    Erro_50_m(sc) = percentileValue(erros, 50);
    Erro_90_m(sc) = percentileValue(erros, 90);
    Erro_95_m(sc) = percentileValue(erros, 95);
end

TCDF_ref = table( ...
    string(scenarioNames(:)), ...
    CN0_ref_CDF * ones(numScenarios,1), ...
    Erro_50_m, ...
    Erro_90_m, ...
    Erro_95_m, ...
    'VariableNames', {'Arquitetura','CN0_ref_dBHz', ...
    'Erro_50_m','Erro_90_m','Erro_95_m'} ...
);

disp('===== Percentis da CDF em C/N0 de referência =====');
disp(TCDF_ref);

writetable(TCDF_ref, fullfile(outDir, ...
    'tabela_percentis_cdf_cn0_50_dbhz.csv'));

%% ##################### Salvar variáveis ###############################====

save(fullfile(outDir, 'scenario_01_results.mat'), ...
    'CNosim', ...
    'scenarioNames', ...
    'RMSE_WLS_all', ...
    'BCRB_all', ...
    'ratio_all', ...
    'PDOP_all', ...
    'numSV_all', ...
    'offset_all', ...
    'PosErr_MC_all', ...
    'Tparam', ...
    'Tpdop', ...
    'Tsummary', ...
    'T_eff', ...
    'TCDF_ref');

%% ##################### FUNÇÕES LOCAIS ###############################======

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

function pos_est_wls = conv2stepsPVT_WLS(wls_input)

    SatPosition   = wls_input.SatPosition;
    rho_meas      = wls_input.rho_meas(:);
    sigma_rho_vec = wls_input.sigma_rho_vec(:);
    numSV         = wls_input.numSV;
    numIts        = wls_input.numIterations;

    xhat = [wls_input.x0; 0];
    W = diag(1 ./ (sigma_rho_vec.^2));

    for it = 1:numIts
        H        = zeros(numSV, 4);
        rho_pred = zeros(numSV, 1);

        for k = 1:numSV
            vec = SatPosition(k,:).' - xhat(1:3);
            dist = norm(vec);
            if dist < eps
                dist = eps;
            end

            u = vec / dist;
            H(k,1:3) = -u.';
            H(k,4)   = 1;
            rho_pred(k) = dist + xhat(4);
        end

        res = rho_meas - rho_pred;
        delta = (H.' * W * H) \ (H.' * W * res);

        xhat = xhat + delta;

        if norm(delta(1:3)) < 1e-3 && abs(delta(4)) < 1e-3
            break;
        end
    end

    pos_est_wls = xhat(1:3).';
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

function val = percentileValue(x, p)
%  Calcula percentil sem depender de toolboxes adicionais.

    x = sort(x(:));
    n = numel(x);

    if n == 0
        val = NaN;
        return;
    end

    if p <= 0
        val = x(1);
        return;
    end

    if p >= 100
        val = x(end);
        return;
    end

    pos = 1 + (n - 1) * p / 100;
    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        val = x(lo);
    else
        w = pos - lo;
        val = (1 - w) * x(lo) + w * x(hi);
    end
end