%% Scenario 4a
% Modelo de qualidade das pseudodistancias + limite posicional + selecao
%
% Entrada:
%   results_scenario4_candidate_pool/scenario4_candidate_pool_26.mat
%
% Modelo controlado (consistente com os cenarios anteriores):
%   C/N0 de referencia = 50 dB-Hz
%   offsets por arquitetura:
%       HAPS +8 dB
%       LEO  +4 dB
%       MEO   0 dB
%       GEO  -4 dB
%
%   beta  = 1.023 MHz
%   Tcoh  = 20 ms
%
%   sigma_rho_i = c/(2*pi*beta*sqrt((C/N0)_i*Tcoh))
%
% Selecao:
%   - inicializacao exaustiva com o melhor subconjunto de 4 transmissores
%   - adicao gulosa do candidato que minimiza o limite posicional
%   - parada quando limite posicional <= 1 m
%
% Observacao:
%   Nesta etapa NAO ha link budget completo nem dependencia adicional de
%   distancia no C/N0. Os offsets representam diferencas relativas
%   controladas entre arquiteturas.

close all;
clear;
clc;
format long;

%% ##################### Entrada ############################################

inputMat = fullfile( ...
    "results_scenario4_candidate_pool", ...
    "scenario4_candidate_pool_26.mat");

if ~isfile(inputMat)
    error("Arquivo MAT nao encontrado: %s",inputMat);
end

S = load(inputMat);

requiredVars = [ ...
    "UserPositionECEF", ...
    "PositionECEF", ...
    "Label", ...
    "Architecture", ...
    "NORAD", ...
    "ObjectName", ...
    "Elevation_deg", ...
    "SlantRange_m" ...
];

for k = 1:numel(requiredVars)
    if ~isfield(S,requiredVars(k))
        error("Variavel ausente no MAT: %s",requiredVars(k));
    end
end

UserPositionECEF = S.UserPositionECEF;
PositionECEF     = S.PositionECEF;
Label            = string(S.Label);
Architecture     = string(S.Architecture);
NORAD            = string(S.NORAD);
ObjectName       = string(S.ObjectName);
Elevation_deg    = S.Elevation_deg;
SlantRange_m     = S.SlantRange_m;

Npool = size(PositionECEF,1);

if Npool ~= 26
    error("Esperados 26 candidatos, obtidos %d.",Npool);
end

outDir = "results_scenario4_selection";

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ##################### Parametros estatisticos ############################

c = 299792458;      % [m/s]
beta = 1.023e6;     % [Hz]
Tcoh = 20e-3;       % [s]

CN0_ref_dBHz = 50;

deltaCN0_dB = zeros(Npool,1);

deltaCN0_dB(Architecture=="HAPS") = +8;
deltaCN0_dB(Architecture=="LEO")  = +4;
deltaCN0_dB(Architecture=="MEO")  =  0;
deltaCN0_dB(Architecture=="GEO")  = -4;

CN0_dBHz = CN0_ref_dBHz + deltaCN0_dB;
CN0_linear = 10.^(CN0_dBHz/10);

sigmaRho_m = c ./ ...
    (2*pi*beta.*sqrt(CN0_linear*Tcoh));

VarianceRho_m2 = sigmaRho_m.^2;

%% ##################### Matriz H completa ##################################

r = PositionECEF - UserPositionECEF(:).';
d = sqrt(sum(r.^2,2));
u = r./d;

Hfull = [-u,ones(Npool,1)];

%% ##################### Tabela de qualidade ################################

PoolIndex = (1:Npool).';

Tquality = table( ...
    PoolIndex,Label,Architecture,NORAD,ObjectName, ...
    Elevation_deg,SlantRange_m/1e3, ...
    deltaCN0_dB,CN0_dBHz,sigmaRho_m,VarianceRho_m2, ...
    'VariableNames',{ ...
    'PoolIndex','Label','Architecture','NORAD_CAT_ID','ObjectName', ...
    'Elevation_deg','SlantRange_km', ...
    'DeltaCN0_dB','CN0_dBHz','SigmaRho_m','VarianceRho_m2'});

disp(' ');
disp('===== Qualidade nominal das pseudodistancias =====');
disp(Tquality);

%% ##################### Diagnostico ponderado ##############################

archList = ["HAPS","LEO","MEO","GEO","ALL"];
nArch = numel(archList);

N = zeros(nArch,1);
RankH = zeros(nArch,1);
PDOP = nan(nArch,1);
PositionalBound_m = nan(nArch,1);

for aidx = 1:nArch

    if archList(aidx) == "ALL"
        idx = (1:Npool).';
    else
        idx = find(Architecture == archList(aidx));
    end

    N(aidx) = numel(idx);

    [PositionalBound_m(aidx),PDOP(aidx),RankH(aidx)] = ...
        subsetMetrics(idx,Hfull,sigmaRho_m);
end

TweightedSummary = table( ...
    archList.',N,RankH,PDOP,PositionalBound_m, ...
    'VariableNames',{ ...
    'Architecture','N','RankH','PDOP','PositionalBound_m'});

disp(' ');
disp('===== Resumo ponderado =====');
disp(TweightedSummary);

%% ##################### Meta ################################################

target_m = 1.0;

fprintf('\nMeta de precisao para selecao: %.3f m\n',target_m);

%% ##################### Melhor subconjunto inicial de 4 ####################

fprintf('\n===== Busca exaustiva do melhor subconjunto de 4 =====\n');

comb4 = nchoosek(1:Npool,4);
nComb4 = size(comb4,1);

best4Bound = inf;
best4 = [];

for k = 1:nComb4

    idx = comb4(k,:);

    [b,~,rankH] = subsetMetrics(idx,Hfull,sigmaRho_m);

    if rankH == 4 && b < best4Bound
        best4Bound = b;
        best4 = idx;
    end
end

if isempty(best4)
    error('Nenhum subconjunto de 4 com posto completo foi encontrado.');
end

fprintf('Combinacoes avaliadas : %d\n',nComb4);
fprintf('Melhor limite posicional com 4 transmissores: %.12f m\n',best4Bound);
fprintf('Indices: ');
fprintf('%d ',best4);
fprintf('\nLabels : ');
fprintf('%s ',Label(best4));
fprintf('\n');

%% ##################### Forward greedy #####################################

selected = best4(:).';
currentBound = best4Bound;

greedyK = numel(selected);
greedyAddedIndex = NaN;
greedyAddedLabel = "";
greedyBound = currentBound;

traceK = greedyK;
traceAddedIndex = greedyAddedIndex;
traceAddedLabel = greedyAddedLabel;
traceBound = greedyBound;

fprintf('\n===== Selecao gulosa =====\n');
fprintf('K=%d | limite posicional = %.12f m | inicializacao\n', ...
    numel(selected),currentBound);

while currentBound > target_m && numel(selected) < Npool

    remaining = setdiff(1:Npool,selected,'stable');

    bestCandidate = NaN;
    bestCandidateBound = inf;

    for k = 1:numel(remaining)

        candidate = remaining(k);
        idxTrial = [selected,candidate];

        [b,~,rankH] = subsetMetrics(idxTrial,Hfull,sigmaRho_m);

        if rankH == 4 && b < bestCandidateBound
            bestCandidateBound = b;
            bestCandidate = candidate;
        end
    end

    if isnan(bestCandidate)
        warning('Nenhum candidato adicional valido encontrado.');
        break;
    end

    selected(end+1) = bestCandidate;
    currentBound = bestCandidateBound;

    traceK(end+1,1) = numel(selected); %#ok<SAGROW>
    traceAddedIndex(end+1,1) = bestCandidate; %#ok<SAGROW>
    traceAddedLabel(end+1,1) = Label(bestCandidate); %#ok<SAGROW>
    traceBound(end+1,1) = currentBound; %#ok<SAGROW>

    fprintf('K=%d | adiciona %-6s | limite posicional = %.12f m\n', ...
        numel(selected),Label(bestCandidate),currentBound);
end

[finalBound,finalPDOP,finalRankH] = ...
    subsetMetrics(selected,Hfull,sigmaRho_m);

targetReached = finalBound <= target_m;

%% ##################### Resultado final ####################################

Tselected = Tquality(selected,:);
Tselected.SelectionOrder = (1:height(Tselected)).';

% Reordenar para deixar SelectionOrder primeiro.
Tselected = movevars(Tselected,'SelectionOrder','Before','PoolIndex');

disp(' ');
disp('===== Subconjunto selecionado =====');
disp(Tselected);

fprintf('\n===== Resultado da selecao =====\n');
fprintf('K final             : %d\n',numel(selected));
fprintf('Limite posicional   : %.12f m\n',finalBound);
fprintf('PDOP                : %.12f\n',finalPDOP);
fprintf('rank(H)             : %d\n',finalRankH);
fprintf('Meta <= %.3f m      : %s\n',target_m,string(targetReached));

fprintf('Composicao: HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
    sum(Architecture(selected)=="HAPS"), ...
    sum(Architecture(selected)=="LEO"), ...
    sum(Architecture(selected)=="MEO"), ...
    sum(Architecture(selected)=="GEO"));

%% ##################### Baselines mesmo K ##################################
% Baseline 1:
%   maior C/N0; em empates, maior elevacao.
%
% Baseline 2:
%   maior elevacao.

Kfinal = numel(selected);

idxAll = (1:Npool).';

TsortCN0 = table( ...
    idxAll,CN0_dBHz,Elevation_deg, ...
    'VariableNames',{'Index','CN0','Elevation'});

TsortCN0 = sortrows(TsortCN0, ...
    {'CN0','Elevation'}, ...
    {'descend','descend'});

baselineCN0 = TsortCN0.Index(1:Kfinal).';

TsortElev = table( ...
    idxAll,Elevation_deg,CN0_dBHz, ...
    'VariableNames',{'Index','Elevation','CN0'});

TsortElev = sortrows(TsortElev, ...
    {'Elevation','CN0'}, ...
    {'descend','descend'});

baselineElev = TsortElev.Index(1:Kfinal).';

[bCN0,pdopCN0,~] = subsetMetrics( ...
    baselineCN0,Hfull,sigmaRho_m);

[bElev,pdopElev,~] = subsetMetrics( ...
    baselineElev,Hfull,sigmaRho_m);

fprintf('\n===== Baselines com K=%d =====\n',Kfinal);

fprintf('Maior C/N0  : limite = %.12f m | PDOP = %.12f\n', ...
    bCN0,pdopCN0);
fprintf('Labels       : ');
fprintf('%s ',Label(baselineCN0));
fprintf('\n');

fprintf('Maior elev.  : limite = %.12f m | PDOP = %.12f\n', ...
    bElev,pdopElev);
fprintf('Labels       : ');
fprintf('%s ',Label(baselineElev));
fprintf('\n');

reductionVsCN0_pct = 100*(bCN0-finalBound)/bCN0;
reductionVsElev_pct = 100*(bElev-finalBound)/bElev;

fprintf('Reducao vs maior C/N0 : %.6f %%\n',reductionVsCN0_pct);
fprintf('Reducao vs maior elev.: %.6f %%\n',reductionVsElev_pct);

%% ##################### Trace ##############################################

Ttrace = table( ...
    traceK,traceAddedIndex,traceAddedLabel,traceBound, ...
    'VariableNames',{ ...
    'K','AddedPoolIndex','AddedLabel','PositionalBound_m'});

disp(' ');
disp('===== Trace guloso =====');
disp(Ttrace);

%% ##################### Salvar #############################################

writetable(Tquality, ...
    fullfile(outDir,'scenario4_measurement_quality.csv'));

writetable(TweightedSummary, ...
    fullfile(outDir,'scenario4_weighted_summary.csv'));

writetable(Tselected, ...
    fullfile(outDir,'scenario4_selected_subset.csv'));

writetable(Ttrace, ...
    fullfile(outDir,'scenario4_greedy_trace.csv'));

save(fullfile(outDir,'scenario4_selection_results.mat'));

%% ##################### FUNCOES LOCAIS #####################################

function [bound_m,PDOP,rankH] = subsetMetrics( ...
    idx,Hfull,sigmaRho_m)

    idx = idx(:);

    H = Hfull(idx,:);
    sigma = sigmaRho_m(idx);

    rankH = rank(H);
    bound_m = inf;
    PDOP = inf;

    if rankH < 4
        return;
    end

    % PDOP geometrico
    A = H.'*H;

    if rcond(A) > 1e-14
        Cgeom = A\eye(4);
        PDOP = sqrt(trace(Cgeom(1:3,1:3)));
    end

    % Limite posicional ponderado
    R = diag(sigma.^2);
    J = H.'*(R\H);

    if rcond(J) <= 1e-14
        return;
    end

    C = J\eye(4);

    trPos = trace(C(1:3,1:3));

    if trPos > 0 && isfinite(trPos)
        bound_m = sqrt(trPos);
    end
end
