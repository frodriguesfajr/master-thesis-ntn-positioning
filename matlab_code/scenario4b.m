%% Scenario 4b
% Restricao de disponibilidade de HAPS
%
% Casos:
%   A) maxHAPS = 4  -> referencia sem restricao adicional
%   B) maxHAPS = 2  -> disponibilidade limitada
%   C) maxHAPS = 0  -> ausencia de HAPS
%
% Mantem exatamente:
%   - os mesmos 26 candidatos
%   - a mesma epoch
%   - o mesmo modelo de C/N0
%   - beta = 1.023 MHz
%   - Tcoh = 20 ms
%   - meta de limite posicional <= 1 m
%
% A restricao de HAPS e aplicada:
%   1) na busca exaustiva do melhor subconjunto inicial de 4;
%   2) em cada passo da selecao gulosa.
%
% Entrada:
%   results_scenario4_candidate_pool/scenario4_candidate_pool_26.mat

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

outDir = "results_scenario4_haps_restrictions";

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ##################### Modelo estatistico #################################

c = 299792458;
beta = 1.023e6;
Tcoh = 20e-3;

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

%% ##################### H completa #########################################

r = PositionECEF - UserPositionECEF(:).';
d = sqrt(sum(r.^2,2));
u = r./d;

Hfull = [-u,ones(Npool,1)];

%% ##################### Configuracoes ######################################

target_m = 1.0;

maxHAPS_cases = [4 2 0];
caseNames = ["Reference_max4";"Limited_max2";"No_HAPS"];

nCases = numel(maxHAPS_cases);

fprintf('\n============================================================\n');
fprintf('SCENARIO 4b - RESTRICOES DE HAPS\n');
fprintf('Meta de limite posicional <= %.3f m\n',target_m);
fprintf('Casos: maxHAPS = 4, 2 e 0\n');
fprintf('============================================================\n');

%% ##################### Resultados agregados ###############################

MaxHAPS = maxHAPS_cases(:);
TargetReached = false(nCases,1);
Kfinal = nan(nCases,1);
FinalBound_m = nan(nCases,1);
FinalPDOP = nan(nCases,1);
CountHAPS = zeros(nCases,1);
CountLEO  = zeros(nCases,1);
CountMEO  = zeros(nCases,1);
CountGEO  = zeros(nCases,1);
Best4Bound_m = nan(nCases,1);
CandidatePoolCount = repmat(Npool,nCases,1);
MaxFeasibleSetSize = nan(nCases,1);
BestMaxFeasibleBound_m = nan(nCases,1);

selectedPerCase = cell(nCases,1);

%% ##################### Loop dos casos #####################################

for caseIdx = 1:nCases

    maxHAPS = maxHAPS_cases(caseIdx);
    caseName = caseNames(caseIdx);

    fprintf('\n\n############################################################\n');
    fprintf('CASO %s | maxHAPS = %d\n',caseName,maxHAPS);
    fprintf('############################################################\n');

    %% ---------- Pool e maior conjunto viavel ----------------------------
    % O pool candidato permanece com 26 objetos para maxHAPS=4 e 2.
    % Entretanto, para maxHAPS=2, o conjunto com todos os 26 NAO e viavel,
    % pois conteria 4 HAPS. Para diagnosticar a factibilidade corretamente,
    % calculamos o melhor conjunto de cardinalidade maxima que respeita a
    % restricao: todos os nao-HAPS + o melhor subconjunto permitido de HAPS.

    allowed = (1:Npool).';

    idxNonHAPS = find(Architecture ~= "HAPS");
    idxHAPS = find(Architecture == "HAPS");

    nHAPSuse = min(maxHAPS,numel(idxHAPS));
    MaxFeasibleSetSize(caseIdx) = numel(idxNonHAPS) + nHAPSuse;

    if nHAPSuse == 0

        bestMaxSet = idxNonHAPS(:).';
        [bestMaxBound,bestMaxPDOP,bestMaxRank] = ...
            subsetMetrics(bestMaxSet,Hfull,sigmaRho_m);

    elseif nHAPSuse == numel(idxHAPS)

        bestMaxSet = [idxNonHAPS(:).',idxHAPS(:).'];
        [bestMaxBound,bestMaxPDOP,bestMaxRank] = ...
            subsetMetrics(bestMaxSet,Hfull,sigmaRho_m);

    else

        combHAPS = nchoosek(idxHAPS,nHAPSuse);

        bestMaxBound = inf;
        bestMaxPDOP = inf;
        bestMaxRank = 0;
        bestMaxSet = [];

        for kk = 1:size(combHAPS,1)

            idxTrial = [idxNonHAPS(:).',combHAPS(kk,:)];

            [b,p,rk] = subsetMetrics(idxTrial,Hfull,sigmaRho_m);

            if rk == 4 && b < bestMaxBound
                bestMaxBound = b;
                bestMaxPDOP = p;
                bestMaxRank = rk;
                bestMaxSet = idxTrial;
            end
        end
    end

    BestMaxFeasibleBound_m(caseIdx) = bestMaxBound;

    fprintf('\nMaior conjunto viavel sob a restricao:\n');
    fprintf('Pool candidato         : %d\n',Npool);
    fprintf('N maximo viavel        : %d\n',MaxFeasibleSetSize(caseIdx));
    fprintf('HAPS no conjunto       : %d\n',sum(Architecture(bestMaxSet)=="HAPS"));
    fprintf('Limite posicional      : %.12f m\n',bestMaxBound);
    fprintf('PDOP                   : %.12f\n',bestMaxPDOP);
    fprintf('rank(H)                : %d\n',bestMaxRank);

    if bestMaxBound > target_m
        fprintf(['ATENCAO: mesmo no melhor conjunto de cardinalidade maxima ' ...
                 'que respeita a restricao, a meta de %.3f m nao e atingida.\n'], ...
                 target_m);
    end

    % Para maxHAPS=0, HAPS sao removidos do pool de selecao.
    if maxHAPS == 0
        allowed = idxNonHAPS;
    end

    %% ---------- Busca exaustiva do melhor 4 ------------------------------

    fprintf('\n===== Busca exaustiva inicial | maxHAPS=%d =====\n',maxHAPS);

    comb4 = nchoosek(allowed,4);

    best4Bound = inf;
    best4 = [];
    nEvaluatedFeasible = 0;

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

        if sum(Architecture(idx)=="HAPS") > maxHAPS
            continue;
        end

        nEvaluatedFeasible = nEvaluatedFeasible + 1;

        [b,~,rankH] = subsetMetrics(idx,Hfull,sigmaRho_m);

        if rankH == 4 && b < best4Bound
            best4Bound = b;
            best4 = idx;
        end
    end

    if isempty(best4)
        warning('Nenhum subconjunto inicial de 4 valido para maxHAPS=%d.',maxHAPS);
        continue;
    end

    Best4Bound_m(caseIdx) = best4Bound;

    fprintf('Combinacoes viaveis avaliadas : %d\n',nEvaluatedFeasible);
    fprintf('Melhor limite com K=4         : %.12f m\n',best4Bound);
    fprintf('Indices                       : ');
    fprintf('%d ',best4);
    fprintf('\nLabels                        : ');
    fprintf('%s ',Label(best4));
    fprintf('\nComposicao K=4                : HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
        sum(Architecture(best4)=="HAPS"), ...
        sum(Architecture(best4)=="LEO"), ...
        sum(Architecture(best4)=="MEO"), ...
        sum(Architecture(best4)=="GEO"));

    %% ---------- Forward greedy com restricao -----------------------------

    selected = best4(:).';
    currentBound = best4Bound;

    traceK = numel(selected);
    traceAddedIndex = NaN;
    traceAddedLabel = "";
    traceBound = currentBound;
    traceHAPS = sum(Architecture(selected)=="HAPS");

    fprintf('\n===== Selecao gulosa | maxHAPS=%d =====\n',maxHAPS);
    fprintf('K=%d | HAPS=%d | limite = %.12f m | inicializacao\n', ...
        numel(selected),traceHAPS,currentBound);

    while currentBound > target_m && numel(selected) < numel(allowed)

        remaining = setdiff(allowed,selected,'stable');

        % Se o limite de HAPS ja foi atingido, remover HAPS restantes.
        currentHAPS = sum(Architecture(selected)=="HAPS");

        if currentHAPS >= maxHAPS
            remaining = remaining(Architecture(remaining) ~= "HAPS");
        end

        if isempty(remaining)
            break;
        end

        bestCandidate = NaN;
        bestCandidateBound = inf;

        for k = 1:numel(remaining)

            candidate = remaining(k);
            idxTrial = [selected,candidate];

            if sum(Architecture(idxTrial)=="HAPS") > maxHAPS
                continue;
            end

            [b,~,rankH] = subsetMetrics(idxTrial,Hfull,sigmaRho_m);

            if rankH == 4 && b < bestCandidateBound
                bestCandidateBound = b;
                bestCandidate = candidate;
            end
        end

        if isnan(bestCandidate)
            break;
        end

        selected(end+1) = bestCandidate;
        currentBound = bestCandidateBound;

        traceK(end+1,1) = numel(selected); %#ok<SAGROW>
        traceAddedIndex(end+1,1) = bestCandidate; %#ok<SAGROW>
        traceAddedLabel(end+1,1) = Label(bestCandidate); %#ok<SAGROW>
        traceBound(end+1,1) = currentBound; %#ok<SAGROW>
        traceHAPS(end+1,1) = sum(Architecture(selected)=="HAPS"); %#ok<SAGROW>

        fprintf('K=%d | HAPS=%d | adiciona %-6s | limite = %.12f m\n', ...
            numel(selected), ...
            sum(Architecture(selected)=="HAPS"), ...
            Label(bestCandidate), ...
            currentBound);
    end

    %% ---------- Resultado ------------------------------------------------

    [finalBound,finalPDOP,finalRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    reached = finalBound <= target_m;

    TargetReached(caseIdx) = reached;
    Kfinal(caseIdx) = numel(selected);
    FinalBound_m(caseIdx) = finalBound;
    FinalPDOP(caseIdx) = finalPDOP;

    CountHAPS(caseIdx) = sum(Architecture(selected)=="HAPS");
    CountLEO(caseIdx)  = sum(Architecture(selected)=="LEO");
    CountMEO(caseIdx)  = sum(Architecture(selected)=="MEO");
    CountGEO(caseIdx)  = sum(Architecture(selected)=="GEO");

    selectedPerCase{caseIdx} = selected;

    fprintf('\n===== RESULTADO | maxHAPS=%d =====\n',maxHAPS);
    fprintf('K final            : %d\n',numel(selected));
    fprintf('Limite posicional  : %.12f m\n',finalBound);
    fprintf('PDOP               : %.12f\n',finalPDOP);
    fprintf('rank(H)            : %d\n',finalRank);
    fprintf('Meta <= %.3f m     : %s\n',target_m,string(reached));
    fprintf('Composicao         : HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
        CountHAPS(caseIdx),CountLEO(caseIdx), ...
        CountMEO(caseIdx),CountGEO(caseIdx));

    fprintf('Labels selecionados: ');
    fprintf('%s ',Label(selected));
    fprintf('\n');

    %% ---------- Salvar trace e subset ------------------------------------

    Ttrace = table( ...
        traceK,traceAddedIndex,traceAddedLabel,traceHAPS,traceBound, ...
        'VariableNames',{ ...
        'K','AddedPoolIndex','AddedLabel','HAPSCount','PositionalBound_m'});

    Tselected = table( ...
        (1:numel(selected)).',selected(:), ...
        Label(selected(:)),Architecture(selected(:)), ...
        NORAD(selected(:)),ObjectName(selected(:)), ...
        Elevation_deg(selected(:)),SlantRange_m(selected(:))/1e3, ...
        CN0_dBHz(selected(:)),sigmaRho_m(selected(:)), ...
        'VariableNames',{ ...
        'SelectionOrder','PoolIndex','Label','Architecture', ...
        'NORAD_CAT_ID','ObjectName','Elevation_deg','SlantRange_km', ...
        'CN0_dBHz','SigmaRho_m'});

    writetable(Ttrace,fullfile(outDir, ...
        sprintf('scenario4_trace_maxHAPS_%d.csv',maxHAPS)));

    writetable(Tselected,fullfile(outDir, ...
        sprintf('scenario4_selected_maxHAPS_%d.csv',maxHAPS)));

end

%% ##################### Comparacao final ###################################

Tcomparison = table( ...
    caseNames,MaxHAPS,CandidatePoolCount,MaxFeasibleSetSize, ...
    BestMaxFeasibleBound_m,Best4Bound_m,TargetReached,Kfinal, ...
    FinalBound_m,FinalPDOP,CountHAPS,CountLEO,CountMEO,CountGEO, ...
    'VariableNames',{ ...
    'Case','MaxHAPS','CandidatePoolCount','MaxFeasibleSetSize', ...
    'BestMaxFeasibleBound_m','Best4Bound_m','TargetReached','Kfinal', ...
    'FinalBound_m','FinalPDOP','HAPS','LEO','MEO','GEO'});

disp(' ');
disp('============================================================');
disp('COMPARACAO FINAL DAS RESTRICOES DE HAPS');
disp('============================================================');
disp(Tcomparison);

writetable(Tcomparison, ...
    fullfile(outDir,'scenario4_haps_restrictions_comparison.csv'));

save(fullfile(outDir,'scenario4_haps_restrictions_results.mat'));

fprintf('\nResultados salvos em: %s\n',outDir);

%% ##################### FUNCAO LOCAL #######################################

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
