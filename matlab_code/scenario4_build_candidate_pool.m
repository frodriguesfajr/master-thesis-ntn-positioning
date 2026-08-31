%% Scenario 4
% Montagem do conjunto candidato final de 26 transmissores em t0
%
% Conjunto:
%   4 HAPS sinteticos
%   8 LEO Starlink reais (TLE historico + SGP4)
%   7 MEO O3B reais      (TLE historico + SGP4)
%   7 GEO reais          (TLE historico + SGP4)
%
% Epoch comum:
%   01-Aug-2026 12:00:00 UTC
%
% Esta etapa e SOMENTE geometrica.
% Ainda nao aplica C/N0, sigma_rho, R ou algoritmo guloso.
%
% Saidas:
%   - tabela completa dos 26 candidatos
%   - verificacao das mascaras
%   - rank(H) e PDOP por arquitetura
%   - rank(H) e PDOP do conjunto completo
%   - diagnostico de pares geometricamente redundantes
%   - CSV e MAT para a proxima etapa

close all;
clear;
clc;
format long;

%% ##################### Arquivos ###########################################

tleLEO = "scenario4_LEO8_eval_2026-08-01T120000Z_combined.tle";
tleMEO = "scenario4_MEO7_eval_2026-08-01T120000Z_combined.tle";
tleGEO = "scenario4_GEO7_v2_eval_2026-08-01T120000Z_combined.tle";

filesRequired = [tleLEO,tleMEO,tleGEO];

for k = 1:numel(filesRequired)
    if ~isfile(filesRequired(k))
        error("Arquivo nao encontrado: %s",filesRequired(k));
    end
end

outDir = "results_scenario4_candidate_pool";

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ##################### Epoch / usuario ####################################

analysisTimeUTC = datetime(2026,8,1,12,0,0,'TimeZone','UTC');

userLat_deg = -22.8596582;
userLon_deg = -43.2303236;
userH_m     = 10;

UserPositionECEF = llh2ecef(userLat_deg,userLon_deg,userH_m);

%% ##################### Mascaras ###########################################

maskHAPS_deg = 15;
maskLEO_deg  = 5;
maskMEO_deg  = 5;
maskGEO_deg  = 5;

%% ##################### 1. HAPS4 sinteticos ###############################

hapsAltitude_m = 20000;

hapsAzDesign_deg = [45;135;225;315];
hapsOffset_km = [15;25;35;45];
hapsOffset_m = hapsOffset_km*1e3;

a = 6378137.0;
f = 1/298.257223563;
e2 = f*(2-f);

lat0 = deg2rad(userLat_deg);
lon0 = deg2rad(userLon_deg);

N0 = a/sqrt(1-e2*sin(lat0)^2);
M0 = a*(1-e2)/(1-e2*sin(lat0)^2)^(3/2);

east_m  = hapsOffset_m .* sind(hapsAzDesign_deg);
north_m = hapsOffset_m .* cosd(hapsAzDesign_deg);

hapsLat_deg = rad2deg(lat0 + north_m/M0);
hapsLon_deg = rad2deg(lon0 + east_m/(N0*cos(lat0)));

posHAPS = zeros(4,3);

for k = 1:4
    posHAPS(k,:) = llh2ecef( ...
        hapsLat_deg(k),hapsLon_deg(k),hapsAltitude_m);
end

%% ##################### 2. LEO / MEO / GEO SGP4 ###########################

[posLEO,nameLEO] = propagateTLEAtEpoch(tleLEO,analysisTimeUTC);
[posMEO,nameMEO] = propagateTLEAtEpoch(tleMEO,analysisTimeUTC);
[posGEO,nameGEO] = propagateTLEAtEpoch(tleGEO,analysisTimeUTC);

if size(posLEO,1) ~= 8
    error('Esperados 8 LEO, obtidos %d.',size(posLEO,1));
end

if size(posMEO,1) ~= 7
    error('Esperados 7 MEO, obtidos %d.',size(posMEO,1));
end

if size(posGEO,1) ~= 7
    error('Esperados 7 GEO, obtidos %d.',size(posGEO,1));
end

%% ##################### 3. Labels e metadados #############################

labelHAPS = ["HAPS1";"HAPS2";"HAPS3";"HAPS4"];
labelLEO  = compose("LEO%d",(1:8).');
labelMEO  = compose("MEO%d",(1:7).');
labelGEO  = compose("GEO%d",(1:7).');

archHAPS = repmat("HAPS",4,1);
archLEO  = repmat("LEO",8,1);
archMEO  = repmat("MEO",7,1);
archGEO  = repmat("GEO",7,1);

nameHAPS = labelHAPS;

% NORADs congelados
noradHAPS = repmat("",4,1);
noradLEO  = ["59340";"59960";"66267";"65388";"64047";"58736";"53155";"66290"];
noradMEO  = ["62362";"40081";"44115";"43231";"43234";"40080";"40079"];
noradGEO  = ["41589";"41904";"42692";"43562";"43175";"43228";"38087"];

%% ##################### 4. Pool final ######################################

PositionECEF = [
    posHAPS
    posLEO
    posMEO
    posGEO
];

Label = [
    labelHAPS
    labelLEO
    labelMEO
    labelGEO
];

Architecture = [
    archHAPS
    archLEO
    archMEO
    archGEO
];

ObjectName = [
    nameHAPS
    nameLEO
    nameMEO
    nameGEO
];

NORAD = [
    noradHAPS
    noradLEO
    noradMEO
    noradGEO
];

Mask_deg = [
    repmat(maskHAPS_deg,4,1)
    repmat(maskLEO_deg,8,1)
    repmat(maskMEO_deg,7,1)
    repmat(maskGEO_deg,7,1)
];

Npool = size(PositionECEF,1);

if Npool ~= 26
    error('O conjunto final deveria conter 26 candidatos, mas contem %d.',Npool);
end

%% ##################### 5. Geometria observada #############################

[Azimuth_deg,Elevation_deg,SlantRange_m,uLOS] = ecefAzElRange( ...
    UserPositionECEF,PositionECEF,userLat_deg,userLon_deg);

Visible = Elevation_deg >= Mask_deg;

PoolIndex = (1:Npool).';

Tpool = table( ...
    PoolIndex,Label,Architecture,NORAD,ObjectName, ...
    PositionECEF(:,1),PositionECEF(:,2),PositionECEF(:,3), ...
    Azimuth_deg,Elevation_deg,SlantRange_m/1e3,Mask_deg,Visible, ...
    'VariableNames',{ ...
    'PoolIndex','Label','Architecture','NORAD_CAT_ID','ObjectName', ...
    'X_ECEF_m','Y_ECEF_m','Z_ECEF_m', ...
    'Azimuth_deg','Elevation_deg','SlantRange_km','Mask_deg','Visible'});

fprintf('\n============================================================\n');
fprintf('SCENARIO 4 - CONJUNTO CANDIDATO FINAL\n');
fprintf('Epoch: %s\n',string(analysisTimeUTC));
fprintf('Total: %d candidatos\n',Npool);
fprintf('HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
    sum(Architecture=="HAPS"),sum(Architecture=="LEO"), ...
    sum(Architecture=="MEO"),sum(Architecture=="GEO"));
fprintf('============================================================\n');

disp(' ');
disp('===== Pool final =====');
disp(Tpool);

if all(Visible)
    fprintf('\nSUCESSO: todos os 26 candidatos satisfazem suas mascaras.\n');
else
    warning('%d candidatos nao satisfazem a mascara.',sum(~Visible));
    disp(Tpool(~Visible,:));
end

%% ##################### 6. Diagnostico por arquitetura #####################

archList = ["HAPS","LEO","MEO","GEO","ALL"];
nArch = numel(archList);

N = zeros(nArch,1);
RankH = zeros(nArch,1);
PDOP = nan(nArch,1);
RcondHtH = nan(nArch,1);

for aidx = 1:nArch

    if archList(aidx) == "ALL"
        idx = true(Npool,1);
    else
        idx = Architecture == archList(aidx);
    end

    [~,PDOP(aidx),RankH(aidx),RcondHtH(aidx)] = ...
        geometryDiagnostics(UserPositionECEF,PositionECEF(idx,:));

    N(aidx) = sum(idx);
end

Tsummary = table(archList.',N,RankH,RcondHtH,PDOP, ...
    'VariableNames',{'Architecture','N','RankH','Rcond_HtH','PDOP'});

disp(' ');
disp('===== Resumo geometrico =====');
disp(Tsummary);

%% ##################### 7. Redundancia angular ############################
% Angulo entre linhas de visada. Pares com separacao < 0.1 deg sao
% reportados para diagnostico, mas NAO sao removidos.

pairI = [];
pairJ = [];
pairAngle_deg = [];

for i = 1:Npool-1
    for j = i+1:Npool
        c = dot(uLOS(i,:),uLOS(j,:));
        c = max(-1,min(1,c));
        ang = acosd(c);

        if ang < 0.1
            pairI(end+1,1) = i; %#ok<SAGROW>
            pairJ(end+1,1) = j; %#ok<SAGROW>
            pairAngle_deg(end+1,1) = ang; %#ok<SAGROW>
        end
    end
end

if isempty(pairI)

    fprintf('\nNenhum par com separacao angular < 0.1 deg.\n');

    Tredundant = table();

else

    Tredundant = table( ...
        pairI,pairJ,Label(pairI),Label(pairJ), ...
        Architecture(pairI),Architecture(pairJ),pairAngle_deg, ...
        'VariableNames',{ ...
        'PoolIndex1','PoolIndex2','Label1','Label2', ...
        'Architecture1','Architecture2','AngularSeparation_deg'});

    disp(' ');
    disp('===== Pares com separacao angular < 0.1 deg =====');
    disp(Tredundant);
end

%% ##################### 8. Salvar ##########################################

writetable(Tpool,fullfile(outDir,'scenario4_candidate_pool_26.csv'));
writetable(Tsummary,fullfile(outDir,'scenario4_geometry_summary.csv'));

if ~isempty(Tredundant)
    writetable(Tredundant,fullfile(outDir,'scenario4_redundant_pairs.csv'));
end

save(fullfile(outDir,'scenario4_candidate_pool_26.mat'), ...
    'analysisTimeUTC', ...
    'userLat_deg','userLon_deg','userH_m','UserPositionECEF', ...
    'PositionECEF','Label','Architecture','NORAD','ObjectName', ...
    'Mask_deg','Azimuth_deg','Elevation_deg','SlantRange_m','Visible','uLOS', ...
    'Tpool','Tsummary','Tredundant');

fprintf('\nArquivos salvos em: %s\n',outDir);
fprintf('Proxima etapa: modelo de C/N0, sigma_rho, R e limite posicional.\n');

%% ##################### FUNCOES LOCAIS #####################################

function [posECEF,name] = propagateTLEAtEpoch(tleFile,epochUTC)

    sc = satelliteScenario(epochUTC,epochUTC + seconds(1),1);
    sat = satellite(sc,tleFile,'OrbitPropagator','sgp4');

    [posRaw,~] = states(sat,epochUTC,'CoordinateFrame','ecef');

    posECEF = reshape(posRaw,3,[]).';

    % Nao usar {sat.Name}: em algumas versoes do MATLAB, sat.Name
    % ja retorna uma colecao de nomes e o cell wrapping cria um elemento
    % nao escalar, causando erro na conversao para string.
    %
    % Os arquivos sao 3LE (nome + linha 1 + linha 2),
    % le os nomes diretamente do arquivo TLE.
    name = read3LENames(tleFile);

    if numel(name) ~= size(posECEF,1)
        error(['Numero de nomes (%d) diferente do numero de satelites ' ...
               'propagados (%d) em %s.'], ...
               numel(name),size(posECEF,1),tleFile);
    end
end

function name = read3LENames(tleFile)

    txt = fileread(tleFile);
    lines = splitlines(string(txt));
    lines = lines(strlength(strtrim(lines)) > 0);

    if mod(numel(lines),3) ~= 0
        error('Arquivo %s nao esta no formato 3LE esperado.',tleFile);
    end

    name = strtrim(lines(1:3:end));
    name = name(:);
end

function ecef = llh2ecef(lat_deg,lon_deg,h_m)

    a = 6378137.0;
    f = 1/298.257223563;
    e2 = f*(2-f);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    N = a/sqrt(1-e2*sin(lat)^2);

    ecef = [ ...
        (N+h_m)*cos(lat)*cos(lon), ...
        (N+h_m)*cos(lat)*sin(lon), ...
        (N*(1-e2)+h_m)*sin(lat)];
end

function [azDeg,elDeg,slantRange_m,uLOS] = ecefAzElRange( ...
    UserPosition,TxPositionECEF,lat_deg,lon_deg)

    UserPosition = UserPosition(:).';

    r = TxPositionECEF - UserPosition;
    slantRange_m = sqrt(sum(r.^2,2));
    uLOS = r./slantRange_m;

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    e_hat = [-sin(lon), cos(lon), 0];
    n_hat = [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)];
    u_hat = [ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)];

    east  = uLOS*e_hat.';
    north = uLOS*n_hat.';
    up    = uLOS*u_hat.';

    up = max(-1,min(1,up));

    elDeg = rad2deg(asin(up));
    azDeg = mod(rad2deg(atan2(east,north)),360);
end

function [H,PDOP,rankH,rcondHtH] = geometryDiagnostics( ...
    UserPosition,TxPositionECEF)

    UserPosition = UserPosition(:).';

    r = TxPositionECEF - UserPosition;
    d = sqrt(sum(r.^2,2));
    u = r./d;

    H = [-u,ones(size(u,1),1)];

    rankH = rank(H);
    A = H.'*H;
    rcondHtH = rcond(A);

    if rankH == 4 && rcondHtH > 1e-14
        C = A\eye(4);
        PDOP = sqrt(trace(C(1:3,1:3)));
    else
        PDOP = NaN;
    end
end
