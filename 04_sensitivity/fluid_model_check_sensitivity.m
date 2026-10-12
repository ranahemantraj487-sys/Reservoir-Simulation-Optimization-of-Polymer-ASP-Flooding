%% polymer_concentration_sensitivity_03.m
% Cairn-Rajasthan inspired polymer flood - polymer concentration sensitivity (v03)
%


clearvars; close all; clc;

mrstModule add ad-core ad-props ad-blackoil ad-eor;

%% 0. SETTINGS
runMode   = 'quick';      % 'quick' -> [0 1.00]  |  'full' -> [0 0.25 0.50 0.75 1.00]
nYears    = 10;           % simulation length; raise (15, 20, 30) until control water cut rises
screening = false;        % true -> 182.5 d steps after year 1 (fast, coarser); false -> 91.25 d
tightTol  = false;        % true -> tighter solver tolerances (slower), test balance error
expectedOOIP = NaN;       % paste the OOIP printed by fluid_model_03.m to guard the baseline (NaN = skip)

qInjection = 50*meter^3/day;     % per injector
nInjectors = 4;

%% 1. PATHS
projectRoot  = 'D:\PROJECT';
baselineFile = fullfile(projectRoot, '01_Baseline', 'fluid_model_03.m');
resultDir    = fullfile(projectRoot, '06_results');
figureDir    = fullfile(projectRoot, '07_Figures', 'v03');

if ~isfile(baselineFile), error('Baseline file not found:\n%s', baselineFile); end
if ~isfolder(resultDir),  mkdir(resultDir);  end
if ~isfolder(figureDir),  mkdir(figureDir);  end

%% 2. CASES AND TIME STEPS
switch runMode
    case 'quick', polymerConcentrations = [0, 1.00];
    case 'full',  polymerConcentrations = [0, 0.25, 0.50, 0.75, 1.00];   % kg/m3
    otherwise,    error('runMode must be ''quick'' or ''full''.');
end
nCases = numel(polymerConcentrations);

dtMainDays = 91.25;
if screening, dtMainDays = 182.5; end

dtA = [repmat(1*day,10,1); repmat(5*day,4,1)];               % first 30 days
dtB = repmat((365-30)*day/12, 12, 1);                         % rest of year 1
if nYears > 1
    nQ  = max(1, round((nYears-1)*365/dtMainDays));
    dtC = repmat((nYears-1)*365*day/nQ, nQ, 1);
else
    dtC = [];
end
dt = [dtA; dtB; dtC];
nSteps = numel(dt);
simulationDays = sum(convertTo(dt, day));

fprintf('\n============================================================\n');
fprintf(' POLYMER CONCENTRATION SENSITIVITY (v03) - mode %s\n', runMode);
fprintf('============================================================\n');
fprintf('Cases: %d | Simulation: %.0f days (%.1f y) | Steps: %d\n', ...
    nCases, simulationDays, simulationDays/365, nSteps);

%% 3. LOAD BASELINE
base = loadBaselineModel(baselineFile);
G = base.G;  rock = base.rock;  baseFluid = base.fluid;
Nx = G.cartDims(1); Ny = G.cartDims(2); Nz = G.cartDims(3);

if Nx ~= 50 || Ny ~= 50 || Nz ~= 3
    error('Baseline grid is %d x %d x %d; well positions in buildWells assume 50 x 50 x 3.', Nx, Ny, Nz);
end
if ~isfield(rock,'poro') || ~isfield(rock,'perm')
    error('Baseline rock must contain poro and perm.');
end
fprintf('Baseline grid check : PASS\n');

%% 4. POLYMER FLUID + MODEL
fluid = makePolymerFluid(baseFluid);
model = OilWaterPolymerModel(G, rock, fluid);

if tightTol
    if isprop(model,'toleranceCNV'), model.toleranceCNV = 1e-6; end
    if isprop(model,'toleranceMB'),  model.toleranceMB  = 1e-9; end
    fprintf('Tight solver tolerances requested.\n');
end

%% 5. OOIP / PORE VOLUME
stateInit = base.state0;
pv   = poreVolume(G, rock);
PV   = sum(pv);
OOIP = sum(pv .* stateInit.s(:,2));
fprintf('Pore volume %.4e m3 | OOIP %.2f m3\n', PV, OOIP);
if ~isnan(expectedOOIP)
    ooipErr = abs(OOIP - expectedOOIP)/expectedOOIP;
    fprintf('OOIP vs expected: rel.err %.3e\n', ooipErr);
    if ooipErr >= 1e-3
        error('OOIP mismatch too large - baseline changed? Sensitivity NOT run.');
    end
end

%% 6. WELLS
W = buildWells(G, rock, 150*barsa, qInjection);
qTotal = sum([W(2:end).val]);                       % m3/s
PVIend = qTotal*sum(dt)/PV;
fprintf('Wells: P1 (25,25) + 4 corner injectors, %.0f m3/day each | water injected: %.2f PV\n', ...
    convertTo(qInjection, meter^3/day), PVIend);

%% 7. RUN ALL CASES
allResults = cell(nCases,1);
partialFile = fullfile(resultDir, 'polymer_sensitivity_v03_partial.mat');

for caseIdx = 1:nCases
    cpInj = polymerConcentrations(caseIdx);
    fprintf('\n============================================================\n');
    fprintf(' CASE %d/%d : %.2f kg/m3\n', caseIdx, nCases, cpInj);
    fprintf('============================================================\n');

    allResults{caseIdx} = runCase(G, model, W, dt, stateInit, cpInj, nInjectors, OOIP);

    r = allResults{caseIdx};
    fprintf('\n--- FINAL RESULTS: %.2f kg/m3 ---\n', cpInj);
    fprintf('Cumulative oil     : %.2f m3\n',    r.cumOil);
    fprintf('Recovery factor    : %.4f %%\n',    r.RF);
    fprintf('Final oil rate     : %.4f m3/day\n', r.oilRate(end));
    fprintf('Final water cut    : %.2f %%\n',    100*r.waterCut(end));
    fprintf('Average pressure   : %.4f bar\n',   r.avgPressure(end));
    fprintf('Max injector BHP   : %.4f bar\n',   max(r.injectorBHP(:)));
    if cpInj > 0
        fprintf('--- POLYMER MASS BALANCE (kg) ---\n');
        fprintf('Injected                  : %.2f\n', r.bal.injected(end));
        fprintf('Produced                  : %.2f\n', r.bal.produced(end));
        fprintf('Stored in reservoir       : %.2f\n', r.bal.stored(end));
        fprintf('Balance error             : %.4f %% of injected\n', r.bal.errPct(end));
        fprintf('Produced implied by closure (injected - stored): %.2f kg\n', ...
            r.bal.injected(end) - r.bal.stored(end));
        if abs(r.bal.errPct(end)) > 2
            warning('Polymer balance error > 2%% in case %.2f kg/m3.', cpInj);
        end
    end

    save(partialFile, 'allResults', 'polymerConcentrations', 'OOIP', 'PV', ...
        'simulationDays', 'nSteps', '-v7.3');
end

%% 8. SUMMARY TABLE
pick = @(f) cellfun(@(r) r.(f), allResults);

cumulativeOil   = pick('cumOil');
recoveryFactor  = pick('RF');
finalOilRate    = cellfun(@(r) r.oilRate(end),     allResults);
finalWaterRate  = cellfun(@(r) r.waterRate(end),   allResults);
finalTotalRate  = cellfun(@(r) r.totalRate(end),   allResults);
finalWaterCut   = cellfun(@(r) r.waterCut(end),    allResults);
finalAvgPress   = cellfun(@(r) r.avgPressure(end), allResults);
finalProdBHP    = cellfun(@(r) r.producerBHP(end), allResults);
maxInjBHP       = cellfun(@(r) max(r.injectorBHP(:)), allResults);
finalAvgSw      = cellfun(@(r) r.avgSw(end),       allResults);
finalAvgSo      = cellfun(@(r) r.avgSo(end),       allResults);
finalAvgCp      = cellfun(@(r) r.avgCp(end),       allResults);
finalMaxCp      = cellfun(@(r) r.maxCp(end),       allResults);
polyInjected    = cellfun(@(r) r.bal.injected(end), allResults);
polyProduced    = cellfun(@(r) r.bal.produced(end), allResults);
polyStored      = cellfun(@(r) r.bal.stored(end),   allResults);
polyErrPct      = cellfun(@(r) r.bal.errPct(end),   allResults);

% incremental vs waterflood control (cp = 0)
iCtl   = find(polymerConcentrations == 0, 1);
incOil = nan(nCases,1);  incRF = nan(nCases,1);  PUF = nan(nCases,1);
if ~isempty(iCtl)
    incOil = cumulativeOil - cumulativeOil(iCtl);
    incRF  = recoveryFactor - recoveryFactor(iCtl);
    ok = incOil > 0;
    PUF(ok) = polyInjected(ok) ./ incOil(ok);       % kg polymer per incremental m3 oil
end

summaryTable = table( ...
    polymerConcentrations(:), cumulativeOil, recoveryFactor, incOil, incRF, PUF, ...
    finalOilRate, finalWaterRate, finalTotalRate, 100*finalWaterCut, finalAvgPress, ...
    finalProdBHP, maxInjBHP, finalAvgSw, finalAvgSo, finalAvgCp, finalMaxCp, ...
    polyInjected, polyProduced, polyStored, polyErrPct, ...
    'VariableNames', { ...
    'Polymer_kg_m3','CumulativeOil_m3','RecoveryFactor_pct', ...
    'IncOilVsControl_m3','IncRF_vsControl_pct','PUF_kg_per_incm3', ...
    'FinalOilRate_m3_day','FinalWaterRate_m3_day','FinalTotalRate_m3_day', ...
    'FinalWaterCut_pct','FinalAvgPressure_bar','ProducerBHP_bar', ...
    'MaxInjectorBHP_bar','FinalAvgSw','FinalAvgSo','FinalAvgCp_kg_m3', ...
    'FinalMaxCp_kg_m3','PolymerInjected_kg','PolymerProduced_kg', ...
    'PolymerStored_kg','PolymerBalanceError_pct'});

fprintf('\n============================================================\n');
fprintf(' FINAL SENSITIVITY SUMMARY\n');
fprintf('============================================================\n');
disp(summaryTable);
fprintf('Water injected: %.2f PV | produced polymer shows whether the bank reached the producer.\n', PVIend);

%% 9. SAVE
matFile = fullfile(resultDir, 'polymer_concentration_sensitivity_v03.mat');
csvFile = fullfile(resultDir, 'polymer_concentration_sensitivity_v03.csv');
save(matFile, 'summaryTable','allResults','polymerConcentrations', ...
    'OOIP','PV','simulationDays','nSteps','-v7.3');
writetable(summaryTable, csvFile);

%% 10. PLOTS
C  = polymerConcentrations;
xl = 'Polymer concentration (kg/m^3)';
plotSensitivity(C, cumulativeOil,     xl, 'Cumulative oil (m^3)',        'Cumulative Oil',       fullfile(figureDir,'17_Sens_Cumulative_Oil.png'));
plotSensitivity(C, recoveryFactor,    xl, 'Recovery factor (%)',         'Recovery Factor',      fullfile(figureDir,'18_Sens_Recovery.png'));
plotSensitivity(C, finalOilRate,      xl, 'Final oil rate (m^3/day)',    'Final Oil Rate',       fullfile(figureDir,'19_Sens_Final_Oil_Rate.png'));
plotSensitivity(C, 100*finalWaterCut, xl, 'Final water cut (%)',         'Water Cut',            fullfile(figureDir,'20_Sens_Water_Cut.png'));
plotSensitivity(C, finalAvgPress,     xl, 'Final average pressure (bar)','Average Pressure',     fullfile(figureDir,'21_Sens_Pressure.png'));
plotSensitivity(C, maxInjBHP,         xl, 'Max injector BHP (bar)',      'Maximum Injector BHP', fullfile(figureDir,'22_Sens_Max_Inj_BHP.png'));
plotSensitivity(C, finalAvgCp,        xl, 'Final average Cp (kg/m^3)',   'Average Cp',           fullfile(figureDir,'23_Sens_Avg_Cp.png'));
if ~isempty(iCtl)
    plotSensitivity(C, incOil, xl, 'Incremental oil vs waterflood (m^3)', 'Incremental Oil', fullfile(figureDir,'34_Sens_Incremental_Oil.png'));
end

fig = figure('Color','w'); hold on;
for k = 1:nCases
    plot(allResults{k}.timeDays, allResults{k}.oilRate, 'LineWidth',1.5, ...
        'DisplayName', sprintf('%.2f kg/m^3', C(k)));
end
grid on; xlabel('Time (days)'); ylabel('Oil rate (m^3/day)');
title('Oil Rate History'); legend('Location','best');
saveas(fig, fullfile(figureDir,'24_Oil_Rate_History.png'));

fig = figure('Color','w'); hold on;
for k = 1:nCases
    plot(allResults{k}.timeDays, 100*allResults{k}.waterCut, 'LineWidth',1.5, ...
        'DisplayName', sprintf('%.2f kg/m^3', C(k)));
end
grid on; xlabel('Time (days)'); ylabel('Water cut (%)');
title('Water Cut History'); legend('Location','best');
saveas(fig, fullfile(figureDir,'35_Water_Cut_History.png'));

fig = figure('Color','w'); hold on;
for k = 1:nCases
    plot(allResults{k}.timeDays, allResults{k}.cumOilHist/1000, 'LineWidth',1.5, ...
        'DisplayName', sprintf('%.2f kg/m^3', C(k)));
end
grid on; xlabel('Time (days)'); ylabel('Cumulative oil (10^3 m^3)');
title('Cumulative Oil History'); legend('Location','best');
saveas(fig, fullfile(figureDir,'36_Cumulative_Oil_History.png'));

fig = figure('Color','w');
plot(maxInjBHP, recoveryFactor, '-o', 'LineWidth',1.5, 'MarkerSize',7); grid on;
xlabel('Maximum injector BHP (bar)'); ylabel('Recovery factor (%)');
title('Recovery Factor vs Maximum Injector BHP');
for k = 1:nCases
    text(maxInjBHP(k), recoveryFactor(k), sprintf('  %.2f kg/m^3', C(k)), 'VerticalAlignment','bottom');
end
saveas(fig, fullfile(figureDir,'25_Recovery_vs_Max_Inj_BHP.png'));

fig = figure('Color','w'); hold on;
for k = 1:nCases
    if C(k) == 0, continue; end
    b = allResults{k}.bal;
    plot(allResults{k}.timeDays, b.injected/1000, '-',  'LineWidth',1.5, 'DisplayName', sprintf('Injected %.2f', C(k)));
    plot(allResults{k}.timeDays, b.stored/1000,   '--', 'LineWidth',1.5, 'DisplayName', sprintf('Stored %.2f',   C(k)));
    plot(allResults{k}.timeDays, b.produced/1000, ':',  'LineWidth',1.5, 'DisplayName', sprintf('Produced %.2f', C(k)));
end
grid on; xlabel('Time (days)'); ylabel('Polymer (tonnes)');
title('Polymer Injected / Stored / Produced'); legend('Location','northwest');
saveas(fig, fullfile(figureDir,'26_Polymer_Injected_vs_Stored.png'));

fprintf('\nSENSITIVITY STUDY COMPLETE\nMAT: %s\nCSV: %s\nFigures: %s\n', matFile, csvFile, figureDir);


%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function base = loadBaselineModel(baselineFile)
% Load the exact baseline model (no new random fields generated).
    run(baselineFile);
    required = {'G','rock','fluid','state0','p0'};
    for k = 1:numel(required)
        if ~exist(required{k}, 'var')
            error('Baseline file did not create variable "%s".', required{k});
        end
    end
    base.G = G;  base.rock = rock;  base.fluid = fluid;
    base.state0 = state0;  base.p0 = p0;
    fprintf('Baseline model loaded successfully.\n');
end


function fluid = makePolymerFluid(baseFluid)
% Add polymer properties to the EXACT baseline fluid.
    fluid = baseFluid;

    cpmax = 1.0;
    fluid.cpmax  = cpmax;
    fluid.mixPar = 0.85;

    % viscosity multiplier (keep in sync with muWMultFun in fluid_model_03.m)
    fluid.muWMult = @(c) 1 + 4*(c./cpmax);

    % no adsorption (0*c keeps AD derivatives valid)
    fluid.adsInx = 1;
    fluid.adsMax = 1.0;
    fluid.ads    = @(c) 0*c;

    fluid.rrf  = 1.0;     % no residual resistance factor
    fluid.dps  = 0.0;     % no inaccessible pore volume
    fluid.rhoR = 2650*kilogram/meter^3;
end


function W = buildWells(G, rock, pProducer, qInjection)
% 5-spot: producer in centre, 4 injectors in corners, all layers.
    nz = G.cartDims(3);

    prodCells = sub2ind(G.cartDims, 25*ones(1,nz), 25*ones(1,nz), 1:nz);
    W = addWell([], G, rock, prodCells, ...
        'Type','bhp', 'Val',pProducer, 'Radius',0.1*meter, ...
        'Comp_i',[0 1], 'Name','P1');

    injIJ = [4 4; 4 47; 47 4; 47 47];
    names = {'I1','I2','I3','I4'};
    for k = 1:4
        cells = sub2ind(G.cartDims, injIJ(k,1)*ones(1,nz), injIJ(k,2)*ones(1,nz), 1:nz);
        W = addWell(W, G, rock, cells, ...
            'Type','rate', 'Val',qInjection, 'Radius',0.1*meter, ...
            'Comp_i',[1 0], 'Name',names{k});
    end

    W = setPolymerConc(W, 1.0);
end


function W = setPolymerConc(W, cpInj)
% Producer: cp = 0, injectors: cp = cpInj. Producer is W(1).
    for n = 1:numel(W)
        if n == 1
            W(n).compi = [0 1];
            W(n).cp    = 0.0;
        else
            W(n).compi = [1 0];
            W(n).cp    = cpInj;
        end
    end
end


function r = runCase(G, model, W, dt, stateInit, cpInj, nInj, OOIP)
% Run one polymer concentration case and return all histories + mass balance.

    nSteps = numel(dt);

    % ---- initial state (polymer fields are REQUIRED) ----
    state0 = stateInit;
    state0.cp    = zeros(G.cells.num, 1);
    state0.cpmax = zeros(G.cells.num, 1);

    % ---- wells / schedule (control must be a STRUCT ARRAY, not a cell) ----
    Wcase = setPolymerConc(W, cpInj);
    schedule.step.val     = dt;
    schedule.step.control = ones(nSteps, 1);
    schedule.control      = struct('W', Wcase);

    fprintf('Starting simulation...\n');
    [wellSols, states] = simulateScheduleAD(state0, model, schedule, 'Verbose', false);
    fprintf('Simulation completed.\n');

    % ---- histories ----
    timeDays    = cumsum(convertTo(dt, day));
    oilRate     = zeros(nSteps,1);  waterRate = zeros(nSteps,1);
    totalRate   = zeros(nSteps,1);  waterCut  = zeros(nSteps,1);
    cumOilHist  = zeros(nSteps,1);  cumWaterHist = zeros(nSteps,1);
    avgPressure = zeros(nSteps,1);
    avgSw = zeros(nSteps,1);  avgSo = zeros(nSteps,1);
    avgCp = zeros(nSteps,1);  maxCp = zeros(nSteps,1);
    producerBHP = zeros(nSteps,1);
    injectorBHP = zeros(nSteps,nInj);

    cumOil = 0;  cumWater = 0;
    for t = 1:nSteps
        ws = wellSols{t};
        st = states{t};

        % production rates are negative in MRST: convert first, then flip sign
        qOil   = max(-convertTo(ws(1).qOs, meter^3/day), 0);
        qWater = max(-convertTo(ws(1).qWs, meter^3/day), 0);

        oilRate(t)   = qOil;
        waterRate(t) = qWater;
        totalRate(t) = qOil + qWater;
        if totalRate(t) > 0, waterCut(t) = qWater/totalRate(t); end

        dtDays = convertTo(dt(t), day);
        cumOil   = cumOil   + qOil*dtDays;
        cumWater = cumWater + qWater*dtDays;
        cumOilHist(t)   = cumOil;
        cumWaterHist(t) = cumWater;

        avgPressure(t) = convertTo(mean(st.pressure), barsa);
        avgSw(t) = mean(st.s(:,1));
        avgSo(t) = mean(st.s(:,2));

        if isfield(st,'cp')
            avgCp(t) = mean(double(st.cp(:)));
            maxCp(t) = max(double(st.cp(:)));
        end

        producerBHP(t) = convertTo(ws(1).bhp, barsa);
        for n = 1:nInj
            injectorBHP(t,n) = convertTo(ws(n+1).bhp, barsa);
        end

        if t == 1 || mod(t,10) == 0 || t == nSteps
            fprintf('Step %3d/%3d | t %7.1f d | Oil %9.3f | Water %9.3f | WC %6.2f%% | P %8.2f bar | Cpmax %6.3f\n', ...
                t, nSteps, timeDays(t), oilRate(t), waterRate(t), ...
                100*waterCut(t), avgPressure(t), maxCp(t));
        end
    end

    % ---- polymer mass balance ----
    bal = polymerBalance(model, states, wellSols, dt, Wcase, cpInj, nInj);

    % ---- pack ----
    r.polymerConcentration = cpInj;
    r.timeDays    = timeDays;
    r.oilRate     = oilRate;      r.waterRate = waterRate;
    r.totalRate   = totalRate;    r.waterCut  = waterCut;
    r.cumOilHist  = cumOilHist;   r.cumWaterHist = cumWaterHist;
    r.cumOil      = cumOil;       r.cumWater  = cumWater;
    r.RF          = 100*cumOil/OOIP;
    r.avgPressure = avgPressure;
    r.avgSw = avgSw;  r.avgSo = avgSo;
    r.avgCp = avgCp;  r.maxCp = maxCp;
    r.producerBHP = producerBHP;
    r.injectorBHP = injectorBHP;
    r.bal         = bal;
    r.finalState  = states{end};
    r.finalWellSol = wellSols{end};
    r.wellSols     = wellSols;      % kept for post-hoc checks (fields, per-perforation rates)
end


function bal = polymerBalance(model, states, wellSols, dt, Wcase, cpInj, nInj) %#ok<INUSD>
% Polymer inventory using the same quantities the simulator uses.
%   stored   = sum( PV(p) * bW(p) * Sw * (1-dps) * cp )   (adsorption = 0 here)
%   injected = integral of injector polymer rate (qWPoly, or qWs*cp)
%   produced = integral of producer polymer rate
% NOTE: if qWPoly is not in wellSol, 'injected' falls back to qWs*cp and is
% therefore not independent of the nominal q*c*t.

    f      = model.fluid;
    nSteps = numel(states);
    pv0    = model.operators.pv;

    dps = 0;
    if isfield(f,'dps'), dps = f.dps; end

    stored   = zeros(nSteps,1);
    injected = zeros(nSteps,1);
    produced = zeros(nSteps,1);
    nominal  = zeros(nSteps,1);

    cumInj = 0;  cumProd = 0;  cumNom = 0;
    usedMean = false;
    nominalRate = 0;
    for n = 2:numel(Wcase)
        nominalRate = nominalRate + Wcase(n).val;   % m3/s (SI)
    end

    for t = 1:nSteps
        st = states{t};
        ws = wellSols{t};
        p  = st.pressure;

        bW = f.bW(p);
        if isfield(f,'pvMultR')
            pv = pv0 .* f.pvMultR(p);
        else
            pv = pv0;
        end
        sW = st.s(:,1);

        stored(t) = sum(pv .* bW .* sW .* (1 - dps) .* st.cp);

        % ---- well polymer rates [kg/s] ----
        hasPoly = isfield(ws,'qWPoly') && ~any(cellfun(@isempty, {ws.qWPoly}));
        if hasPoly
            q    = [ws.qWPoly];
            qInj = sum(q(q > 0));
            qPrd = -sum(q(q < 0));
        else
            qInj = 0;  qPrd = 0;
            for n = 1:numel(ws)
                qw = ws(n).qWs;                       % m3/s, + injection, - production
                if qw > 0
                    qInj = qInj + qw * Wcase(n).cp;
                else
                    % flow-weighted: sum over perforations of (water rate * cell cp).
                    % With layered permeability the arithmetic mean of cp is WRONG
                    % (the high-perm layer carries most of the water and the polymer).
                    cpc   = double(st.cp(Wcase(n).cells));
                    qPerf = [];
                    if isfield(ws,'cqs') && ~isempty(ws(n).cqs)
                        qPerf = ws(n).cqs(:,1);       % water surface rate per perforation [m3/s]
                    end
                    if numel(qPerf) == numel(cpc) && abs(sum(qPerf) - qw) <= 1e-3*abs(qw)
                        qPrd = qPrd - sum(qPerf .* cpc(:));
                    else
                        qPrd = qPrd - qw * mean(cpc);  % fallback (less accurate)
                        usedMean = true;
                    end
                end
            end
        end

        cumInj  = cumInj  + qInj*dt(t);
        cumProd = cumProd + qPrd*dt(t);
        cumNom  = cumNom  + nominalRate*cpInj*dt(t);

        injected(t) = cumInj;
        produced(t) = cumProd;
        nominal(t)  = cumNom;
    end

    if usedMean
        warning(['Producer polymer rate used the arithmetic mean of cp (wellSol.cqs not usable). ' ...
                 'Produced polymer may be underestimated with layered permeability.']);
    end

    bal.stored           = stored;
    bal.injected         = injected;
    bal.produced         = produced;
    bal.injectedNominal  = nominal;
    bal.errPct           = 100*(stored + produced - injected) ./ max(injected, eps);
end


function plotSensitivity(x, y, xLabelText, yLabelText, titleText, fileName)
    fig = figure('Color','w');
    plot(x, y, '-o', 'LineWidth',1.6, 'MarkerSize',7);
    grid on;
    xlabel(xLabelText); ylabel(yLabelText);
    title(['Polymer Concentration Sensitivity - ' titleText]);
    saveas(fig, fileName);
end
