%% polymerflood_08_10year.m
% Cairn Rajasthan-Inspired Polymer EOR Study
%
% 10-year polymer-flood simulation using MRST OilWaterPolymerModel.
%
% IMPORTANT:
% Polymer is transported as a concentration in the aqueous phase.
% It is NOT treated as a third fluid phase.
%
% Current base-case assumptions:
%   - 50 x 50 x 3 layered heterogeneous reservoir (fluid_model_03.m)
%   - 5-spot pattern
%   - 4 injectors + 1 producer
%   - 200 m3/day total injection
%   - 1.0 kg/m3 polymer concentration
%   - No adsorption
%   - No inaccessible pore-space / RRF effect
%   - 10-year simulation
%
% MASS BALANCE:
%   Injected = Produced + Stored + Numerical Residual
%
% For the current fluid model:
%   adsorption = 0
%   dps       = 0
%
% Therefore reservoir polymer storage is calculated from:
%
%   Mstored = sum(phi .* pore-volume .* Sw .* Cp)
%
% where Cp is polymer concentration in the aqueous phase.
%
% -------------------------------------------------------------------------

clear;
clc;
close all;

%% ============================================================
% 1. PROJECT PATHS
% =============================================================

projectRoot = 'D:\PROJECT';

baselineFile = fullfile( ...
    projectRoot, '01_Baseline', 'fluid_model_03.m');

resultsFolder = fullfile( ...
    projectRoot, '06_results');

figuresFolder = fullfile( ...
    projectRoot, '07_Figures', 'v03');

if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

if ~exist(figuresFolder, 'dir')
    mkdir(figuresFolder);
end

%% ============================================================
% 2. LOAD EXISTING BASELINE RESERVOIR
% ============================================================

run(baselineFile);

close all;   % fluid_model_03.m opens plot windows

% Well positions below assume the 50 x 50 x 3 grid of fluid_model_03.m
if G.cartDims(1) ~= 50 || G.cartDims(2) ~= 50 || G.cartDims(3) ~= 3
    error('Baseline grid is %d x %d x %d; well positions assume 50 x 50 x 3.', G.cartDims);
end

% Restore paths in case baseline script changes workspace variables.
projectRoot = 'D:\PROJECT';

resultsFolder = fullfile( ...
    projectRoot, '06_results');

figuresFolder = fullfile( ...
    projectRoot, '07_Figures', 'v03');

%% ============================================================
% 3. LOAD MRST MODULES
% ============================================================

mrstModule add ad-core ad-props ad-eor;

%% ============================================================
% 4. MODEL INFORMATION
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  OPTIMIZED 10-YEAR POLYMER FLOOD - MRST\n');
fprintf('============================================================\n');

fprintf('Grid: %d x %d x %d = %d cells\n', ...
    G.cartDims(1), ...
    G.cartDims(2), ...
    G.cartDims(3), ...
    G.cells.num);

fprintf('Reservoir dimensions: %.0f x %.0f x %.0f m\n', ...
    max(G.nodes.coords(:,1)) - min(G.nodes.coords(:,1)), ...
    max(G.nodes.coords(:,2)) - min(G.nodes.coords(:,2)), ...
    max(G.nodes.coords(:,3)) - min(G.nodes.coords(:,3)));
fprintf('Initial pressure: %.1f bar\n', convertTo(mean(state0.pressure), barsa));
fprintf('Initial Sw: %.2f\n', mean(state0.s(:,1)));
fprintf('Initial So: %.2f\n', mean(state0.s(:,2)));

%% ============================================================
% 5. INITIAL RESERVOIR STATE
% ============================================================

% Initial pressure and saturations come from the baseline model
% (fluid_model_03.m); state0 is NOT rebuilt here.
pInit  = mean(state0.pressure);
SwInit = mean(state0.s(:,1));
SoInit = mean(state0.s(:,2));

% Initially there is no polymer in the reservoir.
state0.cp    = zeros(G.cells.num, 1);
state0.cpmax = zeros(G.cells.num, 1);

%% ============================================================
% 6. POLYMER FLUID MODEL
% ============================================================

fluid.cpmax = 1.0;

fluid.mixPar = 0.85;

% Polymer viscosity relationship:
%
% Cp = 0 kg/m3  -> muW multiplier = 1
% Cp = 1 kg/m3  -> muW multiplier = 5

fluid.muWMult = @(c) ...
    1 + 4*(c/fluid.cpmax);

% ------------------------------------------------------------
% Adsorption
% ------------------------------------------------------------

fluid.adsInx = 1;
fluid.adsMax = 1.0;

% Zero adsorption for the base case.
fluid.ads = @(c) zeros(size(c));

% ------------------------------------------------------------
% Residual resistance / inaccessible pore space
% ------------------------------------------------------------

fluid.rrf = 1.0;
fluid.dps = 0.0;

% ------------------------------------------------------------
% Rock density
% ------------------------------------------------------------

fluid.rhoR = 2650*kilogram/meter^3;

%% ============================================================
% 7. POLYMER MODEL
% ============================================================

modelPolymer = OilWaterPolymerModel( ...
    G, rock, fluid);

modelPolymer.verbose = false;

%% ============================================================
% 8. WELL CONFIGURATION
% ============================================================

nInjectors = 4;

qInjection = 50*meter^3/day;

pProducer = 150*barsa;

producerIJ = [25, 25];

injectorIJ = [ ...
     4,  4;
     4, 47;
    47,  4;
    47, 47];

W = [];

nz = G.cartDims(3);

%% ------------------------------------------------------------
% Producer P1
% ------------------------------------------------------------

producerCells = zeros(nz,1);

for k = 1:nz

    producerCells(k) = ...
        sub2ind( ...
            G.cartDims, ...
            producerIJ(1), ...
            producerIJ(2), ...
            k);

end

W = addWell( ...
    W, ...
    G, ...
    rock, ...
    producerCells, ...
    'Type', 'bhp', ...
    'Val', pProducer, ...
    'Radius', 0.1*meter, ...
    'Name', 'P1', ...
    'Comp_i', [0 1]);

%% ------------------------------------------------------------
% Injectors
% ------------------------------------------------------------

for n = 1:nInjectors

    injCells = zeros(nz,1);

    for k = 1:nz

        injCells(k) = ...
            sub2ind( ...
                G.cartDims, ...
                injectorIJ(n,1), ...
                injectorIJ(n,2), ...
                k);

    end

    W = addWell( ...
        W, ...
        G, ...
        rock, ...
        injCells, ...
        'Type', 'rate', ...
        'Val', qInjection, ...
        'Radius', 0.1*meter, ...
        'Name', sprintf('I%d',n), ...
        'Comp_i', [1 0]);

end

%% ============================================================
% 9. WELL OIL-WATER COMPOSITION
% ============================================================

% P1 = producer
W(1).compi = [0 1];

% I1-I4 = water injectors
for n = 2:numel(W)

    W(n).compi = [1 0];

end

%% ============================================================
% 10. POLYMER CONCENTRATION AT WELLS
% ============================================================

% Producer:
% no polymer injected at producer.
W(1).cp = 0.0;

% Injectors:
% 1 kg/m3 polymer.
for n = 2:numel(W)

    W(n).cp = 1.0;

end

%% ============================================================
% 11. VERIFY WELL POLYMER MAPPING
% ============================================================

fprintf('\n');
fprintf('Well polymer concentrations:\n');

for n = 1:numel(W)

    cpCheck = ...
        modelPolymer.getProp(W(n), 'polymer');

    fprintf( ...
        '  %-3s : %.3f kg/m3\n', ...
        W(n).name, ...
        cpCheck);

end

%% ============================================================
% 12. SIMULATION SCHEDULE
% ============================================================

simulationDays = 3650;

%% First 30 days: daily

dt1 = repmat(1*day, 30, 1);

%% Remaining Year 1: approximately monthly

dt2 = [ ...
    repmat(28*day, 11, 1);
    27*day];

assert( ...
    abs(convertTo(sum(dt1)+sum(dt2),day)-365) < 1e-8, ...
    'First-year timestep construction is incorrect.');

%% Years 2-10: approximately quarterly

dt3 = [ ...
    repmat(91*day, 36, 1);
    9*day];

%% Complete schedule

dt = [dt1; dt2; dt3];

totalScheduleDays = ...
    convertTo(sum(dt),day);

assert( ...
    abs(totalScheduleDays-simulationDays) < 1e-8, ...
    '10-year timestep construction is incorrect.');

schedule = simpleSchedule( ...
    dt, ...
    'W', W);

nSteps = numel(dt);

fprintf('\n');
fprintf('Simulation schedule:\n');
fprintf('  Total simulation time : %.0f days\n', ...
    totalScheduleDays);

fprintf('  Control timesteps     : %d\n', ...
    nSteps);

fprintf('  First 30 days         : daily\n');
fprintf('  Rest of Year 1        : ~monthly\n');
fprintf('  Years 2-10            : ~quarterly\n');

%% ============================================================
% 13. GRID GEOMETRY AND PORE VOLUME
% ============================================================

if ~isfield(G.cells,'volumes')

    G = computeGeometry(G);

end

cellVolumes = G.cells.volumes(:);

poro = rock.poro(:);

if numel(poro) ~= G.cells.num

    error( ...
        'rock.poro must contain one value per grid cell.');

end

%% ============================================================
% 14. OOIP
% ============================================================

poreVolume = sum( ...
    poro .* cellVolumes);

OOIP = ...
    sum(poro .* cellVolumes .* state0.s(:,2));

totalInjectionRate_m3day = ...
    nInjectors * ...
    convertTo(qInjection,meter^3/day);

polymerInjectionConcentration = ...
    1.0;

plannedPolymerInjection = ...
    totalInjectionRate_m3day * ...
    polymerInjectionConcentration * ...
    simulationDays;

fprintf('\n');
fprintf('OOIP = %.2f m3\n',OOIP);

fprintf( ...
    'Total planned polymer injection = %.2f kg\n', ...
    plannedPolymerInjection);

%% ============================================================
% 15. INITIAL CONDITIONS
% ============================================================

state = state0;

initialAvgPressure_bar = ...
    convertTo(mean(state.pressure),barsa);

initialAvgSw = ...
    SwInit;

initialAvgSo = ...
    SoInit;

fprintf( ...
    'Initial average pressure = %.4f bar\n', ...
    initialAvgPressure_bar);

%% ============================================================
% 16. PREALLOCATE HISTORY
% ============================================================

timeDays = ...
    cumsum(convertTo(dt,day));

dtDays = ...
    convertTo(dt,day);

oilRateHistory = zeros(nSteps,1);

waterRateHistory = zeros(nSteps,1);

totalRateHistory = zeros(nSteps,1);

waterCutHistory = zeros(nSteps,1);

pressureHistory = zeros(nSteps,1);

avgSwHistory = zeros(nSteps,1);

avgSoHistory = zeros(nSteps,1);

avgCpHistory = zeros(nSteps,1);

maxCpHistory = zeros(nSteps,1);

producerBHPHistory = zeros(nSteps,1);

injectorBHPHistory = ...
    zeros(nSteps,nInjectors);

cumulativeOil = zeros(nSteps,1);

cumulativeWater = zeros(nSteps,1);

%% Polymer accounting

polymerInjected = zeros(nSteps,1);

polymerProduced = zeros(nSteps,1);

polymerStored = zeros(nSteps,1);

polymerResidual = zeros(nSteps,1);

polymerBalanceErrorPercent = ...
    zeros(nSteps,1);

polymerBalanceClosurePercent = ...
    zeros(nSteps,1);

recoveryFactor = zeros(nSteps,1);

%% ============================================================
% 17. START SIMULATION
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  STARTING 10-YEAR POLYMER SIMULATION\n');
fprintf('============================================================\n');

fprintf( ...
    'This version uses %d control steps.\n\n', ...
    nSteps);

tic;

[wellSols, states, reports] = ...
    simulateScheduleAD( ...
        state, ...
        modelPolymer, ...
        schedule, ...
        'OutputMinisteps', false);

elapsedTime = toc;

fprintf('\n');
fprintf('============================================================\n');
fprintf('  SIMULATION COMPLETE\n');
fprintf('============================================================\n');

fprintf( ...
    'Elapsed wall-clock time: %.2f minutes\n', ...
    elapsedTime/60);

fprintf( ...
    'Control steps solved: %d\n', ...
    nSteps);

%% ============================================================
% 18. PROCESS SIMULATION RESULTS
% ============================================================

for t = 1:nSteps

    currentState = states{t};

    ws = wellSols{t};

    %% --------------------------------------------------------
    % Producer rates
    % ---------------------------------------------------------

    qOil = ...
        max( ...
            -convertTo( ...
                ws(1).qOs, ...
                meter^3/day), ...
            0);

    qWater = ...
        max( ...
            -convertTo( ...
                ws(1).qWs, ...
                meter^3/day), ...
            0);

    oilRateHistory(t) = qOil;

    waterRateHistory(t) = qWater;

    totalRateHistory(t) = ...
        qOil + qWater;

    if totalRateHistory(t) > 0

        waterCutHistory(t) = ...
            qWater / totalRateHistory(t);

    else

        waterCutHistory(t) = 0;

    end

    %% --------------------------------------------------------
    % BHP
    % ---------------------------------------------------------

    producerBHPHistory(t) = ...
        convertTo(ws(1).bhp,barsa);

    for j = 1:nInjectors

        injectorBHPHistory(t,j) = ...
            convertTo( ...
                ws(j+1).bhp, ...
                barsa);

    end

    %% --------------------------------------------------------
    % Average pressure
    % ---------------------------------------------------------

    pressureHistory(t) = ...
        convertTo( ...
            mean(currentState.pressure), ...
            barsa);

    %% --------------------------------------------------------
    % Saturations
    % ---------------------------------------------------------

    sNow = currentState.s;

    if size(sNow,2) >= 2

        SwNow = sNow(:,1);

        SoNow = sNow(:,2);

    else

        error( ...
            'Unexpected saturation structure at timestep %d.', ...
            t);

    end

    SwNow = SwNow(:);

    SoNow = SoNow(:);

    %% --------------------------------------------------------
    % Average saturations
    % ---------------------------------------------------------

    avgSwHistory(t) = ...
        sum(SwNow .* cellVolumes) / ...
        sum(cellVolumes);

    avgSoHistory(t) = ...
        sum(SoNow .* cellVolumes) / ...
        sum(cellVolumes);

    %% --------------------------------------------------------
    % Polymer concentration
    % ---------------------------------------------------------

    cpNow = currentState.cp(:);

    avgCpHistory(t) = ...
        sum(cpNow .* cellVolumes) / ...
        sum(cellVolumes);

    maxCpHistory(t) = ...
        max(cpNow);

    %% ========================================================
    % CUMULATIVE OIL/WATER
    % ========================================================

    if t == 1

        % We assume zero production at t = 0.
        cumulativeOil(t) = ...
            0.5 * qOil * dtDays(t);

        cumulativeWater(t) = ...
            0.5 * qWater * dtDays(t);

    else

        cumulativeOil(t) = ...
            cumulativeOil(t-1) + ...
            0.5 * ...
            (oilRateHistory(t-1)+qOil) * ...
            dtDays(t);

        cumulativeWater(t) = ...
            cumulativeWater(t-1) + ...
            0.5 * ...
            (waterRateHistory(t-1)+qWater) * ...
            dtDays(t);

    end

    %% ========================================================
    % POLYMER INJECTION MASS
    % ========================================================

    % Constant total injection:
    %
    % 4 injectors x 50 m3/day x 1 kg/m3
    % = 200 kg/day

    polymerInjectionRate = ...
        totalInjectionRate_m3day * ...
        polymerInjectionConcentration;

    if t == 1

        polymerInjected(t) = ...
            polymerInjectionRate * ...
            dtDays(t);

    else

        polymerInjected(t) = ...
            polymerInjected(t-1) + ...
            polymerInjectionRate * ...
            dtDays(t);

    end

    %% ========================================================
    % POLYMER PRODUCTION
    % ========================================================

    % Polymer produced with the water at P1.
    %
    % q_water [m3/day] x Cp [kg/m3]
    % = polymer production rate [kg/day]

    cpProducer = ...
        modelPolymer.getProp( ...
            ws(1), ...
            'polymer');

    cpProducer = ...
        double(cpProducer);

    cpProducer = ...
        cpProducer(1);

    polymerProducedRate = ...
        qWater * cpProducer;

    if t == 1

        polymerProduced(t) = ...
            0.5 * ...
            polymerProducedRate * ...
            dtDays(t);

    else

        previousWellSol = ...
            wellSols{t-1};

        cpProducerPrevious = ...
            modelPolymer.getProp( ...
                previousWellSol(1), ...
                'polymer');

        cpProducerPrevious = ...
            double(cpProducerPrevious);

        cpProducerPrevious = ...
            cpProducerPrevious(1);

        previousWaterRate = ...
            waterRateHistory(t-1);

        previousPolymerRate = ...
            previousWaterRate * ...
            cpProducerPrevious;

        polymerProduced(t) = ...
            polymerProduced(t-1) + ...
            0.5 * ...
            (previousPolymerRate + ...
             polymerProducedRate) * ...
            dtDays(t);

    end

    %% ========================================================
    % POLYMER STORAGE IN RESERVOIR
    % ========================================================

    % Current model assumptions:
    %
    % adsorption = zero
    % dps = zero
    %
    % Therefore the polymer inventory in the reservoir aqueous
    % phase is:
    %
    % Mstored = sum(phi * V * Sw * Cp)
    %
    % Cp is kg/m3 of aqueous phase.
    %
    % phi * V * Sw = water volume in the pore space.

    polymerStored(t) = ...
        sum( ...
            poro .* ...
            cellVolumes .* ...
            SwNow .* ...
            cpNow);

    %% ========================================================
    % POLYMER MASS BALANCE
    % ========================================================

    % Physical balance:
    %
    % Injected = Produced + Stored + Residual
    %
    % Therefore:
    %
    % Residual = Injected - Produced - Stored

    polymerResidual(t) = ...
        polymerInjected(t) - ...
        polymerProduced(t) - ...
        polymerStored(t);

    if polymerInjected(t) > 0

        polymerBalanceClosurePercent(t) = ...
            100 * ...
            (polymerProduced(t) + ...
             polymerStored(t)) / ...
            polymerInjected(t);

        polymerBalanceErrorPercent(t) = ...
            100 * ...
            abs(polymerResidual(t)) / ...
            polymerInjected(t);

    else

        polymerBalanceClosurePercent(t) = 0;

        polymerBalanceErrorPercent(t) = 0;

    end

    %% ========================================================
    % RECOVERY FACTOR
    % ========================================================

    recoveryFactor(t) = ...
        100 * ...
        cumulativeOil(t) / ...
        OOIP;

    %% ========================================================
    % PROGRESS OUTPUT
    % ========================================================

    if t == 1 || ...
       mod(t,10) == 0 || ...
       t == nSteps

        fprintf( ...
            ['Step %3d/%3d | Time %7.1f d | ' ...
             'Oil %9.3f | Water %9.3f | ' ...
             'WC %6.2f%% | Pavg %7.2f bar | ' ...
             'Cp,max %6.3f\n'], ...
            t, ...
            nSteps, ...
            timeDays(t), ...
            oilRateHistory(t), ...
            waterRateHistory(t), ...
            100*waterCutHistory(t), ...
            pressureHistory(t), ...
            maxCpHistory(t));

    end

end

%% ============================================================
% 19. FINAL STATE
% ============================================================

finalState = states{end};

finalWellSol = wellSols{end};

%% ============================================================
% 20. FINAL PRODUCTION RESULTS
% ============================================================

finalOilRate = ...
    oilRateHistory(end);

finalWaterRate = ...
    waterRateHistory(end);

finalTotalRate = ...
    totalRateHistory(end);

finalWaterCut = ...
    waterCutHistory(end);

finalCumulativeOil = ...
    cumulativeOil(end);

finalCumulativeWater = ...
    cumulativeWater(end);

finalRecovery = ...
    recoveryFactor(end);

%% ============================================================
% 21. FINAL PRESSURE
% ============================================================

finalAvgPressure = ...
    pressureHistory(end);

pressureChange = ...
    finalAvgPressure - ...
    initialAvgPressure_bar;

finalProducerBHP = ...
    producerBHPHistory(end);

finalInjectorBHP = ...
    injectorBHPHistory(end,:);

maxInjectorBHP = ...
    max(finalInjectorBHP);

%% ============================================================
% 22. FINAL SATURATION
% ============================================================

finalSw = ...
    finalState.s(:,1);

finalSo = ...
    finalState.s(:,2);

finalAvgSw = ...
    avgSwHistory(end);

finalAvgSo = ...
    avgSoHistory(end);

%% ============================================================
% 23. FINAL POLYMER
% ============================================================

finalCp = ...
    finalState.cp(:);

finalAvgCp = ...
    avgCpHistory(end);

finalMaxCp = ...
    maxCpHistory(end);

finalPolymerInjected = ...
    polymerInjected(end);

finalPolymerProduced = ...
    polymerProduced(end);

finalPolymerStored = ...
    polymerStored(end);

finalPolymerResidual = ...
    polymerResidual(end);

finalPolymerClosure = ...
    polymerBalanceClosurePercent(end);

finalPolymerBalanceError = ...
    polymerBalanceErrorPercent(end);

%% ============================================================
% 24. FINAL RESULTS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  10-YEAR POLYMER FLOOD RESULTS\n');
fprintf('============================================================\n');

fprintf( ...
    'Simulation time              : %.0f days (%.2f years)\n', ...
    timeDays(end), ...
    timeDays(end)/365);

fprintf('\n');
fprintf('--- Production ---\n');

fprintf( ...
    'Final oil rate               : %.4f m3/day\n', ...
    finalOilRate);

fprintf( ...
    'Final water rate             : %.4f m3/day\n', ...
    finalWaterRate);

fprintf( ...
    'Final total liquid rate      : %.4f m3/day\n', ...
    finalTotalRate);

fprintf( ...
    'Final water cut              : %.2f %%\n', ...
    100*finalWaterCut);

fprintf( ...
    'Cumulative oil               : %.2f m3\n', ...
    finalCumulativeOil);

fprintf( ...
    'Cumulative water             : %.2f m3\n', ...
    finalCumulativeWater);

fprintf( ...
    'Recovery factor              : %.4f %%\n', ...
    finalRecovery);

fprintf('\n');
fprintf('--- Pressure ---\n');

fprintf( ...
    'Initial average pressure     : %.4f bar\n', ...
    initialAvgPressure_bar);

fprintf( ...
    'Final average pressure       : %.4f bar\n', ...
    finalAvgPressure);

fprintf( ...
    'Pressure change              : %+0.4f bar\n', ...
    pressureChange);

fprintf( ...
    'Final producer BHP           : %.2f bar\n', ...
    finalProducerBHP);

fprintf('\n');
fprintf('--- Saturation ---\n');

fprintf( ...
    'Initial average Sw           : %.4f\n', ...
    initialAvgSw);

fprintf( ...
    'Final average Sw             : %.4f\n', ...
    finalAvgSw);

fprintf( ...
    'Initial average So           : %.4f\n', ...
    initialAvgSo);

fprintf( ...
    'Final average So             : %.4f\n', ...
    finalAvgSo);

fprintf( ...
    'Final Sw range               : %.4f to %.4f\n', ...
    min(finalSw), ...
    max(finalSw));

fprintf( ...
    'Final So range               : %.4f to %.4f\n', ...
    min(finalSo), ...
    max(finalSo));

fprintf('\n');
fprintf('--- Polymer Mass Balance ---\n');

fprintf( ...
    'Polymer injected             : %.2f kg\n', ...
    finalPolymerInjected);

fprintf( ...
    'Polymer produced             : %.2f kg\n', ...
    finalPolymerProduced);

fprintf( ...
    'Polymer stored in reservoir  : %.2f kg\n', ...
    finalPolymerStored);

fprintf( ...
    'Mass-balance residual        : %.2f kg\n', ...
    finalPolymerResidual);

fprintf( ...
    'Mass-balance closure         : %.4f %%\n', ...
    finalPolymerClosure);

fprintf( ...
    'Mass-balance error           : %.4f %%\n', ...
    finalPolymerBalanceError);

fprintf('\n');
fprintf('--- Polymer Concentration ---\n');

fprintf( ...
    'Final average Cp             : %.6f kg/m3\n', ...
    finalAvgCp);

fprintf( ...
    'Final maximum Cp             : %.6f kg/m3\n', ...
    finalMaxCp);

fprintf('\n');
fprintf('--- Injector BHP ---\n');

for n = 1:nInjectors

    fprintf( ...
        'I%d BHP                      : %.2f bar\n', ...
        n, ...
        finalInjectorBHP(n));

end

fprintf( ...
    'Maximum injector BHP         : %.2f bar\n', ...
    maxInjectorBHP);

fprintf('\n');
fprintf('--- Reference fluid values ---\n');

fprintf('Oil viscosity                : 10 cP\n');
fprintf('Water viscosity              : 1 cP\n');
fprintf('Oil density                  : 800 kg/m3\n');
fprintf('Water density                : 1000 kg/m3\n');

%% ============================================================
% 25. VALIDATION CHECKS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  VALIDATION CHECKS\n');
fprintf('============================================================\n');

%% Polymer balance

if finalPolymerBalanceError <= 1.0

    fprintf( ...
        'PASS: Polymer balance error is <= 1%%.\n');

else

    fprintf( ...
        'CHECK: Polymer balance error is > 1%%.\n');

end

%% Explicitly warn if stored polymer exceeds injected polymer

if finalPolymerStored > ...
   finalPolymerInjected + 1e-8

    fprintf( ...
        ['WARNING: Calculated polymer storage exceeds injected ' ...
         'polymer mass by %.2f kg (%.4f%%).\n'], ...
        finalPolymerStored-finalPolymerInjected, ...
        100*(finalPolymerStored-finalPolymerInjected)/ ...
        finalPolymerInjected);

    fprintf( ...
        ['The simulation is numerically close, but the polymer ' ...
         'inventory should be treated as a diagnostic rather ' ...
         'than a strict conservation proof.\n']);

else

    fprintf( ...
        'PASS: Stored polymer does not exceed injected polymer.\n');

end

%% Maximum concentration

if finalMaxCp <= ...
   fluid.cpmax + 1e-8

    fprintf( ...
        'PASS: Maximum Cp does not exceed Cpmax.\n');

else

    fprintf( ...
        'CHECK: Maximum Cp exceeds specified Cpmax.\n');

end

%% Saturation

if abs( ...
       finalAvgSw + ...
       finalAvgSo - 1) < 1e-6

    fprintf( ...
        'PASS: Average oil + water saturation = 1.\n');

else

    fprintf( ...
        'CHECK: Average saturation sum differs from 1.\n');

end

fprintf('\n');
fprintf( ...
    'NOTE: Injector BHP is reported for engineering screening only.\n');

fprintf( ...
    ['A safety/fracture-limit conclusion requires a specified ' ...
     'allowable injection pressure or fracture-pressure criterion.\n']);

%% ============================================================
% 26. RESULTS STRUCTURE
% ============================================================

resultsPolymer10Year = struct();

resultsPolymer10Year.projectName = ...
    'Cairn Rajasthan-Inspired Polymer EOR Study';

resultsPolymer10Year.grid = ...
    G.cartDims;

resultsPolymer10Year.numCells = ...
    G.cells.num;

resultsPolymer10Year.simulationDays = ...
    simulationDays;

resultsPolymer10Year.numControlSteps = ...
    nSteps;

resultsPolymer10Year.elapsedSeconds = ...
    elapsedTime;

%% OOIP

resultsPolymer10Year.OOIP_m3 = ...
    OOIP;

%% Time

resultsPolymer10Year.timeDays = ...
    timeDays;

resultsPolymer10Year.dtDays = ...
    dtDays;

%% Production

resultsPolymer10Year.oilRate_m3day = ...
    oilRateHistory;

resultsPolymer10Year.waterRate_m3day = ...
    waterRateHistory;

resultsPolymer10Year.totalRate_m3day = ...
    totalRateHistory;

resultsPolymer10Year.waterCut = ...
    waterCutHistory;

resultsPolymer10Year.cumulativeOil_m3 = ...
    cumulativeOil;

resultsPolymer10Year.cumulativeWater_m3 = ...
    cumulativeWater;

resultsPolymer10Year.recoveryFactor = ...
    recoveryFactor;

%% Pressure

resultsPolymer10Year.pressure_bar = ...
    pressureHistory;

resultsPolymer10Year.producerBHP_bar = ...
    producerBHPHistory;

resultsPolymer10Year.injectorBHP_bar = ...
    injectorBHPHistory;

%% Saturation

resultsPolymer10Year.avgSw = ...
    avgSwHistory;

resultsPolymer10Year.avgSo = ...
    avgSoHistory;

%% Polymer concentration

resultsPolymer10Year.avgCp_kgm3 = ...
    avgCpHistory;

resultsPolymer10Year.maxCp_kgm3 = ...
    maxCpHistory;

%% Polymer mass balance

resultsPolymer10Year.polymerInjected_kg = ...
    polymerInjected;

resultsPolymer10Year.polymerProduced_kg = ...
    polymerProduced;

resultsPolymer10Year.polymerStored_kg = ...
    polymerStored;

resultsPolymer10Year.polymerResidual_kg = ...
    polymerResidual;

resultsPolymer10Year.polymerClosure_percent = ...
    polymerBalanceClosurePercent;

resultsPolymer10Year.polymerBalanceError_percent = ...
    polymerBalanceErrorPercent;

%% Parameters

resultsPolymer10Year.parameters = struct( ...
    'nInjectors', nInjectors, ...
    'injectionRatePerWell_m3day', ...
        convertTo(qInjection,meter^3/day), ...
    'totalInjectionRate_m3day', ...
        totalInjectionRate_m3day, ...
    'polymerInjectionConcentration_kgm3', ...
        polymerInjectionConcentration, ...
    'cpmax_kgm3', ...
        fluid.cpmax, ...
    'mixPar', ...
        fluid.mixPar, ...
    'muWMultiplierAtCpmax', ...
        fluid.muWMult(fluid.cpmax), ...
    'adsorption', ...
        'zero', ...
    'rrf', ...
        fluid.rrf, ...
    'dps', ...
        fluid.dps);

%% ============================================================
% 27. SAVE RESULTS
% ============================================================

resultsFile = ...
    fullfile( ...
        resultsFolder, ...
        'polymer_10year_results_v03.mat');

save( ...
    resultsFile, ...
    'resultsPolymer10Year', ...
    '-v7.3');

fprintf('\n');
fprintf('Results saved to:\n');
fprintf('%s\n',resultsFile);

%% ============================================================
% 28. PRODUCTION PLOT
% ============================================================

figure( ...
    'Name', ...
    '10-Year Polymer Flood - Production', ...
    'Color', ...
    'w');

plot( ...
    timeDays, ...
    oilRateHistory, ...
    'LineWidth',1.5);

hold on;

plot( ...
    timeDays, ...
    waterRateHistory, ...
    'LineWidth',1.5);

plot( ...
    timeDays, ...
    totalRateHistory, ...
    'LineWidth',1.5);

grid on;

xlabel('Time (days)');
ylabel('Rate (m^3/day)');

title( ...
    '10-Year Polymer Flood Production');

legend( ...
    'Oil', ...
    'Water', ...
    'Total Liquid', ...
    'Location','best');

saveas( ...
    gcf, ...
    fullfile( ...
        figuresFolder, ...
        'polymer_10year_production.png'));

%% ============================================================
% 29. WATER CUT AND PRESSURE
% ============================================================

figure( ...
    'Name', ...
    '10-Year Polymer Flood - Water Cut and Pressure', ...
    'Color', ...
    'w');

yyaxis left;

plot( ...
    timeDays, ...
    100*waterCutHistory, ...
    'LineWidth',1.5);

ylabel('Water Cut (%)');

yyaxis right;

plot( ...
    timeDays, ...
    pressureHistory, ...
    'LineWidth',1.5);

ylabel('Average Pressure (bar)');

grid on;

xlabel('Time (days)');

title( ...
    '10-Year Polymer Flood: Water Cut and Average Pressure');

saveas( ...
    gcf, ...
    fullfile( ...
        figuresFolder, ...
        'polymer_10year_watercut_pressure.png'));

%% ============================================================
% 30. POLYMER CONCENTRATION
% ============================================================

figure( ...
    'Name', ...
    '10-Year Polymer Flood - Polymer Concentration', ...
    'Color', ...
    'w');

plot( ...
    timeDays, ...
    avgCpHistory, ...
    'LineWidth',1.5);

hold on;

plot( ...
    timeDays, ...
    maxCpHistory, ...
    'LineWidth',1.5);

grid on;

xlabel('Time (days)');

ylabel('Polymer concentration (kg/m^3)');

title( ...
    'Polymer Concentration History');

legend( ...
    'Average Cp', ...
    'Maximum Cp', ...
    'Location','best');

saveas( ...
    gcf, ...
    fullfile( ...
        figuresFolder, ...
        'polymer_10year_polymer_concentration.png'));

%% ============================================================
% 31. INJECTOR BHP
% ============================================================

figure( ...
    'Name', ...
    '10-Year Polymer Flood - Injector BHP', ...
    'Color', ...
    'w');

plot( ...
    timeDays, ...
    injectorBHPHistory, ...
    'LineWidth',1.3);

grid on;

xlabel('Time (days)');

ylabel('Injector BHP (bar)');

title( ...
    '10-Year Polymer Flood Injector BHP');

legend( ...
    'I1','I2','I3','I4', ...
    'Location','best');

saveas( ...
    gcf, ...
    fullfile( ...
        figuresFolder, ...
        'polymer_10year_injector_bhp.png'));

%% ============================================================
% 32. CUMULATIVE OIL AND RECOVERY
% ============================================================

figure( ...
    'Name', ...
    '10-Year Polymer Flood - Recovery', ...
    'Color', ...
    'w');

yyaxis left;

plot( ...
    timeDays, ...
    cumulativeOil, ...
    'LineWidth',1.5);

ylabel('Cumulative oil (m^3)');

yyaxis right;

plot( ...
    timeDays, ...
    recoveryFactor, ...
    'LineWidth',1.5);

ylabel('Recovery factor (%)');

grid on;

xlabel('Time (days)');

title( ...
    '10-Year Polymer Flood Recovery');

saveas( ...
    gcf, ...
    fullfile( ...
        figuresFolder, ...
        'polymer_10year_recovery.png'));

%% ============================================================
% 33. POLYMER MASS BALANCE PLOT
% ============================================================

figure( ...
    'Name', ...
    '10-Year Polymer Flood - Mass Balance', ...
    'Color', ...
    'w');

plot( ...
    timeDays, ...
    polymerInjected, ...
    'LineWidth',1.5);

hold on;

plot( ...
    timeDays, ...
    polymerStored, ...
    'LineWidth',1.5);

plot( ...
    timeDays, ...
    polymerProduced, ...
    'LineWidth',1.5);

grid on;

xlabel('Time (days)');

ylabel('Polymer mass (kg)');

title( ...
    'Polymer Mass Accounting');

legend( ...
    'Injected', ...
    'Stored', ...
    'Produced', ...
    'Location','best');

saveas( ...
    gcf, ...
    fullfile( ...
        figuresFolder, ...
        'polymer_10year_mass_balance.png'));

%% ============================================================
% 34. FINAL SUMMARY
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  FINAL SUMMARY\n');
fprintf('============================================================\n');

fprintf( ...
    'Grid                         : %d cells\n', ...
    G.cells.num);

fprintf( ...
    'Control timesteps            : %d\n', ...
    nSteps);

fprintf( ...
    'Simulation time              : %.2f years\n', ...
    simulationDays/365);

fprintf( ...
    'Cumulative oil               : %.2f m3\n', ...
    finalCumulativeOil);

fprintf( ...
    'Recovery factor              : %.4f %%\n', ...
    finalRecovery);

fprintf( ...
    'Final oil rate               : %.4f m3/day\n', ...
    finalOilRate);

fprintf( ...
    'Final water cut              : %.2f %%\n', ...
    100*finalWaterCut);

fprintf( ...
    'Final average pressure       : %.2f bar\n', ...
    finalAvgPressure);

fprintf( ...
    'Maximum injector BHP         : %.2f bar\n', ...
    maxInjectorBHP);

fprintf( ...
    'Polymer injected             : %.2f kg\n', ...
    finalPolymerInjected);

fprintf( ...
    'Polymer stored               : %.2f kg\n', ...
    finalPolymerStored);

fprintf( ...
    'Polymer produced             : %.2f kg\n', ...
    finalPolymerProduced);

fprintf( ...
    'Mass-balance residual        : %.2f kg\n', ...
    finalPolymerResidual);

fprintf( ...
    'Mass-balance closure         : %.4f %%\n', ...
    finalPolymerClosure);

fprintf( ...
    'Mass-balance error           : %.4f %%\n', ...
    finalPolymerBalanceError);

fprintf( ...
    'Elapsed simulation time      : %.2f min\n', ...
    elapsedTime/60);

fprintf( ...
    '============================================================\n');