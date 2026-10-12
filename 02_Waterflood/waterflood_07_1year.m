%% ============================================================
%  1-YEAR WATERFLOOD SIMULATION
%  Cairn Rajasthan Inspired Reservoir Study
%
%  Purpose:
%  Baseline waterflood case for comparison with the
%  1-year polymer flood.
%
%  Reservoir:
%  50 x 50 x 3 grid
%  400 m x 400 m x 30 m  (as defined in fluid_model_03.m)
%
%  Well pattern:
%  5-spot
%  Producer : P1 at center
%  Injectors: I1-I4 at corners
%
%  Injection:
%  50 m3/day per injector
%  200 m3/day total
%
%  Simulation:
%  365 daily timesteps
%
%  ============================================================

clear;
clc;
close all;

fprintf('\n');
fprintf('============================================================\n');
fprintf('              1-YEAR WATERFLOOD SIMULATION\n');
fprintf('============================================================\n\n');

%% ============================================================
%  1. PROJECT PATH
%  ============================================================

projectRoot = 'D:\PROJECT';

baselineFile = fullfile(projectRoot, ...
    '01_Baseline', 'fluid_model_03.m');

resultsFolder = fullfile(projectRoot, '06_results');

if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% ============================================================
%  2. LOAD BASELINE RESERVOIR MODEL
%  ============================================================

fprintf('Loading baseline reservoir model...\n');

run(baselineFile);

close all;   % fluid_model_03.m opens plot windows

% Well positions below assume the 50 x 50 x 3 grid of fluid_model_03.m
if G.cartDims(1) ~= 50 || G.cartDims(2) ~= 50 || G.cartDims(3) ~= 3
    error('Baseline grid is %d x %d x %d; well positions assume 50 x 50 x 3.', G.cartDims);
end

fprintf('Baseline model loaded successfully.\n\n');

%% ============================================================
%  3. RESERVOIR INFORMATION
%  ============================================================

fprintf('============================================================\n');
fprintf('                    RESERVOIR DATA\n');
fprintf('============================================================\n');

% Grid and size taken from the baseline model (fluid_model_03.m)
Nx = G.cartDims(1);
Ny = G.cartDims(2);
Nz = G.cartDims(3);

Lx = max(G.nodes.coords(:,1)) - min(G.nodes.coords(:,1));
Ly = max(G.nodes.coords(:,2)) - min(G.nodes.coords(:,2));
Lz = max(G.nodes.coords(:,3)) - min(G.nodes.coords(:,3));

fprintf('Grid                  : %d x %d x %d\n', Nx, Ny, Nz);
fprintf('Total cells           : %d\n', Nx*Ny*Nz);
fprintf('Reservoir dimensions  : %.0f m x %.0f m x %.0f m\n', ...
    Lx, Ly, Lz);

fprintf('Initial pressure      : %.1f bar\n', convertTo(mean(state0.pressure), barsa));
fprintf('Initial water sat.    : %.2f\n', mean(state0.s(:,1)));
fprintf('Initial oil sat.      : %.2f\n', mean(state0.s(:,2)));

fprintf('Oil viscosity         : %.1f cP\n', convertTo(muO, centi*poise));
fprintf('Water viscosity       : %.1f cP\n', convertTo(muW, centi*poise));
fprintf('Oil density           : %.0f kg/m3\n', convertTo(rhoO, kilogram/meter^3));
fprintf('Water density         : %.0f kg/m3\n', convertTo(rhoW, kilogram/meter^3));

%% ============================================================
%  4. INITIAL STATE
%  ============================================================

state = state0;

%% ============================================================
%  5. CALCULATE OOIP
%  ============================================================

% Pore volume
poreVolume = sum(G.cells.volumes .* rock.poro);

% Original Oil In Place from the baseline initial state (fluid_model_03.m)
OOIP = sum(G.cells.volumes .* rock.poro .* state0.s(:,2));

fprintf('\n');
fprintf('OOIP                  : %.2f m3\n', OOIP);

%% ============================================================
%  6. WELL CONFIGURATION
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                    WELL CONFIGURATION\n');
fprintf('============================================================\n');

% Producer location
iP = 25;
jP = 25;

% Injector locations
iI = [4 4 47 47];
jI = [4 47 4 47];

% Number of injectors
nInjectors = 4;

% Producer BHP
pProducer = 150 * barsa;

% Injection rate per injector
qInjection = 50 * meter^3/day;

% Total injection rate
totalInjectionRate = nInjectors * qInjection;

fprintf('Producer              : P1\n');
fprintf('Producer location     : (%d, %d)\n', iP, jP);
fprintf('Producer BHP          : %.2f bar\n', ...
    convertTo(pProducer, barsa));

fprintf('\nInjectors:\n');

for n = 1:nInjectors
    fprintf('I%d                    : (%d, %d)\n', ...
        n, iI(n), jI(n));
end

fprintf('\nInjection per well    : %.2f m3/day\n', ...
    convertTo(qInjection, meter^3/day));

fprintf('Total injection       : %.2f m3/day\n', ...
    convertTo(totalInjectionRate, meter^3/day));

%% ============================================================
%  7. CREATE WELLS
%  ============================================================

fprintf('\nCreating wells...\n');

W = [];

%% ------------------------------------------------------------
% Producer P1
% ------------------------------------------------------------

cellP = sub2ind([Nx Ny Nz], ...
    iP, jP, 1:Nz);

W = addWell(W, G, rock, cellP, ...
    'Type', 'bhp', ...
    'Val', pProducer, ...
    'Name', 'P1', ...
    'Comp_i', [0 1]);

%% ------------------------------------------------------------
% Injectors
% ------------------------------------------------------------

for n = 1:nInjectors

    cellI = sub2ind([Nx Ny Nz], ...
        iI(n), jI(n), 1:Nz);

    W = addWell(W, G, rock, cellI, ...
        'Type', 'rate', ...
        'Val', qInjection, ...
        'Name', sprintf('I%d', n), ...
        'Comp_i', [1 0]);

end

fprintf('Wells created successfully.\n');

%% ============================================================
%  8. SIMULATION TIME
%  ============================================================

simulationDays = 365;

% One timestep per day
dt = repmat(1*day, simulationDays, 1);

schedule = simpleSchedule(dt, 'W', W);

fprintf('\n');
fprintf('Simulation duration  : %.0f days\n', simulationDays);
fprintf('Number of timesteps  : %d\n', numel(dt));
fprintf('Timestep size        : 1 day\n');

%% ============================================================
%  9. RUN SIMULATION
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('              STARTING WATERFLOOD SIMULATION\n');
fprintf('============================================================\n\n');

tic;

[wellSols, states, report] = ...
    simulateScheduleAD(state, model, schedule);

simulationTime = toc;

fprintf('\nSimulation completed.\n');
fprintf('Simulation runtime    : %.2f seconds\n', simulationTime);

%% ============================================================
%  10. INITIAL PRESSURE AND SATURATION
% ============================================================

initialPressure = mean(state.pressure);

initialSw = mean(state.s(:,1));
initialSo = mean(state.s(:,2));

%% ============================================================
%  11. PREALLOCATE HISTORY ARRAYS
% ============================================================

nSteps = numel(wellSols);

timeDays = zeros(nSteps,1);

oilRateHistory = zeros(nSteps,1);
waterRateHistory = zeros(nSteps,1);
totalRateHistory = zeros(nSteps,1);
waterCutHistory = zeros(nSteps,1);

pressureHistory = zeros(nSteps,1);

avgSwHistory = zeros(nSteps,1);
avgSoHistory = zeros(nSteps,1);

producerBHPHistory = zeros(nSteps,1);

injectorBHPHistory = zeros(nSteps,nInjectors);

cumulativeOil = zeros(nSteps,1);
cumulativeWater = zeros(nSteps,1);

recoveryFactor = zeros(nSteps,1);

%% ============================================================
%  12. EXTRACT PRODUCTION HISTORY
%  ============================================================

fprintf('\nProcessing production history...\n');

for t = 1:nSteps

    %% Time

    timeDays(t) = t;

    %% --------------------------------------------------------
    % Producer
    % --------------------------------------------------------

    wsP = wellSols{t}(1);

    % Oil production rate
    oilRate = abs(convertTo(wsP.qOs, meter^3/day));

    % Water production rate
    waterRate = abs(convertTo(wsP.qWs, meter^3/day));

    % Total production rate
    totalRate = oilRate + waterRate;

    % Water cut
    if totalRate > 0
        waterCut = waterRate / totalRate;
    else
        waterCut = 0;
    end

    oilRateHistory(t) = oilRate;
    waterRateHistory(t) = waterRate;
    totalRateHistory(t) = totalRate;
    waterCutHistory(t) = waterCut * 100;

    %% --------------------------------------------------------
    % Producer BHP
    % --------------------------------------------------------

    if isfield(wsP, 'bhp')
        producerBHPHistory(t) = ...
            convertTo(wsP.bhp, barsa);
    else
        producerBHPHistory(t) = NaN;
    end

    %% --------------------------------------------------------
    % Reservoir pressure
    % --------------------------------------------------------

    pressureHistory(t) = ...
        convertTo(mean(states{t}.pressure), barsa);

    %% --------------------------------------------------------
    % Saturation
    % --------------------------------------------------------

    avgSwHistory(t) = mean(states{t}.s(:,1));
    avgSoHistory(t) = mean(states{t}.s(:,2));

    %% --------------------------------------------------------
    % Injector BHP
    % --------------------------------------------------------

    for n = 1:nInjectors

        wsI = wellSols{t}(n+1);

        if isfield(wsI, 'bhp')
            injectorBHPHistory(t,n) = ...
                convertTo(wsI.bhp, barsa);
        else
            injectorBHPHistory(t,n) = NaN;
        end

    end

    %% --------------------------------------------------------
    % Cumulative production
    % --------------------------------------------------------

    if t == 1

        % Same convention as the 1-year polymer case:
        cumulativeOil(t) = ...
            oilRateHistory(t) * timeDays(t);

        cumulativeWater(t) = ...
            waterRateHistory(t) * timeDays(t);

    else

        dtDays = timeDays(t) - timeDays(t-1);

        % Trapezoidal integration
        cumulativeOil(t) = cumulativeOil(t-1) + ...
            0.5 * ...
            (oilRateHistory(t-1) + oilRateHistory(t)) ...
            * dtDays;

        cumulativeWater(t) = cumulativeWater(t-1) + ...
            0.5 * ...
            (waterRateHistory(t-1) + waterRateHistory(t)) ...
            * dtDays;

    end

    %% --------------------------------------------------------
    % Recovery factor
    % --------------------------------------------------------

    recoveryFactor(t) = ...
        cumulativeOil(t) / OOIP * 100;

end

%% ============================================================
%  13. FINAL RESULTS
%  ============================================================

finalOilRate = oilRateHistory(end);
finalWaterRate = waterRateHistory(end);
finalTotalRate = totalRateHistory(end);
finalWaterCut = waterCutHistory(end);

finalPressure = pressureHistory(end);

finalSw = avgSwHistory(end);
finalSo = avgSoHistory(end);

finalProducerBHP = producerBHPHistory(end);

finalInjectorBHP = injectorBHPHistory(end,:);

maxInjectorBHP = max(injectorBHPHistory(:));

finalCumulativeOil = cumulativeOil(end);
finalCumulativeWater = cumulativeWater(end);

finalRecoveryFactor = recoveryFactor(end);

%% ============================================================
%  14. FINAL RESERVOIR STATISTICS
%  ============================================================

finalState = states{end};

finalSwMin = min(finalState.s(:,1));
finalSwMax = max(finalState.s(:,1));

finalSoMin = min(finalState.s(:,2));
finalSoMax = max(finalState.s(:,2));

finalPressureMin = ...
    convertTo(min(finalState.pressure), barsa);

finalPressureMax = ...
    convertTo(max(finalState.pressure), barsa);

finalPressureMean = ...
    convertTo(mean(finalState.pressure), barsa);

%% ============================================================
%  15. PRESSURE CHANGE
%  ============================================================

pressureChange = finalPressure - ...
                 convertTo(initialPressure, barsa);

%% ============================================================
%  16. RESULTS STRUCTURE
%  ============================================================

resultsWaterflood = struct();

% Simulation information
resultsWaterflood.simulationDays = simulationDays;
resultsWaterflood.simulationRuntime = simulationTime;

% Reservoir
resultsWaterflood.OOIP = OOIP;
resultsWaterflood.grid = [Nx Ny Nz];
resultsWaterflood.dimensions = [Lx Ly Lz];

% Production
resultsWaterflood.timeDays = timeDays;
resultsWaterflood.oilRateHistory = oilRateHistory;
resultsWaterflood.waterRateHistory = waterRateHistory;
resultsWaterflood.totalRateHistory = totalRateHistory;
resultsWaterflood.waterCutHistory = waterCutHistory;

resultsWaterflood.cumulativeOil = cumulativeOil;
resultsWaterflood.cumulativeWater = cumulativeWater;
resultsWaterflood.recoveryFactor = recoveryFactor;

% Pressure
resultsWaterflood.pressureHistory = pressureHistory;
resultsWaterflood.initialPressure = ...
    convertTo(initialPressure, barsa);

resultsWaterflood.finalPressure = finalPressure;
resultsWaterflood.pressureChange = pressureChange;

% Saturation
resultsWaterflood.avgSwHistory = avgSwHistory;
resultsWaterflood.avgSoHistory = avgSoHistory;

resultsWaterflood.initialSw = initialSw;
resultsWaterflood.finalSw = finalSw;

resultsWaterflood.initialSo = initialSo;
resultsWaterflood.finalSo = finalSo;

% Wells
resultsWaterflood.producerBHPHistory = ...
    producerBHPHistory;

resultsWaterflood.injectorBHPHistory = ...
    injectorBHPHistory;

resultsWaterflood.finalInjectorBHP = ...
    finalInjectorBHP;

resultsWaterflood.maxInjectorBHP = ...
    maxInjectorBHP;

% Final saturation statistics
resultsWaterflood.finalSwMin = finalSwMin;
resultsWaterflood.finalSwMax = finalSwMax;

resultsWaterflood.finalSoMin = finalSoMin;
resultsWaterflood.finalSoMax = finalSoMax;

% Final pressure statistics
resultsWaterflood.finalPressureMin = ...
    finalPressureMin;

resultsWaterflood.finalPressureMax = ...
    finalPressureMax;

resultsWaterflood.finalPressureMean = ...
    finalPressureMean;

%% ============================================================
%  17. SAVE RESULTS
%  ============================================================

% Re-define project/results paths to ensure they exist
projectRoot = 'D:\PROJECT';
resultsFolder = fullfile(projectRoot, '06_results');

if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

resultsFile = fullfile( ...
    resultsFolder, ...
    'waterflood_1year_results_v03.mat');

save(resultsFile, 'resultsWaterflood');
%% ============================================================
%  18. DISPLAY FINAL PRODUCTION RESULTS
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                    PRODUCTION\n');
fprintf('============================================================\n');

fprintf('Final oil rate          : %.4f m3/day\n', ...
    finalOilRate);

fprintf('Final water rate        : %.4f m3/day\n', ...
    finalWaterRate);

fprintf('Final total rate        : %.4f m3/day\n', ...
    finalTotalRate);

fprintf('Final water cut         : %.2f %%\n', ...
    finalWaterCut);

fprintf('Final producer BHP      : %.2f bar\n', ...
    finalProducerBHP);

%% ============================================================
%  19. CUMULATIVE PRODUCTION
%  ============================================================

fprintf('\n');
fprintf('Cumulative production:\n');

fprintf('Cumulative oil          : %.2f m3\n', ...
    finalCumulativeOil);

fprintf('Cumulative water        : %.2f m3\n', ...
    finalCumulativeWater);

fprintf('Recovery factor         : %.4f %%\n', ...
    finalRecoveryFactor);

%% ============================================================
%  20. PRESSURE
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                       PRESSURE\n');
fprintf('============================================================\n');

fprintf('Initial avg pressure    : %.2f bar\n', ...
    convertTo(initialPressure, barsa));

fprintf('Final avg pressure      : %.2f bar\n', ...
    finalPressure);

fprintf('Pressure change         : %+0.2f bar\n', ...
    pressureChange);

%% ============================================================
%  21. SATURATION
%  ============================================================

fprintf('\n');
fprintf('Saturation:\n');

fprintf('Initial avg Sw          : %.4f\n', ...
    initialSw);

fprintf('Final avg Sw            : %.4f\n', ...
    finalSw);

fprintf('Initial avg So          : %.4f\n', ...
    initialSo);

fprintf('Final avg So            : %.4f\n', ...
    finalSo);

%% ============================================================
%  22. INJECTOR BHP
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                    INJECTOR BHP\n');
fprintf('============================================================\n');

for n = 1:nInjectors

    fprintf('I%d final BHP           : %.2f bar\n', ...
        n, finalInjectorBHP(n));

end

fprintf('Maximum injector BHP    : %.2f bar\n', ...
    maxInjectorBHP);

%% ============================================================
%  23. WATER INJECTION
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                  INJECTION INFORMATION\n');
fprintf('============================================================\n');

fprintf('Injection per injector : %.2f m3/day\n', ...
    convertTo(qInjection, meter^3/day));

fprintf('Total water injection   : %.2f m3/day\n', ...
    convertTo(totalInjectionRate, meter^3/day));

fprintf('Total injected water    : %.2f m3\n', ...
    convertTo(totalInjectionRate, meter^3/day) ...
    * simulationDays);

%% ============================================================
%  24. FINAL SATURATION STATISTICS
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('             FINAL SATURATION STATISTICS\n');
fprintf('============================================================\n');

fprintf('Minimum Sw : %.4f\n', finalSwMin);
fprintf('Maximum Sw : %.4f\n', finalSwMax);

fprintf('Minimum So : %.4f\n', finalSoMin);
fprintf('Maximum So : %.4f\n', finalSoMax);

%% ============================================================
%  25. FINAL PRESSURE STATISTICS
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('              FINAL PRESSURE STATISTICS\n');
fprintf('============================================================\n');

fprintf('Minimum P : %.2f bar\n', ...
    finalPressureMin);

fprintf('Maximum P : %.2f bar\n', ...
    finalPressureMax);

fprintf('Mean P    : %.2f bar\n', ...
    finalPressureMean);

%% ============================================================
%  26. SAVE CONFIRMATION
%  ============================================================

fprintf('\n');
fprintf('Results saved to:\n');
fprintf('%s\n', resultsFile);

%% ============================================================
%  27. FINAL SUMMARY
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                 FINAL 1-YEAR SUMMARY\n');
fprintf('============================================================\n');

fprintf('Simulation time          : %.1f days\n', ...
    simulationDays);

fprintf('Cumulative oil           : %.2f m3\n', ...
    finalCumulativeOil);

fprintf('Cumulative water         : %.2f m3\n', ...
    finalCumulativeWater);

fprintf('Recovery factor          : %.4f %%\n', ...
    finalRecoveryFactor);

fprintf('Final oil rate           : %.4f m3/day\n', ...
    finalOilRate);

fprintf('Final water cut          : %.2f %%\n', ...
    finalWaterCut);

fprintf('Final average pressure   : %.2f bar\n', ...
    finalPressure);

fprintf('Final average Sw         : %.4f\n', ...
    finalSw);

fprintf('Final average So         : %.4f\n', ...
    finalSo);

fprintf('Total water injected     : %.2f m3/day\n', ...
    convertTo(totalInjectionRate, meter^3/day));

fprintf('Maximum injector BHP     : %.2f bar\n', ...
    maxInjectorBHP);

fprintf('============================================================\n');

fprintf('\n');
fprintf('1-YEAR WATERFLOOD SIMULATION FINISHED.\n');
fprintf('\n');