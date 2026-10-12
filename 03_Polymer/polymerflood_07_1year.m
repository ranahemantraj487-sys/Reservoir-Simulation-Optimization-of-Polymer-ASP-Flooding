%% polymerflood_07_1year.m
% =============================================================
% 1-YEAR POLYMER FLOOD SIMULATION
% Cairn Rajasthan Inspired Polymer EOR Study
%
% MRST 2026a
% MATLAB R2026b
%
% Based on validated 30-day polymer model
%
% Reservoir:
%   50 x 50 x 3 grid
%   400 m x 400 m x 30 m  (as defined in fluid_model_03.m)
%   Layered heterogeneous permeability / random porosity
%
% Well pattern:
%   1 central producer
%   4 corner injectors
%
% Polymer:
%   Injection concentration = 1 kg/m3
%   Todd-Longstaff mixing parameter = 0.85
%   Water viscosity multiplier = 1 + 4*C/Cmax
%   Adsorption = 0
%   RRF = 1
%   DPS = 0
%
% Simulation:
%   365 days
%   365 daily timesteps
%
% Outputs:
%   Oil rate history
%   Water rate history
%   Water cut history
%   Cumulative oil
%   Cumulative water
%   Recovery factor
%   Pressure history
%   Polymer concentration history
%   Injector BHP history
% =============================================================

clear;
clc;
close all;

%% =============================================================
% 1. LOAD MRST MODULES
% =============================================================

mrstModule add ad-core ad-blackoil ad-eor ad-props;

fprintf('\n');
fprintf('============================================================\n');
fprintf('        1-YEAR POLYMER FLOOD SIMULATION - MRST\n');
fprintf('============================================================\n\n');


%% =============================================================
% 2. PROJECT PATH
% =============================================================

projectRoot = fileparts(fileparts(mfilename('fullpath')));

baselineFile = fullfile(projectRoot, ...
    '01_Baseline', 'fluid_model_03.m');

fprintf('Project root:\n%s\n\n', projectRoot);


%% =============================================================
% 3. LOAD BASELINE MODEL
% =============================================================

fprintf('Loading baseline reservoir model...\n');

run(baselineFile);

close all;   % fluid_model_03.m opens plot windows

% Well positions below assume the 50 x 50 x 3 grid of fluid_model_03.m
if G.cartDims(1) ~= 50 || G.cartDims(2) ~= 50 || G.cartDims(3) ~= 3
    error('Baseline grid is %d x %d x %d; well positions assume 50 x 50 x 3.', G.cartDims);
end

fprintf('Baseline model loaded successfully.\n\n');


%% =============================================================
% 4. RESERVOIR INFORMATION
% =============================================================

fprintf('--- Reservoir Information ---\n');

fprintf('Grid cells          : %d\n', G.cells.num);

fprintf('Initial pressure    : %.2f bar\n', ...
    convertTo(mean(state0.pressure), barsa));

fprintf('Initial Sw          : %.4f\n', ...
    mean(state0.s(:,1)));

fprintf('Initial So          : %.4f\n', ...
    mean(state0.s(:,2)));

fprintf('Oil viscosity       : %.1f cP\n', convertTo(muO, centi*poise));
fprintf('Water viscosity     : %.1f cP\n', convertTo(muW, centi*poise));
fprintf('Oil density         : %.0f kg/m3\n', convertTo(rhoO, kilogram/meter^3));
fprintf('Water density       : %.0f kg/m3\n\n', convertTo(rhoW, kilogram/meter^3));


%% =============================================================
% 5. POLYMER FLUID PROPERTIES
% =============================================================

fprintf('Setting polymer properties...\n');

% Maximum polymer concentration
fluid.cpmax = 1.0;

% Todd-Longstaff mixing parameter
fluid.mixPar = 0.85;

% Polymer-induced water viscosity multiplier
fluid.muWMult = @(c) 1 + 4*(c/fluid.cpmax);

% Adsorption index
fluid.adsInx = 1;

% Maximum adsorption
fluid.adsMax = 1.0;

% Zero adsorption for current base case
fluid.ads = @(c) zeros(size(c));

% Resistance reduction factor
fluid.rrf = 1.0;

% DPS / inaccessible pore-space parameter
fluid.dps = 0.0;

% Rock density
fluid.rhoR = 2650*kilogram/meter^3;


fprintf('Maximum polymer concentration : %.2f kg/m3\n', ...
    fluid.cpmax);

fprintf('Mixing parameter              : %.2f\n', ...
    fluid.mixPar);

fprintf('Water viscosity multiplier    : 1 + 4*C/Cmax\n');

fprintf('Adsorption                    : 0\n');

fprintf('RRF                           : %.2f\n', ...
    fluid.rrf);

fprintf('DPS                           : %.2f\n\n', ...
    fluid.dps);


%% =============================================================
% 6. CREATE POLYMER MODEL
% =============================================================

fprintf('Creating OilWaterPolymerModel...\n');

modelPolymer = OilWaterPolymerModel(G, rock, fluid);

fprintf('Polymer model created successfully.\n\n');


%% =============================================================
% 7. INITIAL POLYMER STATE
% =============================================================

state0Polymer = state0;

state0Polymer.cp = zeros(G.cells.num,1);

state0Polymer.cpmax = zeros(G.cells.num,1);


fprintf('Initial polymer concentration:\n');

fprintf('Average Cp = %.6f kg/m3\n', ...
    mean(state0Polymer.cp));

fprintf('Maximum Cp = %.6f kg/m3\n\n', ...
    max(state0Polymer.cp));


%% =============================================================
% 8. WELL CONFIGURATION
% =============================================================

fprintf('Creating 5-spot well pattern...\n');

% Grid size taken from the baseline model (fluid_model_03.m)
Nx = G.cartDims(1);
Ny = G.cartDims(2);
Nz = G.cartDims(3);


% -------------------------------------------------------------
% Producer
% -------------------------------------------------------------

iP = 25;
jP = 25;


% -------------------------------------------------------------
% Four injectors
% -------------------------------------------------------------

iI = [4 4 47 47];
jI = [4 47 4 47];


W = [];


% -------------------------------------------------------------
% Producer P1
% -------------------------------------------------------------

cellP = sub2ind([Nx Ny Nz], ...
    iP, jP, 1:Nz);

W = addWell(W, G, rock, cellP, ...
    'Type', 'bhp', ...
    'Val', 150*barsa, ...
    'Name', 'P1', ...
    'Comp_i', [0 1]);


% -------------------------------------------------------------
% Injectors
% -------------------------------------------------------------

qInjection = 50*meter^3/day;

for n = 1:4

    cellI = sub2ind([Nx Ny Nz], ...
        iI(n), jI(n), 1:Nz);

    W = addWell(W, G, rock, cellI, ...
        'Type', 'rate', ...
        'Val', qInjection, ...
        'Name', sprintf('I%d', n), ...
        'Comp_i', [1 0]);

end


%% =============================================================
% 9. POLYMER WELL CONCENTRATIONS
% =============================================================

fprintf('Assigning polymer concentrations to wells...\n');

% Producer
W(1).cp = 0.0;


% Injectors
for n = 2:numel(W)

    W(n).cp = 1.0;

end


% -------------------------------------------------------------
% Verify polymer well-property mapping
% -------------------------------------------------------------

fprintf('\nMRST polymer well-property check:\n');

for n = 1:numel(W)

    cpCheck = modelPolymer.getProp(W(n), 'polymer');

    fprintf('%s cp = %.2f kg/m3\n', ...
        W(n).name, cpCheck);

end

fprintf('\n');


%% =============================================================
% 10. TOTAL INJECTION RATE
% =============================================================

nInjectors = 4;

totalInjectionRate = nInjectors * qInjection;

totalInjectionRate_m3day = ...
    convertTo(totalInjectionRate, meter^3/day);


fprintf('Total water injection rate : %.2f m3/day\n', ...
    totalInjectionRate_m3day);

fprintf('Polymer concentration      : %.2f kg/m3\n', ...
    W(2).cp);

fprintf('Polymer injection rate     : %.2f kg/day\n\n', ...
    totalInjectionRate_m3day * W(2).cp);


%% =============================================================
% 11. OOIP
% =============================================================

poreVolume = G.cells.volumes .* rock.poro;

initialOilVolume = ...
    sum(poreVolume .* state0Polymer.s(:,2));


fprintf('Initial oil volume (OOIP)  : %.2f m3\n\n', ...
    initialOilVolume);


%% =============================================================
% 12. 1-YEAR SIMULATION SCHEDULE
% =============================================================

simulationDays = 365;

dt = repmat(1*day, simulationDays, 1);

schedule = simpleSchedule(dt, 'W', W);


fprintf('Simulation schedule:\n');

fprintf('Number of timesteps : %d\n', ...
    numel(dt));

fprintf('Simulation time     : %.1f days\n\n', ...
    sum(convertTo(dt, day)));


%% =============================================================
% 13. RUN SIMULATION
% =============================================================

fprintf('============================================================\n');
fprintf('             STARTING 1-YEAR SIMULATION\n');
fprintf('============================================================\n\n');


[wellSols, states, report] = ...
    simulateScheduleAD( ...
    state0Polymer, ...
    modelPolymer, ...
    schedule);


fprintf('\n');
fprintf('============================================================\n');
fprintf('            1-YEAR SIMULATION COMPLETE\n');
fprintf('============================================================\n\n');


%% =============================================================
% 14. PREALLOCATE HISTORY ARRAYS
% =============================================================

nSteps = numel(states);

timeDays = zeros(nSteps,1);

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

injectorBHPHistory = zeros(nSteps,4);


%% =============================================================
% 15. EXTRACT PRODUCTION HISTORY
% =============================================================

fprintf('Extracting production history...\n');


for t = 1:nSteps

    % Current state
    currentState = states{t};

    % Current well solution
    currentWellSol = wellSols{t};


    % Time
    timeDays(t) = sum(convertTo(dt(1:t), day));


    % ---------------------------------------------------------
    % Producer rates
    % ---------------------------------------------------------

    oilRate = abs(currentWellSol(1).qOs);

    waterRate = abs(currentWellSol(1).qWs);

    totalRate = oilRate + waterRate;


    oilRateHistory(t) = ...
        convertTo(oilRate, meter^3/day);

    waterRateHistory(t) = ...
        convertTo(waterRate, meter^3/day);

    totalRateHistory(t) = ...
        convertTo(totalRate, meter^3/day);


    % Water cut
    waterCutHistory(t) = ...
        100 * waterRate / max(totalRate, eps);


    % Producer BHP
    producerBHPHistory(t) = ...
        convertTo(currentWellSol(1).bhp, barsa);


    % ---------------------------------------------------------
    % Reservoir pressure
    % ---------------------------------------------------------

    pressureHistory(t) = ...
        convertTo(mean(currentState.pressure), barsa);


    % ---------------------------------------------------------
    % Saturation
    % ---------------------------------------------------------

    avgSwHistory(t) = ...
        mean(currentState.s(:,1));

    avgSoHistory(t) = ...
        mean(currentState.s(:,2));


    % ---------------------------------------------------------
    % Polymer
    % ---------------------------------------------------------

    avgCpHistory(t) = ...
        mean(currentState.cp);

    maxCpHistory(t) = ...
        max(currentState.cp);


    % ---------------------------------------------------------
    % Injector BHP
    % ---------------------------------------------------------

    for n = 1:4

        injectorBHPHistory(t,n) = ...
            convertTo(currentWellSol(n+1).bhp, barsa);

    end

end


fprintf('Production history extracted successfully.\n\n');


%% =============================================================
% 16. CUMULATIVE PRODUCTION
% =============================================================

% Time step in days
dtDays = convertTo(dt, day);


% -------------------------------------------------------------
% Cumulative oil
% -------------------------------------------------------------

cumulativeOil = zeros(nSteps,1);

cumulativeWater = zeros(nSteps,1);

cumulativeOil(1) = ...
    oilRateHistory(1) * dtDays(1);

cumulativeWater(1) = ...
    waterRateHistory(1) * dtDays(1);


for t = 2:nSteps

    cumulativeOil(t) = ...
        cumulativeOil(t-1) + ...
        0.5 * ...
        (oilRateHistory(t-1) + oilRateHistory(t)) * ...
        dtDays(t);


    cumulativeWater(t) = ...
        cumulativeWater(t-1) + ...
        0.5 * ...
        (waterRateHistory(t-1) + waterRateHistory(t)) * ...
        dtDays(t);

end


%% =============================================================
% 17. RECOVERY FACTOR
% =============================================================

recoveryFactor = ...
    100 * cumulativeOil / initialOilVolume;


finalCumulativeOil = cumulativeOil(end);

finalCumulativeWater = cumulativeWater(end);

finalRecoveryFactor = recoveryFactor(end);


%% =============================================================
% 18. FINAL STATE
% =============================================================

finalState = states{end};

finalWellSol = wellSols{end};


%% =============================================================
% 19. FINAL RESULTS
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('             1-YEAR POLYMER FLOOD RESULTS\n');
fprintf('============================================================\n');


fprintf('\nProduction:\n');

fprintf('Final oil rate          : %.4f m3/day\n', ...
    oilRateHistory(end));

fprintf('Final water rate        : %.4f m3/day\n', ...
    waterRateHistory(end));

fprintf('Final total rate        : %.4f m3/day\n', ...
    totalRateHistory(end));

fprintf('Final water cut         : %.2f %%\n', ...
    waterCutHistory(end));

fprintf('Final producer BHP      : %.2f bar\n', ...
    producerBHPHistory(end));


fprintf('\nCumulative production:\n');

fprintf('Cumulative oil          : %.2f m3\n', ...
    finalCumulativeOil);

fprintf('Cumulative water        : %.2f m3\n', ...
    finalCumulativeWater);

fprintf('Recovery factor         : %.4f %%\n', ...
    finalRecoveryFactor);


fprintf('\nPressure:\n');

fprintf('Initial avg pressure    : %.2f bar\n', ...
    pressureHistory(1));

fprintf('Final avg pressure      : %.2f bar\n', ...
    pressureHistory(end));

fprintf('Pressure change         : %+0.2f bar\n', ...
    pressureHistory(end) - pressureHistory(1));


fprintf('\nSaturation:\n');

fprintf('Initial avg Sw          : %.4f\n', ...
    avgSwHistory(1));

fprintf('Final avg Sw            : %.4f\n', ...
    avgSwHistory(end));

fprintf('Initial avg So          : %.4f\n', ...
    avgSoHistory(1));

fprintf('Final avg So            : %.4f\n', ...
    avgSoHistory(end));


fprintf('\nPolymer:\n');

fprintf('Final average Cp        : %.6f kg/m3\n', ...
    avgCpHistory(end));

fprintf('Final maximum Cp        : %.6f kg/m3\n', ...
    maxCpHistory(end));


fprintf('\nInjector BHP:\n');

for n = 1:4

    fprintf('I%d final BHP           : %.2f bar\n', ...
        n, injectorBHPHistory(end,n));

end

fprintf('Maximum injector BHP    : %.2f bar\n', ...
    max(injectorBHPHistory(:)));


%% =============================================================
% 20. POLYMER MASS BALANCE
% =============================================================

fprintf('\nPolymer mass balance:\n');


totalSimulationDays = ...
    sum(dtDays);


injectionPolymerConcentration = 1.0;


% Expected polymer injected
polymerInjectedExpected = ...
    totalInjectionRate_m3day * ...
    injectionPolymerConcentration * ...
    totalSimulationDays;


% Polymer currently dissolved in water
finalWaterSaturation = finalState.s(:,1);

polymerMassInWater = ...
    sum(poreVolume .* ...
    finalWaterSaturation .* ...
    finalState.cp);


massBalanceDifference = ...
    polymerInjectedExpected - polymerMassInWater;


massBalanceClosure = ...
    100 * polymerMassInWater / ...
    polymerInjectedExpected;


fprintf('Simulation time             : %.1f days\n', ...
    totalSimulationDays);

fprintf('Total water injection       : %.2f m3/day\n', ...
    totalInjectionRate_m3day);

fprintf('Expected polymer injected   : %.2f kg\n', ...
    polymerInjectedExpected);

fprintf('Polymer in water            : %.2f kg\n', ...
    polymerMassInWater);

fprintf('Mass difference              : %.2f kg\n', ...
    massBalanceDifference);

fprintf('Mass-balance closure         : %.4f %%\n', ...
    massBalanceClosure);


%% =============================================================
% 21. FINAL SATURATION STATISTICS
% =============================================================

fprintf('\nFinal saturation statistics:\n');

fprintf('Minimum Sw : %.4f\n', ...
    min(finalState.s(:,1)));

fprintf('Maximum Sw : %.4f\n', ...
    max(finalState.s(:,1)));

fprintf('Minimum So : %.4f\n', ...
    min(finalState.s(:,2)));

fprintf('Maximum So : %.4f\n', ...
    max(finalState.s(:,2)));


%% =============================================================
% 22. FINAL POLYMER STATISTICS
% =============================================================

fprintf('\nFinal polymer statistics:\n');

fprintf('Minimum Cp : %.6f kg/m3\n', ...
    min(finalState.cp));

fprintf('Maximum Cp : %.6f kg/m3\n', ...
    max(finalState.cp));

fprintf('Mean Cp    : %.6f kg/m3\n', ...
    mean(finalState.cp));


%% =============================================================
% 23. FINAL PRESSURE STATISTICS
% =============================================================

fprintf('\nFinal pressure statistics:\n');

fprintf('Minimum P : %.2f bar\n', ...
    convertTo(min(finalState.pressure), barsa));

fprintf('Maximum P : %.2f bar\n', ...
    convertTo(max(finalState.pressure), barsa));

fprintf('Mean P    : %.2f bar\n', ...
    convertTo(mean(finalState.pressure), barsa));


%% =============================================================
% 24. PLOT OIL RATE
% =============================================================

figure('Name', '1-Year Polymer Flood - Oil Rate');

plot(timeDays, oilRateHistory, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Oil Rate (m3/day)');

title('Polymer Flood - Oil Production Rate');


%% =============================================================
% 25. PLOT WATER RATE
% =============================================================

figure('Name', '1-Year Polymer Flood - Water Rate');

plot(timeDays, waterRateHistory, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Water Rate (m3/day)');

title('Polymer Flood - Water Production Rate');


%% =============================================================
% 26. PLOT WATER CUT
% =============================================================

figure('Name', '1-Year Polymer Flood - Water Cut');

plot(timeDays, waterCutHistory, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Water Cut (%)');

title('Polymer Flood - Water Cut');


%% =============================================================
% 27. PLOT CUMULATIVE OIL
% =============================================================

figure('Name', '1-Year Polymer Flood - Cumulative Oil');

plot(timeDays, cumulativeOil, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Cumulative Oil (m3)');

title('Polymer Flood - Cumulative Oil Production');


%% =============================================================
% 28. PLOT RECOVERY FACTOR
% =============================================================

figure('Name', '1-Year Polymer Flood - Recovery');

plot(timeDays, recoveryFactor, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Recovery Factor (%)');

title('Polymer Flood - Recovery Factor');


%% =============================================================
% 29. PLOT RESERVOIR PRESSURE
% =============================================================

figure('Name', '1-Year Polymer Flood - Pressure');

plot(timeDays, pressureHistory, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Average Reservoir Pressure (bar)');

title('Polymer Flood - Average Reservoir Pressure');


%% =============================================================
% 30. PLOT POLYMER CONCENTRATION
% =============================================================

figure('Name', '1-Year Polymer Flood - Polymer');

plot(timeDays, avgCpHistory, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Average Polymer Concentration (kg/m3)');

title('Polymer Flood - Average Polymer Concentration');


%% =============================================================
% 31. PLOT MAXIMUM POLYMER
% =============================================================

figure('Name', '1-Year Polymer Flood - Maximum Cp');

plot(timeDays, maxCpHistory, 'LineWidth', 1.5);

grid on;

xlabel('Time (days)');
ylabel('Maximum Cp (kg/m3)');

title('Polymer Flood - Maximum Polymer Concentration');


%% =============================================================
% 32. PLOT INJECTOR BHP
% =============================================================

figure('Name', '1-Year Polymer Flood - Injector BHP');

plot(timeDays, injectorBHPHistory, 'LineWidth', 1.2);

grid on;

xlabel('Time (days)');
ylabel('Injector BHP (bar)');

title('Polymer Flood - Injector BHP History');

legend('I1','I2','I3','I4', ...
    'Location','best');


%% =============================================================
% 33. FINAL RESERVOIR WATER SATURATION
% =============================================================

figure('Name', '1-Year Polymer Flood - Final Sw');

plotCellData(G, finalState.s(:,1));

view(3);
axis tight;
colorbar;

title('1-Year Polymer Flood - Final Water Saturation');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% =============================================================
% 34. FINAL POLYMER DISTRIBUTION
% =============================================================

figure('Name', '1-Year Polymer Flood - Final Cp');

plotCellData(G, finalState.cp);

view(3);
axis tight;
colorbar;

title('1-Year Polymer Flood - Final Polymer Concentration');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% =============================================================
% 35. FINAL PRESSURE DISTRIBUTION
% =============================================================

figure('Name', '1-Year Polymer Flood - Final Pressure');

plotCellData(G, ...
    convertTo(finalState.pressure, barsa));

view(3);
axis tight;
colorbar;

title('1-Year Polymer Flood - Final Pressure');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% =============================================================
% 36. FINAL OIL SATURATION
% =============================================================

figure('Name', '1-Year Polymer Flood - Final So');

plotCellData(G, finalState.s(:,2));

view(3);
axis tight;
colorbar;

title('1-Year Polymer Flood - Final Oil Saturation');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% =============================================================
% 37. SAVE RESULTS
% =============================================================

% Re-establish project root for saving results
projectRoot = 'D:\PROJECT';

results1Year.timeDays = timeDays;

results1Year.oilRate = oilRateHistory;

results1Year.waterRate = waterRateHistory;

results1Year.totalRate = totalRateHistory;

results1Year.waterCut = waterCutHistory;

results1Year.cumulativeOil = cumulativeOil;

results1Year.cumulativeWater = cumulativeWater;

results1Year.recoveryFactor = recoveryFactor;

results1Year.pressure = pressureHistory;

results1Year.avgSw = avgSwHistory;

results1Year.avgSo = avgSoHistory;

results1Year.avgCp = avgCpHistory;

results1Year.maxCp = maxCpHistory;

results1Year.producerBHP = producerBHPHistory;

results1Year.injectorBHP = injectorBHPHistory;

results1Year.OOIP = initialOilVolume;

results1Year.polymerInjected = polymerInjectedExpected;

results1Year.polymerInWater = polymerMassInWater;

results1Year.massBalanceClosure = massBalanceClosure;


save(fullfile(projectRoot, ...
    '06_results', ...
    'polymer_1year_results_v03.mat'), ...
    'results1Year');


fprintf('\nResults saved to:\n');

fprintf('%s\n', ...
    fullfile(projectRoot, ...
    '06_results', ...
    'polymer_1year_results_v03.mat'));


%% =============================================================
% 38. FINAL SUMMARY
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                 FINAL 1-YEAR SUMMARY\n');
fprintf('============================================================\n');

fprintf('Simulation time          : %.1f days\n', ...
    totalSimulationDays);

fprintf('Cumulative oil           : %.2f m3\n', ...
    finalCumulativeOil);

fprintf('Cumulative water         : %.2f m3\n', ...
    finalCumulativeWater);

fprintf('Recovery factor          : %.4f %%\n', ...
    finalRecoveryFactor);

fprintf('Final oil rate           : %.4f m3/day\n', ...
    oilRateHistory(end));

fprintf('Final water cut          : %.2f %%\n', ...
    waterCutHistory(end));

fprintf('Final average pressure   : %.2f bar\n', ...
    pressureHistory(end));

fprintf('Final average Sw         : %.4f\n', ...
    avgSwHistory(end));

fprintf('Final average Cp         : %.6f kg/m3\n', ...
    avgCpHistory(end));

fprintf('Final maximum Cp         : %.6f kg/m3\n', ...
    maxCpHistory(end));

fprintf('Polymer injected         : %.2f kg\n', ...
    polymerInjectedExpected);

fprintf('Polymer in water         : %.2f kg\n', ...
    polymerMassInWater);

fprintf('Mass-balance closure     : %.2f %%\n', ...
    massBalanceClosure);

fprintf('Maximum injector BHP     : %.2f bar\n', ...
    max(injectorBHPHistory(:)));

fprintf('============================================================\n');

fprintf('\n1-YEAR POLYMER FLOOD SIMULATION FINISHED.\n\n');