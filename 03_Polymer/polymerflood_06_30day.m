%% polymerflood_06_30day.m
% 30-Day Polymer Flood Validation
% Cairn Rajasthan Inspired Polymer EOR Study
%
% MRST 2026a
% MATLAB R2026b
%
% Model:
%   50 x 50 x 3 layered heterogeneous reservoir (fluid_model_03.m)
%   5-spot water/polymer flood
%   Producer P1 at center
%   Four corner injectors
%
% Polymer:
%   Injection concentration = 1 kg/m^3
%   Todd-Longstaff mixing parameter = 0.85
%   Water viscosity multiplier = 1 + 4*C/Cmax
%   Adsorption = 0
%   RRF = 1
%   DPS = 0
%
% Simulation:
%   30 days
%   30 daily timesteps

clear;
clc;
close all;

%% ============================================================
% 1. LOAD MRST MODULES
% =============================================================

mrstModule add ad-core ad-blackoil ad-eor ad-props;

fprintf('\n');
fprintf('============================================================\n');
fprintf('      30-DAY POLYMER FLOOD VALIDATION - MRST\n');
fprintf('============================================================\n\n');


%% ============================================================
% 2. PROJECT PATH
% =============================================================

projectRoot = fileparts(fileparts(mfilename('fullpath')));

baselineFile = fullfile(projectRoot, ...
    '01_Baseline', 'fluid_model_03.m');

fprintf('Project root:\n%s\n\n', projectRoot);


%% ============================================================
% 3. LOAD BASELINE RESERVOIR MODEL
% =============================================================

fprintf('Loading baseline reservoir model...\n');

run(baselineFile);

close all;   % fluid_model_03.m opens plot windows

% Well positions below assume the 50 x 50 x 3 grid of fluid_model_03.m
if G.cartDims(1) ~= 50 || G.cartDims(2) ~= 50 || G.cartDims(3) ~= 3
    error('Baseline grid is %d x %d x %d; well positions assume 50 x 50 x 3.', G.cartDims);
end

fprintf('Baseline model loaded successfully.\n\n');


%% ============================================================
% 4. BASIC RESERVOIR INFORMATION
% =============================================================

fprintf('Reservoir information:\n');

fprintf('Grid cells          : %d\n', G.cells.num);

fprintf('Initial pressure    : %.2f bar\n', ...
    convertTo(mean(state0.pressure), barsa));

fprintf('Initial Sw          : %.4f\n', ...
    mean(state0.s(:,1)));

fprintf('Initial So          : %.4f\n', ...
    mean(state0.s(:,2)));

% Fluid properties defined in the baseline model
fprintf('Oil viscosity       : %.1f cP\n', convertTo(muO, centi*poise));
fprintf('Water viscosity     : %.1f cP\n', convertTo(muW, centi*poise));
fprintf('Oil density         : %.0f kg/m3\n', convertTo(rhoO, kilogram/meter^3));
fprintf('Water density       : %.0f kg/m3\n\n', convertTo(rhoW, kilogram/meter^3));


%% ============================================================
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

% Zero adsorption for current validation
fluid.ads = @(c) zeros(size(c));

% Resistance Reduction Factor
fluid.rrf = 1.0;

% Dispersion / inaccessible pore-space parameter
fluid.dps = 0.0;

% Rock density
fluid.rhoR = 2650*kilogram/meter^3;

fprintf('Maximum polymer concentration : %.2f kg/m3\n', ...
    fluid.cpmax);

fprintf('Mixing parameter              : %.2f\n', ...
    fluid.mixPar);

fprintf('Polymer viscosity multiplier   : 1 + 4*C/Cmax\n');

fprintf('Adsorption                     : 0\n');
fprintf('RRF                            : %.2f\n', ...
    fluid.rrf);

fprintf('DPS                            : %.2f\n\n', ...
    fluid.dps);


%% ============================================================
% 6. CREATE POLYMER MODEL
% =============================================================

fprintf('Creating OilWaterPolymerModel...\n');

modelPolymer = OilWaterPolymerModel(G, rock, fluid);

fprintf('Polymer model created successfully.\n\n');


%% ============================================================
% 7. INITIAL POLYMER STATE
% =============================================================

state0Polymer = state0;

% Polymer concentration
state0Polymer.cp = zeros(G.cells.num,1);

% Maximum polymer concentration state variable
state0Polymer.cpmax = zeros(G.cells.num,1);

fprintf('Initial polymer concentration:\n');
fprintf('Average Cp = %.6f kg/m3\n', ...
    mean(state0Polymer.cp));

fprintf('Maximum Cp = %.6f kg/m3\n\n', ...
    max(state0Polymer.cp));


%% ============================================================
% 8. WELL CONFIGURATION
% =============================================================

fprintf('Creating 5-spot well pattern...\n');

% Grid size taken from the baseline model (fluid_model_03.m)
Nx = G.cartDims(1);
Ny = G.cartDims(2);
Nz = G.cartDims(3);

% ------------------------------------------------------------
% Producer
% ------------------------------------------------------------

iP = 25;
jP = 25;

% ------------------------------------------------------------
% Four injectors
% ------------------------------------------------------------

iI = [4 4 47 47];
jI = [4 47 4 47];

W = [];

% ------------------------------------------------------------
% Producer P1
% ------------------------------------------------------------

cellP = sub2ind([Nx Ny Nz], ...
    iP, jP, 1:Nz);

W = addWell(W, G, rock, cellP, ...
    'Type', 'bhp', ...
    'Val', 150*barsa, ...
    'Name', 'P1', ...
    'Comp_i', [0 1]);


% ------------------------------------------------------------
% Injectors
% ------------------------------------------------------------

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


%% ============================================================
% 9. POLYMER WELL CONCENTRATIONS
% ============================================================

fprintf('Assigning polymer concentrations to wells...\n');

% Producer receives zero polymer
W(1).cp = 0.0;

% Injectors inject polymer at 1 kg/m3
for n = 2:numel(W)

    W(n).cp = 1.0;

end


% ------------------------------------------------------------
% Verify MRST polymer well-property mapping
% ------------------------------------------------------------

fprintf('\nMRST polymer well-property check:\n');

for n = 1:numel(W)

    cpCheck = modelPolymer.getProp(W(n), 'polymer');

    fprintf('%s cp = %.2f kg/m3\n', ...
        W(n).name, cpCheck);

end

fprintf('\n');


%% ============================================================
% 10. TOTAL INJECTION RATE
% =============================================================

% Four injectors, each injecting 50 m3/day
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

%% ============================================================
% 11. 30-DAY SIMULATION SCHEDULE
% =============================================================

fprintf('Creating 30-day simulation schedule...\n');

% One timestep per day
dt = repmat(1*day, 30, 1);

schedule = simpleSchedule(dt, 'W', W);

fprintf('Number of timesteps : %d\n', numel(dt));

fprintf('Total simulation time: %.1f days\n\n', ...
    sum(convertTo(dt, day)));


%% ============================================================
% 12. RUN POLYMER FLOOD SIMULATION
% =============================================================

fprintf('============================================================\n');
fprintf('              STARTING 30-DAY SIMULATION\n');
fprintf('============================================================\n\n');

[wellSols, states, report] = simulateScheduleAD( ...
    state0Polymer, modelPolymer, schedule);

fprintf('\n');
fprintf('============================================================\n');
fprintf('             30-DAY SIMULATION COMPLETE\n');
fprintf('============================================================\n\n');


%% ============================================================
% 13. FINAL STATE
% =============================================================

finalState = states{end};

finalWellSol = wellSols{end};


%% ============================================================
% 14. FINAL PRODUCER RESULTS
% =============================================================

fprintf('30-DAY POLYMER FLOOD RESULTS\n');
fprintf('------------------------------------------------------------\n');

producerIndex = 1;

% Surface rates
oilRate = abs(finalWellSol(producerIndex).qOs);
waterRate = abs(finalWellSol(producerIndex).qWs);

totalRate = oilRate + waterRate;

waterCut = waterRate / max(totalRate, eps) * 100;

producerBHP = finalWellSol(producerIndex).bhp;


fprintf('\nProducer P1:\n');

fprintf('Oil rate             : %.4f m3/day\n', ...
    convertTo(oilRate, meter^3/day));

fprintf('Water rate           : %.4f m3/day\n', ...
    convertTo(waterRate, meter^3/day));

fprintf('Total rate           : %.4f m3/day\n', ...
    convertTo(totalRate, meter^3/day));

fprintf('Water cut            : %.2f %%\n', ...
    waterCut);

fprintf('BHP                  : %.2f bar\n\n', ...
    convertTo(producerBHP, barsa));


%% ============================================================
% 15. FINAL POLYMER DISTRIBUTION
% =============================================================

cpFinal = finalState.cp;

avgPolymer = mean(cpFinal);

maxPolymer = max(cpFinal);

fprintf('Polymer:\n');

fprintf('Average polymer      : %.6f kg/m3\n', ...
    avgPolymer);

fprintf('Maximum polymer      : %.6f kg/m3\n\n', ...
    maxPolymer);


%% ============================================================
% 16. PRESSURE RESULTS
% =============================================================

initialPressure = mean(state0Polymer.pressure);

finalPressure = mean(finalState.pressure);

pressureChange = finalPressure - initialPressure;

fprintf('Pressure:\n');

fprintf('Initial avg P        : %.4f bar\n', ...
    convertTo(initialPressure, barsa));

fprintf('Final avg P          : %.4f bar\n', ...
    convertTo(finalPressure, barsa));

fprintf('Pressure change      : %+0.4f bar\n\n', ...
    convertTo(pressureChange, barsa));


%% ============================================================
% 17. SATURATION RESULTS
% =============================================================

initialSw = mean(state0Polymer.s(:,1));
finalSw   = mean(finalState.s(:,1));

initialSo = mean(state0Polymer.s(:,2));
finalSo   = mean(finalState.s(:,2));

fprintf('Saturation:\n');

fprintf('Initial Sw           : %.4f\n', initialSw);
fprintf('Final Sw             : %.4f\n', finalSw);

fprintf('Initial So           : %.4f\n', initialSo);
fprintf('Final So             : %.4f\n\n', finalSo);


%% ============================================================
% 18. INJECTOR RESULTS
% =============================================================

fprintf('Injector performance:\n');

for n = 2:numel(finalWellSol)

    injectorName = finalWellSol(n).name;

    injectorRate = abs(finalWellSol(n).qWs);

    injectorBHP = finalWellSol(n).bhp;

    fprintf('%s: Rate = %.2f m3/day, BHP = %.2f bar\n', ...
        injectorName, ...
        convertTo(injectorRate, meter^3/day), ...
        convertTo(injectorBHP, barsa));

end

fprintf('\n');


%% ============================================================
% 19. POLYMER MASS BALANCE
% =============================================================

fprintf('Polymer mass balance:\n');

simulationDays = sum(convertTo(dt, day));

injectionPolymerConcentration = 1.0;

% Expected polymer injected
polymerInjectedExpected = ...
    totalInjectionRate_m3day * ...
    injectionPolymerConcentration * ...
    simulationDays;


% ------------------------------------------------------------
% Polymer mass currently dissolved in water
% ------------------------------------------------------------

poreVolume = G.cells.volumes .* rock.poro;

waterSaturation = finalState.s(:,1);

polymerMassInWater = ...
    sum(poreVolume .* waterSaturation .* cpFinal);


% ------------------------------------------------------------
% Mass balance difference
% ------------------------------------------------------------

massBalanceDifference = ...
    polymerInjectedExpected - polymerMassInWater;


massBalanceClosure = ...
    100 * polymerMassInWater / ...
    polymerInjectedExpected;


fprintf('Simulation time             : %.1f days\n', ...
    simulationDays);

fprintf('Total water injection       : %.4f m3/day\n', ...
    totalInjectionRate_m3day);

fprintf('Expected polymer injected   : %.4f kg\n', ...
    polymerInjectedExpected);

fprintf('Polymer in water            : %.4f kg\n', ...
    polymerMassInWater);

fprintf('Mass difference              : %.4f kg\n', ...
    massBalanceDifference);

fprintf('Mass-balance closure         : %.4f %%\n\n', ...
    massBalanceClosure);


%% ============================================================
% 20. OOIP AND RECOVERY
% =============================================================

% Initial oil pore volume
initialOilVolume = ...
    sum(G.cells.volumes .* rock.poro .* ...
    state0Polymer.s(:,2));

fprintf('Reservoir oil volume:\n');

fprintf('Initial oil volume : %.2f m3\n', ...
    initialOilVolume);


%% ============================================================
% 21. FINAL SATURATION STATISTICS
% =============================================================

fprintf('\nFinal saturation statistics:\n');

fprintf('Minimum Sw : %.4f\n', min(finalState.s(:,1)));
fprintf('Maximum Sw : %.4f\n', max(finalState.s(:,1)));

fprintf('Minimum So : %.4f\n', min(finalState.s(:,2)));
fprintf('Maximum So : %.4f\n', max(finalState.s(:,2)));


%% ============================================================
% 22. POLYMER STATISTICS
% =============================================================

fprintf('\nFinal polymer statistics:\n');

fprintf('Minimum Cp : %.6f kg/m3\n', min(cpFinal));
fprintf('Maximum Cp : %.6f kg/m3\n', max(cpFinal));
fprintf('Mean Cp    : %.6f kg/m3\n', mean(cpFinal));


%% ============================================================
% 23. FINAL PRESSURE STATISTICS
% =============================================================

fprintf('\nFinal pressure statistics:\n');

fprintf('Minimum P : %.2f bar\n', ...
    convertTo(min(finalState.pressure), barsa));

fprintf('Maximum P : %.2f bar\n', ...
    convertTo(max(finalState.pressure), barsa));

fprintf('Mean P    : %.2f bar\n', ...
    convertTo(mean(finalState.pressure), barsa));


%% ============================================================
% 24. PLOT WATER SATURATION
% =============================================================

figure('Name', '30-Day Polymer Flood - Water Saturation');

plotCellData(G, finalState.s(:,1));

view(3);
axis tight;
colorbar;

title('30-Day Polymer Flood - Final Water Saturation');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% ============================================================
% 25. PLOT POLYMER CONCENTRATION
% =============================================================

figure('Name', '30-Day Polymer Flood - Polymer Concentration');

plotCellData(G, cpFinal);

view(3);
axis tight;
colorbar;

title('30-Day Polymer Concentration');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% ============================================================
% 26. PLOT PRESSURE
% =============================================================

figure('Name', '30-Day Polymer Flood - Pressure');

plotCellData(G, ...
    convertTo(finalState.pressure, barsa));

view(3);
axis tight;
colorbar;

title('30-Day Polymer Flood - Final Pressure (bar)');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% ============================================================
% 27. PLOT OIL SATURATION
% =============================================================

figure('Name', '30-Day Polymer Flood - Oil Saturation');

plotCellData(G, finalState.s(:,2));

view(3);
axis tight;
colorbar;

title('30-Day Polymer Flood - Final Oil Saturation');

xlabel('X');
ylabel('Y');
zlabel('Z');


%% ============================================================
% 28. FINAL SUMMARY
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                 FINAL 30-DAY SUMMARY\n');
fprintf('============================================================\n');

fprintf('Simulation time       : %.1f days\n', ...
    simulationDays);

fprintf('Producer oil rate     : %.4f m3/day\n', ...
    convertTo(oilRate, meter^3/day));

fprintf('Producer water rate   : %.4f m3/day\n', ...
    convertTo(waterRate, meter^3/day));

fprintf('Water cut             : %.2f %%\n', ...
    waterCut);

fprintf('Producer BHP          : %.2f bar\n', ...
    convertTo(producerBHP, barsa));

fprintf('Average reservoir P   : %.2f bar\n', ...
    convertTo(finalPressure, barsa));

fprintf('Average Sw            : %.4f\n', ...
    finalSw);

fprintf('Average polymer       : %.6f kg/m3\n', ...
    avgPolymer);

fprintf('Maximum polymer       : %.6f kg/m3\n', ...
    maxPolymer);

fprintf('Polymer injected      : %.2f kg\n', ...
    polymerInjectedExpected);

fprintf('Polymer in water      : %.2f kg\n', ...
    polymerMassInWater);

fprintf('Mass balance closure  : %.2f %%\n', ...
    massBalanceClosure);

fprintf('============================================================\n');

fprintf('\n30-DAY POLYMER VALIDATION FINISHED SUCCESSFULLY.\n\n');