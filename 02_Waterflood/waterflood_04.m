%% ============================================================
% STAGE 5 - 5-SPOT WATERFLOOD
% 1-DAY VALIDATION
% =============================================================

clear;
clc;
close all;

mrstModule add ad-core ad-blackoil ad-props

%% ------------------------------------------------------------
% 1. Load the common reservoir model
% ------------------------------------------------------------

projectRoot = fileparts(fileparts(mfilename('fullpath')));

run(fullfile(projectRoot, ...
    '01_Baseline', 'fluid_model_02.m'));

fprintf('\n');
fprintf('====================================================\n');
fprintf('          STAGE 5 - WATERFLOOD TEST\n');
fprintf('====================================================\n');

%% ------------------------------------------------------------
% 2. Grid dimensions
% ------------------------------------------------------------

Nx = G.cartDims(1);
Ny = G.cartDims(2);
Nz = G.cartDims(3);

fprintf('\nGrid: %d x %d x %d\n', Nx, Ny, Nz);

%% ------------------------------------------------------------
% 3. Well locations
% ------------------------------------------------------------

% Central producer
iP = round(Nx/2);
jP = round(Ny/2);

% Offset from corners
offset = 4;

% Four injectors
iI = [offset, offset, Nx-offset+1, Nx-offset+1];
jI = [offset, Ny-offset+1, offset, Ny-offset+1];

%% ------------------------------------------------------------
% 4. Well controls
% ------------------------------------------------------------

producerBHP = 150 * barsa;

injRate = 50 * meter^3/day;

%% ------------------------------------------------------------
% 5. Create producer
% ------------------------------------------------------------

W = verticalWell( ...
    [], G, rock, ...
    iP, jP, 1:Nz, ...
    'Name', 'P1', ...
    'Type', 'bhp', ...
    'Val', producerBHP, ...
    'Radius', 0.1 * meter, ...
    'Comp_i', [0 1]);

%% ------------------------------------------------------------
% 6. Create four water injectors
% ------------------------------------------------------------

for n = 1:4

    W = verticalWell( ...
        W, G, rock, ...
        iI(n), jI(n), 1:Nz, ...
        'Name', ['I', num2str(n)], ...
        'Type', 'rate', ...
        'Val', injRate, ...
        'Radius', 0.1 * meter, ...
        'Comp_i', [1 0]);

end

%% ------------------------------------------------------------
% 7. Display well locations
% ------------------------------------------------------------

fprintf('\nProducer:\n');
fprintf('P1 -> i = %d, j = %d\n', iP, jP);

fprintf('\nInjectors:\n');

for n = 1:4
    fprintf('I%d -> i = %d, j = %d\n', ...
        n, iI(n), jI(n));
end

fprintf('\nInjection rate per injector = %.2f m3/day\n', ...
    convertTo(injRate, meter^3/day));

fprintf('Total injection rate = %.2f m3/day\n', ...
    convertTo(4*injRate, meter^3/day));

%% ------------------------------------------------------------
% 8. One-day schedule
% ------------------------------------------------------------

testTime = 1 * day;

schedule = simpleSchedule( ...
    testTime, ...
    'W', W);

%% ------------------------------------------------------------
% 9. Run simulation
% ------------------------------------------------------------

fprintf('\nRunning 1-day waterflood validation...\n');

[wellSols, states, reports] = ...
    simulateScheduleAD( ...
        state0, model, schedule);

fprintf('Simulation completed.\n');

%% ------------------------------------------------------------
% 10. Extract final well solution
% ------------------------------------------------------------

ws = wellSols{end};

%% ------------------------------------------------------------
% 11. Validation results
% ------------------------------------------------------------

fprintf('\n');
fprintf('--------------------------------------------\n');
fprintf('     1-DAY WATERFLOOD VALIDATION\n');
fprintf('--------------------------------------------\n');

fprintf('Producer BHP       : %.2f bar\n', ...
    convertTo(ws(1).bhp, barsa));

fprintf('Producer oil rate  : %.4f m3/day\n', ...
    -convertTo(ws(1).qOs, meter^3/day));

fprintf('Producer water rate: %.4f m3/day\n', ...
    -convertTo(ws(1).qWs, meter^3/day));

fprintf('Producer total rate: %.4f m3/day\n', ...
    -convertTo(ws(1).qTs, meter^3/day));

for n = 2:5

    fprintf('I%d water rate      : %.4f m3/day\n', ...
        n-1, ...
        convertTo(ws(n).qWs, meter^3/day));

end

%% ------------------------------------------------------------
% 12. Average pressure
% ------------------------------------------------------------

P_initial = convertTo(mean(state0.pressure), barsa);

P_final = convertTo(mean(states{end}.pressure), barsa);

fprintf('\nInitial average pressure : %.2f bar\n', ...
    P_initial);

fprintf('After 1 day              : %.2f bar\n', ...
    P_final);

fprintf('Pressure change          : %.4f bar\n', ...
    P_final - P_initial);

%% ------------------------------------------------------------
% 13. Final water saturation
% ------------------------------------------------------------

figure;

plotCellData(G, states{end}.s(:,1));

view(3);
axis tight equal;
colorbar;

title('Water Saturation After 1 Day');

%% ------------------------------------------------------------
% 14. Final pressure
% ------------------------------------------------------------

figure;

plotCellData(G, ...
    convertTo(states{end}.pressure, barsa));

view(3);
axis tight equal;
colorbar;

title('Pressure After 1 Day');

%% ------------------------------------------------------------
% END
% ------------------------------------------------------------

fprintf('\n');
fprintf('====================================================\n');
fprintf('       1-DAY VALIDATION COMPLETED\n');
fprintf('====================================================\n');
%% ============================================================
%  STAGE 5 - 10-YEAR WATERFLOOD SIMULATION
% =============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('       STAGE 5 - 10-YEAR WATERFLOOD\n');
fprintf('====================================================\n');

%% Simulation time

simTime = 10 * year;

% Use approximately 100 ramp-up timesteps
dt = rampupTimesteps(simTime, simTime/100);

schedule10 = simpleSchedule(dt, 'W', W);

fprintf('\nTotal simulation time : %.2f years\n', ...
    convertTo(sum(schedule10.step.val), year));

fprintf('Number of timesteps   : %d\n', ...
    numel(schedule10.step.val));

fprintf('\nRunning 10-year waterflood...\n');

%% Run simulation

[wellSols10, states10, reports10] = ...
    simulateScheduleAD(state0, model, schedule10);

fprintf('\n10-year simulation completed.\n');


%% ============================================================
%  EXTRACT WELL RESULTS
% =============================================================

nSteps = numel(wellSols10);

time_days = convertTo(cumsum(schedule10.step.val), day);
time_years = convertTo(cumsum(schedule10.step.val), year);

oilRate   = zeros(nSteps,1);
waterRate = zeros(nSteps,1);
totalRate = zeros(nSteps,1);

injRate1 = zeros(nSteps,1);
injRate2 = zeros(nSteps,1);
injRate3 = zeros(nSteps,1);
injRate4 = zeros(nSteps,1);

producerBHP = zeros(nSteps,1);

injBHP = zeros(nSteps,4);


for k = 1:nSteps

    ws = wellSols10{k};

    % --------------------------------------------
    % Producer P1
    % --------------------------------------------

    producerBHP(k) = convertTo(ws(1).bhp, barsa);

    oilRate(k) = ...
        -convertTo(ws(1).qOs, meter^3/day);

    waterRate(k) = ...
        -convertTo(ws(1).qWs, meter^3/day);

    totalRate(k) = ...
        -convertTo(ws(1).qTs, meter^3/day);


    % --------------------------------------------
    % Injectors
    % --------------------------------------------

    injRate1(k) = ...
        convertTo(ws(2).qWs, meter^3/day);

    injRate2(k) = ...
        convertTo(ws(3).qWs, meter^3/day);

    injRate3(k) = ...
        convertTo(ws(4).qWs, meter^3/day);

    injRate4(k) = ...
        convertTo(ws(5).qWs, meter^3/day);


    % Injector BHPs

    for n = 1:4
        injBHP(k,n) = ...
            convertTo(ws(n+1).bhp, barsa);
    end

end


%% ============================================================
%  CUMULATIVE PRODUCTION
% =============================================================

% Numerical integration of production rates

cumOil = cumtrapz( ...
    [0; time_days], ...
    [oilRate(1); oilRate]);

cumWater = cumtrapz( ...
    [0; time_days], ...
    [waterRate(1); waterRate]);


%% ============================================================
%  OOIP AND RECOVERY
% =============================================================

pv = poreVolume(G, rock);

OOIP_m3 = sum(pv .* state0.s(:,2));

finalCumOil = cumOil(end);

recoveryFactor = ...
    100 * finalCumOil / OOIP_m3;


%% ============================================================
%  FINAL RESERVOIR PROPERTIES
% =============================================================

finalState = states10{end};

P_initial = convertTo(mean(state0.pressure), barsa);

P_final = convertTo(mean(finalState.pressure), barsa);

Sw_initial = mean(state0.s(:,1));

Sw_final = mean(finalState.s(:,1));

So_initial = mean(state0.s(:,2));

So_final = mean(finalState.s(:,2));


%% ============================================================
%  WATER CUT
% =============================================================

waterCut = waterRate ./ ...
    max(oilRate + waterRate, eps);

finalWaterCut = waterCut(end) * 100;


%% ============================================================
%  INJECTION RESULTS
% =============================================================

totalInjectionRate = ...
    injRate1 + injRate2 + injRate3 + injRate4;

cumInjection = cumtrapz( ...
    [0; time_days], ...
    [totalInjectionRate(1); totalInjectionRate]);


%% ============================================================
%  PRINT FINAL RESULTS
% =============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('          10-YEAR WATERFLOOD RESULTS\n');
fprintf('====================================================\n');

fprintf('\n--- Production ---\n');

fprintf('Final oil rate       : %.4f m3/day\n', ...
    oilRate(end));

fprintf('Final water rate     : %.4f m3/day\n', ...
    waterRate(end));

fprintf('Final total rate     : %.4f m3/day\n', ...
    totalRate(end));

fprintf('Final water cut      : %.2f %%\n', ...
    finalWaterCut);

fprintf('Cumulative oil       : %.2f m3\n', ...
    finalCumOil);

fprintf('Cumulative water     : %.2f m3\n', ...
    cumWater(end));


fprintf('\n--- Recovery ---\n');

fprintf('OOIP                 : %.2f m3\n', ...
    OOIP_m3);

fprintf('Recovery factor      : %.4f %%\n', ...
    recoveryFactor);


fprintf('\n--- Pressure ---\n');

fprintf('Initial average P    : %.2f bar\n', ...
    P_initial);

fprintf('Final average P      : %.2f bar\n', ...
    P_final);

fprintf('Pressure change      : %.2f bar\n', ...
    P_final - P_initial);


fprintf('\n--- Saturation ---\n');

fprintf('Initial average Sw   : %.4f\n', ...
    Sw_initial);

fprintf('Final average Sw     : %.4f\n', ...
    Sw_final);

fprintf('Initial average So   : %.4f\n', ...
    So_initial);

fprintf('Final average So     : %.4f\n', ...
    So_final);


fprintf('\n--- Injection ---\n');

fprintf('Final total injection rate : %.2f m3/day\n', ...
    totalInjectionRate(end));

fprintf('Cumulative injection       : %.2f m3\n', ...
    cumInjection(end));


fprintf('\n--- Injector BHP ---\n');

fprintf('I1 maximum BHP : %.2f bar\n', max(injBHP(:,1)));

fprintf('I2 maximum BHP : %.2f bar\n', max(injBHP(:,2)));

fprintf('I3 maximum BHP : %.2f bar\n', max(injBHP(:,3)));

fprintf('I4 maximum BHP : %.2f bar\n', max(injBHP(:,4)));


fprintf('\n====================================================\n');
fprintf('       10-YEAR WATERFLOOD COMPLETED\n');
fprintf('====================================================\n');


%% ============================================================
%  PLOT 1 - OIL PRODUCTION RATE
% =============================================================

figure;

plot(time_years, oilRate, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Oil Rate (m^3/day)');
title('Waterflood Oil Production Rate');
grid on;


%% ============================================================
%  PLOT 2 - WATER PRODUCTION RATE
% =============================================================

figure;

plot(time_years, waterRate, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Water Rate (m^3/day)');
title('Waterflood Water Production Rate');
grid on;


%% ============================================================
%  PLOT 3 - WATER CUT
% =============================================================

figure;

plot(time_years, waterCut * 100, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Water Cut (%)');
title('Water Cut During Waterflood');
grid on;


%% ============================================================
%  PLOT 4 - CUMULATIVE OIL
% =============================================================

figure;

plot([0; time_years], cumOil, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Cumulative Oil (m^3)');
title('Cumulative Oil Production');
grid on;

%% ============================================================
%  PLOT 5 - AVERAGE RESERVOIR PRESSURE
% =============================================================

avgPressure = zeros(nSteps,1);

for k = 1:nSteps
    avgPressure(k) = ...
        convertTo(mean(states10{k}.pressure), barsa);
end

figure;

plot(time_years, avgPressure, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Average Reservoir Pressure (bar)');
title('Average Reservoir Pressure During Waterflood');
grid on;


%% ============================================================
%  PLOT 6 - INJECTOR BHP
% =============================================================

figure;

plot(time_years, injBHP, 'LineWidth', 1.2);

xlabel('Time (years)');
ylabel('Injector BHP (bar)');
title('Injector Bottomhole Pressure');

legend('I1','I2','I3','I4', ...
    'Location','best');

grid on;


%% ============================================================
%  PLOT 7 - FINAL WATER SATURATION
% =============================================================

figure;

plotCellData(G, finalState.s(:,1));

colorbar;

title('Final Water Saturation - 10 Year Waterflood');

view(3);
axis tight;

%% ============================================================
%  TOP-LAYER WATER SATURATION MAP
% =============================================================

Nx = G.cartDims(1);
Ny = G.cartDims(2);
Nz = G.cartDims(3);

Sw_final = finalState.s(:,1);

% Reshape cell data into reservoir dimensions
Sw_3D = reshape(Sw_final, [Nx, Ny, Nz]);

figure;

imagesc(Sw_3D(:,:,1)');

axis equal tight;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('X Grid Index');
ylabel('Y Grid Index');

title('Final Water Saturation - Top Layer');

hold on;

% Producer location
plot(iP, jP, 'ko', 'MarkerSize', 8, 'LineWidth', 2);

% Injector locations
plot(iI, jI, 'ks', 'MarkerSize', 8, 'LineWidth', 2);

text(iP+1, jP, 'P1');
text(iI(1)+1, jI(1), 'I1');
text(iI(2)+1, jI(2), 'I2');
text(iI(3)+1, jI(3), 'I3');
text(iI(4)+1, jI(4), 'I4');

hold off;
%% ============================================================
%  FINAL PRESSURE - TOP LAYER
% =============================================================

P_final_cells = convertTo(finalState.pressure, barsa);

P_3D = reshape(P_final_cells, [Nx, Ny, Nz]);

figure;

imagesc(P_3D(:,:,1)');

axis equal tight;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('X Grid Index');
ylabel('Y Grid Index');

title('Final Pressure Distribution - Top Layer');

hold on;

plot(iP, jP, 'ko', 'MarkerSize', 8, 'LineWidth', 2);

plot(iI, jI, 'ks', 'MarkerSize', 8, 'LineWidth', 2);

text(iP+1, jP, 'P1');
text(iI(1)+1, jI(1), 'I1');
text(iI(2)+1, jI(2), 'I2');
text(iI(3)+1, jI(3), 'I3');
text(iI(4)+1, jI(4), 'I4');

hold off;