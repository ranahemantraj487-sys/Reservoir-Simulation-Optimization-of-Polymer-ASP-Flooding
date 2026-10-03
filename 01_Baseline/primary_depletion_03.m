%% Cairn-Rajasthan Inspired Polymer EOR Study
% Stage 4: Primary Depletion

clear;
clc;
close all;

mrstModule add ad-core ad-blackoil ad-props

%% 1. Load Reservoir and Fluid Model

run('fluid_model_02.m');

%% 2. Producer Location

iP = round(Nx/2);
jP = round(Ny/2);

W = verticalWell( ...
    [], G, rock, ...
    iP, jP, 1:Nz, ...
    'Name', 'P1', ...
    'Type', 'bhp', ...
    'Val', 150 * barsa, ...
    'Radius', 0.1 * meter, ...
    'Comp_i', [0 1]);

%% 3. Simulation Time

simTime = 10 * year;

dt = rampupTimesteps( ...
    simTime, ...
    simTime/100);

schedule = simpleSchedule(dt, 'W', W);

%% 4. Run Simulation

fprintf('\nRunning primary depletion simulation...\n');

[wellSols, states, reports] = ...
    simulateScheduleAD( ...
    state0, model, schedule);

fprintf('Simulation completed.\n');

%% 5. Extract Well Production Data

nSteps = numel(wellSols);

time_years = zeros(nSteps,1);
oilRate = zeros(nSteps,1);
waterRate = zeros(nSteps,1);
bhp = zeros(nSteps,1);

for k = 1:nSteps

    time_years(k) = ...
        convertTo(schedule.step.val(k), year);

    oilRate(k) = ...
        -convertTo(wellSols{k}(1).qOs, meter^3/day);

    waterRate(k) = ...
        -convertTo(wellSols{k}(1).qWs, meter^3/day);

    bhp(k) = ...
        convertTo(wellSols{k}(1).bhp, barsa);

end

%% 6. Cumulative Oil Production

dt_days = convertTo(schedule.step.val, day);
dt_days = dt_days(:);

cumOil = zeros(nSteps,1);

for k = 2:nSteps

    cumOil(k) = cumOil(k-1) + ...
        0.5 * (oilRate(k) + oilRate(k-1)) * ...
        (dt_days(k) - dt_days(k-1));

end

%% 7. OOIP

pv = poreVolume(G, rock);

OOIP_m3 = sum(pv .* state0.s(:,2));

recoveryFactor = ...
    100 * cumOil(end) / OOIP_m3;

%% 8. Results

fprintf('\n--- Primary Depletion Results ---\n');

fprintf('OOIP: %.2f m3\n', OOIP_m3);

fprintf('Cumulative oil production: %.2f m3\n', ...
    cumOil(end));

fprintf('Primary depletion recovery: %.2f %%\n', ...
    recoveryFactor);

fprintf('Initial average pressure: %.2f bar\n', ...
    mean(state0.pressure)/barsa);

fprintf('Final average pressure: %.2f bar\n', ...
    mean(states{end}.pressure)/barsa);

fprintf('Final oil rate: %.2f m3/day\n', ...
    oilRate(end));

fprintf('Final water rate: %.2f m3/day\n', ...
    waterRate(end));

%% 9. Final Pressure Distribution

figure;

plotCellData(G, ...
    convertTo(states{end}.pressure, barsa));

view(3);
axis tight equal;
colorbar;

title('Final Reservoir Pressure - Primary Depletion');

%% 10. Oil Production Rate

figure;

plot(time_years, oilRate, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Oil rate (m^3/day)');
title('Oil Production Rate - Primary Depletion');

grid on;

%% 11. Water Production Rate

figure;

plot(time_years, waterRate, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Water rate (m^3/day)');
title('Water Production Rate - Primary Depletion');

grid on;

%% 12. BHP

figure;

plot(time_years, bhp, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('BHP (bar)');
title('Bottomhole Pressure - Primary Depletion');

grid on;

%% 13. Cumulative Oil

figure;

plot(time_years, cumOil, 'LineWidth', 1.5);

xlabel('Time (years)');
ylabel('Cumulative oil (m^3)');
title('Cumulative Oil Production - Primary Depletion');

grid on;