%% ============================================================
% Cairn Rajasthan Inspired Polymer EOR Study
% Stage 6: Polymer Flood - 1 Day Validation
%
% MRST 2026a
% MATLAB R2026b
%
% Model:
%   50 x 50 x 3 layered heterogeneous reservoir (fluid_model_03.m)
%   Initial pressure / saturations taken from fluid_model_03.m
%   5-spot water/polymer flood
%
% Polymer:
%   1 kg/m3 injected concentration
%   Todd-Longstaff mixing parameter = 0.85
%   Polymer viscosity multiplier:
%       muW_mult = 1 + 4*C/Cmax
%
% Therefore:
%       C = 0       -> 1x
%       C = 1 kg/m3 -> 5x
%
% IMPORTANT:
%   This file performs ONLY a 1-day validation.
%   Do not start the long polymer simulation until
%   this validation completes successfully.
%% ============================================================

clear;
clc;
close all;


%% ------------------------------------------------------------
% 1. LOAD MRST MODULES
% ------------------------------------------------------------

mrstModule add ad-core ad-blackoil ad-eor ad-props;


%% ------------------------------------------------------------
% 2. LOAD COMMON RESERVOIR MODEL
% ------------------------------------------------------------

projectRoot = fileparts(fileparts(mfilename('fullpath')));

run(fullfile(projectRoot, ...
             '01_Baseline', ...
             'fluid_model_03.m'));

close all;   % fluid_model_03.m opens plot windows

% Well positions below assume the 50 x 50 x 3 grid of fluid_model_03.m
if G.cartDims(1) ~= 50 || G.cartDims(2) ~= 50 || G.cartDims(3) ~= 3
    error('Baseline grid is %d x %d x %d; well positions assume 50 x 50 x 3.', G.cartDims);
end


%% ------------------------------------------------------------
% 3. DISPLAY RESERVOIR INFORMATION
% ------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' POLYMER FLOOD - RESERVOIR MODEL\n');
fprintf('============================================================\n');

fprintf('Grid cells          : %d\n', G.cells.num);

fprintf('Grid dimensions     : %d x %d x %d\n', ...
        Nx, Ny, Nz);

fprintf('Initial pressure    : %.2f bar\n', ...
        convertTo(p0, barsa));

fprintf('Initial Sw          : %.3f\n', Sw0);
fprintf('Initial So          : %.3f\n', So0);

fprintf('Oil viscosity       : %.2f cP\n', ...
        convertTo(muO, centi*poise));

fprintf('Water viscosity     : %.2f cP\n', ...
        convertTo(muW, centi*poise));

fprintf('============================================================\n');


%% ------------------------------------------------------------
% 4. POLYMER FLUID PROPERTIES
% ------------------------------------------------------------

% Maximum polymer concentration
fluid.cpmax = 1.0;

% Todd-Longstaff mixing parameter
fluid.mixPar = 0.85;

% Polymer-induced water viscosity multiplier
%
% C = 0 kg/m3 -> multiplier = 1
% C = 1 kg/m3 -> multiplier = 5

fluid.muWMult = @(c) ...
    1 + 4*(c/fluid.cpmax);


%% ------------------------------------------------------------
% 5. POLYMER ADSORPTION / RETENTION PARAMETERS
% ------------------------------------------------------------

% adsInx = 1:
% adsorption is evaluated as fluid.ads(cp)

fluid.adsInx = 1;

% Maximum adsorption capacity
fluid.adsMax = 1.0;

% Zero adsorption for initial validation
fluid.ads = @(c) zeros(size(c));

% Resistance reduction factor
fluid.rrf = 1.0;

% Inaccessible pore space fraction
fluid.dps = 0.0;

% Rock grain density
fluid.rhoR = 2650*kilogram/meter^3;


%% ------------------------------------------------------------
% 6. CREATE POLYMER MODEL
% ------------------------------------------------------------

modelPolymer = OilWaterPolymerModel( ...
    G, ...
    rock, ...
    fluid);

fprintf('\nPolymer model created successfully.\n');


%% ------------------------------------------------------------
% 7. INITIAL POLYMER STATE
% ------------------------------------------------------------

% Start with the same reservoir state as baseline
state0Polymer = state0;

% Initial polymer concentration in reservoir cells
state0Polymer.cp = zeros(G.cells.num, 1);

% Maximum polymer concentration experienced by each cell
state0Polymer.cpmax = zeros(G.cells.num, 1);


%% ------------------------------------------------------------
% 8. WELL LOCATIONS
% ------------------------------------------------------------

% Producer at reservoir centre
iP = 25;
jP = 25;

% Four corner injectors
iI = [4  4  47 47];
jI = [4 47   4 47];


%% ------------------------------------------------------------
% 9. CREATE WELLS
% ------------------------------------------------------------

% Start with an empty well list
W = [];


%% ------------------------------------------------------------
% Producer P1
% BHP-controlled producer
% ------------------------------------------------------------

cellP = sub2ind( ...
    [Nx, Ny, Nz], ...
    iP, ...
    jP, ...
    1:Nz);

W = addWell( ...
    W, ...
    G, ...
    rock, ...
    cellP, ...
    'Type', 'bhp', ...
    'Val', 150*barsa, ...
    'Name', 'P1', ...
    'Comp_i', [0 1]);


%% ------------------------------------------------------------
% Four polymer injectors
% ------------------------------------------------------------

for n = 1:4

    cellI = sub2ind( ...
        [Nx, Ny, Nz], ...
        iI(n), ...
        jI(n), ...
        1:Nz);

    W = addWell( ...
        W, ...
        G, ...
        rock, ...
        cellI, ...
        'Type', 'rate', ...
        'Val', 50*meter^3/day, ...
        'Name', sprintf('I%d', n), ...
        'Comp_i', [1 0]);

end


%% ------------------------------------------------------------
% 10. POLYMER CONCENTRATION AT WELLS
% ------------------------------------------------------------
%
% IMPORTANT:
%
% MRST's OilWaterPolymerModel maps the polymer variable
% 'polymer' to the state/well field 'cp'.
%
% Therefore use:
%
%       W(i).cp
%
% NOT:
%
%       W(i).polymer
%
% Producer:
%       cp = 0 kg/m3
%
% Polymer injectors:
%       cp = 1 kg/m3
% ------------------------------------------------------------

% Producer
W(1).cp = 0.0;

% Four polymer injectors
for n = 2:numel(W)
    W(n).cp = 1.0;
end


%% ------------------------------------------------------------
% 11. VERIFY WELL POLYMER CONCENTRATIONS
% ------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' WELL CONFIGURATION\n');
fprintf('============================================================\n');

for n = 1:numel(W)

    fprintf('\nWell: %s\n', W(n).name);

    fprintf('  Type       : %s\n', W(n).type);

    fprintf('  Polymer cp : %.2f kg/m3\n', W(n).cp);

    if strcmpi(W(n).type, 'bhp')

        fprintf('  BHP        : %.2f bar\n', ...
                convertTo(W(n).val, barsa));

    else

        fprintf('  Rate       : %.2f m3/day\n', ...
                convertTo(W(n).val, meter^3/day));

    end

end

fprintf('\n============================================================\n');


%% ------------------------------------------------------------
% 12. VERIFY MRST POLYMER WELL PROPERTY
% ------------------------------------------------------------

fprintf('\nChecking MRST polymer well properties...\n');

for n = 1:numel(W)

    cpCheck = modelPolymer.getProp(W(n), 'polymer');

    fprintf('  %s -> polymer = %.2f kg/m3\n', ...
            W(n).name, cpCheck);

end

fprintf('MRST polymer well-property check passed.\n');


%% ------------------------------------------------------------
% 13. CREATE 1-DAY VALIDATION SCHEDULE
% ------------------------------------------------------------

dt = 1*day;

schedule = simpleSchedule( ...
    dt, ...
    'W', W);


%% ------------------------------------------------------------
% 14. RUN 1-DAY POLYMER VALIDATION
% ------------------------------------------------------------

fprintf('\n');
fprintf('============================================\n');
fprintf(' Running 1-day polymer validation...\n');
fprintf('============================================\n\n');

[wellSols, states, reports] = ...
    simulateScheduleAD( ...
        state0Polymer, ...
        modelPolymer, ...
        schedule);


%% ------------------------------------------------------------
% 15. EXTRACT FINAL STATE
% ------------------------------------------------------------

finalState = states{end};
finalWellSol = wellSols{end};


%% ------------------------------------------------------------
% 16. PRODUCER RESULTS
% ------------------------------------------------------------

P1 = finalWellSol(1);

qOil = -P1.qOs;
qWater = -P1.qWs;

qTotal = qOil + qWater;


%% ------------------------------------------------------------
% 17. INJECTOR RESULTS
% ------------------------------------------------------------

qInjection = zeros(4,1);
injectorBHP = zeros(4,1);

for n = 1:4

    WI = finalWellSol(n+1);

    qInjection(n) = WI.qWs;

    injectorBHP(n) = WI.bhp;

end


%% ------------------------------------------------------------
% 18. POLYMER CONCENTRATION RESULTS
% ------------------------------------------------------------

cpFinal = finalState.cp;

averageCp = mean(cpFinal);

maximumCp = max(cpFinal);


%% ------------------------------------------------------------
% 19. PRESSURE RESULTS
% ------------------------------------------------------------

initialPressure = mean(state0Polymer.pressure);

finalPressure = mean(finalState.pressure);

pressureChange = finalPressure - initialPressure;


%% ------------------------------------------------------------
% 20. SATURATION RESULTS
% ------------------------------------------------------------

initialSw = mean(state0Polymer.s(:,1));

finalSw = mean(finalState.s(:,1));

initialSo = mean(state0Polymer.s(:,2));

finalSo = mean(finalState.s(:,2));


%% ------------------------------------------------------------
% 21. WATER CUT
% ------------------------------------------------------------

if qTotal > 0

    waterCut = qWater/qTotal;

else

    waterCut = NaN;

end


%% ------------------------------------------------------------
% 22. DISPLAY RESULTS
% ------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' 1-DAY POLYMER VALIDATION RESULTS\n');
fprintf('============================================================\n');


fprintf('\n--- Producer P1 ---\n');

fprintf('Oil rate             : %.4f m3/day\n', qOil);

fprintf('Water rate           : %.4f m3/day\n', qWater);

fprintf('Total rate           : %.4f m3/day\n', qTotal);

fprintf('Water cut            : %.2f %%\n', waterCut*100);

fprintf('BHP                  : %.2f bar\n', ...
        convertTo(P1.bhp, barsa));


fprintf('\n--- Polymer ---\n');

fprintf('Average polymer      : %.6f kg/m3\n', averageCp);

fprintf('Maximum polymer      : %.6f kg/m3\n', maximumCp);


fprintf('\n--- Reservoir Pressure ---\n');

fprintf('Initial average P    : %.4f bar\n', ...
        convertTo(initialPressure, barsa));

fprintf('Final average P      : %.4f bar\n', ...
        convertTo(finalPressure, barsa));

fprintf('Pressure change      : %.4f bar\n', ...
        convertTo(pressureChange, barsa));


fprintf('\n--- Saturation ---\n');

fprintf('Initial average Sw   : %.4f\n', initialSw);

fprintf('Final average Sw     : %.4f\n', finalSw);

fprintf('Initial average So   : %.4f\n', initialSo);

fprintf('Final average So     : %.4f\n', finalSo);


fprintf('\n--- Polymer Injection ---\n');

for n = 1:4

    fprintf('I%d rate              : %.4f m3/day\n', ...
            n, ...
            convertTo(qInjection(n), meter^3/day));

    fprintf('I%d BHP               : %.2f bar\n', ...
            n, ...
            convertTo(injectorBHP(n), barsa));

end


fprintf('\n============================================================\n');
fprintf(' Validation simulation completed successfully.\n');
fprintf('============================================================\n');


%% ------------------------------------------------------------
% 23. PLOT FINAL WATER SATURATION
% ------------------------------------------------------------

Sw_final = finalState.s(:,1);

Sw_3D = reshape( ...
    Sw_final, ...
    [Nx, Ny, Nz]);


figure;

imagesc(Sw_3D(:,:,1)');

axis equal tight;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('X Grid Index');
ylabel('Y Grid Index');

title('1-Day Polymer Flood - Water Saturation (Top Layer)');

hold on;


% Producer
plot( ...
    iP, ...
    jP, ...
    'ko', ...
    'MarkerSize', 9, ...
    'LineWidth', 2);

text( ...
    iP+1, ...
    jP, ...
    'P1', ...
    'FontWeight', 'bold');


% Injectors
plot( ...
    iI, ...
    jI, ...
    'ks', ...
    'MarkerSize', 8, ...
    'LineWidth', 2);

for n = 1:4

    text( ...
        iI(n)+1, ...
        jI(n), ...
        sprintf('I%d', n), ...
        'FontWeight', 'bold');

end

hold off;


%% ------------------------------------------------------------
% 24. PLOT POLYMER CONCENTRATION
% ------------------------------------------------------------

cp_3D = reshape( ...
    cpFinal, ...
    [Nx, Ny, Nz]);


figure;

imagesc(cp_3D(:,:,1)');

axis equal tight;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('X Grid Index');
ylabel('Y Grid Index');

title('1-Day Polymer Concentration - Top Layer (kg/m3)');

hold on;


% Producer
plot( ...
    iP, ...
    jP, ...
    'ko', ...
    'MarkerSize', 9, ...
    'LineWidth', 2);


% Injectors
plot( ...
    iI, ...
    jI, ...
    'ks', ...
    'MarkerSize', 8, ...
    'LineWidth', 2);


text( ...
    iP+1, ...
    jP, ...
    'P1', ...
    'FontWeight', 'bold');


for n = 1:4

    text( ...
        iI(n)+1, ...
        jI(n), ...
        sprintf('I%d', n), ...
        'FontWeight', 'bold');

end

hold off;


%% ------------------------------------------------------------
% 25. PLOT FINAL PRESSURE
% ------------------------------------------------------------

P_final_cells = convertTo( ...
    finalState.pressure, ...
    barsa);

P_3D = reshape( ...
    P_final_cells, ...
    [Nx, Ny, Nz]);


figure;

imagesc(P_3D(:,:,1)');

axis equal tight;

set(gca, 'YDir', 'normal');

colorbar;

xlabel('X Grid Index');
ylabel('Y Grid Index');

title('1-Day Polymer Flood - Pressure Distribution (bar)');

hold on;


% Producer
plot( ...
    iP, ...
    jP, ...
    'ko', ...
    'MarkerSize', 9, ...
    'LineWidth', 2);


% Injectors
plot( ...
    iI, ...
    jI, ...
    'ks', ...
    'MarkerSize', 8, ...
    'LineWidth', 2);


text( ...
    iP+1, ...
    jP, ...
    'P1', ...
    'FontWeight', 'bold');


for n = 1:4

    text( ...
        iI(n)+1, ...
        jI(n), ...
        sprintf('I%d', n), ...
        'FontWeight', 'bold');

end

hold off;

%% ------------------------------------------------------------
% 26. POLYMER MASS-BALANCE CHECK
% ------------------------------------------------------------

totalInjectionRate = sum(qInjection);

% Injected polymer concentration = 1 kg/m3
polymerInjectedExpected = ...
    convertTo(totalInjectionRate, meter^3/day) * 1.0 * 1.0;

% Approximate polymer mass currently dissolved in reservoir water
poreVolume = G.cells.volumes .* rock.poro;

polymerMassInWater = sum( ...
    poreVolume .* finalState.s(:,1) .* cpFinal);

fprintf('\n');
fprintf('============================================================\n');
fprintf(' POLYMER MASS-BALANCE CHECK\n');
fprintf('============================================================\n');

fprintf('Total water injection : %.4f m3/day\n', ...
    convertTo(totalInjectionRate, meter^3/day));

fprintf('Expected polymer      : %.4f kg\n', ...
    polymerInjectedExpected);

fprintf('Polymer in water phase: %.4f kg\n', ...
    polymerMassInWater);

fprintf('============================================================\n');


%% ============================================================
% END OF STAGE 6 - 1 DAY POLYMER VALIDATION
% ============================================================
%% ------------------------------------------------------------