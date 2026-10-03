%% Cairn-Rajasthan Inspired Polymer EOR Study
% Stage 3: Fluid Properties and Initial Reservoir State

clear;
clc;
close all;

%% Load MRST modules
mrstModule add ad-core ad-blackoil ad-props

%% -------------------------------
% 1. Reservoir Grid
% -------------------------------

Nx = 50;
Ny = 50;
Nz = 3;

Lx = 1200 * meter;
Ly = 1200 * meter;
Lz = 30 * meter;

G = cartGrid([Nx, Ny, Nz], [Lx, Ly, Lz]);
G = computeGeometry(G);

%% -------------------------------
% 2. Rock Properties
% -------------------------------

rng(42);

phi = 0.20 + 0.02 * randn(G.cells.num, 1);
phi = min(max(phi, 0.15), 0.25);

kBase = 100;

perm_mD = exp(log(kBase) + 0.8 * randn(G.cells.num, 1));
perm_mD = min(max(perm_mD, 10), 500);

rock = makeRock(G, ...
    perm_mD * milli * darcy, ...
    phi);

rock.cr = 4.0e-10/Pascal;

%% -------------------------------
% 3. Fluid Properties
% -------------------------------
muW = 1*centi*poise;
muO = 10*centi*poise;

rhoW = 1000*kilogram/meter^3;
rhoO = 800*kilogram/meter^3;

p0 = 250*barsa;

Sw0 = 0.25;
So0 = 0.75;

fluid = initSimpleADIFluid('phases','WO', ...
    'mu',[1 10]*centi*poise, ...
    'rho',[1000 800]*kilogram/meter^3, ...
    'pRef',250*barsa, ...
    'c',[4.5e-10 1.0e-9]/Pascal, ...
    'cR',4.0e-10/Pascal, ...
    'n',[2 2]);
%% -------------------------------
% 4. Reservoir Model
% -------------------------------

model = TwoPhaseOilWaterModel(G,rock,fluid);

state0 = initResSol(G,p0,[Sw0 So0]);
%% -------------------------------
% 5. Initial Reservoir Conditions
% -------------------------------

p0 = 250 * barsa;

Sw0 = 0.25;
So0 = 0.75;

state0 = initResSol(G, p0, [Sw0, So0]);

%% -------------------------------
% 6. Display Initial Conditions
% -------------------------------

fprintf('\n--- Reservoir Model ---\n');

fprintf('Grid cells: %d\n', G.cells.num);

fprintf('Initial pressure: %.1f bar\n', ...
    convertTo(p0, barsa));

fprintf('Initial water saturation: %.2f\n', Sw0);

fprintf('Initial oil saturation: %.2f\n', So0);

fprintf('Oil viscosity: %.1f cP\n', ...
    convertTo(muO, centi*poise));

fprintf('Water viscosity: %.1f cP\n', ...
    convertTo(muW, centi*poise));

fprintf('Oil density: %.0f kg/m3\n', ...
    convertTo(rhoO, kilogram/meter^3));

fprintf('Water density: %.0f kg/m3\n', ...
    convertTo(rhoW, kilogram/meter^3));

%% -------------------------------
% 7. Initial Saturation Plot
% -------------------------------

figure;

plotCellData(G, state0.s(:,1));

view(3);
axis tight equal;
colorbar;

title('Initial Water Saturation');

%% -------------------------------
% 8. Initial Pressure Plot
% -------------------------------

figure;

plotCellData(G, convertTo(state0.pressure, barsa));

view(3);
axis tight equal;
colorbar;

title('Initial Reservoir Pressure (bar)');
