%% Cairn-Rajasthan Inspired Polymer EOR Study
% Stage 2: Synthetic Reservoir Model

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

fprintf('Number of cells: %d\n', G.cells.num);

%% -------------------------------
% 2. Porosity
% -------------------------------

rng(42);

phi = 0.20 + 0.02 * randn(G.cells.num, 1);

% Keep porosity within realistic synthetic limits
phi = min(max(phi, 0.15), 0.25);

%% -------------------------------
% 3. Permeability
% -------------------------------

kBase = 100;       % mD

perm_mD = exp(log(kBase) + 0.8 * randn(G.cells.num, 1));

% Limit permeability range
perm_mD = min(max(perm_mD, 10), 500);

%% -------------------------------
% 4. Create MRST Rock Object
% -------------------------------

rock = makeRock( ...
    G, ...
    perm_mD * milli * darcy, ...
    phi);

%% -------------------------------
% 5. Display Reservoir Properties
% -------------------------------

fprintf('Average porosity: %.3f\n', mean(phi));
fprintf('Minimum porosity: %.3f\n', min(phi));
fprintf('Maximum porosity: %.3f\n', max(phi));

fprintf('Average permeability: %.2f mD\n', mean(perm_mD));
fprintf('Minimum permeability: %.2f mD\n', min(perm_mD));
fprintf('Maximum permeability: %.2f mD\n', max(perm_mD));

%% -------------------------------
% 6. Plot Grid
% -------------------------------

figure;

plotGrid(G);

view(3);
axis tight equal;
title('Synthetic Reservoir Grid');

%% -------------------------------
% 7. Plot Porosity
% -------------------------------

figure;

plotCellData(G, phi);

view(3);
axis tight equal;
colorbar;

title('Synthetic Reservoir Porosity');

%% -------------------------------
% 8. Plot Permeability
% -------------------------------

figure;

plotCellData(G, log10(perm_mD));

view(3);
axis tight equal;
colorbar;

title('Log_{10} Permeability Distribution (mD)');
