%% Cairn-Rajasthan Inspired Polymer EOR Study
% Stage 3 (v03): Fluid Properties and Initial Reservoir State
% MODEL ASSUMPTIONS (polymer physics simplified)
%  - No polymer adsorption (ads = 0)
%  - No inaccessible pore volume (dps = 0)
%  - No residual resistance factor (rrf = 1)
%  - Linear viscosity law: muWMult = 1 + 4*c/cpmax (no Flory-Huggins curve)
%  - No shear thinning, polymer degradation, or salinity/temperature effects
%  - No injector BHP limit (fracture pressure not enforced)
%  - Results (incremental oil, PUF, injector BHP) are therefore optimistic
%    upper-bound estimates, suited to screening and not to forecasting.
close all;

%% Load MRST modules
mrstModule add ad-core ad-blackoil ad-props

%% -------------------------------
% 0. USER SETTINGS
% -------------------------------
% grid / pattern
Nx = 50;  Ny = 50;  Nz = 3;
Lx = 400*meter;  Ly = 400*meter;  Lz = 30*meter;   % 400 m -> 16 ha, try 600 m for a larger pattern

% rock
phiMean = 0.20;  phiStd = 0.02;  phiMin = 0.15;  phiMax = 0.25;
kLayer_mD = [400 100 30];      % horizontal perm of layers 1..Nz (mD)  <- set from your data
lnSigma   = 0.3;               % random noise (natural-log std) inside each layer
kvkh      = 0.1;               % vertical / horizontal permeability
cR        = 4.0e-10;           % rock compressibility (1/Pa)

% fluids
muW  = 1*centi*poise;   muO  = 10*centi*poise;
rhoW = 1000*kilogram/meter^3;   rhoO = 800*kilogram/meter^3;
cW   = 4.5e-10;   cO = 1.0e-9;           % compressibilities (1/Pa)

% relative permeability (Corey with end-points)
Swc = 0.25;   Sor = 0.25;
krwMax = 0.30;   kroMax = 0.80;
nw = 2;   no = 2;

% initial state
p0  = 250*barsa;
Sw0 = Swc;                 % initial water immobile
So0 = 1 - Sw0;

% polymer viscosity law used in the sensitivity script (only for the check below)
muWMultFun = @(c) 1 + 4*(c./1.0);

%% -------------------------------
% 1. Reservoir Grid
% -------------------------------
G = cartGrid([Nx, Ny, Nz], [Lx, Ly, Lz]);
G = computeGeometry(G);

%% -------------------------------
% 2. Rock Properties
% -------------------------------
rng(42);

phi = phiMean + phiStd*randn(G.cells.num, 1);
phi = min(max(phi, phiMin), phiMax);

if numel(kLayer_mD) ~= Nz
    error('kLayer_mD must have %d entries (one per layer).', Nz);
end
[~, ~, kIdx] = ind2sub(G.cartDims, (1:G.cells.num)');
kLay = kLayer_mD(:);
kh   = kLay(kIdx) .* exp(lnSigma*randn(G.cells.num, 1));
kh   = min(max(kh, 1), 2000);              % guard against extreme values (mD)
kv   = kvkh * kh;

rock = makeRock(G, [kh, kh, kv]*milli*darcy, phi);
rock.cr = cR/Pascal;

%% -------------------------------
% 3. Fluid Properties
% -------------------------------
fluid = initSimpleADIFluid('phases','WO', ...
    'mu',  [muW muO], ...
    'rho', [rhoW rhoO], ...
    'pRef', p0, ...
    'c',   [cW cO]/Pascal, ...
    'cR',  cR/Pascal);

% Override the default relperms with Corey curves WITH end-points.
% Check that your MRST version uses these field names: fieldnames(fluid)
if ~isfield(fluid,'krW') || ~isfield(fluid,'krO')
    disp(fieldnames(fluid));
    error('fluid has no krW/krO fields in this MRST version - adapt the override below.');
end
den = 1 - Swc - Sor;
fluid.krW = @(sw) krwMax * (max(sw - Swc, 0)/den).^nw;
fluid.krO = @(so) kroMax * (max(so - Sor, 0)/den).^no;

%% -------------------------------
% 4. Reservoir Model and Initial State
% -------------------------------
model  = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, p0, [Sw0, So0]);

%% -------------------------------
% 5. Sanity Checks
% -------------------------------
pv   = poreVolume(G, rock);
PV   = sum(pv);
OOIP = sum(pv .* state0.s(:,2));

dx = Lx/Nx;  dy = Ly/Ny;
dInjProd = hypot((25-4)*dx, (25-4)*dy);          % corner injector (4,4) to producer (25,25)

krw0 = fluid.krW(Sw0);   kro0 = fluid.krO(So0);
fw0  = (krw0/muW) / (krw0/muW + kro0/muO);

qTot = 4*50;                                      % m3/day, 4 injectors x 50 (as in the sensitivity script)
PVI10 = qTot*3650 / PV;

fprintf('\n--- Reservoir Model (v03) ---\n');
fprintf('Grid cells            : %d (%d x %d x %d)\n', G.cells.num, Nx, Ny, Nz);
fprintf('Cell size             : %.1f x %.1f x %.1f m\n', dx, dy, Lz/Nz);
fprintf('Pattern area          : %.1f ha | injector-producer distance %.0f m\n', Lx*Ly/1e4, dInjProd);
fprintf('Pore volume           : %.4e m3\n', PV);
fprintf('OOIP                  : %.2f m3   <- set expectedOOIP in the sensitivity script\n', OOIP);
fprintf('Initial pressure      : %.1f bar\n', convertTo(p0, barsa));
fprintf('Swc / Sor / Sw0       : %.2f / %.2f / %.2f\n', Swc, Sor, Sw0);
fprintf('Layer perm (mD)       : %s | kv/kh %.2f | mean kh %.1f mD\n', mat2str(kLayer_mD), kvkh, mean(kh));
fprintf('Oil / water viscosity : %.1f / %.1f cP\n', convertTo(muO, centi*poise), convertTo(muW, centi*poise));
fprintf('Oil / water density   : %.0f / %.0f kg/m3\n', convertTo(rhoO, kilogram/meter^3), convertTo(rhoW, kilogram/meter^3));
fprintf('Initial water cut     : %.2f %% (should be ~0)\n', 100*fw0);
fprintf('Water injected, 10 yr : %.2f PV at %d m3/day total\n', PVI10, qTot);
fprintf('End-point mobility ratio M = (krwMax/muW_eff)/(kroMax/muO):\n');
for c = [0 0.25 0.50 0.75 1.00]
    muEff = muW*muWMultFun(c);
    M = (krwMax/muEff) / (kroMax/muO);
    fprintf('   polymer %.2f kg/m3 -> muW_eff %.1f cP, M = %.2f\n', c, convertTo(muEff, centi*poise), M);
end
if fw0 > 0.01
    warning('Initial water cut > 1%%: check Sw0 vs Swc.');
end

%% -------------------------------
% 6. Plots
% -------------------------------
figure; plotCellData(G, state0.s(:,1)); view(3); axis tight equal; colorbar;
title('Initial Water Saturation');

figure; plotCellData(G, convertTo(state0.pressure, barsa)); view(3); axis tight equal; colorbar;
title('Initial Reservoir Pressure (bar)');

figure; plotCellData(G, log10(rock.perm(:,1)/(milli*darcy))); view(3); axis tight equal; colorbar;
title('log_{10} horizontal permeability (mD)');

sw = linspace(0, 1, 201)';
figure; hold on;
plot(sw, fluid.krW(sw), 'b', 'LineWidth', 1.6);
plot(sw, fluid.krO(1 - sw), 'r', 'LineWidth', 1.6);
grid on; xlabel('S_w'); ylabel('Relative permeability'); legend('k_{rw}','k_{ro}');
title('Relative permeability curves');

figure; hold on;
for c = [0 0.25 0.50 0.75 1.00]
    muEff = muW*muWMultFun(c);
    fw = (fluid.krW(sw)/muEff) ./ (fluid.krW(sw)/muEff + fluid.krO(1 - sw)/muO + eps);
    plot(sw, fw, 'LineWidth', 1.5, 'DisplayName', sprintf('%.2f kg/m^3', c));
end
grid on; xlabel('S_w'); ylabel('Fractional flow f_w'); legend('Location','southeast');
title('Fractional flow vs polymer concentration');