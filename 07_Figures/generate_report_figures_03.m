%% generate_report_figures_03.m
% Cairn-Rajasthan inspired polymer flood - REPORT FIGURES (v03 model)
% Figures
%   R01  permeability maps (3 layers)        R08  final Sw: waterflood vs polymer
%   R02  porosity / permeability histograms  R09  final polymer concentration
%   R03  relperm + fractional flow curves    R10  final pressure: waterflood vs polymer
%   R04  5-spot well pattern                 R11  oil saturation gain (polymer - waterflood)
%   R05  production history (4 panels)       R12  polymer inventory + balance error
%   R06  recovery factor vs PV injected      R13  summary vs polymer concentration
%   R07  injector BHP history

clearvars; close all; clc;

%% 0. PATHS ---------------------------------------------------------------
projectRoot  = 'D:\PROJECT';
baselineFile = fullfile(projectRoot, '01_Baseline', 'fluid_model_03.m');
resultDir    = fullfile(projectRoot, '06_results');
figDir       = fullfile(projectRoot, '07_Figures', 'v03', 'report');
readmeDir    = fullfile(projectRoot, '08_Readme');
sensFile     = fullfile(resultDir, 'polymer_concentration_sensitivity_v03.mat');
ecoFile      = fullfile(resultDir, 'polymer_economics_v03.mat');

if ~isfile(baselineFile), error('Baseline file not found:\n%s', baselineFile); end
if ~isfile(sensFile),     error('Sensitivity results not found:\n%s', sensFile); end
if ~isfolder(figDir),     mkdir(figDir);    end
if ~isfolder(readmeDir),  mkdir(readmeDir); end

%% 1. LOAD BASELINE MODEL (same rock as the simulations: rng(42)) ----------
run(baselineFile);
close all;                                   % baseline script opens plot windows

Nx = G.cartDims(1);  Ny = G.cartDims(2);  Nz = G.cartDims(3);
dims = [Nx Ny Nz];
xm = ((1:Nx) - 0.5) * (Lx/Nx);               % cell-centre coordinates (m)
ym = ((1:Ny) - 0.5) * (Ly/Ny);

iP = 25;  jP = 25;
iI = [4 4 47 47];  jI = [4 47 4 47];
wells.px = xm(iP);  wells.py = ym(jP);
wells.ix = xm(iI);  wells.iy = ym(jI);

kmD = rock.perm(:,1) / (milli*darcy);
[~, ~, layerIdx] = ind2sub(G.cartDims, (1:G.cells.num)');

%% 2. LOAD RESULTS ---------------------------------------------------------
S      = load(sensFile);
Cc     = S.polymerConcentrations(:);
Rs     = S.allResults;
nCases = numel(Cc);
iCtl   = find(Cc == 0, 1);
if isempty(iCtl), error('0 kg/m3 waterflood control missing from the results.'); end

cand = find(Cc > 0);
if isempty(cand), error('No polymer case in the results.'); end
bestIdx = NaN;
if isfile(ecoFile)
    Ee = load(ecoFile, 'bestIdx');
    if isfield(Ee, 'bestIdx'), bestIdx = Ee.bestIdx; end
end
if isnan(bestIdx)
    [~, j] = max(cellfun(@(r) r.RF, Rs(cand)));
    bestIdx = cand(j);
end
fprintf('Best polymer case used for the maps: %.2f kg/m3\n', Cc(bestIdx));

%% 3. STYLE ----------------------------------------------------------------
set(groot, 'defaultAxesFontSize',12, 'defaultLineLineWidth',1.6, 'defaultFigureColor','w');
cols = lines(nCases);
cols(iCtl,:) = [0 0 0];                       % waterflood control in black
labs = arrayfun(@(c) sprintf('%.2f kg/m^3', c), Cc, 'UniformOutput', false);
labs{iCtl} = 'Waterflood (0)';
tYr = @(r) r.timeDays/365;

figList = cell(0,2);

%% R01 - permeability maps -------------------------------------------------
K3 = reshape(log10(kmD), dims);
fig = figure('Position',[60 60 1400 430]);
tl = tiledlayout(1, Nz, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Horizontal permeability by layer');
clim = [min(K3(:)) max(K3(:))];
for k = 1:Nz
    nexttile;  imagesc(xm, ym, K3(:,:,k)');  set(gca,'YDir','normal');  axis equal tight;
    caxis(clim);  hold on;  drawWells(wells);
    title(sprintf('Layer %d (mean %.0f mD)', k, mean(kmD(layerIdx == k))));
    xlabel('x (m)');  if k == 1, ylabel('y (m)'); end
end
cb = colorbar;  cb.Layout.Tile = 'east';  cb.Label.String = 'log_{10} k_h (mD)';
figList = reg(figList, fig, figDir, 'R01_Permeability_Layers', ...
    'Horizontal permeability (log10 mD) of the three layers with well positions (triangles = injectors, circle = producer).');

%% R02 - rock property distributions ---------------------------------------
fig = figure('Position',[60 60 1200 420]);
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');
nexttile;  histogram(rock.poro, 30, 'FaceColor',[0.2 0.4 0.8]);
grid on;  xlabel('Porosity (-)');  ylabel('Number of cells');  title('Porosity distribution');
nexttile;  hold on;
for k = 1:Nz
    histogram(log10(kmD(layerIdx == k)), 20, 'DisplayName', sprintf('Layer %d', k));
end
grid on;  xlabel('log_{10} k_h (mD)');  ylabel('Number of cells');
title('Permeability distribution by layer');  legend('Location','best');
figList = reg(figList, fig, figDir, 'R02_Rock_Distributions', ...
    'Histograms of porosity and of horizontal permeability for each layer.');

%% R03 - relative permeability and fractional flow -------------------------
swv = linspace(0, 1, 201)';
fig = figure('Position',[60 60 1200 440]);
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');
nexttile;  hold on;
plot(swv, fluid.krW(swv), 'b', 'DisplayName','k_{rw}');
plot(swv, fluid.krO(1 - swv), 'r', 'DisplayName','k_{ro}');
grid on;  xlabel('S_w');  ylabel('Relative permeability');  title('Corey relative permeability');
legend('Location','best');
nexttile;  hold on;
for cc = [0 0.25 0.50 0.75 1.00]
    muEff = muW * muWMultFun(cc);
    fwc = (fluid.krW(swv)/muEff) ./ (fluid.krW(swv)/muEff + fluid.krO(1 - swv)/muO + eps);
    plot(swv, fwc, 'DisplayName', sprintf('%.2f kg/m^3', cc));
end
grid on;  xlabel('S_w');  ylabel('Fractional flow f_w');
title('Fractional flow vs polymer concentration');  legend('Location','southeast');
figList = reg(figList, fig, figDir, 'R03_RelPerm_FractionalFlow', ...
    'Corey relative permeability curves and the fractional-flow curves obtained for each polymer concentration (linear viscosity law).');

%% R04 - well pattern ------------------------------------------------------
fig = figure('Position',[60 60 760 700]);
imagesc(xm, ym, K3(:,:,1)');  set(gca,'YDir','normal');  axis equal tight;  hold on;
for n = 1:4
    plot([wells.ix(n) wells.px], [wells.iy(n) wells.py], 'w--');
    text(wells.ix(n)+8, wells.iy(n)+8, sprintf('I%d', n), 'Color','w', 'FontWeight','bold');
end
drawWells(wells);
text(wells.px+8, wells.py+8, 'P1', 'Color','w', 'FontWeight','bold');
cb = colorbar;  cb.Label.String = 'log_{10} k_h (mD), layer 1';
xlabel('x (m)');  ylabel('y (m)');
title(sprintf('5-spot pattern (injector-producer distance %.0f m)', ...
    hypot((iP - iI(1))*Lx/Nx, (jP - jI(1))*Ly/Ny)));
figList = reg(figList, fig, figDir, 'R04_Well_Pattern', ...
    '5-spot well pattern: four corner injectors (50 m3/day each) and a central BHP-controlled producer, over the layer-1 permeability map.');

%% R05 - production history ------------------------------------------------
fig = figure('Position',[60 60 1300 850]);
tl = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Production performance vs polymer concentration');
nexttile;  hold on;
for k = 1:nCases, plot(tYr(Rs{k}), Rs{k}.oilRate, 'Color',cols(k,:), 'DisplayName',labs{k}); end
grid on;  xlabel('Time (years)');  ylabel('Oil rate (m^3/day)');  title('Oil rate');  legend('Location','best');
nexttile;  hold on;
for k = 1:nCases, plot(tYr(Rs{k}), 100*Rs{k}.waterCut, 'Color',cols(k,:)); end
grid on;  xlabel('Time (years)');  ylabel('Water cut (%)');  title('Water cut');
nexttile;  hold on;
for k = 1:nCases, plot(tYr(Rs{k}), Rs{k}.cumOilHist/1000, 'Color',cols(k,:)); end
grid on;  xlabel('Time (years)');  ylabel('Cumulative oil (10^3 m^3)');  title('Cumulative oil');
nexttile;  hold on;
for k = 1:nCases, plot(tYr(Rs{k}), Rs{k}.avgPressure, 'Color',cols(k,:)); end
grid on;  xlabel('Time (years)');  ylabel('Average pressure (bar)');  title('Average reservoir pressure');
figList = reg(figList, fig, figDir, 'R05_Production_History', ...
    'Oil rate, water cut, cumulative oil and average reservoir pressure for the waterflood control and each polymer concentration.');

%% R06 - recovery factor vs PV injected ------------------------------------
fig = figure('Position',[60 60 900 650]);  hold on;
for k = 1:nCases
    pvi = cumInjected(Rs{k}) / S.PV;
    plot(pvi, 100*Rs{k}.cumOilHist/S.OOIP, 'Color',cols(k,:), 'DisplayName',labs{k});
end
grid on;  xlabel('Pore volumes of water injected (PV)');  ylabel('Recovery factor (% OOIP)');
title('Recovery factor vs injected pore volumes');  legend('Location','northwest');
figList = reg(figList, fig, figDir, 'R06_Recovery_vs_PVI', ...
    'Recovery factor against pore volumes of water injected - the clearest like-for-like comparison of polymer and waterflood.');

%% R07 - injector BHP -------------------------------------------------------
fig = figure('Position',[60 60 900 650]);  hold on;
for k = 1:nCases
    plot(tYr(Rs{k}), max(Rs{k}.injectorBHP, [], 2), 'Color',cols(k,:), 'DisplayName',labs{k});
end
yline(convertTo(p0, barsa), 'k--', 'Initial pressure');
grid on;  xlabel('Time (years)');  ylabel('Maximum injector BHP (bar)');
title('Injection pressure requirement');  legend('Location','best');
figList = reg(figList, fig, figDir, 'R07_Injector_BHP', ...
    'Maximum injector bottom-hole pressure over time (no BHP limit is imposed in the model).');

%% R08 - final water saturation maps ---------------------------------------
cmapW = flipud(parula(256));
fig = mapRows({Rs{iCtl}.finalState.s(:,1), Rs{bestIdx}.finalState.s(:,1)}, ...
    {'Waterflood', sprintf('Polymer %.2f kg/m^3', Cc(bestIdx))}, dims, xm, ym, wells, ...
    [S_min(Rs, 's') 1], 'S_w', 'Final water saturation (end of simulation)', cmapW);
figList = reg(figList, fig, figDir, 'R08_Final_Sw_Maps', ...
    'Final water saturation by layer: waterflood (top row) vs polymer flood (bottom row), same colour scale.');

%% R09 - final polymer concentration --------------------------------------
fig = mapRows({Rs{bestIdx}.finalState.cp(:)}, {sprintf('Polymer %.2f kg/m^3', Cc(bestIdx))}, ...
    dims, xm, ym, wells, [0 max(1e-6, max(Rs{bestIdx}.finalState.cp(:)))], ...
    'c_p (kg/m^3)', 'Final polymer concentration', parula(256));
figList = reg(figList, fig, figDir, 'R09_Final_Polymer_Concentration', ...
    'Final polymer concentration by layer for the best polymer case: shows how far the polymer bank has swept each layer.');

%% R10 - final pressure maps ------------------------------------------------
pC = convertTo(Rs{iCtl}.finalState.pressure,   barsa);
pB = convertTo(Rs{bestIdx}.finalState.pressure, barsa);
fig = mapRows({pC, pB}, {'Waterflood', sprintf('Polymer %.2f kg/m^3', Cc(bestIdx))}, ...
    dims, xm, ym, wells, [min([pC; pB]) max([pC; pB])], 'Pressure (bar)', ...
    'Final reservoir pressure', jet(256));
figList = reg(figList, fig, figDir, 'R10_Final_Pressure_Maps', ...
    'Final pressure by layer for the waterflood and the best polymer case, same colour scale.');

%% R11 - oil saturation gain -------------------------------------------------
dSo = Rs{iCtl}.finalState.s(:,2) - Rs{bestIdx}.finalState.s(:,2);   % positive = oil removed by polymer
lim = max(1e-3, max(abs(dSo)));
cmapDiv = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
fig = mapRows({dSo}, {'S_o(waterflood) - S_o(polymer)'}, dims, xm, ym, wells, ...
    [-lim lim], '\Delta S_o', 'Additional oil displaced by polymer (end of simulation)', cmapDiv);
figList = reg(figList, fig, figDir, 'R11_Oil_Saturation_Gain', ...
    'Difference in final oil saturation (waterflood minus polymer). Red = oil recovered thanks to polymer.');

%% R12 - polymer inventory and balance --------------------------------------
b  = Rs{bestIdx}.bal;  tt = tYr(Rs{bestIdx});
fig = figure('Position',[60 60 950 650]);
yyaxis left;  hold on;
plot(tt, b.injected/1000, '-',  'DisplayName','Injected');
plot(tt, b.stored/1000,   '--', 'DisplayName','Stored in reservoir');
plot(tt, b.produced/1000, ':',  'DisplayName','Produced');
ylabel('Polymer (tonnes)');
yyaxis right;  plot(tt, b.errPct, 'k-.', 'DisplayName','Balance error');
ylabel('Balance error (% of injected)');
grid on;  xlabel('Time (years)');
title(sprintf('Polymer mass balance - %.2f kg/m^3', Cc(bestIdx)));  legend('Location','northwest');
figList = reg(figList, fig, figDir, 'R12_Polymer_Mass_Balance', ...
    'Polymer injected, stored and produced (tonnes) and the closure error of the balance.');

%% R13 - summary vs concentration ---------------------------------------------
RFv  = cellfun(@(r) r.RF, Rs);
wcv  = cellfun(@(r) 100*r.waterCut(end), Rs);
cov  = cellfun(@(r) r.cumOil, Rs);
bhpv = cellfun(@(r) max(r.injectorBHP(:)), Rs);
fig = figure('Position',[60 60 1200 850]);
tl = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Effect of polymer concentration');
nexttile;  bar(Cc, RFv);  grid on;  ylabel('Recovery factor (%)');  xlabel('Polymer (kg/m^3)');
for k = 1:nCases, text(Cc(k), RFv(k), sprintf('%.1f', RFv(k)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom'); end
nexttile;  bar(Cc, cov - cov(iCtl));  grid on;  ylabel('Incremental oil vs waterflood (m^3)');  xlabel('Polymer (kg/m^3)');
nexttile;  bar(Cc, wcv);  grid on;  ylabel('Final water cut (%)');  xlabel('Polymer (kg/m^3)');
nexttile;  bar(Cc, bhpv);  grid on;  ylabel('Max injector BHP (bar)');  xlabel('Polymer (kg/m^3)');
figList = reg(figList, fig, figDir, 'R13_Summary_vs_Concentration', ...
    'Recovery factor, incremental oil, final water cut and maximum injector BHP against polymer concentration.');

%% INDEX FILE ---------------------------------------------------------------
idxFile = fullfile(readmeDir, 'FIGURE_INDEX_v03.txt');
fid = fopen(idxFile, 'w');
fprintf(fid, 'REPORT FIGURES (v03)  -  folder: %s\n\n', figDir);
for k = 1:size(figList,1)
    fprintf(fid, '%s.png\n    %s\n\n', figList{k,1}, figList{k,2});
end
fprintf(fid, 'Economics figures (40-45) and the dashboard (50) are in 07_Figures\\v03.\n');
fclose(fid);

fprintf('\n%d figures saved in:\n%s\nIndex: %s\n', size(figList,1), figDir, idxFile);


%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function drawWells(w)
    plot(w.px, w.py, 'ko', 'MarkerFaceColor','y', 'MarkerSize',7);
    plot(w.ix, w.iy, 'k^', 'MarkerFaceColor','w', 'MarkerSize',7);
end


function fig = mapRows(dataCells, rowLabels, dims, xm, ym, wells, clim, cbLabel, ttl, cmap)
% One row per data vector, one column per layer, shared colour bar.
    nR = numel(dataCells);  nz = dims(3);
    fig = figure('Position', [60 60 1400 330*nR + 90]);
    tl = tiledlayout(nR, nz, 'TileSpacing','compact', 'Padding','compact');
    title(tl, ttl);
    for r = 1:nR
        A = reshape(dataCells{r}, dims);
        for k = 1:nz
            nexttile;  imagesc(xm, ym, A(:,:,k)');  set(gca, 'YDir','normal');  axis equal tight;
            caxis(clim);  colormap(gca, cmap);  hold on;  drawWells(wells);
            if r == 1, title(sprintf('Layer %d', k)); end
            if k == 1, ylabel(rowLabels{r}); end
            if r == nR, xlabel('x (m)'); end
        end
    end
    cb = colorbar;  cb.Layout.Tile = 'east';  cb.Label.String = cbLabel;
end


function m = S_min(Rs, fld)
% Smallest water saturation over all cases (colour-scale lower limit).
    m = 1;
    for k = 1:numel(Rs)
        m = min(m, min(Rs{k}.finalState.(fld)(:,1)));
    end
end


function cumInj = cumInjected(r)
% Cumulative injected water [m3] at the end of every step.
    t   = r.timeDays(:);
    dtd = diff([0; t]);
    qIn = 200*ones(numel(t),1);                         % fallback (m3/day)
    if isfield(r,'wellSols') && ~isempty(r.wellSols)
        for i = 1:numel(t)
            ws = r.wellSols{i};
            qIn(i) = sum(max([ws(2:end).qWs], 0))*86400;
        end
    end
    cumInj = cumsum(qIn .* dtd);
end


function figList = reg(figList, fig, folder, name, caption)
% Save a figure as 300-dpi PNG, register its caption, close it.
    f = fullfile(folder, [name '.png']);
    try
        exportgraphics(fig, f, 'Resolution', 300);
    catch
        saveas(fig, f);
    end
    figList(end+1,:) = {name, caption};
    close(fig);
end