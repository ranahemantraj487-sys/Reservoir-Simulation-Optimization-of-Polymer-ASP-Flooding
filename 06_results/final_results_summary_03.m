%% final_results_summary_03.m
% Cairn-Rajasthan inspired polymer flood - STAGE 6: FINAL RESULTS (v03 model)


clearvars; close all; clc;

%% 0. PATHS ---------------------------------------------------------------
projectRoot = 'D:\PROJECT';
resultDir   = fullfile(projectRoot, '06_results');
figureDir   = fullfile(projectRoot, '07_Figures', 'v03');
readmeDir   = fullfile(projectRoot, '08_Readme');
sensFile    = fullfile(resultDir, 'polymer_concentration_sensitivity_v03.mat');
ecoFile     = fullfile(resultDir, 'polymer_economics_v03.mat');

if ~isfile(sensFile), error('Not found: %s\nRun the sensitivity script first.', sensFile); end
if ~isfile(ecoFile),  error('Not found: %s\nRun polymer_economics_03.m first.', ecoFile); end
if ~isfolder(figureDir), mkdir(figureDir); end
if ~isfolder(readmeDir), mkdir(readmeDir); end

%% 1. LOAD ----------------------------------------------------------------
S = load(sensFile);
E = load(ecoFile);

C      = S.polymerConcentrations(:);
R      = S.allResults;
nCases = numel(C);
iCtl   = find(C == 0, 1);
iBest  = E.bestIdx;
hasBest = ~isnan(iBest);

%% 2. MASTER WORKBOOK -----------------------------------------------------
xlsxFile = fullfile(resultDir, 'FINAL_RESULTS_v03.xlsx');

fn   = fieldnames(E.eco);
vals = cellfun(@(f) E.eco.(f), fn);
ecoInputs = table(fn, vals, 'VariableNames', {'Parameter','Value'});

modelAssumptions = table({ ...
    'No polymer adsorption (ads = 0)'; ...
    'No inaccessible pore volume (dps = 0)'; ...
    'No residual resistance factor (rrf = 1)'; ...
    'Linear viscosity law: muWMult = 1 + 4*c/cpmax'; ...
    'No shear thinning, degradation, salinity or temperature effects'; ...
    'No injector BHP limit (fracture pressure not enforced)'; ...
    'Economic inputs are placeholders - replace with project values'; ...
    'Results are optimistic screening estimates, not forecasts'}, ...
    'VariableNames', {'Assumption'});

try
    if isfile(xlsxFile), delete(xlsxFile); end
    writetable(S.summaryTable,   xlsxFile, 'Sheet', 'Sensitivity');
    writetable(E.economicsTable, xlsxFile, 'Sheet', 'Economics');
    writetable(ecoInputs,        xlsxFile, 'Sheet', 'EconomicInputs');
    writetable(modelAssumptions, xlsxFile, 'Sheet', 'ModelAssumptions');
    fprintf('Master workbook written: %s\n', xlsxFile);
catch ME
    warning('Excel export failed (%s). Writing CSV files instead.', ME.message);
    writetable(S.summaryTable,   fullfile(resultDir, 'FINAL_sensitivity_v03.csv'));
    writetable(E.economicsTable, fullfile(resultDir, 'FINAL_economics_v03.csv'));
    writetable(ecoInputs,        fullfile(resultDir, 'FINAL_economic_inputs_v03.csv'));
end

%% 3. DASHBOARD -----------------------------------------------------------
fig = figure('Color','w', 'Position',[50 50 1500 850]);
tl  = tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Polymer flood - final results dashboard (v03 model)');

% (1) oil rate
nexttile;  hold on;
for k = 1:nCases
    plot(R{k}.timeDays/365, R{k}.oilRate, 'LineWidth',1.5, 'DisplayName', sprintf('%.2f kg/m^3', C(k)));
end
grid on;  xlabel('Time (years)');  ylabel('Oil rate (m^3/day)');  title('Oil rate');  legend('Location','best');

% (2) water cut
nexttile;  hold on;
for k = 1:nCases
    plot(R{k}.timeDays/365, 100*R{k}.waterCut, 'LineWidth',1.5);
end
grid on;  xlabel('Time (years)');  ylabel('Water cut (%)');  title('Water cut');

% (3) cumulative oil
nexttile;  hold on;
for k = 1:nCases
    plot(R{k}.timeDays/365, R{k}.cumOilHist/1000, 'LineWidth',1.5);
end
grid on;  xlabel('Time (years)');  ylabel('Cumulative oil (10^3 m^3)');  title('Cumulative oil');

% (4) recovery factor
nexttile;
RF = cellfun(@(r) r.RF, R);
bar(C, RF);
grid on;  xlabel('Polymer concentration (kg/m^3)');  ylabel('Recovery factor (% OOIP)');  title('Recovery factor');
for k = 1:nCases
    text(C(k), RF(k), sprintf('%.1f', RF(k)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end

% (5) incremental NPV
nexttile;
bar(C, E.incNPV/1e6);
grid on;  xlabel('Polymer concentration (kg/m^3)');  ylabel('Incremental NPV (MUSD)');  title('Incremental NPV vs waterflood');

% (6) polymer inventory of the best case
nexttile;  hold on;
if hasBest
    b = R{iBest}.bal;  tt = R{iBest}.timeDays/365;
    plot(tt, b.injected/1000, '-',  'LineWidth',1.5, 'DisplayName','Injected');
    plot(tt, b.stored/1000,   '--', 'LineWidth',1.5, 'DisplayName','Stored');
    plot(tt, b.produced/1000, ':',  'LineWidth',1.5, 'DisplayName','Produced');
    legend('Location','northwest');
    title(sprintf('Polymer inventory - %.2f kg/m^3', C(iBest)));
else
    title('Polymer inventory');
end
grid on;  xlabel('Time (years)');  ylabel('Polymer (tonnes)');

saveas(fig, fullfile(figureDir, '50_Final_Dashboard.png'));

%% 4. TEXT SUMMARY FOR THE REPORT -----------------------------------------
txtFile = fullfile(readmeDir, 'RESULTS_SUMMARY_v03.txt');
fid = fopen(txtFile, 'w');
if fid < 0, error('Cannot write %s', txtFile); end

fprintf(fid, 'CAIRN-RAJASTHAN INSPIRED POLYMER FLOOD - RESULTS SUMMARY (v03 model)\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM'));

fprintf(fid, 'MODEL\n');
fprintf(fid, '  Pore volume            : %.4e m3\n', S.PV);
fprintf(fid, '  OOIP                   : %.0f m3\n', S.OOIP);
fprintf(fid, '  Simulation length      : %.0f days (%.1f years)\n', S.simulationDays, S.simulationDays/365);
fprintf(fid, '  Polymer cases (kg/m3)  : %s\n\n', mat2str(C', 3));

fprintf(fid, 'RESERVOIR PERFORMANCE\n');
for k = 1:nCases
    fprintf(fid, '  %.2f kg/m3 : cum oil %9.0f m3 | RF %5.2f %% | final water cut %5.1f %% | max inj. BHP %6.1f bar\n', ...
        C(k), R{k}.cumOil, R{k}.RF, 100*R{k}.waterCut(end), max(R{k}.injectorBHP(:)));
end

if hasBest
    T = E.economicsTable;
    fprintf(fid, '\nBEST POLYMER CASE (highest incremental NPV): %.2f kg/m3\n', C(iBest));
    fprintf(fid, '  Waterflood RF          : %.2f %%\n', R{iCtl}.RF);
    fprintf(fid, '  Polymer RF             : %.2f %%\n', R{iBest}.RF);
    fprintf(fid, '  Incremental oil        : %.0f m3 (%.0f kbbl, +%.2f RF points)\n', ...
        R{iBest}.cumOil - R{iCtl}.cumOil, T.IncOil_kbbl(iBest), T.IncRF_pct(iBest));
    fprintf(fid, '  Polymer injected       : %.0f t | PUF %.2f kg per incremental bbl\n', ...
        T.PolymerInjected_t(iBest), T.PUF_kg_per_bbl(iBest));
    fprintf(fid, '  Incremental NPV        : %.2f MUSD (discount %.0f %%)\n', ...
        T.IncNPV_MUSD(iBest), 100*E.eco.discount);
    fprintf(fid, '  Incremental IRR        : %.1f %%\n', T.IncIRR_pct(iBest));
    fprintf(fid, '  Payback                : %.2f years\n', T.Payback_years(iBest));
    fprintf(fid, '  Break-even oil price   : %.1f USD/bbl (base case %.0f USD/bbl)\n', ...
        T.BreakEvenOil_USD_bbl(iBest), E.eco.oilPrice);
end

fprintf(fid, '\nASSUMPTIONS / LIMITATIONS\n');
for k = 1:height(modelAssumptions)
    fprintf(fid, '  - %s\n', modelAssumptions.Assumption{k});
end

fprintf(fid, '\nOUTPUT FILES\n');
fprintf(fid, '  06_results\\polymer_concentration_sensitivity_v03.mat / .csv\n');
fprintf(fid, '  06_results\\polymer_economics_v03.mat / .csv\n');
fprintf(fid, '  06_results\\FINAL_RESULTS_v03.xlsx\n');
fprintf(fid, '  07_Figures\\v03\\ (sensitivity 17-36, economics 40-45, dashboard 50)\n');
fclose(fid);

fprintf('Summary written: %s\n', txtFile);
fprintf('Dashboard saved: %s\n', fullfile(figureDir, '50_Final_Dashboard.png'));
fprintf('\nFINAL RESULTS COMPLETE\n');