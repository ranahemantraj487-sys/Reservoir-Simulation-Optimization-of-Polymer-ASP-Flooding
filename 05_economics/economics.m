%% polymer_economics_03.m
% Cairn-Rajasthan inspired polymer flood - STAGE 5: ECONOMICS (v03 model)

clearvars; close all; clc;

%% 0. SETTINGS ------------------------------------------------------------
projectRoot = 'D:\PROJECT';
sensFile    = fullfile(projectRoot, '06_results', 'polymer_concentration_sensitivity_v03.mat');
resultDir   = fullfile(projectRoot, '06_results');
figureDir   = fullfile(projectRoot, '07_Figures', 'v03');

% ---- economic inputs (PLACEHOLDERS - edit) --------------------------------
eco.oilPrice             = 70;      % USD/bbl
eco.royalty              = 0.20;    % fraction of gross revenue (set to your PSC / royalty terms)
eco.lifting              = 8;       % USD per bbl of oil produced
eco.waterProdCost        = 1.0;     % USD per m3 of water produced (lift + treat + dispose)
eco.waterInjCost         = 0.5;     % USD per m3 of water injected
eco.nWells               = 5;       % 1 producer + 4 injectors
eco.wellCapex            = 1.5e6;   % USD per well (set 0 if the wells already exist; it cancels in the incremental result)
eco.opexPerWell          = 1.0e5;   % USD per well per year (fixed)
eco.polymerPrice         = 3.0;     % USD per kg of polymer (delivered)
eco.polymerHandling      = 0.2;     % USD per kg (mixing, pumping, QC)
eco.polymerFacilityCapex = 1.5e6;   % USD, one-off polymer mixing / injection skid (year 0)
eco.polymerFixedOpex     = 2.0e5;   % USD per year while polymer is being injected
eco.discount             = 0.10;    % discount rate (fraction)
eco.taxRate              = 0.30;    % income tax on (EBITDA - depreciation), no loss carry-forward
eco.bblPerM3             = 6.28981; % unit conversion

fallbackInjRate_m3day = 200;        % only used if wellSols were not saved in the results file

if ~isfile(sensFile)
    error('Sensitivity results not found:\n%s\nRun polymer_concentration_sensitivity_03.m first.', sensFile);
end
if ~isfolder(resultDir), mkdir(resultDir); end
if ~isfolder(figureDir), mkdir(figureDir); end

%% 1. LOAD SENSITIVITY RESULTS --------------------------------------------
S = load(sensFile);
need = {'allResults','polymerConcentrations','simulationDays','OOIP'};
for k = 1:numel(need)
    if ~isfield(S, need{k}), error('Variable "%s" missing in %s', need{k}, sensFile); end
end

C       = S.polymerConcentrations(:);
nCases  = numel(C);
nYears  = floor(S.simulationDays/365 + 1e-9);
OOIP    = S.OOIP;
iCtl    = find(C == 0, 1);
if isempty(iCtl)
    error('The 0 kg/m3 waterflood control is not in the results; the economics are incremental vs the control.');
end

fprintf('\n============================================================\n');
fprintf(' POLYMER FLOOD ECONOMICS (v03)\n');
fprintf('============================================================\n');
fprintf('Cases: %d | Horizon: %d years | Discount rate %.1f %% | Oil %.0f USD/bbl\n', ...
    nCases, nYears, 100*eco.discount, eco.oilPrice);
fprintf('NOTE: prices/costs are placeholders - edit the "economic inputs" block.\n');

%% 2. YEARLY VOLUMES -------------------------------------------------------
Prod = cell(nCases,1);
for k = 1:nCases
    Prod{k} = annualProduction(S.allResults{k}, nYears, fallbackInjRate_m3day, eco.bblPerM3);
end

%% 3. CASH FLOW PER CASE ---------------------------------------------------
Eco = cell(nCases,1);
for k = 1:nCases
    Eco{k} = evaluateCase(Prod{k}, C(k) > 0, eco);
end

incOil_bbl = nan(nCases,1);  incRF = nan(nCases,1);   incNPV = nan(nCases,1);
incCash    = nan(nCases,1);  polyKg = nan(nCases,1);  polyCost = nan(nCases,1);
PUF        = nan(nCases,1);  unitCost = nan(nCases,1);
IRR        = nan(nCases,1);  payback = nan(nCases,1); breakEven = nan(nCases,1);
incCF      = cell(nCases,1);

for k = 1:nCases
    polyKg(k) = sum(Prod{k}.polymer_kg);
    incOil_bbl(k) = sum(Prod{k}.oil_bbl) - sum(Prod{iCtl}.oil_bbl);
    incRF(k)      = 100*(incOil_bbl(k)/eco.bblPerM3)/OOIP;
    incCF{k}      = Eco{k}.cf - Eco{iCtl}.cf;
    incNPV(k)     = Eco{k}.npv - Eco{iCtl}.npv;
    incCash(k)    = sum(incCF{k});
    if k == iCtl, continue; end

    polyCost(k) = Eco{k}.polyCostTotal + Eco{k}.polyCapex;          % undiscounted, USD
    if incOil_bbl(k) > 0
        PUF(k)      = polyKg(k)/incOil_bbl(k);                      % kg per incremental bbl
        unitCost(k) = polyCost(k)/incOil_bbl(k);                    % USD per incremental bbl
    end
    IRR(k)     = calcIRR(incCF{k});
    payback(k) = calcPayback(incCF{k});
    breakEven(k) = calcBreakEven(Prod{k}, Prod{iCtl}, eco);
end

economicsTable = table( ...
    C, incOil_bbl/1e3, incRF, polyKg/1e3, PUF, polyCost/1e6, unitCost, ...
    cellfun(@(e) e.npv, Eco)/1e6, incNPV/1e6, incCash/1e6, IRR*100, payback, breakEven, ...
    'VariableNames', {'Polymer_kg_m3','IncOil_kbbl','IncRF_pct','PolymerInjected_t', ...
    'PUF_kg_per_bbl','PolymerRelatedCost_MUSD','PolymerCost_USD_per_incbbl', ...
    'NPV_MUSD','IncNPV_MUSD','IncNetCash_MUSD','IncIRR_pct','Payback_years','BreakEvenOil_USD_bbl'});

fprintf('\n--- ECONOMIC SUMMARY (incremental vs waterflood control) ---\n');
disp(economicsTable);

cand = find(C > 0);
bestIdx = NaN;
if ~isempty(cand)
    [~, j] = max(incNPV(cand));
    bestIdx = cand(j);
    fprintf('Best case by incremental NPV: %.2f kg/m3  (IncNPV %.2f MUSD, %.0f kbbl extra oil, PUF %.2f kg/bbl)\n', ...
        C(bestIdx), incNPV(bestIdx)/1e6, incOil_bbl(bestIdx)/1e3, PUF(bestIdx));
    if numel(cand) == 1
        fprintf('Only one polymer concentration was simulated (runMode = ''quick'').\n');
        fprintf('Run the sensitivity in ''full'' mode to find the optimum concentration.\n');
    end
end

%% 4. SAVE TABLE -----------------------------------------------------------
writetable(economicsTable, fullfile(resultDir, 'polymer_economics_v03.csv'));
save(fullfile(resultDir, 'polymer_economics_v03.mat'), 'economicsTable', 'eco', 'C', ...
    'bestIdx', 'nYears', 'incNPV', 'incCF', 'Prod', '-v7.3');

%% 5. SENSITIVITY OF THE BEST CASE ----------------------------------------
labelOf = @(c) sprintf('%.2f kg/m^3', c);
yrs = (0:nYears)';

if ~isnan(bestIdx)
    Pk = Prod{bestIdx};  Pc = Prod{iCtl};

    % 5a. oil price x polymer price map
    oilGrid = 30:10:120;   polGrid = 1:1:10;
    M = zeros(numel(oilGrid), numel(polGrid));
    for a = 1:numel(oilGrid)
        for b = 1:numel(polGrid)
            pp = eco;  pp.oilPrice = oilGrid(a);  pp.polymerPrice = polGrid(b);
            M(a,b) = incNPVfun(Pk, Pc, pp)/1e6;
        end
    end

    % 5b. tornado (+/-25 %)
    fn  = {'oilPrice','polymerPrice','discount','lifting','waterProdCost', ...
           'waterInjCost','polymerFacilityCapex','royalty','taxRate'};
    lab = {'Oil price','Polymer price','Discount rate','Lifting cost','Water production cost', ...
           'Water injection cost','Polymer facility capex','Royalty','Tax rate'};
    nT = numel(fn);
    baseNPV = incNPVfun(Pk, Pc, eco);
    vLo = zeros(nT,1);  vHi = zeros(nT,1);
    for i = 1:nT
        pLo = eco;  pLo.(fn{i}) = 0.75*eco.(fn{i});
        pHi = eco;  pHi.(fn{i}) = 1.25*eco.(fn{i});
        vLo(i) = (incNPVfun(Pk, Pc, pLo) - baseNPV)/1e6;
        vHi(i) = (incNPVfun(Pk, Pc, pHi) - baseNPV)/1e6;
    end
    [~, ord] = sort(max(abs([vLo vHi]), [], 2), 'ascend');
end

%% 6. FIGURES --------------------------------------------------------------
% 40 - incremental NPV (and PUF) vs concentration
fig = figure('Color','w');
yyaxis left;
bar(C, incNPV/1e6);  ylabel('Incremental NPV vs waterflood (MUSD)');
yyaxis right;
plot(C, PUF, '-o', 'LineWidth',1.6, 'MarkerSize',7);  ylabel('PUF (kg polymer / incremental bbl)');
grid on;  xlabel('Polymer concentration (kg/m^3)');
title('Incremental NPV and Polymer Utilization Factor');
saveas(fig, fullfile(figureDir, '40_Econ_NPV_PUF_vs_Concentration.png'));

% 41 - cumulative discounted incremental cash flow
fig = figure('Color','w');  hold on;
for k = cand(:)'
    d = (1+eco.discount).^(0:nYears)';
    plot(yrs, cumsum(incCF{k}./d)/1e6, '-o', 'LineWidth',1.6, 'DisplayName', labelOf(C(k)));
end
hz = yline(0, 'k--');  hz.HandleVisibility = 'off';
grid on;  xlabel('Year');  ylabel('Cumulative discounted incremental cash flow (MUSD)');
title('Incremental Cash Flow vs Waterflood');  legend('Location','best');
saveas(fig, fullfile(figureDir, '41_Econ_Cumulative_Discounted_CF.png'));

if ~isnan(bestIdx)
    % 42 - annual net cash flow, control vs best case
    fig = figure('Color','w');
    bar((1:nYears)', [Eco{iCtl}.cf(2:end), Eco{bestIdx}.cf(2:end)]/1e6, 'grouped');
    grid on;  xlabel('Year');  ylabel('Annual net cash flow after tax (MUSD)');
    title('Annual Cash Flow: Waterflood vs Polymer');
    legend('Waterflood (0 kg/m^3)', ['Polymer ' labelOf(C(bestIdx))], 'Location','best');
    saveas(fig, fullfile(figureDir, '42_Econ_Annual_Cash_Flow.png'));

    % 43 - oil price x polymer price map
    fig = figure('Color','w');
    imagesc(polGrid, oilGrid, M);  set(gca, 'YDir','normal');
    cb = colorbar;  ylabel(cb, 'Incremental NPV (MUSD)');
    hold on;
    try
        contour(polGrid, oilGrid, M, [0 0], 'k', 'LineWidth',2);
    catch
    end
    xlabel('Polymer price (USD/kg)');  ylabel('Oil price (USD/bbl)');
    title(sprintf('Incremental NPV map - %s (black line = break-even)', labelOf(C(bestIdx))));
    saveas(fig, fullfile(figureDir, '43_Econ_NPV_Map_Oil_vs_Polymer_Price.png'));

    % 44 - tornado
    fig = figure('Color','w');
    barh(1:nT, [vLo(ord) vHi(ord)]);
    set(gca, 'YTick',1:nT, 'YTickLabel',lab(ord));
    grid on;  xlabel('Change in incremental NPV (MUSD)');
    title(sprintf('Tornado (+/-25 %%) - %s', labelOf(C(bestIdx))));
    legend('-25 %','+25 %', 'Location','best');
    saveas(fig, fullfile(figureDir, '44_Econ_Tornado.png'));
end

% 45 - polymer cost per incremental barrel vs oil price
fig = figure('Color','w');
plot(C(cand), unitCost(cand), '-o', 'LineWidth',1.6, 'MarkerSize',7);  hold on;
yline(eco.oilPrice*(1-eco.royalty) - eco.lifting, 'r--', 'Net oil margin (before water costs)');
grid on;  xlabel('Polymer concentration (kg/m^3)');  ylabel('Polymer-related cost (USD / incremental bbl)');
title('Unit Cost of Incremental Oil');
saveas(fig, fullfile(figureDir, '45_Econ_Unit_Cost.png'));

fprintf('\nECONOMICS COMPLETE\nCSV: %s\nMAT: %s\nFigures: %s\n', ...
    fullfile(resultDir,'polymer_economics_v03.csv'), ...
    fullfile(resultDir,'polymer_economics_v03.mat'), figureDir);


%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function P = annualProduction(r, nYears, fallbackInj, bblPerM3)
% Yearly oil / water production, water injection and polymer injection of one case.
    t   = [0; r.timeDays(:)];                 % days (end of each step, plus t = 0)
    n   = numel(r.timeDays);
    dtd = diff(t);

    cumOil = [0; r.cumOilHist(:)];            % m3
    cumWat = [0; r.cumWaterHist(:)];          % m3

    % injected water (sum of injector surface water rates, m3/day)
    qIn = fallbackInj*ones(n,1);
    if isfield(r,'wellSols') && ~isempty(r.wellSols)
        for i = 1:n
            ws = r.wellSols{i};
            qIn(i) = sum(max([ws(2:end).qWs], 0))*86400;
        end
    end
    cumInj = [0; cumsum(qIn.*dtd)];           % m3

    cumPol = zeros(n+1,1);                    % kg
    if r.polymerConcentration > 0 && isfield(r,'bal') && isfield(r.bal,'injected')
        cumPol = [0; r.bal.injected(:)];
    end

    ye = min((0:nYears)*365, t(end));
    at = @(c) interp1(t, c, ye(:), 'linear');
    P.oil_m3       = diff(at(cumOil));
    P.oil_bbl      = P.oil_m3*bblPerM3;
    P.waterProd_m3 = diff(at(cumWat));
    P.waterInj_m3  = diff(at(cumInj));
    P.polymer_kg   = diff(at(cumPol));
end


function out = evaluateCase(P, hasPolymer, p)
% Yearly cash flow of one case. cf(1) = year 0 (capex), cf(2:end) = years 1..N (after tax).
    N = numel(P.oil_bbl);

    revenue = P.oil_bbl*p.oilPrice*(1 - p.royalty);
    opex    = p.lifting*P.oil_bbl + p.waterProdCost*P.waterProd_m3 ...
            + p.waterInjCost*P.waterInj_m3 + p.opexPerWell*p.nWells*ones(N,1);

    polyCost = zeros(N,1);
    if hasPolymer
        polyCost = P.polymer_kg*(p.polymerPrice + p.polymerHandling) ...
                 + p.polymerFixedOpex*(P.polymer_kg > 0);
    end

    capex0 = p.nWells*p.wellCapex + hasPolymer*p.polymerFacilityCapex;
    dep    = capex0/N*ones(N,1);                         % straight-line over the horizon

    ebitda = revenue - opex - polyCost;
    tax    = p.taxRate*max(0, ebitda - dep);

    out.cf            = [-capex0; ebitda - tax];
    out.npv           = sum(out.cf ./ (1 + p.discount).^(0:N)');
    out.polyCostTotal = sum(polyCost);
    out.polyCapex     = hasPolymer*p.polymerFacilityCapex;
end


function v = incNPVfun(Pk, Pc, p)
% Incremental NPV of a polymer case vs the waterflood control.
    ek = evaluateCase(Pk, true,  p);
    ec = evaluateCase(Pc, false, p);
    v  = ek.npv - ec.npv;
end


function irr = calcIRR(cf)
    irr = NaN;
    if ~(any(cf < 0) && any(cf > 0)), return; end
    f = @(r) sum(cf(:) ./ (1 + r).^(0:numel(cf)-1)');
    try
        irr = fzero(f, [-0.9 100]);
    catch
        irr = NaN;
    end
end


function pb = calcPayback(cf)
% Undiscounted payback in years (linear interpolation inside the year).
    cum = cumsum(cf(:));
    if cum(1) >= 0, pb = 0; return; end
    idx = find(cum >= 0 & (1:numel(cum))' > 1, 1);
    if isempty(idx), pb = NaN; return; end
    pb = (idx-2) + (-cum(idx-1))/(cum(idx) - cum(idx-1));
end


function be = calcBreakEven(Pk, Pc, p)
% Oil price (USD/bbl) at which the incremental NPV is zero.
    f = @(x) incNPVfun(Pk, Pc, setPrice(p, x));
    be = NaN;
    try
        if f(0) >= 0
            be = 0;
        else
            be = fzero(f, [0 1000]);
        end
    catch
        be = NaN;
    end
end


function p = setPrice(p, x)
    p.oilPrice = x;
end