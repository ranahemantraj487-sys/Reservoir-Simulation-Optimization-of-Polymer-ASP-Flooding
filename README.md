# Synthetic Reservoir Polymer-Flood Simulation and Economic Screening (MRST)

A MATLAB Reservoir Simulation Toolbox (MRST) study comparing conventional waterflooding with polymer flooding in a synthetic, heterogeneous, three-layer reservoir. The project includes reservoir modelling, production simulation, polymer transport, injection-pressure analysis, sensitivity analysis, and preliminary economic screening.

> **Scope:** The current model implements polymer flooding, not a complete alkaline–surfactant–polymer (ASP) process. It does not explicitly model alkali or surfactant transport or chemistry. The reservoir is synthetic and inspired by Rajasthan field conditions; it is not a history-matched field model.

> **Sensitivity status:** The currently saved v03 results compare waterflooding at 0 kg/m³ polymer concentration with polymer flooding at 1.00 kg/m³. The five-case concentration sweep remains to be completed before an optimal concentration can be established.

## Project Highlights

- Developed a 3D, heterogeneous reservoir model using MATLAB and MRST.
- Simulated a 10-year waterflood and polymer-flood scenario.
- Compared cumulative oil production, recovery factor, water cut, and injector bottom-hole pressure (BHP).
- Implemented polymer transport and concentration-dependent water viscosity.
- Prepared a preliminary economic comparison using incremental cash flow and screening metrics.
- Exported simulation results to CSV, MAT, and Excel formats and generated comparison figures.

## Key Results

The saved v03 simulation reports the following results over a 10-year period.

| Metric | Waterflood Control | Polymer Flood |
|---|---:|---:|
| Polymer injection concentration | 0 kg/m³ | 1.00 kg/m³ |
| Cumulative oil production | 238,394.96 m³ | 323,801.18 m³ |
| Recovery factor | 33.14% OOIP | 45.01% OOIP |
| Incremental oil versus waterflood | — | 85,406.22 m³ |
| Incremental recovery | — | 11.87 percentage points |
| Final oil production rate | 21.03 m³/day | 29.18 m³/day |
| Final water cut | 89.50% | 85.42% |
| Maximum injector BHP | 269.72 bar | 302.42 bar |

The polymer-flood case produces approximately **85,406 m³ of additional oil** relative to the waterflood control in the saved simulation. It also produces a lower final water cut, while requiring a higher maximum injector BHP.

These are simulated results from a synthetic model. The model does not currently impose an injector fracture-pressure constraint, so operational pressure feasibility has not been established.

## Reservoir Model

| Parameter | Value |
|---|---|
| Grid | 50 × 50 × 3 cells |
| Total grid cells | 7,500 |
| Reservoir dimensions | 400 m × 400 m × 30 m |
| Pattern area | 16 hectares |
| Porosity range | 0.15–0.25 |
| Nominal horizontal permeability by layer | 400, 100, and 30 mD |
| Vertical-to-horizontal permeability ratio | 0.1 |
| Initial reservoir pressure | 250 bar |
| Initial water saturation | 0.25 |
| Initial oil saturation | 0.75 |
| Water viscosity | 1 cP |
| Oil viscosity | 10 cP |
| Approximate initial oil in place (OOIP) | 719,320 m³ |
| Simulation period | 10 years (3,650 days) |

### Well Configuration

The model uses a five-spot injection pattern:

- **Producer:** One central producer controlled at a BHP of 150 bar.
- **Injectors:** Four corner injectors.
- **Injection rate:** 50 m³/day per injector.
- **Total injection rate:** 200 m³/day.

The same reservoir configuration and production-control approach are used for the waterflood and polymer-flood comparison.

## Polymer Flood Model

The polymer-flood simulation uses MRST's `OilWaterPolymerModel`.

The polymer viscosity multiplier is represented by:

`muWMult = 1 + 4*c/cpmax`

where:

- `c` is the local polymer concentration.
- `cpmax = 1 kg/m³` is the reference polymer concentration.

The model uses a Todd–Longstaff mixing parameter of `0.85`. The reported polymer-flood case injects polymer at a concentration of `1.00 kg/m³`.

### Current Polymer Assumptions

- Zero polymer adsorption.
- Zero inaccessible pore volume.
- Residual resistance factor equal to 1.
- Simplified concentration-dependent viscosity.
- No explicit polymer degradation, shear-thinning, or salinity/temperature effects.
- No explicit alkali or surfactant mechanisms.

These assumptions simplify the study and should be considered when interpreting the results.

## Preliminary Economic Screening

The saved economic analysis reports the following for the 1.00 kg/m³ polymer case relative to the waterflood control.

| Economic Metric | Reported Result |
|---|---:|
| Incremental NPV | $8.87 million |
| Incremental IRR | 93.4% |
| Undiscounted payback period | 1.60 years |
| Break-even oil price | $25.9/bbl |
| Polymer injected | 730 tonnes |
| Polymer usage factor | Approximately 1.36 kg per incremental barrel |

The current economic script uses illustrative assumptions, including an oil price of $70/bbl, polymer cost of $3/kg, a 10% discount rate, and assumed royalty and operating costs.

**These economics are preliminary screening estimates, not a validated field-development forecast.** The results should be reassessed using defensible cost assumptions and a broader sensitivity analysis.

## Polymer Mass-Balance Check

The saved diagnostic for the 1.00 kg/m³ case reports approximately:

- Polymer injected: 730,000 kg
- Polymer stored: 438,294 kg
- Polymer produced: 301,616 kg
- Mass-balance discrepancy: approximately 1.36%

The calculated stored and produced quantities exceed the injected quantity. This discrepancy remains unresolved and should not be interpreted as exact polymer mass conservation.

A further validation step is required to reconcile the inventory calculation with the transport formulation and numerical simulation.

## Repository Structure

```text
Reservoir-Simulation-Optimization-of-Polymer-ASP-Flooding/
│
├── 01_Baseline/       # Reservoir and fluid setup; primary depletion
├── 02_Waterflood/     # Waterflood simulations and checks
├── 03_Polymer/       # Polymer-flood simulations and validation
├── 04_sensitivity/   # Polymer concentration sensitivity script
├── 05_economics/     # Preliminary economic analysis
├── 06_results/       # CSV, MAT, and Excel results
├── 07_Figures/       # Results plots and report figures
└── 08_Readme/        # Supporting documentation and figure index
```

The current concentration-comparison script is located at:

`04_sensitivity/fluid_model_check_sensitivity.m`

## Software and Reproducibility

### Requirements

- MATLAB
- MATLAB Reservoir Simulation Toolbox (MRST)
- MRST modules: `ad-core`, `ad-props`, `ad-blackoil`, and `ad-eor`

### Reproducing the Saved Results

1. Install a compatible MATLAB release and MRST installation with the required modules.
2. Clone or download this repository.
3. Update hard-coded project paths, including `projectRoot`, to match your local directory.
4. Open `04_sensitivity/fluid_model_check_sensitivity.m`.
5. The current `runMode = 'quick'` setting compares 0 and 1.00 kg/m³ polymer concentration.
6. To run the five planned cases, set `runMode = 'full'` and execute the sensitivity script.
7. Execute `05_economics/polymer_economics_03.m` after the required sensitivity outputs have been generated.
8. Run `06_results/final_results_summary_03.m` and `07_Figures/generate_report_figures_03.m` to update the summary and figures.
9. Review the generated CSV files, MAT files, workbook, and figures to ensure they all correspond to the same simulation run.

**Important:** The currently saved results represent the two-case quick comparison, not a completed five-concentration optimization.

## Limitations and Future Work

1. Complete the five-case polymer concentration sweep at 0, 0.25, 0.50, 0.75, and 1.00 kg/m³.
2. Compare recovery, polymer consumption, injection pressure, water cut, and economics across the full concentration range.
3. Resolve the polymer mass-balance discrepancy and document the calculation methodology.
4. Introduce a defensible maximum injector BHP or fracture-pressure constraint.
5. Test economic sensitivity to oil price, polymer price, operating costs, and discount rate.
6. Improve the polymer model by considering adsorption, inaccessible pore volume, retention, degradation, and other relevant mechanisms where suitable.
7. If full ASP flooding is the intended scope, implement and validate the alkali and surfactant mechanisms. Otherwise, retain the narrower polymer-flooding project description.

## Conclusion

The saved v03 results demonstrate the potential for improved simulated oil recovery from polymer flooding compared with conventional waterflooding in a synthetic heterogeneous reservoir. The polymer case also requires higher injector pressure, highlighting the importance of balancing recovery improvement against operational constraints and chemical cost.

The work provides a foundation for further reservoir simulation, concentration optimization, pressure-constrained evaluation, and economic screening. Its results should be interpreted within the stated model assumptions and validation limitations.

---

**Technologies:** MATLAB · MRST · Reservoir Simulation · Waterflooding · Polymer EOR · Production Analysis · Sensitivity Analysis · Economic Screening
