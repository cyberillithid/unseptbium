# Mathematical model

As the first estimation, the block model of the supermatter reactor (its pneumohydraulic scheme mostly) could be presented as a... semi-directed graph?

Typical timescale in all cases given here is the 
<!-- atmospheric module tick of `SSprocessing` (which works on pipe networks) $\Delta \tau = 1$ second. -->
`SSmachinery` tick (for SM core) $\Delta\tau = 2$ seconds.

## Blocks

- [Radiator pipes](./blocks/radiator_pipes.md)
- [Supermatter core](./blocks/supermatter.md)
- [x] Gas  (space: 5 x 3)
- [x] Emitter   `/obj/machinery/emitter`
    - TODO: check the power draw (seems faulty)
- [x] Canister     `machinery/portable_atmospherics/canister` (`empty`, `nitrogen/engine_setup`)
- [x] Pipe `machinery/atmospherics/pipe/simple      `
- [x] Heat exchanger `atmospherics/unary/heat_exchanger`
- [ ] Turbine                     `binary/circulator`
- [ ] TEG      
- [ ] Pump           `atmospherics/binary/pump`
    - [ ] Hi power                     `./pump/high_power`
- [ ] Injector       `atmospherics/unary/outlet_injector`
- [ ] Vent                             `./vent_pump/engine`
- [ ] Filter                      `omni/filter`


## Terminology
XGM/Atmos terms are particularly confounding for those acquainted with the established computational fluid dynamics argot.

Thus, a few notes are to be of use:
- `gas_specific_heat` is $c_{V,m}=c_{p,m}-R$ -- molar specific isochoric heat;
- gases are perfect (ideal, $pV=\nu RT$, and calorically perfect, $c_{v,m}=\mathrm{const}$) unless implied otherwise;
- the natural state of pipenetwork is obviously isochoric ($V=\mathrm{const}$), and independent variables are $T$,$V$ and $\nu_j$ (for species #j);
- XGM usually uses kilomoles, kilopascals, Kelvins and meters;
- something about `group_multiplier`s usage?

## Designations
- $C$ -- heat capacity, J/K;
- $c$ -- specific heat, J/(kg K);
- $c_{,m}$ -- molar specific heat, J/(kg mol);
- $p$ -- pressure, Pa;
- $P$ -- power, W;
- $Q$ -- heat, J;
- $q$ -- heat flux density, W/m²;
- $R$ -- ideal gas constant, ≈8.314 J/(mol K);
<!-- - $S$ -- entropy, -->
- $T$ -- temperature, K;
- $V_m$ -- molar volume, m³/mol;
- $\nu$ -- amount of substance (total), mol;
- $\tau$ -- time, sec;