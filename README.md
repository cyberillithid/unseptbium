# Unseptbium (Supermatter Core Mathmodel)

The goal of this project is twofold:
- to automatically parse the structure of the supermatter reactor from raw .dmm files (and optionally to render it in a traditional engineering way as well);
- and to evaluate a mathematical model of the aforementioned reactor.

The latter as well allows for solving power balance (or, here -- heat balance) of the system, approximating the stationary (or almost-stationary)
mode of operation.

The project (currently) consists of:
- Rust package, which implements (using [Spacemaniac's DMM tools](https://github.com/SpaceManiac/SpacemanDMM), ):
  - binary CLI executable for extraction of relevant map segment (status: WIP, dirty hack)
  - dynamic library with Python bindings for map parsing (inspired by [StrongDMM](https://github.com/SpaiR/StrongDMM) SpacemanDMM bindings) (status: not implemented),
- Python realisation of blocks of supermatter core mathematical model (status: to be reimplemented),
- a plethora of Jupyter Notebook scripts with tests.

Planned extra features:
- rendering of block model as draw.io-diagram (or maybe plain SVG? or CAD drawing?)
- power/heat balance solution.

## Building blocks of math model

- [x] Gas  (space: 5 x 3)
- [x] Core 
- [x] Emitter   `/obj/machinery/emitter`
    - TODO: check the power draw (seems faulty)
- [x] Canister     `machinery/portable_atmospherics/canister` (`empty`, `nitrogen/engine_setup`)
- [x] Pipe `machinery/atmospherics/pipe/simple      `
- [x] Radiator pipes            `./pipe/simple/heat_exchanging`
- [x] Heat exchanger `atmospherics/unary/heat_exchanger`
- [ ] Pump           `atmospherics/binary/pump`
    - [ ] Hi power                     `./pump/high_power`
- [ ] Injector       `atmospherics/unary/outlet_injector`
- [ ] Vent                             `./vent_pump/engine`
- [ ] Filter                      `omni/filter`
- [ ] Turbine                     `binary/circulator`
- [ ] TEG      `/obj/machinery/generator`
