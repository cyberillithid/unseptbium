# Mathematical model

As the first estimation, the block model of the supermatter reactor (its pneumohydraulic scheme mostly) could be presented as a... semi-directed graph?

Typical timescale in all cases given here is the atmospheric module tick of $\Delta \tau = 2$ seconds.

## Blocks

### Radiators
`/obj/machinery/atmospherics/pipe/simple/heat_exchanging`

Brief stats: $P_{\max{}}=360,\,P_{fatigue}=300$ (atm); $\overline{\lambda} = 0.4$

Radiator pipes, as implemented in the Nebula codebase, have two drastically different modes of operation. In both cases, though, each plate successively modifies state of the whole pipe-network.

### Radiator pipes in space
This actually works **almost** like radiation cooling proper.

(hints: $\alpha$ is W/(m²×K); $\lambda$ is W/(m×K); ...)

#### Physical basis

Stefan-Boltzmann law defines *radiant emittance* (heat flux) of blackbody as $M^\circ = \dfrac{P}{S}=\sigma_0 T^4$, where $\sigma$ is Stefan-Boltzman constant, $S$ is area and $P$ is power ($dQ/d\tau$); then, using *gray body emissivity* $\varepsilon$, we can write elementary heat emitted as $dQ = \varepsilon M^\circ S d\tau$.
Considering the heat backflux from environment with temperature $T_0$, we can write discretized equation for heat difference as $\Delta Q = \varepsilon \sigma T^4 S \Delta\tau$.

Moreover, we usually consider gray body emissivity to be equal to the *absorptance* $A$ -- `ratio of the absorbed to the incident radiant power`. Thus, assuming the body receives a heat flux of $q'$ (W/m²) onto a surface of $S'$, the actual heat gained would be $\Delta Q' = \varepsilon q' S' \Delta\tau$.

*Technically*, all of this usually applies to a single (rigid) body or a gas in thermodynamical equilibrium state. To apply this for a radiator, we could either model the internal convective and conductive heat transfer extensively... or apply for a number of known approximations, e.g., using Nusseldt etc numbers. For example, Ivliev suggests using $q/q_0 \propto (p/p_0)^{0.85}$ (with a similar cross section area and other limitations).

#### How it works in Nebula

Initial gray-body-emissivity of radiator pipe (as given by `Process()` to `radiate_heat_to_space()` as `thermal_conductivity` argument) is $\varepsilon_0=1$ (which makes sense, since they *are* pretty black).

Every pipe segment calculates the `gas_density` $\nu/V = 1/V_m = p/RT = \rho/\mu$ (mol/m³) (actually *molar concentration*, $c_{[gas]}$), which is used to determine effective multiplicative coefficient to `thermal_conductivity` (yet again, implied to be emissivity): 
$\varepsilon=\varepsilon_0 \cdot \min\left(1, \dfrac{c_{[gas]}}{c_{[ref]}}\right)$, where $c_{[ref]} = \dfrac{p_{crit, air}}{R T_{crit, air}}$ -- probably maximum concentration for gaseous air.

## you are located here

calculates thermal-conductivity multiplier $\lambda$ as $\min(1, \dfrac{\nu}{V} / \dfrac{p_{c, air}}{RT_{c, air}})$
= $\min(1, \dfrac{1/V_m}{1/V_{m,c,air}}) = \min(1, \dfrac{V_{m,c,air}}{V_m}) $

$S = 2\,\mathrm{m}^2$; $f = 0.04$ => $S_{eff} = 0.08$?

$\Delta Q = \Delta Q_\odot - \Delta Q_r$; 

solar radiation surplus is $\Delta Q_\odot = E_\odot fS \lambda$, $E_\odot = 200 \mathrm{W}/\mathrm{m}^2$;

radiation loss is $\Delta Q_r = S \sigma \lambda (T - T_0)^4$, $\sigma$ -- Стефана-Больцмана, $T_0 = 3.15$ K

### Radiator pipes on turfs
Whenever the radiator pipe is located on a turf that is `/turf/simulated` or `t/turf/exterior`, the conductivity cooling is used instead (with a code limitation on $\Delta T>20$ K to prevent excessive calculations).

Fourier equation of thermal conductivity usually looks like $dQ = -\lambda \dfrac{\partial T}{\partial \vec{n}} dF d\tau$