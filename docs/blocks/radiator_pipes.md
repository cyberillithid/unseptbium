# Radiators
`/obj/machinery/atmospherics/pipe/simple/heat_exchanging` (`code/module/atmospherics/he_pipes.dm` and `./datum_pipeline.dm`)

Is technically a kind of [pipe](./pipe.md), thus preserving default volume per tile of $V_1$ = 70 L and implied cross-section of 0.07 m² (diameter ≈0.3m).

Brief stats:
- $P_{\max{}}=360,\,P_{fatigue}=300$ (atm) -- for pressure damage?
- $\overline{h} = 0.4$ -- heat transfer coefficient?
- $S = 2$ m² -- surface area
- $\overline{F} = 0.04$ -- "exposed surface area ratio": $\dfrac{3\text{ cm}+100\text{ cm}\cdot\sin 3\degree}{2\cdot(3\text{ cm}+100\text{ cm})}$ per comment

Radiator pipes, as implemented in the Nebula codebase, have two drastically different modes of operation. In both cases, though, each plate successively modifies state of the whole pipe-network.

## Physical basis

(hints: $h$ is; $\lambda$ is ; ...)
Stefan-Boltzmann law defines *radiant emittance* (heat flux) of blackbody as $M^\circ = \dfrac{P}{S}=\sigma_0 T^4$, where $\sigma$ is Stefan-Boltzman constant, $S$ is area and $P$ is power ($dQ/d\tau$); then, using *gray body emissivity* $\varepsilon$, we can write elementary heat emitted as $dQ = \varepsilon M^\circ S d\tau$.
Considering the heat backflux from environment with temperature $T_0$, we can write discretized equation for heat difference as $\Delta Q = \varepsilon \sigma T^4 S \Delta\tau$.

Moreover, we usually consider gray body emissivity to be equal to the *absorptance* $A$ -- `ratio of the absorbed to the incident radiant power`. Thus, assuming the body receives a heat flux of $q'$ (W/m²) onto a surface of $S'$, the actual heat gained would be $\Delta Q' = \varepsilon q' S' \Delta\tau$.

*Technically*, all of this usually applies to a single (rigid) body or a gas in thermodynamical equilibrium state. To apply this for a radiator, we could either model the internal convective and conductive heat transfer extensively... or apply for a number of known approximations, e.g., using various values or empirical correlations for Nusselt number values to determine heat transfer coefficient (W/(m²×K), $h$ per Newton's law of cooling of $q=h\Delta T$). For example, Ivliev suggests using $q/q_0 \propto (p/p_0)^{0.85}$ (with a similar cross section area and other limitations).

Fourier equation of thermal conductivity usually looks like $dQ = -\lambda \dfrac{\partial T}{\partial \vec{n}} dS d\tau$, $\lambda$ is thermal conductivity coefficient (W/(m×K)), but it is not used here.
Though if we assume the flow inside pipes to be 'fully developed laminar flow', then we can use value of $\mathrm{Nu} = \dfrac{hD}{\lambda}$ of 3.66 or 4.36; thus giving us heat transfer coefficient between gas and pipe surface of 0.03 .. 0.36 (per air at 25°C). This is actually close enough to given value of $\overline{h}=0.4$.

## Radiator pipes in space
This actually works **almost** like radiation calculations proper.

Initial gray-body-emissivity of radiator pipe (as given by `Process()` to `radiate_heat_to_space()` as `thermal_conductivity` argument) is $\varepsilon=1$ (which makes sense, since they *are* pretty black).

Every pipe segment calculates the `gas_density` $\nu/V = 1/V_m = p/RT = \rho/\mu$ (mol/m³) (not actually *molar concentration* -- just inverse of molar volume, $V_m^{-1}$), which is used to determine effective multiplicative coefficient to `thermal_conductivity` (yet again, implied to be emissivity): 
$\eta= min\left(1, \dfrac{V_m^{-1}}{V_{m,ref}^{-1}}\right)$, where $V_{m,ref}^{-1} = \dfrac{p_{crit, air}}{R T_{crit, air}}$ -- probably minimal volume for gaseous air (before going supercritical).

Then we calculate total heat gained as $\Delta Q = \Delta Q_\odot - \Delta Q_r$, where 
- incident solar radiation heat $\Delta Q_\odot = \varepsilon \eta q_\odot (S \overline{F})$, $q_\odot = 200$ W/m² -- solar flux density (equal to one of Sol at distance somewhere in the middle asteroid belt, 2.6 au, ca. Ceres);
- radiation emitted is $\Delta Q_r = \eta \varepsilon S \sigma (T - T_{CMB})^4$, with $T_{CMB} = 3.15$ K -- `Cosmic Microwave Background`.

We may note that the error of $\Delta Q_r$ calculation ($(T-T_{CMB})^4$ instead of $T^4-T_0^4$) is ≈12% at 100K and exponentially recedes to sub-1% by 1250K.

Obviously, the given equation stabilizes with $\Delta Q=0$ at $T=\sqrt[4]{q_\odot\overline{F}/\sigma} + T_{CMB} \approx 111$ K (or at ~109K in proper formulation).

## Radiator pipes on turfs
Whenever the radiator pipe is located on a turf that is `/turf/simulated` or `/turf/exterior`, the heat transfer is used **instead** of radiation (with a code limitation on $\Delta T>20$ K to prevent excessive calculations). And in the previous iteration of ministation some of the radiator pipes were indeed located upon `/turf/simulated/floor/plating`.

Given turf temperature of $T_t$ and pipenet temperature of $T$, we first determine the $\Delta T=T_t-T$ and heat capacity of the current pipe segment $C_1 = c_{v,m} V_1/V_m$; $C_t$ of `/turf/simulated/floor` is 10 kJ/K (which would mean, if we deem `hull plating` to be steel, that the floor thickness is 25mm).

Then, easily enough, $\Delta Q = \overline{h} \Delta T \dfrac{C_tC_1}{C_t+C_1}$, which... looks like a scaled solution to a infinitely long heat transfer: with $\overline{h}=1$, $\Delta Q$ would be the heat required to equilibrate both gas in pipe and turf to an equal $T'$ such that $C_t (T_t-T')=C_1(T'-T)$. 

Thus effective Newton's heat transfer coefficient $h=\overline{h} / (C_t^{-1}+C_1^{-1})$.

## Model evaluation
Technically, atmos module processes all the pipes in succession. Thus, for a single pipe network that is cooled by $n$ spaced and $m$ turfed radiator pipes, $n$ heat radiations and $m$ heat transfers are performed isochorically iteratively.

Since neither amount of substance $\nu$ nor volume $V$ change with isochoric heat radiation, molar volume $V_m$ stays constant and efficiency $\eta$ does not change; nor does change $h$, since gases are calorically perfect. This means that the only 'moving parts' here are the temperature (singular $T$ for the pipe network and $T_i,\,i\in[1;m]$ for turf temperatures).

Using $C_v=c_{v,m}\nu=C_1V$ as total pipe network gas heat capacity, we can write
$\Delta Q_{1space} = C_v \Delta T_{1space} = \eta\varepsilon S (\overline{F}q_\odot - \sigma(T-T_{CMB})^4)$ and $\Delta Q_{turf\#i} = C_v \Delta T_{turf\#i}=h(T-T_i)$.

Numerical experiments show that coefficient $\eta$ actually limits cooling.
For a given $V$ volume of pipe network, we can define $\nu_{\eta} = V/V_{m,ref}$ (where $V_{ref}$ is molar air volume at air critical point). While amount of gas inside pipe network is $\nu<\nu_{\eta}$, it does not affect the cooling process in any way whatsoever; with $\nu>\nu_{\eta}$, the cooling exponentially slows to none.

Let's rephrase the equation considering this (and applying $\varepsilon=1$ to unclutter it, as well as using $q' = q_\odot \overline{F} - \sigma (T-T_{CMB})^4$ for effective heat flux density, W/m²):

$\Delta Q_{1space} = c_{v,m}\nu \Delta T_{1space} = \min(\dfrac{\nu}{\nu_\eta}, 1) S_1 q' \Delta \tau$,

$\Delta T_{1space} = \dfrac{S_1q'\Delta\tau}{\nu_\eta c_{v,m}}$ for $\nu < \nu_\eta$, $ = \dfrac{S_1q'\Delta\tau}{\nu c_{v,m}}$ for $\nu > \nu_\eta$.

This can be perceived as limiting $C_{v,eff}=\max(\nu,\nu_\eta)c_{v,m}$
to get $\Delta T_{1s} = Sq'\Delta\tau/C_{v,eff}$. The physical implications of this are unclear, even if setting $\eta=1$ lets temperature drop instantly close to 0K for small $\nu$.

Now let us consider the effective heat flux $q'$. A way to simplify the iterative process would be to reinterpret it as continuous: $\dfrac{dy}{dx} = a - by^4$, where $y=T-T_{CMB}$, $x=n_{space}$, $a=\Delta T_\odot= \overline{F}Sq_\odot/C_{ve}$, $b=\sigma S/C_{ve}$ (alternatively, with $y=T$, $a'=a+\sigma T_{CMB}^4/C_{ve}$ we could get a bit more physically correct equation).

Assuming that since $\Delta T_\odot << \Delta T_{rad}$ (for $T>300$ K), we can ignore inhomogeneous part to solve $y'+by^4=0$ to $y=(3bx+k)^{-1/3}$, with $y(0)=T$ giving us $k=(T-T_{CMB})^{-3}$. Then we can approximate $n$ runs of radiative cooling by calculating $T_n = y(x=n)+T_{CMB} +n \Delta T_\odot$.

For temperature diapason that may interest us in this task, relative error is $<0.01 T_0$ for $T_0\in[300;3200]$ K or for $T_0<600, N>30$. Technically this could have indicated that the step of 1 meter (1 pipe tile) is too rough for numerical integration?

And, overall, heat loss would be $\Delta Q_\forall = \nu c_{v,m} (T - T_n)$

TODO: later: check $\Delta\tau$ / $x_n$ conditions by $\Delta T_{1s}$ limits.

### Summary of numerical experiments
(TODO: reexplicate this properly)

Linear application of solar flux is good enough for T>200K

turn downs? with temp too high cooling is unnaturally fast? analyse it!


alternatively: calc this as a single radiating body, compare with physicalish results

## Final equation form

For $n$ tiles of space-radiating pipes belonging to a network of total volume $V$ with $\nu$ moles of gas mixture with average specific heat of $c_{v,m}$, this heat exchange drops temperature $T$ by
$$\dfrac{T - T_{sp}}{\Delta \tau} = T_d - (3n\sigma S_1/C_{ve}+T_d^{-3})^{-1/3}  - n\overline{F}S_1q_\odot/C_{ve},$$

$$\dfrac{\Delta Q_{sp}}{\Delta \tau} = \dfrac{T - T_{sp}}{\Delta \tau} c_{v,m} \nu,$$
where $T_d={T-T_{CMB}}$, $C_{ve} = c_{v,m} \max(\nu, \dfrac{V}{V_{m,ref}})$,

or, alternatively, $(T-T_{sp})/\Delta\tau = T_d(1 - \dfrac{1}{\sqrt[3]{k_1T_d^3+1}})-k_2$ -- and this looks unintegratable. On the other hand, we may easily substitute time with space with the approach described above; then for $\tau=\tau_1$ temperature would be $T(n\tau_1)$? Still, this would work only for standalone cooling.