# Supermatter core
`/obj/machinery/power/supermatter` (`code/modules/supermatter/supermatter.dm`)

Tick is $\Delta\tau=2$ sec.

Units used: 
- {EER} is a placeholder unit; currently displayed as MeV/cm3; looks unphysical-ish, but may require some more handwavium anyway;
- `hp` is hit points (of own supermatter damage counter) ; as well as (unpun?) projectile damage (used by emitter calcs)

Brief stats:
- $T_c=5000$ K -- critical temperature (damage start);
- $\lambda = 700$ -- decay factor, {EER}^{2/3};
- $\eta_r = 1.1$ -- `reaction power modifier` (??);
- $k_Q = 10^4$ -- J / {EER} coeff, `thermal_release_modifier`?
- $k_{O_2} = 15\cdot 10^3$ -- {EER}/mol, `oxygen_release_modifier`
- $k_{Ex} = 1.5\cdot 10^3$ -- {EER}/mol, `product_release_modifier`
- $\mathcal{H}_0 = 1000$ `hp` -- default maximum damage (`explosion point`)
- $d_{\max} = 4.5$ `hp` -- reference maximum damage per tick (at 300 {EER} and default $\mathcal{H}_0$) (`damage_rate_limit`)
- $\eta_g = 0.25$ -- `gas_efficiency`
- $\eta_{N_2} = 0.15$ -- `nitrogen_retardation_factor`

Useless stats:
- $\eta_P = 1.0$ -- `power factor`; just ignoring that;
- `w_class` STRUCTURE 20; `matter_amount_modifier` = 50
- => matter amout: 10000 cm³ EXOTIC, 5000 cm³ STEEL -- may  look like much, but it means only 0.6% of 2.5 m³ (tile volume) is taken by 15 L of crystal. What?  

Hardcoded constants:
- $D_T = 150$ K/`hp` -- damage rate from temperature
- $D_P = 10$ {EER}/`hp` -- damage rate from own power in vacuum
- $\mathcal{H}_d = 1000$ `hp` -- reference SM health for maximum per-tick damage calc
- $\mathcal{P}_d = 300$ {EER} -- reference power for maximum per-tick damage calc
- $T_r = 800$ K -- temperature factor for power reaction increase

Variables:
- $\mathcal{P}$ -- current power, {EER}
- $\mathcal{H}$ -- current integrity (1-damage) (?), `hp`
- $\tau$ -- server time (seconds)
- $t$ -- subsystem time (in iterations;  $\tau = t \cdot \Delta \tau$)
- $T$, $p$, $X_{[?]}$ -- properties of surrounding gases (temperature, pressure, mole fractions).

### Terminology

The originating commits to Baystation (made by the same `atlantiscz`) of the abbreviations EER and EPR used two different phrases for EER. 

In [February 2016 pull request](https://github.com/Baystation12/Baystation12/pull/12104) `Relative EER` is deciphered as `Relative Energy Emission Ratio`.

In [February 2017 pull request](https://github.com/Baystation12/Baystation12/pull/16197) `EER` is re-introduced as `Energy Emission Rate` (dimensionalized as `MeV/cm3`), and `Chamber EPR (Engine Pressure Ratio)` is defined as well.

EER is just variable `power` ($\mathcal{P}$ in these designations) displayed.

EPR is claimed to `show real amount of coolant in the core (in standard canisters worth)`.

EPR is calculated from surrounding gas mixture properties  as $\nu_1 / (V/V_1) / 23.1$, where $V_1=2.5$ m³, $\nu_1$ is amount of substance in $V_1$. Since usually $V=n_tV_1$, where $n_t$ is number of tiles in core room, we can rewrite this in terms of molar volume $V_m = V/\nu$  (m³/mol):
$\nu_1 = V_1/V_m$; EPR = $\dfrac{V_1}{V_m}\dfrac{V_1}{V}\dfrac{1}{23.1} = \dfrac{V_1^2}{V^2}\dfrac{\nu}{23.1\text{ mol}}$



## Processing flow
`/Process()` (which is called once in $\Delta\tau$) describes most of valuable interactions here.

## Damage

SM core tracks damage incurred, which is limited by $\mathcal{H}=1000$. Effectively we may just consider that SM has 1000 HP.

SM core takes damage in two ways. When it is in vacuum, it takes $\max(D_P(\mathcal{P} - 15\eta_P), 0)$ damage per processing tick. This means that it does not heal itself; and if it's energized more than $15\eta_P$ (in MeV/cm3), it starts taking even more damage (and when it's `exploding`, it consumes all gas in the zone, `grav_pulling` it).

Otherwise (when there is some gas to pressure SM core), damage depends on the current gas temperature $T$ in the pressure chamber.

First, some intermediate values are calculated:
- declared maximum increment damage (`damage_inc_limit`) $d_i = f_i d_{\max}$, $f_i = \dfrac{\mathcal{P}}{\mathcal{P}_d} \dfrac{\mathcal{H}_0}{\mathcal{H}_d}$  -- linearly scaled to `ensure that damage doesn't increase too quickly due to super high temperatures resulting from no coolant, for example. We dont want the SM exploding before anyone can react.` (original comment)
- temperature-affected damage $d_t = (T-T_c)/D_T$.

<!-- Then damage increment is calculating by clamping $-d_i$ between $d_T$ and $\delta d_{\max}$, so,
 effectively (since "forgiving" `clamp` reorders Low and Hi if needed), low-bound is $d_L=\min(d_T, \delta d_{\max})$, hi-bound is $d_H=\max(d_T, \delta d_{\max})$, and resulting $\Delta d  =\min(\max(-d_i, d_L), d_H)$.

This kinda looks like that the intent was to clamp $d_T$ to $[-\delta d_{\max}; f_l \delta d_{\max}]$, to linearly scale damage *taken*, not *healed*; we shall denote the assumed intended result of this calculation as $\tilde\Delta d$. 

Let us analyze the possible branches in detail. Note that $\delta d_{\max}>0$; we should expect $f_m>0$, thus $d_i$>0. ($\Delta d > 0$ means damage taken, <0 -- healed). Denoting $T_\delta = \delta d_{\max} D_T$ (=675K on defaults), we get that:

- If $T<T_c$, actual temperature is lower than critical, $d_L=d_T<0$, $d_H =\delta d_{\max}$:
  - $\Delta d = \min(-d_i,d_T)$; $\tilde\Delta d = \min(-\delta d_{\max}, d_T)$ (damage healed)
- If $T>T_c$, temperature is over critical, yet $0<d_L = d_T < \delta d_{\max} = d_H$:
  - $\Delta d = \tilde\Delta d = d_T$ (damage taken)
- If $T>(T_c + T_\delta)$, temperature is sufficiently over critical, thus $d_L = \delta d_{\max} < d_T = d_H$:
  - $\Delta d = \delta d_{\max}$; $\tilde\Delta d = d_i$.

So it seems like this clamping actually hinders healing of underpowered cooled crystals and does not affect -->

Then damage increment should be calculated by clamping $d_T$ between $-d_{\max}$ and $d_i$, thus effectively making *damaging* speed linearly proportional to power and maximal health, and limiting *healing* speed to the fixed value.

### Delamination
Original code notes implied that within that 30-second interval between going into negative integrity and actually exploding something can be done. It doesn't seem so from the first glance at `/explode()` here; this also needs some investigation and documentation, even if it goes a bit beyond the purpose of the current mathmodel.


## Reaction

After determining damage incurred (or healed) on current processing tick, the reaction conditions are considered.

Technically the process of changing ownership of an amount of substance is done before calculating damage, but that does not affect the logics here.

So, given the environment, SM core removes $\nu_r$ moles of it to process. 
$\nu_r = \eta_g \nu_1$, where $\eta_g$ is gas efficiency, and $\nu_1$ is amount of moles in a single `gas_mixture` datum. For most of SM Core uses, the environment should be a single zone, so the group multiplier is equal to number of tiles, and $\nu_1$ represents amount of substance in the gas content of a single tile volume $V_1 = 2.5$ m³. So technically we may just say that $\nu_r$ is amount of substance in volume of $\eta_g V_1$, which may as well describe the envelope of gas around SM crystal or something.
(If it would have been porous, that would have worked as well, but it is usually described as `"Superdense crystalline structure - appears to have been shaped or hewn, lattice is approximately 20 times denser than should be possible."`.)

In any case: given this amount of substance, then we calculate effective oxidizers fraction by summing mole fractions of oxidizers and subtracting nitrogen mole fraction multiplied by retardation coefficient:
$f_O = X_{[O_2]} + X_{[N_2O]} + X_{[NO]} + X_{[NO_2]} - \eta_{N_2} X_{[N_2]}$.

The value of $f_O$ determines target equilibrium power $\mathcal{P}_=$ (400 EER if $f_O>0.8$, 250 EER otherwise), which is then applied to calculate 
$f_T = \left(\dfrac{\mathcal{P}_=}{\lambda}\right)^3 \dfrac{1}{T_r}$ and the power increase $\Delta \mathcal{P}_{rxn} = T f_T f_O = \left(\dfrac{\mathcal{P}_=}{\lambda}\right)^3 \dfrac{T}{T_r} f_O$.

Let us denote intermediate value of core power $\mathcal{P}_1 = \mathcal{P} + \Delta \mathcal{P}_{rxn}$.

Then, with another intermediate values of $E_P = \eta_r \mathcal{P}_1$ and $E_T = (T - 0\degree\mathrm{C})\cdot 1\mathrm{\dfrac{\{EER\}}{K}}$, we determine amount of substance generated by the reaction:
- $\nu_{Ex} = \max(E_P / k_{Ex}, 0)$
- $\nu_{O_2} = \max((E_P+E_T) / k_{O_2}, 0)$

These exhaust gases are added to the consumed mixture at temperature $T$, effectively releasing extra $c_{v,Ex}\nu_{Ex} + c_{v,O_2}\nu_{O_2}$ J of thermal energy. Though it's obvious that under default parameters this heat is negligibly small.

Then the resulting mixture receives additional $Q_r = k_Q E_P$ J of thermal energy; then the resulting temperature is clamped to $[0; 10^4]$ K interval 
and then it is released back to the core chamber.

Thus, since gases are equalized every second, the pressure by itself does not affect the SM core operation in any way whatsoever.

... something about fake time with double-seconds ...



In strict power decay conditions (i.e., $f_O=0$),
the internal power is stably decaying:
$\dfrac{\Delta \mathcal{P}_{loss}}{\Delta\tau} = -(\mathcal{P}/\lambda)^3$;
easily enough we can write 
$\mathcal{P}(\tau) = \mathcal{P}_0 / \sqrt{\dfrac{\mathcal{P}_0^2}{\lambda^3}t+1}$.

According to original code comment, this is `power losses due to radiation`. 

We may note that to keep total $\Delta \mathcal{P} = \Delta \mathcal{P}_{loss} + \Delta \mathcal{P}_{rxn} = 0$, $\dfrac{T}{T_r} f_O=\left(\dfrac{\mathcal{P}}{\mathcal{P}_=}\right)^3$; for example, gas mixture of 65% N2, 35% O2 ($f_O \approx 0.25$) constant $T=4T_r=3200$ K will keep the power stable for $\mathcal{P}=\mathcal{P}_= = 250$ EER.
 
<!-- do I need math here? --->

Heat emitted by SMCore reaction is 
$\Delta Q = k_Q \eta_r \mathcal{P}$ (if we ignore the waste gases).
Thus rough estimate of SMcore thermal power is ~5.5 kW/(MeV/cm3), i.e., you should multiply EER by 5.5 to get kW.
