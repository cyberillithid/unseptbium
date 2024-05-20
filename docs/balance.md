# Balance eqns

Let's designate the volumes in a following way:
- in hot loop:
    - $V_{h0}$ -- SM core room;
    - $V_{h1}$ -- hot turbine (T1) input;
    - $V_{h2}$ -- hot turbine (T1) output;
- in cold loop:
    - $V_{c1}$ -- cool turbine (T2) input;
    - $V_{c2}$ -- cool turbine (T2) output.

All of these volumes are constants (parameters loaded from the map), and each of those has their own $p$, $T$ and $\nu$.

There is a mass limitation of $\sum \nu_{c,i}=\nu_c$, $\sum \nu_{h,i}=\nu_h$ -- total amount of substance loaded into a relevant contour.

Basically the question of whether there is a stationary state means that, given $\mathcal{E}_{SM}$ -- EER of the core (ignoring its decay or runaway states) and $\nu_c$, $\nu_h$ -- amount of substance pumped into contours -- we must find all pressures and temperatures for dynamic equilibrium. And if $T_{h0}<5000$ K -- SM core does not take damage -- all is good. 

## Block model
Without graphing diags or drawing loop lines, linearized model of a working contour looks somewhat like:
```
h0{+SM} -[PI1]-> h1 -[T1]-> h2 -[F1,2]--[PO1] -> h0;
                     ||  TEG
c2{-HE} -[HP1]-> c1 -[T2]-> c2;
```
Wherein supermatter core `SM` inserts heat (and waste gases, filtered by `F1,F2`), `TEG` exchanges it between contours, and heat exchanging pipes `HE` remove it.

- Heats:
    - [+SM](./blocks/supermatter.md): SM Core emits thermic power of $N_\Theta = k_Q \eta_r \mathcal{E}_{SM}$.
    - [TEG](./blocks/teg.md) removes $\Delta E_{TEG}$ heat from hot loop and puts $(1-\eta_T)\Delta E_{TEG}$ into cold loop.
    - [-HE](./blocks/radiator_pipes.md) pipes remove $\Delta Q_{sp}$ heat from cold loop.
- Amounts of substance:
    - [PI1, PO1, HP1](./blocks/pump.md) -- `vent-pump`, `outlet-injector` and `pump/high-power` respectively -- all work like pumps (with different power ratings and internal volumes); they move $\Delta\nu_{?P?1}$ moles of matter through themselves.
    - [T1,T2](./blocks/teg.md) are turbines; and $\nu_{T1}$,$\nu_{T2}$ moles are not just moved, but also heat-exchanged in TEG.

## Heat balance

To reach the stationary state,the total heat of the system must be unchanging. Thus, we may assume that all heat emitted by SM core is processed by TEG, and all TEG waste heat is successfully radiated into space. That is,

$$\begin{equation}
\begin{cases}
\Delta E_{TEG} = N_\Theta, \\
\Delta Q_{sp} = (1-\eta_T) \Delta E_{TEG} = (1-\eta_T) N_\Theta.
\end{cases}
\end{equation}$$

### Roughest model

Ignoring the pumping mechanics, we may just assume that both loops are saturated enough that turbines are consuming all their intake (of $V_t$=400L). Assuming both loops are pure nitrogen, the TEG heat exchange rate is determined by $\nu_h$ and $\nu_c$, and the system (1) with two unknowns -- $T_h$ and $T_c$ -- is cleanly separable: first we solve $\Delta Q_{sp} = (1-\eta_T)N_\Theta$ for cold loop temperature, and then find out hot loop temperature from the TEG heat exchange rate.

Unfortunately, this simplified system gives us for new map and 3 canisters of N2 in each loop results of 565+ kW TEG power, 600K cold loop and 1000K hot loop -- which is obviously too high of a discrepancy from the actual runs.

## Mass balance

For the stationary state the mass distribution in the loop must also be static. Thus we can write that

$$\begin{equation}\begin{cases}
\dot\nu_{c1} = \Delta\nu_{HP1} - \Delta\nu_{T2} = 0\\
\dot\nu_{c2} = \Delta\nu_{T2} - \Delta\nu_{HP1} = 0\\
\dot\nu_{h1} = \Delta\nu_{PI1} - \Delta\nu_{T1} = 0\\
\dot\nu_{h2} = \Delta\nu_{T1} - \Delta\nu_{PO1} = 0\\
\dot\nu_{h0} = \Delta\nu_{PO1} - \Delta\nu_{PI1} = 0\\
\end{cases}\end{equation}$$

Removing equations for $\dot\nu_{c2}$ and $\dot\nu_{h2}$ as extraneous (linear combination of others), we get three effective mass balance equations.

Expanding heat balance equations of different volumes -- using  $\Delta\nu_{T2} = \Delta\nu_{HP1} = \dot\nu_c$ and
$\Delta\nu_{PI1} = \Delta\nu_{PO1} = \Delta\nu_{T1} = \dot\nu_h$ -- assuming that all loops are pure nitrogen, yet again -- we get

$$\begin{equation}\begin{cases}
\nu_{c1} c_v T_{c1} = (\nu_{c1} -\dot\nu_c) c_v T_{c1} + \dot\nu_c c_v T_{c2}\\
\nu_{c2} c_v T_{c2} = (\nu_{c2} - \dot\nu_c) c_v T_{c2} + \dot\nu_c c_v T_{c1} + (1-\eta_T)\Delta E_{TEG} - \Delta Q_{sp}\\
\nu_{h1} c_v T_{h1} = (\nu_{h1} - \dot\nu_h) c_v T_{h1} + \dot\nu_h c_v T_{h0}\\
\nu_{h2} c_v T_{h2} = (\nu_{h2} - \dot\nu_h) c_v T_{h2} + (\dot\nu_h c_v T_{h1} - \Delta E_{TEG})\\
\nu_{h0} c_v T_{h0} = (\nu_{h0} - \dot\nu_h) c_v T_{h0} + \dot\nu_h c_v T_{h2} + N_\Theta\\
\end{cases}\end{equation}$$

which may be resolved to

$$\begin{equation}\begin{cases}
T_{c1}=T_{c2}=T_c\\
(1-\eta_T)\Delta E_{TEG} = \Delta Q_{sp}\\
T_{h1} = T_{h0} = T_h\\
\dot\nu_h c_v (T_h - T_{h2}) = \Delta E_{TEG}\\
\dot\nu_h c_v (T_h - T_{h2}) = N_\Theta\\
\end{cases}\end{equation}$$

which supplement system (1) with temperature designations, as well as equation for `warm` part of hot loop: 
 $T_{h2} = T_w = T_h - \dfrac{\Delta E_{TEG}}{\dot\nu_h c_v}$.

Thus, we get 5 unknown pressures (used by $\Delta\nu$ functions), 3 unknown temperatures, for total of 8 variables; combining 3 mass balance equations with 3 heat balance equations and 2 mass conservation equations ($\Sigma\nu$), we get a system that is possibly solvable.

Unfortunately, the `min/max`-less purely mathematical closure looks not very analytically solvable (with two Lambert functions at  least due to enthropy calculation algorithm), and the programmatic iterative realisation stabilizes with results similar to mentioned above: with cold loop at 315K, core temp at 613K and warm loop at 554K (with total power of whopping 6MW for some unfathomable reason). This looks like an error in either assumptions or my implementation.
