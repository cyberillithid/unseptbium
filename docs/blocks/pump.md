# Pump

Tick is $\Delta\tau=2$ sec (`SSmachinery`)

High_pump:
- $\mathcal{P}_i=450$ W -- idle power cons (150 for normal)
- $N_{\max} = 45$ kW -- power rating (30 kW for normal)
- $P_{\max} = 15$ MPa  (for both)
- $V_i = V_o = 0.2$ m3

Also: vent_pump has $N_{\max}$=30kW; injector has rating of 45 kW and 0.7 m3

$P^*$ -- target pressure;

$\Delta P = P^* - P_o$

$\Delta \nu = F(?_i,?_o,\Delta P, V_o^*)$
### calculate transfer moles

$\nu_{in} = \dfrac{P_iV_i}{RT_i}$ -- `source.moles`

$\nu_a = \Delta P V_o^* / RT_o $ -- approx moles

$C_o = c_{v,m,o} P_oV_o/RT_o$; $C_1 = c_{v,m,i} \nu_a$

$T_{1} = (T_oC_o+T_iC_1)/(C_o+C_1)$ -- approx temo

$\Delta\nu_1 =  \Delta P V_o^* / RT_1$

### pump gas
$\Delta\nu_2 = \min(\nu_in,\Delta\nu_1)$

$s_j = R (15 + \log (1+(\mu_j c_{vj} T')^{2/3}\dfrac{C_S V}{T' \nu_j}))$

$S = \sum X_j s_j$; 


$S_s = S_o - S_i$; $N_s' = -S_sT_o$ if $S_s<0$ else 0

$N_s = N_s'/2.5$ -- specific power ??...

$\Delta\nu = \min(\Delta\nu_2, N_{\max}/N_s)$ if N_s > 0.

Assuming pure N2 and T>3K, $X_{[N_2]}=1$,

$S = s_{[N_2]}=15R + R \log(1  + (cT)^{2/3} \dfrac{C_S V}{\nu T})$,
$c = \mu_{N_2}c_{v,m}$ (J/(mol K) * (kg/mol)), (cT: J kg /mol2 K);

$S_s = R \log \dfrac{1 + (cT_o)^{2/3} R/p_o}{1 + (cT_i)^{2/3} R/p_i}$

$c_1 = R (\mu_{N_2} c_{v,m})^{2/3}$;

$N_s = -T_o\dfrac{R}{2.5} \log\dfrac{1+c_1T_o^{2/3}/p_o}{1+c_1T_i^{2/3}/p_i}$
