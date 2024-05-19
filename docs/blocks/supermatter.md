# Supermatter core
`/obj/machinery/power/supermatter` (`code/modules/supermatter/supermatter.dm`)

Tick is $\Delta\tau=2$ sec.

Brief stats:
- $T_c=5000$ K -- critical temperature (damage start);
- $\lambda = 700$ -- decay factor (??);
- $\eta_r = 1.1$ -- `reaction power modified` (??);
- $k_Q = 10^4$ -- J / (MeV/cm3) coeff `k_thermal`?

## States

In non-runaway conditions (i.e., when O2 effective is =0),
the internal power is stably decaying:
$\Delta P /\Delta\tau = -(P/\lambda)^3$;
easily enough we can write 
$P(\tau) = P_0 / \sqrt{\dfrac{P_0^2}{\lambda^3}t+1}$

Heat emitted by SMCore reaction is 
$\Delta Q = k_Q \eta_r P$ (if we ignore the waste gases).
Thus rough estimate of SMcore thermal power is ~5.5 kW/(MeV/cm3), i.e., you should multiply EER by 5.5 to get kW.

