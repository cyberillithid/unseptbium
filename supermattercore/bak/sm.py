import cantera as ct
import numpy as np
import scipy.constants as C
from typing import NamedTuple, Optional
from dataclasses import dataclass

MACHINERY_TICK = 2 # seconds

POWERED_EMITTER_DMG = 142.22
EMITTER_POWER_CONS = 100e3 # W
EMITTER_FULL_TIME = 6.4 # sec
EMITTER_BURST = 4 # or 3 tho?

def gas_qt_ixr(q: ct.Quantity, mul: float=1) -> dict[str, float]:
    return dict(zip(q.species_names, q.X * mul))

def gas_qt_epr(q: ct.Quantity) -> float:
    if q.moles == 0: return float('nan')
    return q.moles / q.volume**2 * 5**4 / 2.31

class SupermatterCore:
    gas_eff = 0.25 # 1
    decay = 700   # Pa^[2/3]??? Used in radiation loss 
    reaction = 1.1 # `reaction_power_modified`
    k_thermal = 10e3 # J / (MeV/cm3)
    k_exotic = 15e2  
    k_o2 = 15e3
    k_n2_ret = 0.15 # 1; N2 retardation coef.
    k_damage = 4.5e-3 # scaled by 1000
    T_crit = 5000 # K

    damage: float # scaled by 1000 to 0..1 interval; integrity = 1-damage
    power: float # EER, MeV / cm3 -- meaningless!

    def __init__(self, damage = 0, power = 0):
        self.damage = damage
        self.power = power

    def powerup(self, dmg: float):
        self.power += dmg/10.
    
    def tick_damage(self, N: float, T: float) -> bool:
        """N -- moles; T -- temperature"""
        if N == 0:
            self.damage += max((self.power - 15)/10e3, 0)
            return True
        max_dmg = SupermatterCore.k_damage * self.power/300
        techn_dmg = (T - SupermatterCore.T_crit) / 150e3
        self.damage = max(0, self.damage + np.clip(techn_dmg, -SupermatterCore.k_damage, max_dmg))
        return False

    def react(self, o2: float, n2: float, T: float) -> float:
        o2_eff = max(o2 - SupermatterCore.k_n2_ret * n2, 0)
        P_eq = 250 if o2_eff <= 0.8 else 400 # Target equilibrium power?
        fT = T/800 * (P_eq / SupermatterCore.decay)**3
        return o2_eff * fT + self.power

    def tick(self, air: ct.Quantity) -> ct.Quantity:
        if self.damage >= 1:
            air.moles = 0 # exploded
            air.phase.TP = 2.7, 1 # vacuum
            return air
        if self.tick_damage(air.moles, air.T):
            return air
        X = gas_qt_ixr(air)
        P = self.react(X['O2'], X['N2'], air.T)
        
        E = P * SupermatterCore.reaction # ????
        frac_taken = 2.5 * SupermatterCore.gas_eff / air.volume
        new_N = air.moles * frac_taken
        new_X = gas_qt_ixr(air, new_N) # removing a quarter...
        new_Ph = E / SupermatterCore.k_exotic
        new_O2 = max((E + air.T - 273.15) / SupermatterCore.k_o2, 0)
        new_X['Ph'] += new_Ph
        new_X['O2'] += new_O2
        new_N += new_Ph + new_O2

        dQ = E * SupermatterCore.k_thermal

        reag = ct.Solution(air.phase.source)
        reag.TPX = air.T, air.P, new_X
        reag.UV = reag.u + dQ, air.phase.volume_mass
        if reag.T > 10e3: reag.TP = 10e3, reag.P

        reag_q = ct.Quantity(reag, moles=new_N)
        reag_q.constant = 'UV'
        air *= (1 - frac_taken)
        air.constant = 'UV'
        air += reag_q
        
        self.power = P - (P/SupermatterCore.decay)**3
        return air        
    

PUMP_DEFAULT_VOLUME = 0.2 # m^3 == 200 L
SIMPLE_PIPE_VOLUME  = 0.07 # m^3 = 70 L

PIPE_HE_SURFACE = 2 # m^2
SPACE_THC = 1

CANISTER_VOLUME = 1 # m^3 = 1000 L
ENG_CANISTER_PRESSURE = 90*101325 # Pa

AIR_PCRIT = 3771e3 # Pa
AIR_TCRIT = 132.65 # K, 

R_IDEAL_GAS = 8.31
AIR_VMOLAR = R_IDEAL_GAS * AIR_TCRIT / AIR_PCRIT

T_COSMIC_MICROWAVE_BACKGROUND = 3.15
T_MINIMAL = 2.7

class HeatExchanger:
    """
    
    Heat-exchangers work over the total pipeline, piece by piece.
It takes `density` = $\nu/V = 1 /V_m = p/RT$;

calculates thermal-conductivity multiplier $\lambda$ as $\min(1, \dfrac{\nu}{V} / \dfrac{p_{c, air}}{RT_{c, air}})$
= $\min(1, \dfrac{1/V_m}{1/V_{m,c,air}}) = \min(1, \dfrac{V_{m,c,air}}{V_m}) $

$S = 2\,\mathrm{m}^2$; $f = 0.04$ => $S_{eff} = 0.08$?

$\Delta Q = \Delta Q_\odot - \Delta Q_r$; 

solar radiation surplus is $\Delta Q_\odot = E_\odot fS \lambda$, $E_\odot = 200 \mathrm{W}/\mathrm{m}^2$;

radiation loss is $\Delta Q_r = S \sigma \lambda (T - T_0)^4$, $\sigma$ -- Стефана-Больцмана, $T_0 = 3.15$ K

    """
    space_no: int
    plate_no: int
    T_plates: list[float]

    def __init__(self, space: int, plate: int, T_plate: float = 293.15):
        self.space_no = space
        self.plate_no = plate
        self.T_plates = [T_plate]*plate

    def cool_iter(self, q: ct.Quantity):
        """Practically always swappable to cool-lin (even for 3K inital discrepancy <30K, <7 kPa)"""
        s = q.phase # ct.Solution(q.phase.source)
        u1 = s.u
        for j in range(self.space_no):
            lamda = min( AIR_VMOLAR * 1e3 / s.volume_mole , 1 ) 
            dQ1 = 200 * 0.08 - 2 * C.sigma * (s.T - T_COSMIC_MICROWAVE_BACKGROUND)**4
            dQ = (dQ1 * lamda) / q.mass
            s.UV = s.u + dQ, s.v
        u2 = s.u
        for j in range(self.plate_no):
            dT = s.T - self.T_plates[j]
            cp_segm = (0.07 / s.volume_mole * s.cv_mole)
            Q2 = 0.4 * dT * (cp_segm * 10e3 / (cp_segm + 10e3))
            self.T_plates[j] += Q2 / 10e3
            s.UV = s.u - (Q2 / q.mass), s.v
        u3 = s.u
        return ct.Quantity(s, mass = q.mass), (u2-u1), (u3-u2)
    
    def cool_lin(self, q: ct.Quantity):
        s = q.phase
        t_plate = np.sum(self.T_plates) / self.plate_no
        lamda = min( AIR_VMOLAR * 1e3 / q.volume_mole , 1 ) 
        dQ1 = 200 * 0.08 - 2 * C.sigma * (q.T - T_COSMIC_MICROWAVE_BACKGROUND)**4
        dQ = (dQ1 * lamda) / q.mass * self.space_no
        dT = s.T - t_plate
        cp_segm = (0.07 / s.volume_mole * s.cv_mole)
        Q2 = 0.4 * dT * (cp_segm * 10e3 / (cp_segm + 10e3))
        t_plate += Q2 / 10e3
        self.T_plates = [t_plate] * self.plate_no
        dQ2 = - (Q2 / q.mass) * self.plate_no
        s.UV = s.u + dQ + dQ2, s.v
        return ct.Quantity(s, mass = q.mass), dQ, dQ2


IDEAL_GAS_ENTROPY_CONSTANT = 1164    # (mol^3 * s^3) / (kg^3 * L)
@dataclass
class PipeNet: # (NamedTuple):
    """Simulates pipe network; qty == None if there is no gas in it.
    Also simulates segments of it, if V != qty.volume."""
    V: float
    qty: Optional[ct.Quantity]
    @property
    def T(self) -> float:
        return 0 if self.qty is None else self.qty.phase.T
    @property
    def Nu(self) -> float:
        """kmole"""
        return 0 if self.qty is None else self.qty.moles * self.qty.volume / self.V
    @property
    def Cp(self) -> float:
        if self.qty is None: return 0
        return self.Nu * self.qty.phase.cp_mole
    @property
    def P(self) -> float:
        return 0 if self.qty is None else self.qty.phase.P
    @property
    def S_13(self) -> float:
        """Uses the approximate entropy from xgm-gas-mixture"""
        if self.qty is None: return 150e3
        ret = 0
        T = max(self.T, T_MINIMAL)
        nu = self.Nu * 1e3
        for (i, x) in enumerate(self.qty.phase.X):
            sp = self.qty.phase.species(i)
            scv = sp.thermo.cp(self.T)/1000 - 8.314
            mu = (self.qty.phase.atomic_weights @ np.array([sp.composition.get(x, 0) for x in self.qty.phase.element_names])) * 1e-3
            snu = nu * x
            if snu > 0:
                ss = R_IDEAL_GAS * (np.log( 1 +
                    (mu * scv * T)**(2/3) * IDEAL_GAS_ENTROPY_CONSTANT * self.V * 1e3 / T / snu
                ) + 15)
                ret += snu * ss
        return ret / nu
    def sub(self, moles: float) -> Optional[ct.Quantity]:
        if self.qty is None:
            return None
        q = ct.Quantity(self.qty.phase, moles=moles)
        self.qty.moles -= moles
        return q
    def __add__(self, q: Optional[ct.Quantity]):
        if self.qty is None:
            self.qty = q
            return self
        if q is None:
            return self
        




def calc_transfer_moles(inlet: PipeNet, outlet: PipeNet, deltaP: float, defV: float) -> float:
    """Imitates calc_transfer_moles of original code"""
    if inlet.T == 0 or inlet.Nu == 0: 
        return 0 # Nothing to pump
    outV = outlet.V + defV
    src_tot_Nu = inlet.Nu/inlet.V*defV # Scaled for 200L or something
    airT = inlet.T
    if outlet.T > 0 and outlet.Nu > 0:
        # approximate mixture temp to re-calc
        eNu = deltaP*outV/outlet.T/R_IDEAL_GAS/1e3 # kMoles
        sinkCp = outlet.Cp
        srcCp = inlet.Cp * eNu / src_tot_Nu
        airT = (outlet.T * sinkCp + inlet.T*srcCp) / (sinkCp + srcCp)
    return deltaP*outV/airT/R_IDEAL_GAS/1e3 # kMoles

MIN_KMOLES_PUMP = 1e-5 # 0.01 moles
ATMOS_PUMP_EFFICIENCY = 2.5 # KPD 40%

def pump_gas(inlet: PipeNet, outlet: PipeNet, transferNu: float = 0, availablePower: float = 0):
    """Imitates pump_gas of original code; modifies nets, returns power draw in W"""
    if inlet.Nu < MIN_KMOLES_PUMP: return 0
    transferNu = min(transferNu, inlet.Nu) if transferNu else inlet.Nu
    specific_power = specific_power()

class Pump:
    idle_power = 150 # W?
    power_rating = 30e3 # W
    target_P = 15e6 # Pa
    own_V = PUMP_DEFAULT_VOLUME # 0.2 m^3

    def process(self, inlet: PipeNet, outlet: PipeNet) -> float:
        """Modifies inlet and outlet, returns power draw"""
        deltaP = self.target_P - outlet.phase.P 
        if (deltaP > 10) and (inlet.phase.T > 0 ): # or `inlet.mass > 0`?
            deltaNu = calc_transfer_moles(inlet, outlet, deltaP, self.own_V)

    pass