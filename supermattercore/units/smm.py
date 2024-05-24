"""An attempt to rewrite model to be cantera-less and with proper physical units inlined everywhere"""

import astropy.units as u
from astropy.units import Quantity
import math
# import cantera as ct
import numpy as np
# import scipy.constants as C
import astropy.constants as C
from typing import NamedTuple, Optional, Literal
from typing_extensions import Self
from dataclasses import dataclass

import logging
SmLogger   = logging.getLogger('ct.sm')
PumpLogger = logging.getLogger('ct.sm.Pump')
# SmLogger.setLevel(logging.DEBUG)

hp =  u.def_unit('hitpts')
atm = u.def_unit('atm', 101325*u.Pa)
EER = u.def_unit('EER')

MACHINERY_TICK = 2 * u.s # seconds
PIPENET_TICK = 1 * u.s # seconds

EMITTER_BURST = 4 # or 3 tho?
POWERED_EMITTER_DMG = 142.22 * hp

EMITTER_POWER_CONS = 100 * u.kW # W
EMITTER_FULL_TIME = 6.4 * u.s # sec

PUMP_DEFAULT_VOLUME = 200 * u.L
SIMPLE_PIPE_VOLUME  = 70 * u.L

MIN_MOLES_PUMP = 0.01 * u.mol

PIPE_HE_SURFACE = 2 * u.m**2

CANISTER_VOLUME = 1000 * u.L
ENG_CANISTER_PRESSURE = 90 * atm

AIR_PCRIT = 3771 * u.kPa 
AIR_TCRIT = 132.65 * u.K

R_IDEAL_GAS = 8.31 * u.J / u.K / u.mol
AIR_VMOLAR: Quantity[u.L/u.mol] = R_IDEAL_GAS * AIR_TCRIT / AIR_PCRIT

T_COSMIC_MICROWAVE_BACKGROUND = 3.15 * u.K
T_MINIMAL = 2.7 * u.K

SPECIFIC_ENTROPY_VACUUM = 150000 * u.J / u.mol / u.K
IDEAL_GAS_ENTROPY_CONSTANT = 1164 * u.K*u.mol/u.l / (u.J*u.kg/u.mol**2)**(2/3) # * (u.mol * u.s / u.kg)**3 / (u.L) # -- looks like shit

CELL_VOLUME = 2500 * u.L

STEFAN_BOLTZMANN_CONSTANT = 5.6704e-8 * u.W / u.m**2 / u.K**4

ATMOS_PUMP_EFFICIENCY = 2.5

KnownGas = Literal['N2'] | Literal['O2'] | Literal['Ph']

@dataclass
class GasSpecies:
    name: KnownGas
    mu: Quantity[u.kg/u.mol]
    cv: Quantity[u.J/u.mol/u.K] = 20 * u.J/u.mol/u.K
    @classmethod
    def get(cls, n: KnownGas):
        match n:
            case 'O2': return GasO2
            case 'N2': return GasN2
            case 'Ph': return GasPh
            case _: raise NotImplementedError

GasO2 = GasSpecies('O2', 0.032 *u.kg/u.mol)
GasN2 = GasSpecies('N2', 0.028 *u.kg/u.mol)
GasPh = GasSpecies('Ph', 0.405 *u.kg/u.mol, 200*u.J/u.mol/u.K)

@dataclass
class GasMixture:
    """closely modeled after /xgm/gas_mixture"""
    T: Quantity[u.K]
    V: Quantity[u.L] = CELL_VOLUME
    nus: dict[KnownGas, Quantity[u.mol]] = None
    _V0 = None
    @classmethod
    def pure(cls, gas: KnownGas, T: Quantity[u.K], V: Quantity[u.L], P: Quantity[u.Pa]):
        nu = (P*V/R_IDEAL_GAS/T).to(u.mol)
        return GasMixture(
            T.to(u.K),
            V.to(u.L),
            {gas: nu}
        )
    @property
    def Nu(self) -> Quantity[u.mol]:
        if self.nus is None or len(self.nus) == 0: return 0*u.mol
        return sum(self.nus.values()).to(u.mol)
    @property
    def p(self) -> Quantity[u.Pa]:
        return (self.Nu*R_IDEAL_GAS*self.T/self.V).to(u.Pa)
    @property
    def C(self) -> Quantity[u.J/u.K]:
        c = 0*u.J/u.K
        if self.nus is None: return c
        for k in self.nus:
            c += (GasSpecies.get(k).cv * self.nus[k])
        return c.to(u.J/u.K)
    @property
    def cv_mole(self) -> Quantity[u.J/u.K/u.mol]:
        nu = self.Nu
        if nu == 0: return 0*u.J/u.K/u.mol
        return self.C/nu
    @property
    def S(self) -> Quantity[u.J/u.K/u.mol]:
        """specific entropy as per SS13 code"""
        if self.nus is None or len(self.nus) == 0:
            return SPECIFIC_ENTROPY_VACUUM
        Ss = sum([self.s_gas(k)*self.nus[k] for k in self.nus])
        return (Ss/self.Nu).to(u.J/u.K/u.mol)
    def s_gas(self, k: KnownGas) -> Quantity[u.J/u.K/u.mol]:
        """very magical"""
        if k not in self.nus or self.nus[k] == 0:
            return SPECIFIC_ENTROPY_VACUUM
        gas = GasSpecies.get(k)
        mu, cv = gas.mu, gas.cv
        safe_T = max(self.T, T_MINIMAL) #TCMB in code
        # partial_P = self.nus[k] * R_IDEAL_GAS * self.T / self.V
        return R_IDEAL_GAS * (
            math.log(
                IDEAL_GAS_ENTROPY_CONSTANT * self.V / (
                    self.nus[k] * safe_T) * (mu*cv*safe_T)**(2/3) + 1
            ) + 15
        )

    @property
    def thermal_E(self):
        return self.T * self.C
    @thermal_E.setter
    def thermal_E(self, value: Quantity[u.J]):
        self.T = value / self.C

    def __add__(self, other: Self | None) -> Self:
        """Consumes both mixtures and combines their volumes"""
        if other is None: return self
        V1 = self.V + other.V
        ret = GasMixture(self.T, V1, self.nus)
        ret += other
        return ret
    
    def __mul__(self, other: float) -> Self:
        """Returns a fractional fragment (for networks and subs)"""
        ret = GasMixture(
            self.T,
            self.V * other,
            {k: v * other for k,v in self.nus.items()} if self.nus else {}
        )
        ret._V0 = self.V
        return ret
    
    def sub(self, other: Quantity[u.mol]):
        """Removes a number of moles and returns them, 
            keeping the volume constant."""
        if other.value == 0: return None
        frac = other / self.Nu
        ret = self * frac
        self.nus = {k: v * (1-frac) for k,v in self.nus.items()}
        return ret
    
    def __iadd__(self, other: Self | None):
        """Merges matter into current, keeping V=const"""
        if other is None: return self
        Q1 = self.C*self.T + other.C*other.T
        snus = self.nus if self.nus else {}
        onus = other.nus if other.nus else {}
        gases = set(snus.keys()) | set(onus.keys())
        nus = {k: snus.get(k, 0*u.mol) + onus.get(k, 0*u.mol) 
               for k in gases}
        self.nus = nus
        self.T = Q1 / self.C
        return self

@dataclass
class PipeNetwork:
    gases: list[GasMixture]

    def equalize(self):
        if len(self.gases) <= 1: return
        ref = self.gases[0] + self.gases[1]
        for i in range(2, len(self.gases)):
            ref += self.gases[i]
        for i in range(len(self.gases)):
            self.gases[i] = ref * (ref.V / self.gases[i].V)

class HeatExchanger:
    n: int

    area = 2 * u.m**2
    qSol = 200 * u.W / u.m**2
    f_insolated = 0.04

    def __init__(self, n_space_tiles: int):
        self.n = n_space_tiles
    
    @property
    def V(self) -> Quantity[u.L]:
        return self.n * SIMPLE_PIPE_VOLUME

    def cool(self, gas: GasMixture):
        PSolEff = self.qSol * self.area * self.f_insolated * u.s
        for _ in range(self.n):
            eta = min(AIR_VMOLAR/(gas.V/gas.Nu), 1)
            Teff = (gas.T - T_COSMIC_MICROWAVE_BACKGROUND)
            dQ = PSolEff - u.s* self.area * STEFAN_BOLTZMANN_CONSTANT * Teff**4
            dT = dQ * eta / gas.C
            gas.T += dT
        return gas

    def approx_cool(self, gas: GasMixture) -> float:
        """Approximates without modifying"""
        PSolEff = self.qSol * self.area * self.f_insolated * self.n * u.s
        Td = gas.T - T_COSMIC_MICROWAVE_BACKGROUND
        Cve = gas.cv_mole * max(gas.Nu, gas.V / AIR_VMOLAR)
        k1 = 3*self.n*STEFAN_BOLTZMANN_CONSTANT*self.area*u.s/Cve
        Tcool = (k1 + Td**-3)**(-1/3.)
        dTrad = PSolEff/Cve
        dT = Td - Tcool - dTrad
        return dT

def calc_transfer_moles(source: GasMixture, sink: GasMixture, delta_p: Quantity[u.Pa]) -> Quantity[u.mol]:
    """Estimates moles to transfer, as per orig. code"""
    if source.T == 0 or source.Nu == 0: 
        return 0 # Nothing to pump
    outV = sink._V0 # full volume
    src_Nu = source.Nu # just the owned source -- `source total moles`
    airT = source.T
    if sink.T > 0 and sink.Nu > 0:
        # approximate mixture temp to re-calc
        eNu = delta_p*outV/sink.T/R_IDEAL_GAS # should be moles
        sinkC = sink.C # total heat cap
        srcC = source.C * eNu.to(u.mol) / src_Nu
        airT = (sink.T * sinkC + source.T*srcC) / (sinkC + srcC)
    return (delta_p*outV/airT/R_IDEAL_GAS).to(u.mol) # kMoles

def pump_gas(
        source: GasMixture, sink: GasMixture, 
        delta_nu: Quantity[u.mol], power_rating: Quantity[u.W],
        debug_name: str = ''
        ) -> tuple[Quantity[u.W], Quantity[u.mol]]:
    if source.Nu < MIN_MOLES_PUMP: 
        return 0*u.W, 0*u.mol # Nothing to move
    transferNu = min(delta_nu, source.Nu) if delta_nu is None else source.Nu
    # V Calculate specific power
    usedT = sink.T if sink.T > 0 else source.T
    srcS = source.S
    sinkS = sink.S
    specific_entropy = sinkS - srcS
    specific_power = -specific_entropy*usedT if specific_entropy<0 else 0*u.J/u.mol
    specific_power /= ATMOS_PUMP_EFFICIENCY # magical coeff.
    PumpLogger.info(f'{debug_name}: Source entropy {srcS:.6g} -> Sink entropy {sinkS:.6g}')
    PumpLogger.info(f'{debug_name}: Specific entropy change {specific_entropy:.6g}')
    PumpLogger.info(f'{debug_name}: Specific power {specific_power:.6g}')
    transferNu = min(transferNu, power_rating*u.s/specific_power) if specific_power.value>0 else transferNu
    transferNu = transferNu.to(u.mol)
    PumpLogger.info(f'{debug_name}: Transferred {transferNu:.6g} in fact')
    return specific_power*transferNu,transferNu


@dataclass
class Pump:
    name: str
    V_in : Quantity[u.L]= 200*u.L
    V_out: Quantity[u.L] = 200*u.L
    P_max:    Quantity[u.W] = 30*u.kW # .to(u.W)
    target_p: Quantity[u.Pa] = 15*u.MPa #.to(u.Pa)
    inlet : GasMixture|None = None
    outlet: GasMixture|None = None

    def reconnect_in(self, inlet: GasMixture):
        inlet.V += self.V_in
        self.inlet = inlet
    def connect(self, inlet: GasMixture, outlet: GasMixture):
        inlet.V += self.V_in
        self.inlet = inlet
        outlet.V += self.V_out
        self.outlet = outlet

    def pump(self): #, inlet: GasMixture, outlet: GasMixture):
        if self.inlet is None: return 0
        if self.outlet is None: return 0
        inlet,outlet = self.inlet, self.outlet
        delta_p = self.target_p - outlet.p
        if delta_p <= 10*u.Pa or inlet.T <= 0*u.K:
            return 0*u.W
        # Gets a own-volume segment -- or everything for env.
        source = inlet * (self.V_in/inlet.V)  if self.V_in !=0 else inlet *1
        sink = outlet * (self.V_out/outlet.V) if self.V_out!=0 else outlet*1
        # PumpLogger.info(source)
        delta_nu = calc_transfer_moles(source, sink, delta_p)
        w_draw, moles = pump_gas(source, sink, delta_nu, self.P_max, self.name)
        if moles > 0:
            qty = inlet.sub(moles)
            outlet += qty
        return w_draw

@dataclass
class Turbine:
    name: str
    V_in = 400*u.L
    V_out = 200*u.L
    inlet: GasMixture|None = None
    outlet: GasMixture|None = None
    def connect(self, inlet: GasMixture, outlet: GasMixture):
        inlet.V += self.V_in
        self.inlet = inlet
        outlet.V += self.V_out
        self.outlet = outlet
    def take_air(self) -> tuple[GasMixture|None, Quantity[u.J]]:
        dP = max(0*u.Pa, self.inlet.p-self.outlet.p-5*u.kPa)
        if dP < 5*u.kPa:
            return None, 0*u.J
        dnumax = dP * self.inlet.V /3/self.inlet.T/R_IDEAL_GAS
        E_i = min(dP*self.inlet.V, self.inlet.p*self.V_in)
        etaK,fv,kappa = 0.04, 0.2, 0.66
        etaV = etaK / kappa * (1 - fv**kappa)
        E = (etaV*E_i).to(u.kJ)
        vol_cap_frac = min (dP * self.inlet.V/3/self.inlet.p/self.V_in, 1)
        nu_can_use_max = self.V_in/self.inlet.V*self.inlet.Nu
        nu_use = min(dnumax, nu_can_use_max)
        gret = self.inlet.sub(nu_use)
        SmLogger.info(f'{self.name}: inP = {self.inlet.p.to(u.kPa):.4g}, outP={self.outlet.p.to(u.kPa):.4g}')
        SmLogger.info(f'{self.name}: dP={dP:.4g}, {dnumax=}, {E=}')
        SmLogger.info(f'{self.name}: {vol_cap_frac=:.4g}, frac nu{gret.Nu/nu_can_use_max:.4g}')
        return gret, E
    def push_air(self, new: GasMixture) -> GasMixture:
        self.outlet += new
        SmLogger.info(f'{self.name}: merged: {self.inlet.p=}, {self.outlet.p=}')

@dataclass
class TEG:
    hot: Turbine
    cold: Turbine
    eta:float = 0.65
    def process(self):
        cold, E1 = self.cold.take_air()
        hot , E2 = self.hot.take_air()
        if hot is None or cold is None:
            self.cold.push_air(cold)
            self.hot.push_air(hot)
            return E1+E2    
        dT = hot.T - cold.T
        C1 = cold.C*hot.C/(cold.C+hot.C)
        dE = dT*C1
        ET = self.eta*dE if dE>0 else 0
        dQ = dE - ET
        if dT>0:
            cold.T += dQ/cold.C
            hot.T -= dE/hot.C
        self.cold.push_air(cold)
        self.hot.push_air(hot)
        return ET+E1+E2

SM_EER_EMIT = 10 * hp / EER

SM_T_CRIT = 5000 * u.K
SM_LAMDA = 700 * (EER**(2/3))
SM_THERM_K = 1e4 * u.J / EER
SM_OXY_K = 15e3 * EER / u.mol
SM_USB_K = 1.5e3 * EER / u.mol
SM_MAX_DMG = 4.5 * hp
SM_ETA_GAS = 0.25
SM_ETA_N2 = 0.15
SM_ETA_RXN = 1.1
SM_DMG_TEMP = 150 * u.K / hp
SM_DMG_EER_VAC = 10 * EER / hp
# SM_REF_HEALTH = 1000 * hp
# SM_MAX_HEALTH = 1000 * hp
SM_REF_EER = 300 * EER
SM_REF_TEMP = 800 * u.K
SM_EER_OXY = 450 * EER
SM_EER_N2 = 200 * EER

@dataclass
class Supermatter:
    room: GasMixture
    eer: Quantity[EER] = 0*EER # `power`
    health = 1000 * hp # starting
    # --- stats
    def process(self):
        if not self.damage():
            self.react()
            self.decay()
    def damage(self):
        if self.room.p == 0:
            self.health -= max(self.eer-15*EER, 0)/ SM_DMG_EER_VAC
            return True
        d_i = (self.eer / SM_REF_EER) * SM_MAX_DMG
        d_T = (self.room.T - SM_T_CRIT) / SM_DMG_TEMP
        if d_T > d_i :
            self.health -= d_i
            return False
        if d_T < -SM_MAX_DMG: # healing
            self.health += SM_MAX_DMG
        else:
            self.health -= d_T
        self.health = max(self.health, 1000*hp)
        return self.health < 0
    def react(self):
        f_O = self.room.nus.get('O2', 0*u.mol) - SM_ETA_N2 * self.room.nus.get('N2',0*u.mol)
        f_O /= self.room.Nu
        f_O = max(min(f_O, 1), 0)
        P_E = (400 if f_O>0.8 else 250)*EER
        dP_rxn = (P_E/SM_LAMDA)**3 * (self.room.T / SM_REF_TEMP) * f_O
        self.eer += dP_rxn
        EP = SM_ETA_RXN*self.eer
        ET = (self.room.T - 273.15*u.K) * EER / u.K
        nu_Ex = max(EP/SM_USB_K, 0*u.mol)
        nu_O2 = max((EP+ET)/SM_OXY_K, 0*u.mol)
        self.room.nus['Ph'] = self.room.nus.get('Ph', 0*u.mol) + nu_Ex
        self.room.nus['O2'] = self.room.nus.get('O2', 0*u.mol) + nu_O2
        self.room.thermal_E += SM_THERM_K * EP
        if self.room.T > 1e4*u.K:
            print('OVERHEATED')
            self.room.T = 1e4*u.K
    def decay(self):
        self.eer -= (self.eer / SM_LAMDA)**3
