"""Pipe networks and heating pipes."""

from dataclasses import dataclass
import logging 

from .consts import *
from .gases import GasMixture, KnownGas, GasSpecies

PipeLogger = logging.getLogger('ct.sm.Pipe')

SIMPLE_PIPE_VOLUME  = 0.070  #* u.m**3

PIPE_HE_SURFACE = 2 #* u.m**2

AIR_PCRIT = 3771e3 #* u.Pa 
AIR_TCRIT = 132.65 #* u.K
AIR_VMOLAR= R_IDEAL_GAS * AIR_TCRIT / AIR_PCRIT
#: Quantity[u.L/u.mol] 

CANISTER_VOLUME = 1# *u.m**3
NORM_CANISTER_PRESSURE = 45 * 101325
ENG_CANISTER_PRESSURE = 2 * NORM_CANISTER_PRESSURE # Pa  atm

def mk_canister(gas: KnownGas, engine=False, cryo=False):
    p0 = NORM_CANISTER_PRESSURE * (2 if engine else 1) * (20/45 if cryo else 1)
    T0 = 80 if cryo else 293.15
    # pV= nuRT
    return GasMixture.pure(gas, T0, CANISTER_VOLUME, p0)
    # nu = p0 * CANISTER_VOLUME / R_IDEAL_GAS / T0

@dataclass
class PipeNetwork:
    gases: dict[str, GasMixture]

    def equalize(self):
        if len(self.gases) <= 1: return
        ks = [*self.gases.keys()]
        ref = self.gases[ks[0]] + self.gases[ks[1]]
        Q0 = self.gases[ks[0]].thermal_E + self.gases[ks[1]].thermal_E 
        for i in ks[2:]:
            ref = ref + self.gases[i]
            Q0 += self.gases[i].thermal_E
        ref.T = Q0 / ref.C if ref.C > 0 else ref.T
        for i in ks:
            self.gases[i] = ref * ((self.gases[i].V) / ref.V) if self.gases[i] else None
    
    def add(self, name: str, V: float, T0=T_MINIMAL,nus=None):
        #  Quantity[u.L]
        if name in self.gases:
            self.gases[name].V = V
        else:
            if V == 0:
                self.gases[name] = None
            else:
                self.gases[name] = GasMixture(T0, V, nus if nus else {})
    
    @property
    def V(self) -> float: # Quantity[u.L]:
        return sum([k.V for k in self.gases.values() if k])
    def _repr_latex_(self) -> str:
        s = "Pipe network with members {" + ', '.join([
            k + f" (V = {v.V*1e3:.0f} L ) " 
            for k,v in self.gases.items()
            ]) + "}"
        if len(self.gases) > 0:
            s += f", total volume {self.V*1e3:.0f} L"
            k0 = [*self.gases.keys()][0]
            s += ";\n- first member props: " + self.gases[k0]._repr_latex_()[12:]
        return s
class HeatExchanger:
    n: int

    area = 2   #* u.m**2
    qSol = 200 # * u.W / u.m**2
    f_insolated = 0.04

    def __init__(self, n_space_tiles: int):
        self.n = n_space_tiles
    
    @property
    def V(self) -> float: # Quantity[u.L]:
        return self.n * SIMPLE_PIPE_VOLUME

    def cool(self, gas: GasMixture):
        PSolEff = self.qSol * self.area * self.f_insolated # * u.s
        if gas.C == 0: return gas #* u.J/u.K # nothing to do here
        t0 = gas.T * 1
        for _ in range(self.n):
            eta = min(AIR_VMOLAR/(gas.V/gas.Nu), 1)
            Teff4 = (gas.T**4 - T_COSMIC_MICROWAVE_BACKGROUND**4)
            dQ = PSolEff - self.area * STEFAN_BOLTZMANN_CONSTANT * Teff4
            dT = dQ * eta / gas.C
            gas.T = gas.T + dT
            # SmLogger.info(f'{(eta*dQ).to(u.J)}, {gas.T}')
        PipeLogger.info(f'last {(eta*dQ)} J; Per {self.n} pipes {t0:.2f} -> {gas.T:.2f} ({t0-gas.T:.4g} loss)')
        return gas

    def approx_cool(self, gas: GasMixture) -> float:
        """Approximates without modifying"""
        PSolEff = self.qSol * self.area * self.f_insolated * self.n #* u.s
        Td = gas.T - T_COSMIC_MICROWAVE_BACKGROUND
        if gas.C == 0: return 0 #* u.J/u.K: return 0 # nothing to do here
        Cve = gas.cv_mole * max(gas.Nu, gas.V / AIR_VMOLAR)
        k1 = 3*self.n*STEFAN_BOLTZMANN_CONSTANT*self.area/Cve # *u.s
        Tcool = (k1 + Td**-3)**(-1/3.)
        dTrad = PSolEff/Cve
        dT = Td - Tcool - dTrad
        return dT
