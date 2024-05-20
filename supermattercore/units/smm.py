"""An attempt to rewrite model to be cantera-less and with proper physical units inlined everywhere"""

import astropy.units as u
from astropy.units import Quantity
import math
# import cantera as ct
import numpy as np
import scipy.constants as C
from typing import NamedTuple, Optional, Literal
from dataclasses import dataclass

hp =  u.def_unit('hitpts')
atm = u.def_unit('atm', 101325*u.Pa)

MACHINERY_TICK = 2 * u.s # seconds
PIPENET_TICK = 1 * u.s # seconds

EMITTER_BURST = 4 # or 3 tho?
POWERED_EMITTER_DMG = 142.22 * hp

EMITTER_POWER_CONS = 100 * u.kW # W
EMITTER_FULL_TIME = 6.4 * u.s # sec

PUMP_DEFAULT_VOLUME = 200 * u.L
SIMPLE_PIPE_VOLUME  = 70 * u.L

PIPE_HE_SURFACE = 2 * u.m**2

CANISTER_VOLUME = 1000 * u.L
ENG_CANISTER_PRESSURE = 90 * atm

AIR_PCRIT = 3771 * u.kPa 
AIR_TCRIT = 132.65 * u.K

R_IDEAL_GAS = 8.31 * u.J / u.K / u.mol
AIR_VMOLAR: Quantity[u.L/u.mol] = R_IDEAL_GAS * AIR_TCRIT / AIR_PCRIT

T_COSMIC_MICROWAVE_BACKGROUND = 3.15 * u.K
T_MINIMAL = 2.7 * u.K

SPECIFIC_ENTROPY_VACUUM = 15000 * u.J / u.mol / u.K
IDEAL_GAS_ENTROPY_CONSTANT = 1164 * u.mol * u.K**(1/3) / u.L # * (u.mol * u.s / u.kg)**3 / (u.L) # -- looks like shit

CELL_VOLUME = 2500 * u.L

KnownGas = Literal['N2'] | Literal['O2'] | Literal['Ph']

@dataclass
class GasSpecies:
    name: KnownGas
    mu: Quantity[u.kg/u.mol]
    cv: Quantity[u.J/u.mol/u.K] = 20
    @classmethod
    def get(cls, n: KnownGas):
        match n:
            case 'O2': return GasO2
            case 'N2': return GasN2
            case 'Ph': return GasPh
            case _: raise NotImplementedError

GasO2 = GasSpecies('O2', 0.032)
GasN2 = GasSpecies('N2', 0.028)
GasPh = GasSpecies('Ph', 0.405, 200)

@dataclass
class GasMixture:
    T: Quantity[u.K]
    V: Quantity[u.L] = CELL_VOLUME
    nus: dict[KnownGas, Quantity[u.mol]] = None
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
    def P(self) -> Quantity[u.Pa]:
        return (self.Nu*R_IDEAL_GAS*self.T/self.V).to(u.Pa)
    @property
    def C(self) -> Quantity[u.J/u.K]:
        c = 0*u.J/u.K
        if self.nus is None: return c
        for k in self.nus:
            c += GasSpecies.get(k).cv * self.nus[k]
        return c
    @property
    def S(self) -> Quantity[u.J/u.K/u.mol]:
        """specific entropy as per SS13 code"""
        if self.nus is None or len(self.nus) == 0:
            return SPECIFIC_ENTROPY_VACUUM
        Ss = sum([self.s_gas(k)*self.nus[k] for k in self.nus])
        return Ss/self.Nu
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
