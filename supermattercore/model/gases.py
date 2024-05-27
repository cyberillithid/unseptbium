"""Gases"""
from typing import NamedTuple, Optional, Literal
from typing_extensions import Self
from dataclasses import dataclass
import math
import logging 

GasLogger = logging.getLogger('ct.sm.Gas')

from .consts import *

KnownGas = Literal['N2', 'O2', 'Ex*', 'He', 'CO2', 'H2']
"""A typing hint"""

@dataclass
class GasSpecies:
    """/decl/material, basically"""
    name: KnownGas
    mu: float # Quantity[u.kg/u.mol]
    cv: float=20 # Quantity[u.J/u.mol/u.K] = 20 * u.J/u.mol/u.K
    is_fuel: bool = False
    is_oxy: bool  = False
    burn_product: KnownGas|None = None
    @classmethod
    def get(cls, n: KnownGas):
        return KnownGases.get(n)


KnownGases = {
    'N2'  : GasSpecies('N2',  0.028),
    'O2'  : GasSpecies('O2',  0.032, is_oxy=True),
    'Ex*' : GasSpecies('Ex*', 0.405, cv=200, is_fuel=True),
    'He'  : GasSpecies('He',  0.004, cv= 80),
    'CO2' : GasSpecies('CO2', 0.044, cv=30),
    'H2'  : GasSpecies('H2', 0.002, cv=100),
}
"""local repository of known gases"""


FLAMMABLE_GAS_MINIMUM_BURN_TEMPERATURE = (273.15 + 126) #*u.K
FIRE_FUEL_ENERGY_RELEASE = 866000 #* u.J/u.mol
MIN_REACT_FUEL = 0.005 #* u.mol
FUEL_STOICHIO = 2
OXY_STOICHIO = 3

@dataclass
class GasMixture:
    """closely modeled after /xgm/gas_mixture"""
    T: float #Quantity[u.K]
    V: float = CELL_VOLUME #Quantity[u.L] 
    nus: dict[KnownGas, float] = None  # Quantity[u.mol]
    _gr = 1
    _V0 = None
    @classmethod
    def pure(cls, gas: KnownGas, T:float,V:float,P:float):
             # T: Quantity[u.K], V: Quantity[u.L], P: Quantity[u.Pa]):
        nu = (P*V/R_IDEAL_GAS/T) #.to(u.mol)
        return GasMixture(
            T, #.to(u.K),
            V, #.to(u.L),
            {gas: nu}
        )
    @property
    def Nu(self) -> float:# Quantity[u.mol]:
        if self.nus is None or len(self.nus) == 0: return 0 #*u.mol
        return sum(self.nus.values()) #.to(u.mol)
    @property
    def X(self) -> dict[KnownGas, float]:
        if len(self.nus) == 0: return {}
        nu = self.Nu
        if nu == 0: # * u.mol: 
            return {}
        return {k: (v/nu) for k,v in self.nus.items()}

    @property
    def p(self) -> float: #Quantity[u.Pa]:
        return (self.Nu*R_IDEAL_GAS*self.T/self.V) #.to(u.Pa)
    @property
    def C(self) -> float: #Quantity[u.J/u.K]:
        c = 0 #*u.J/u.K
        if self.nus is None: return c
        for k in self.nus:
            c += (GasSpecies.get(k).cv * self.nus[k])
        return c #.to(u.J/u.K)
    @property
    def cv_mole(self) -> float: # Quantity[u.J/u.K/u.mol]:
        nu = self.Nu
        if nu == 0: return 0 #*u.J/u.K/u.mol
        return self.C/nu
    @property
    def S(self) -> float: # Quantity[u.J/u.K/u.mol]:
        """specific entropy as per SS13 code"""
        if self.nus is None or len(self.nus) == 0:
            return SPECIFIC_ENTROPY_VACUUM
        Ss = sum([self.s_gas(k)*self.nus[k] for k in self.nus])
        return (Ss/self.Nu) #.to(u.J/u.K/u.mol)
    def s_gas(self, k: KnownGas) -> float: # Quantity[u.J/u.K/u.mol]:
        """very magical"""
        if k not in self.nus or self.nus[k] == 0:
            return SPECIFIC_ENTROPY_VACUUM
        gas = GasSpecies.get(k)
        mu, cv = gas.mu, gas.cv
        safe_T = max(self.T, T_MINIMAL) #TCMB in code
        # partial_P = self.nus[k] * R_IDEAL_GAS * self.T / self.V
        return R_IDEAL_GAS * (
            math.log(
                IDEAL_GAS_ENTROPY_CONSTANT * self.V*1e3 / (
                    self.nus[k] * safe_T) * (mu*cv*safe_T)**(2/3) + 1
            ) + 15
        )

    @property
    def thermal_E(self):
        return self.T * self.C
    @thermal_E.setter
    def thermal_E(self, value: float): #Quantity[u.J]):
        self.T = value / self.C

    def __add__(self, other: Self | None) -> Self:
        """Consumes both mixtures and combines their volumes"""
        if other is None: return GasMixture(self.T, self.V, self.nus)
        V1 = self.V + other.V #.to(u.L)
        ret = GasMixture(self.T, V1, self.nus)
        ret += other
        return ret
    
    def __mul__(self, other: float) -> Self:
        """Returns a fractional fragment (for networks and subs)"""
        ret = GasMixture(
            self.T * 1,
            self.V * other,
            {k: v * other for k,v in self.nus.items() if v*other>=1e-12} if self.nus else {}
        )
        ret._V0 = self.V
        return ret
    
    def sub(self, other: float):#  Quantity[u.mol]):
        """Removes a number of moles and returns them, 
            keeping the volume constant."""
        if other == 0: return None
        frac = other / self.Nu
        ret = self * frac
        self.nus = {k: v * (1-frac) for k,v in self.nus.items() if v*(1-frac) >= 1e-12}
        return ret
    
    def __iadd__(self, other: Self | None):
        """Merges matter into current, keeping V=const"""
        if other is None: return self
        # Q1 = self.C*self.T + other.C*other.T if abs(self.T - other.T)>0.5*u.K else None
        if abs(self.T - other.T) > 0.5: #*u.K:
            Cs = self.C + other.C
            if Cs != 0: #*u.J/u.K:
                self.T = (self.T*self.C + other.T * other.C)/Cs
        snus = self.nus if self.nus else {}
        onus = other.nus if other.nus else {}
        gases = set(snus.keys()) | set(onus.keys())
        nus = {k: (snus.get(k, 0) + onus.get(k, 0)) # .to(u.mol)
               for k in gases}
        self.nus = nus
        return self

    def react(self):
        if self.T < FLAMMABLE_GAS_MINIMUM_BURN_TEMPERATURE: return
        fuel, oxy = 0,0 #*u.mol, 0*u.mol
        for k,v in self.nus.items():
            g = GasSpecies.get(k)
            if g.is_fuel: fuel += v
            if g.is_oxy:  oxy += v
        if fuel <= MIN_REACT_FUEL: return 
        start_E = self.thermal_E
        rxn_max = min(fuel, oxy*(FUEL_STOICHIO/OXY_STOICHIO))
        firelevel = self.calc_firelevel(fuel, oxy, rxn_max, self.V)
        min_burn = 0.3 * self.V/CELL_VOLUME # *u.mol
        tot_rxn_prog = min(max(min_burn, firelevel*fuel), fuel)
        used_fuel = min(tot_rxn_prog, rxn_max)
        used_oxy = used_fuel*(OXY_STOICHIO/FUEL_STOICHIO)
        used_fuel = min(used_fuel, fuel)
        self.remove(used_oxy, is_oxy=True)
        rm_fuel = self.remove(used_fuel, is_fuel=True).nus
        add_els = [(GasSpecies.get(k).burn_product,v) for k,v in rm_fuel.items()]
        add_nus = {k: v for (k,v) in add_els if k is not None}
        if len(add_nus)>0:
            self += GasMixture(self.T, 0, add_nus) # L
        self.thermal_E = start_E + FIRE_FUEL_ENERGY_RELEASE*used_fuel
        self.nus = {k: v for k,v in self.nus.items() if v>=1e-10}

    def remove(self, qty: float, is_oxy=False, is_fuel=False) -> Self:
        # qty: Quantity[u.mol]
        ref = {}
        # print(qty)
        for k,v in self.nus.items():
            g = GasSpecies.get(k)
            if (is_oxy and g.is_oxy) or (is_fuel and g.is_fuel):
                ref[k] = v #.to(u.mol).value
        refqty = sum(ref.values())
        if refqty < qty: 
            GasLogger.warn('asked to remove by flag too much')
            qty = refqty
        self.nus.update({
            k: max(0,self.nus[k] - (v*qty/refqty)) for k,v in ref.items()
            })
        
        return GasMixture(
            self.T, self.V, 
            nus={k: v*qty/refqty for k,v in ref.items()}
            )
    
    def calc_firelevel(self, fuel: float, oxy: float, rxn_lim: float, V: float):
        #(self, fuel: Quantity[u.mol], oxy: Quantity[u.mol], rxn_lim: Quantity[u.mol], V: Quantity[u.L]):
        """/datum/gas_mixture/proc/calculate_firelevel(total_fuel, total_oxidizers, reaction_limit, gas_volume)"""
        total = fuel+oxy
        active = (1 + OXY_STOICHIO/FUEL_STOICHIO)*rxn_lim
        firelevel = 0
        if total>0: #*u.mol:
            damping = min (1, active / (self.Nu/self._gr))
            damping = damping*(2 - damping)
            mix_mult = 1 / (1 + (5 * ((fuel / total) ** 2)))
            firelevel = mix_mult * damping
        return max(firelevel, 0)

    def _repr_latex_(self) -> str:
        s = f'Gas mixture:\n T = {self.T:.2f} K, '
        s += f' V =  {self.V*1e3:.2f} L;\n'
        s += f' p =  {self.p/1e3:.2f} kPa, '
        s += f' nu = {self.Nu   :.2f} mol\n'
        s += '(Mole fractions:\n'
        s += ', '.join([f'{k}: {100*v:.2f}%' for k,v in self.X.items()])
        s += ')'
        return s

