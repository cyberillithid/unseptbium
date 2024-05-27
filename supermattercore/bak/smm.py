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

KnownGas = Literal['N2'] | Literal['O2'] | Literal['Ex*'] | Literal['He']

@dataclass
class GasSpecies:
    name: KnownGas
    mu: Quantity[u.kg/u.mol]
    cv: Quantity[u.J/u.mol/u.K] = 20 * u.J/u.mol/u.K
    is_fuel: bool = False
    is_oxy: bool  = False
    burn_product: KnownGas|None = None
    @classmethod
    def get(cls, n: KnownGas):
        match n:
            case 'O2': return GasO2
            case 'N2': return GasN2
            case 'Ex*': return GasPh
            case 'He': return GasHe
            case _: raise NotImplementedError

GasO2 = GasSpecies('O2',  0.032 *u.kg/u.mol, is_oxy=True)
GasN2 = GasSpecies('N2',  0.028 *u.kg/u.mol)
GasHe = GasSpecies('He',  0.004 *u.kg/u.mol, cv= 80*u.J/u.mol/u.K)
GasPh = GasSpecies('Ex*', 0.405 *u.kg/u.mol, cv=200*u.J/u.mol/u.K, is_fuel=True)

FLAMMABLE_GAS_MINIMUM_BURN_TEMPERATURE = (273.15 + 126)*u.K
FIRE_FUEL_ENERGY_RELEASE = 866000 * u.J/u.mol
MIN_REACT_FUEL = 0.005 * u.mol

@dataclass
class GasMixture:
    """closely modeled after /xgm/gas_mixture"""
    T: Quantity[u.K]
    V: Quantity[u.L] = CELL_VOLUME
    nus: dict[KnownGas, Quantity[u.mol]] = None
    _gr = 1
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
    def X(self) -> dict[KnownGas, float]:
        if len(self.nus) == 0: return {}
        nu = self.Nu
        if nu == 0 * u.mol: return {}
        return {k: (v/nu).value for k,v in self.nus.items()}

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
        if other is None: return GasMixture(self.T, self.V, self.nus)
        V1 = self.V.to(u.L) + other.V.to(u.L)
        ret = GasMixture(self.T, V1, self.nus)
        ret += other
        return ret
    
    def __mul__(self, other: float) -> Self:
        """Returns a fractional fragment (for networks and subs)"""
        ret = GasMixture(
            self.T * 1,
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
        # Q1 = self.C*self.T + other.C*other.T if abs(self.T - other.T)>0.5*u.K else None
        if abs(self.T - other.T) > 0.5*u.K:
            Cs = self.C + other.C
            if Cs != 0 *u.J/u.K:
                self.T = (self.T*self.C + other.T * other.C)/Cs
        snus = self.nus if self.nus else {}
        onus = other.nus if other.nus else {}
        gases = set(snus.keys()) | set(onus.keys())
        nus = {k: (snus.get(k, 0*u.mol) + onus.get(k, 0*u.mol)).to(u.mol)
               for k in gases}
        self.nus = nus
        return self

    def react(self):
        if self.T < FLAMMABLE_GAS_MINIMUM_BURN_TEMPERATURE: return
        fuel, oxy = 0*u.mol, 0*u.mol
        for k,v in self.nus.items():
            g = GasSpecies.get(k)
            if g.is_fuel: fuel += v
            if g.is_oxy:  oxy += v
        if fuel <= MIN_REACT_FUEL: return 
        start_E = self.thermal_E
        rxn_max = min(fuel, oxy*(FUEL_STOICHIO/OXY_STOICHIO))
        firelevel = self.calc_firelevel(fuel, oxy, rxn_max, self.V)
        min_burn = 0.3*u.mol * self.V/CELL_VOLUME
        tot_rxn_prog = min(max(min_burn, firelevel*fuel), fuel)
        used_fuel = min(tot_rxn_prog, rxn_max)
        used_oxy = used_fuel*(OXY_STOICHIO/FUEL_STOICHIO)
        used_fuel = min(used_fuel, fuel)
        self.remove(used_oxy, is_oxy=True)
        rm_fuel = self.remove(used_fuel, is_fuel=True).nus
        add_els = [(GasSpecies.get(k).burn_product,v) for k,v in rm_fuel.items()]
        add_nus = {k: v for (k,v) in add_els if k is not None}
        if len(add_nus)>0:
            self += GasMixture(self.T, 0*u.L, add_nus)
        self.thermal_E = start_E + FIRE_FUEL_ENERGY_RELEASE*used_fuel

    def remove(self, qty: Quantity[u.mol], is_oxy=False, is_fuel=False) -> Self:
        ref = {}
        # print(qty)
        for k,v in self.nus.items():
            g = GasSpecies.get(k)
            if (is_oxy and g.is_oxy) or (is_fuel and g.is_fuel):
                ref[k] = v.to(u.mol).value
        refqty = sum(ref.values())
        if refqty * u.mol < qty: 
            SmLogger.warn('asked to remove by flag too much')
            qty = refqty
        self.nus.update({
            k: max(0*u.mol,self.nus[k] - (v*qty/refqty)) for k,v in ref.items()
            })
        
        return GasMixture(
            self.T, self.V, 
            nus={k: v*qty/refqty for k,v in ref.items()}
            )
    
    def calc_firelevel(self, fuel: Quantity[u.mol], oxy: Quantity[u.mol], rxn_lim: Quantity[u.mol], V: Quantity[u.L]):
        """/datum/gas_mixture/proc/calculate_firelevel(total_fuel, total_oxidizers, reaction_limit, gas_volume)"""
        total = fuel+oxy
        active = (1 + OXY_STOICHIO/FUEL_STOICHIO)*rxn_lim
        firelevel = 0
        if total>0*u.mol:
            damping = min (1, active / (self.Nu/self._gr))
            damping = damping*(2 - damping)
            mix_mult = 1 / (1 + (5 * ((fuel / total) ** 2)))
            firelevel = mix_mult * damping
        return max(firelevel, 0)

    def _repr_latex_(self) -> str:
        s = 'Gas mixture:\n $T = ' + self.T.to(u.K).round(2)._repr_latex_().strip('$') + '$, '
        s += ' $V = ' + self.V.to(u.L).round(2)._repr_latex_().strip('$') + '$;\n'
        s += ' $p = ' + self.p.to(u.kPa).round(2)._repr_latex_().strip('$') + '$, '
        s += ' $\\nu = ' + self.Nu.to(u.mol).round(2)._repr_latex_().strip('$') + '$\n'
        s += '(Mole fractions:\n'
        s += ', '.join([f'{k}: {100*v:.2f}%' for k,v in self.X.items()])
        s += ')'
        return s

FUEL_STOICHIO = 2
OXY_STOICHIO = 3

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
            self.gases[i] = ref * ((self.gases[i].V.to(u.L)) / ref.V) if self.gases[i] else None
    
    def add(self, name: str, V: Quantity[u.L], T0=T_MINIMAL,nus=None):
        if name in self.gases:
            self.gases[name].V = V
        else:
            if V == 0:
                self.gases[name] = None
            else:
                self.gases[name] = GasMixture(T0, V, nus if nus else {})
    
    @property
    def V(self) -> Quantity[u.L]:
        return sum([k.V for k in self.gases.values() if k])
    def _repr_latex_(self) -> str:
        s = "Pipe network with members {" + ', '.join([
            k + " $(V=" + v.V.to(u.L).round(2)._repr_latex_().strip('$') + "$ ) " 
            for k,v in self.gases.items()
            ]) + "}"
        if len(self.gases) > 0:
            s += ", total volume " + self.V.to(u.L).round(2)._repr_latex_()
            k0 = [*self.gases.keys()][0]
            s += ";\n- first member props: " + self.gases[k0]._repr_latex_()[12:]
        return s
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
        if gas.C == 0 * u.J/u.K: return gas # nothing to do here
        t0 = gas.T * 1
        for _ in range(self.n):
            eta = min(AIR_VMOLAR/(gas.V/gas.Nu), 1)
            Teff4 = (gas.T**4 - T_COSMIC_MICROWAVE_BACKGROUND**4)
            dQ = PSolEff - u.s* self.area * STEFAN_BOLTZMANN_CONSTANT * Teff4
            dT = dQ * eta / gas.C
            gas.T = gas.T + dT
            # SmLogger.info(f'{(eta*dQ).to(u.J)}, {gas.T}')
        SmLogger.info(f'last {(eta*dQ).to(u.J)}; Per {self.n} pipes {t0:.2f} -> {gas.T:.2f} ({t0-gas.T:.4g} loss)')
        return gas

    def approx_cool(self, gas: GasMixture) -> float:
        """Approximates without modifying"""
        PSolEff = self.qSol * self.area * self.f_insolated * self.n * u.s
        Td = gas.T - T_COSMIC_MICROWAVE_BACKGROUND
        if gas.C == 0 * u.J/u.K: return 0 # nothing to do here
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
    outV = sink.V + (sink._V0 if sink._V0 is not None else 0*u.L) # full volume
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
    transferNu = min(delta_nu, source.Nu) if delta_nu is not None else source.Nu
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
    PumpLogger.info(f'{debug_name}: Sink stat = {sink.T:.5g}, {sink.Nu:.6g}, {sink.p.to(u.kPa):.6g}')
    return specific_power*transferNu,transferNu


@dataclass
class Pump:
    name: str
    V_in : Quantity[u.L]= 200*u.L
    V_out: Quantity[u.L] = 200*u.L
    P_max:    Quantity[u.W] = 30*u.kW # .to(u.W)
    target_p: Quantity[u.Pa] = 15*u.MPa #.to(u.Pa)
    inlet : PipeNetwork|None = None
    outlet: PipeNetwork|None = None

    def reconnect_in(self, inlet: PipeNetwork):
        inlet.add(self.name, self.V_in)
        self.inlet = inlet
    def connect(self, inlet: PipeNetwork, outlet: PipeNetwork):
        inlet.add(self.name, self.V_in)
        self.inlet = inlet
        inlet.equalize()
        outlet.add(self.name, self.V_out)
        self.outlet = outlet
        outlet.equalize()

    def pump(self): #, inlet: GasMixture, outlet: GasMixture):
        if self.inlet is None: return 0
        if self.outlet is None: return 0
        inlet,outlet = self.inlet.gases[self.name], self.outlet.gases[self.name]
        if inlet is None: inlet = self.inlet.gases['room']
        if outlet is None: outlet = self.outlet.gases['room']
        delta_p = self.target_p - outlet.p
        if delta_p <= 10*u.Pa or inlet.T <= 0*u.K:
            PumpLogger.info(f'{self.name}: dP {delta_p} skip')
            return 0*u.W
        # Gets a own-volume segment -- or everything for env.
        source = inlet # * (self.V_in/inlet.V)  if self.V_in !=0 else inlet *1
        sink = outlet  # * (self.V_out/outlet.V) if self.V_out!=0 else outlet*1
        # PumpLogger.info(source)
        delta_nu = calc_transfer_moles(source, sink, delta_p)
        w_draw, moles = pump_gas(source, sink, delta_nu, self.P_max, self.name)
        if moles > 0:
            qty = inlet.sub(moles)
            outlet += qty
            PumpLogger.info(f'{self.name} Merged-sink stat: {outlet.T:.7g}, {outlet.Nu:.4g}, {outlet.p.to(u.kPa):.6g}')
        return w_draw


@dataclass
class OutletInjector:
    name: str
    V_in : Quantity[u.L]= 700*u.L
    P_max:    Quantity[u.W] = 45*u.kW # .to(u.W)
    target_V: Quantity[u.L] = 50*u.L # actually L/2s
    inlet : PipeNetwork|None = None
    room: GasMixture|None = None

    def connect(self, inlet: PipeNetwork, room: GasMixture):
        inlet.add(self.name, self.V_in)
        self.inlet = inlet
        inlet.equalize()
        self.room = room
        
    def pump(self): #, inlet: GasMixture, outlet: GasMixture):
        if self.inlet is None: return 0
        if self.room is None: return 0
        inlet,outlet = self.inlet.gases[self.name], self.room
        assert inlet is not None #: inlet = self.inlet.gases['room']
        
        source = inlet # * (self.V_in/inlet.V)  if self.V_in !=0 else inlet *1
        sink = outlet  # * (self.V_out/outlet.V) if self.V_out!=0 else outlet*1
        # PumpLogger.info(source)
        delta_nu = (self.target_V / source.V) * source.Nu
        w_draw, moles = pump_gas(source, sink, delta_nu, self.P_max, self.name)
        if moles > 0:
            qty = inlet.sub(moles)
            outlet += qty
            PumpLogger.info(f'{self.name} Merged-sink stat: {outlet.T:.7g}, {outlet.Nu / (outlet.V/CELL_VOLUME):.4g}, {outlet.p.to(u.kPa):.6g}')
        return w_draw



@dataclass
class VentPump:
    name: str
    V_out : Quantity[u.L]= 200*u.L
    P_max:    Quantity[u.W] = 30*u.kW # .to(u.W)
    ex_P : Quantity[u.Pa] = 1*atm
    outlet : PipeNetwork|None = None
    room: GasMixture|None = None

    def connect(self,  room: GasMixture, outlet: PipeNetwork):
        outlet.add(self.name, self.V_out)
        self.outlet = outlet
        outlet.equalize()
        self.room = room
        
    def pump(self): #, inlet: GasMixture, outlet: GasMixture):
        if self.outlet is None: return 0
        if self.room is None: return 0
        inlet,outlet = self.room, self.outlet.gases[self.name]
        assert outlet is not None #: inlet = self.inlet.gases['room']
        
        source = inlet # * (self.V_in/inlet.V)  if self.V_in !=0 else inlet *1
        sink = outlet  # * (self.V_out/outlet.V) if self.V_out!=0 else outlet*1
        # PumpLogger.info(source)
        delta_p = 10*u.MPa # default pressure delta val
        delta_p = min(delta_p, self.room.p - self.ex_P)
        PumpLogger.info(f'{self.name} target dP={delta_p.to(u.kPa)}')
        if delta_p < 0.5*u.kPa:
            return 0
        delta_nu = calc_transfer_moles(self.room, sink, delta_p) / (self.room.V / CELL_VOLUME)
        w_draw, moles = pump_gas(source, sink, delta_nu, self.P_max, self.name)
        if moles > 0:
            qty = inlet.sub(moles)
            outlet += qty
            PumpLogger.info(f'{self.name} Merged-sink stat: {outlet.T:.7g}, {outlet.Nu:.4g}, {outlet.p.to(u.kPa):.6g}')
        return w_draw


@dataclass
class Turbine:
    name: str
    V_in = 400*u.L
    V_out = 200*u.L
    inlet:  PipeNetwork|None = None
    outlet: PipeNetwork|None = None
    def connect(self, inlet: PipeNetwork, outlet: PipeNetwork):
        inlet.add(self.name, self.V_in)
        self.inlet = inlet
        outlet.add(self.name, self.V_out)
        self.outlet = outlet
    def take_air(self) -> tuple[GasMixture|None, Quantity[u.J]]:
        src = self.inlet.gases[self.name]
        sink = self.outlet.gases[self.name]
        dP = max(0*u.Pa, src.p-sink.p-5*u.kPa)
        if dP < 5*u.kPa:
            return None, 0*u.J
        dnumax = dP * self.inlet.V /3/src.T/R_IDEAL_GAS
        E_i = min(dP*self.inlet.V,    src.p*self.V_in)
        etaK,fv,kappa = 0.04, 0.2, 0.66
        etaV = etaK / kappa * (1 - fv**kappa)
        E = (etaV*E_i).to(u.kJ)
        vol_cap_frac = min (dP * self.inlet.V/3/src.p/self.V_in, 1)
        nu_can_use_max = src.Nu # self.V_in/self.inlet.V*
        nu_use = min(dnumax, nu_can_use_max)
        SmLogger.info(f'{self.name}: inP = {src.p.to(u.kPa):.4g}, outP={sink.p.to(u.kPa):.4g}')
        gret = src.sub(nu_use)
        SmLogger.info(f'{self.name}: inP[red] = {src.p.to(u.kPa):.4g}')
        SmLogger.info(f'{self.name}: dP={dP:.4g}, {dnumax=}, {E=}')
        SmLogger.info(f'{self.name}: {vol_cap_frac=:.4g}, frac nu{gret.Nu/nu_can_use_max:.4g}')
        return gret, E
    def push_air(self, new: GasMixture) -> GasMixture:
        self.outlet.gases[self.name] += new
        SmLogger.info(f'{self.name}: merged: inP={self.inlet.gases[self.name].p:.5g}, outP={self.outlet.gases[self.name].p:.5g}')

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

def calc_epr(q: GasMixture):
    if q.V == 0 : return 0*atm
    # Vm = q.V/q.Nu
    return CELL_VOLUME**2 / q.V**2 * q.Nu / (23.1 * u.mol)

@dataclass
class Supermatter:
    room: GasMixture
    eer: Quantity[EER] = 0*EER # `power`
    health = 1000 * hp # starting
    # --- stats
    was_overheated = False
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
        self.room.nus['Ex*'] = self.room.nus.get('Ex*', 0*u.mol) + nu_Ex
        self.room.nus['O2'] = self.room.nus.get('O2', 0*u.mol) + nu_O2
        self.room.thermal_E += SM_THERM_K * EP
        if self.room.T > 1e4*u.K:
            # print('OVERHEATED')
            if not self.was_overheated:
                SmLogger.warn('OVERHEATED')
            self.room.T = 1e4*u.K
        else:
            self.was_overheated = False
    def decay(self):
        self.eer -= (self.eer / SM_LAMDA)**3
    def _repr_latex_(self) -> str:
        comps = ', '.join([
            'temp ' + self.room.T.to(u.K).round(2)._repr_latex_(),
            'pres ' + self.room.p.to(u.kPa).round(2)._repr_latex_(),
            self.eer.to(EER).round(2)._repr_latex_(),
            'EPR ' + calc_epr(self.room).round(2)._repr_latex_(),
            ', '.join([f'{k}: {100*v:.2f}%' for k,v in self.room.X.items()])
        ])
        return "Supermatter Core: " + comps