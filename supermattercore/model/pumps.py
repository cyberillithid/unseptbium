"""Pumps of all varieties. TODO: refactor to cleanliness!"""
from dataclasses import dataclass

from .gases import GasMixture, KnownGas
from .pipes import PipeNetwork
from .consts import *

import logging
PumpLogger = logging.getLogger('ct.sm.Pump')

MIN_MOLES_PUMP = 0.01 #* u.mol
MIN_MOLES_FILTER = 0.04

ATMOS_PUMP_EFFICIENCY = 2.5
ATMOS_FILTER_EFFICIENCY = 2.5

PUMP_DEFAULT_VOLUME = 0.2 # * u.m**3


def calc_transfer_moles(source: GasMixture, sink: GasMixture, delta_p: float)->float:# Quantity[u.Pa]) -> Quantity[u.mol]:
    """Estimates moles to transfer, as per orig. code"""
    if source.T == 0 or source.Nu == 0: 
        return 0 # Nothing to pump
    outV = sink.V + (sink._V0 if sink._V0 is not None else 0) #*u.L) # full volume
    src_Nu = source.Nu # just the owned source -- `source total moles`
    airT = source.T
    if sink.T > 0 and sink.Nu > 0:
        # approximate mixture temp to re-calc
        eNu = delta_p*outV/sink.T/R_IDEAL_GAS # should be moles
        sinkC = sink.C # total heat cap
        srcC = source.C * eNu / src_Nu
        airT = (sink.T * sinkC + source.T*srcC) / (sinkC + srcC)
    return (delta_p*outV/airT/R_IDEAL_GAS) #

def pump_gas(
        source: GasMixture, sink: GasMixture, 
        delta_nu:    float,#Quantity[u.mol], 
        power_rating:float,# Quantity[u.W],
        debug_name: str = ''
        ) -> tuple[float,float]: #Quantity[u.W], Quantity[u.mol]]:
    if source.Nu < MIN_MOLES_PUMP: 
        return 0,0#*u.W, 0*u.mol # Nothing to move
    transferNu = min(delta_nu, source.Nu ) if delta_nu is not None else source.Nu
    # V Calculate specific power
    usedT = sink.T if sink.T > 0 else source.T
    srcS = source.S
    sinkS = sink.S
    specific_entropy = sinkS - srcS
    specific_power = -specific_entropy*usedT if specific_entropy<0 else 0#*u.J/u.mol
    specific_power /= ATMOS_PUMP_EFFICIENCY # magical coeff.
    PumpLogger.info(f'{debug_name}: Source entropy {srcS:.6g} -> Sink entropy {sinkS:.6g}')
    PumpLogger.info(f'{debug_name}: Specific entropy change {specific_entropy:.6g}')
    PumpLogger.info(f'{debug_name}: Specific power {specific_power:.6g}')
    transferNu = min(transferNu, power_rating/specific_power) if specific_power>0 else transferNu
    # transferNu = transferNu.to(u.mol)
    PumpLogger.info(f'{debug_name}: Transferred {transferNu:.6g} in fact')
    PumpLogger.info(f'{debug_name}: Sink stat = {sink.T:.5g}, {sink.Nu:.6g}, {sink.p/1e3:.6g}')
    return specific_power*transferNu,transferNu


@dataclass
class Pump:
    name: str
    V_in :   float = PUMP_DEFAULT_VOLUME  # Quantity[u.L]= 200*u.L
    V_out:   float = PUMP_DEFAULT_VOLUME  # Quantity[u.L] = 200*u.L
    P_max:   float = 30e3 # Quantity[u.W] = 30*u.kW # .to(u.W)
    target_p:float = 15e6 # Quantity[u.Pa] = 15*u.MPa #.to(u.Pa)
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
        if delta_p <= 10 or inlet.T <= 0: #uPa uK
            PumpLogger.info(f'{self.name}: dP {delta_p} skip')
            return 0 #*u.W
        # Gets a own-volume segment -- or everything for env.
        source = inlet # * (self.V_in/inlet.V)  if self.V_in !=0 else inlet *1
        sink = outlet  # * (self.V_out/outlet.V) if self.V_out!=0 else outlet*1
        # PumpLogger.info(source)
        delta_nu = calc_transfer_moles(source, sink, delta_p)
        w_draw, moles = pump_gas(source, sink, delta_nu, self.P_max, self.name)
        if moles > 0:
            qty = inlet.sub(moles)
            outlet += qty
            PumpLogger.info(f'{self.name} Merged-sink stat: {outlet.T:.7g}, {outlet.Nu:.4g}, {outlet.p/1e3:.6g}')
        return w_draw


@dataclass
class OutletInjector:
    name: str
    V_in :    float = 0.7 #Quantity[u.L]= 700*u.L
    P_max:    float = 45e3 #Quantity[u.W] = 45*u.kW # .to(u.W)
    target_V: float = 0.05 #Quantity[u.L] = 50*u.L # actually L/2s
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
            PumpLogger.info(f'{self.name} Merged-sink stat: {outlet.T:.7g},'
                            f' {outlet.Nu / (outlet.V/CELL_VOLUME):.4g}, {outlet.p/1e3:.6g}')
        return w_draw



@dataclass
class VentPump:
    name: str
    V_out : float = 0.2 #: Quantity[u.L]= 200*u.L
    P_max : float = 30e3 #:    Quantity[u.W] = 30*u.kW # .to(u.W)
    ex_P  : float = 101325 #: Quantity[u.Pa] = 1*atm
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
        delta_p = 10e6  #*u.MPa # default pressure delta val
        delta_p = min(delta_p, self.room.p - self.ex_P)
        PumpLogger.info(f'{self.name} target dP={delta_p/1e3}')
        if delta_p < 500: #*u.kPa:
            return 0
        delta_nu = calc_transfer_moles(self.room, sink, delta_p) / (self.room.V / CELL_VOLUME)
        w_draw, moles = pump_gas(source, sink, delta_nu, self.P_max, self.name)
        if moles > 0:
            qty = inlet.sub(moles)
            outlet += qty
            PumpLogger.info(f'{self.name} Merged-sink stat: {outlet.T:.7g}, {outlet.Nu:.4g}, {outlet.p/1e3:.6g}')
        return w_draw

def calc_Pspec(g: KnownGas, source: GasMixture, sink: GasMixture) -> float:
    airT = sink.T if sink.T>0 else source.T
    dS = sink.s_gas(g) - source.s_gas(g)
    return -dS*airT if dS<0 else 0

def filter_gas_multi(name: str, filts: dict[KnownGas, GasMixture], source: GasMixture, sink: GasMixture, nu: float, P_max: float)->float:
    if source.Nu < MIN_MOLES_FILTER: return 0
    filtering = {k: v for k,v in filts.items() if k in source.nus}
    filter_nu = skip_nu = 0
    tot_Pspec = 0
    gas_Pspec = {}
    for k,v in source.nus.items():
        if v < MIN_MOLES_FILTER: continue
        if k in filtering:
            sink_f = filtering[k]
            gas_Pspec[k] = calc_Pspec(k, source, sink_f)/ATMOS_FILTER_EFFICIENCY
            filter_nu += v
        else:
            gas_Pspec[k] = calc_Pspec(k, source, sink)/ATMOS_FILTER_EFFICIENCY
            skip_nu += v
        tot_Pspec += gas_Pspec[k] * v / source.Nu
    nu = source.Nu if nu else min(source.Nu, nu)
    if tot_Pspec>0 and P_max:
        nu = min(nu, P_max/tot_Pspec)
    if nu < MIN_MOLES_FILTER: return 0
    rm = source.sub(nu)
    tot_pwrfact = 0
    for k,v in rm.nus.items():
        pwr = v * gas_Pspec.get(k, 0)
        if k in filtering:
            filtering[k] += GasMixture(rm.T, 0, {k: v})
            rm.nus[k] = 0
        tot_pwrfact += pwr*v
    rm.nus = {k: v for k,v in rm.nus.items() if v>0}
    sink += rm
    return pwr

@dataclass
class Filter:
    name: str
    inlet : PipeNetwork|None = None
    outlet : PipeNetwork|None = None
    filters : dict[KnownGas,PipeNetwork]|None = None
    max_out_p: float = 15e6 #check!
    P_max: float = 30e3 #check!
    flow_rate: float = 0.5 # m³/"s"
    def connect(self, inlet: PipeNetwork, outlet: PipeNetwork, filters: dict[KnownGas,PipeNetwork]|None):
        inlet.add(self.name+'_i', 0.5)
        self.inlet = inlet
        inlet.equalize()
        outlet.add(self.name+'_o', 0.5)
        self.outlet = outlet
        outlet.equalize()
        for k,v in filters.items():
            v.add(self.name+f'_{k}', 0.5)
            v.equalize()
        self.filters = filters
    def process(self):
        in_air  =  self.inlet.gases[self.name+'_i']  if self.inlet is not None else None
        out_air = self.outlet.gases[self.name+'_o'] if self.outlet is not None else None
        delta = min(max(0, (0 if out_air is None else self.max_out_p-out_air.p)),self.max_out_p)
        dnu_max = calc_transfer_moles(in_air, out_air, delta)
        filters = {k: v.gases[self.name+f'_{k}'] for k,v in self.filters.items()}
        for fn in filters.values():
            delta = min(max(0, (0 if fn is None else self.max_out_p-fn.p)),self.max_out_p)
            dnu_max = min(dnu_max, calc_transfer_moles(in_air, fn, delta))
        dnu_ex = min(max(0, self.flow_rate/in_air.V*in_air.Nu), dnu_max)
        if dnu_ex < MIN_MOLES_PUMP: return 0
        return filter_gas_multi(self.name, filters, in_air, out_air, dnu_max, self.P_max)
    pass