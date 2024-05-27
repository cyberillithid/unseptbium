"""TEG and turbines."""
from dataclasses import dataclass
import logging

from .consts import *
from .pipes import PipeNetwork
from .gases import GasMixture

TegLogger = logging.getLogger('ct.sm.TEG')

@dataclass
class Turbine:
    name: str
    V_in =  0.4#00*u.L
    V_out = 0.2#00*u.L
    inlet:  PipeNetwork|None = None
    outlet: PipeNetwork|None = None
    def connect(self, inlet: PipeNetwork, outlet: PipeNetwork):
        inlet.add(self.name, self.V_in)
        self.inlet = inlet
        outlet.add(self.name, self.V_out)
        self.outlet = outlet
    def take_air(self) -> tuple[GasMixture|None, float]: #  Quantity[u.J]
        src = self.inlet.gases[self.name]
        sink = self.outlet.gases[self.name]
        dP = max(0, src.p-sink.p-5e3) # Pa,*u.kPa)
        if dP < 5e3: #*u.kPa:
            return None, 0#*u.J
        dnumax = dP * self.inlet.V /3/src.T/R_IDEAL_GAS
        E_i = min(dP*self.inlet.V,    src.p*self.V_in)
        etaK,fv,kappa = 0.04, 0.2, 0.66
        etaV = etaK / kappa * (1 - fv**kappa)
        E = (etaV*E_i)#.to(u.kJ)
        vol_cap_frac = min (dP * self.inlet.V/3/src.p/self.V_in, 1)
        nu_can_use_max = src.Nu # self.V_in/self.inlet.V*
        nu_use = min(dnumax, nu_can_use_max)
        TegLogger.info(f'{self.name}: inP = {src.p/1e3:.4g}, outP={sink.p/1e3:.4g}')
        gret = src.sub(nu_use)
        TegLogger.info(f'{self.name}: inP[red] = {src.p/1e3:.4g}')
        TegLogger.info(f'{self.name}: dP={dP:.4g}, {dnumax=}, {E=}')
        TegLogger.info(f'{self.name}: {vol_cap_frac=:.4g}, frac nu{gret.Nu/nu_can_use_max:.4g}')
        return gret, E
    def push_air(self, new: GasMixture) -> GasMixture:
        self.outlet.gases[self.name] += new
        TegLogger.info(f'{self.name}: merged: inP={self.inlet.gases[self.name].p:.5g}, outP={self.outlet.gases[self.name].p:.5g}')

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
