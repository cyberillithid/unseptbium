"""Building blocks of detailed unitless pure Python math model."""
from datetime import timedelta

import logging
from .gases import GasSpecies, GasMixture, KnownGas
from .consts import *
from .supermatter import Supermatter, SM_EER_EMIT
from .pipes import PipeNetwork, HeatExchanger, mk_canister
from .pumps import Pump, OutletInjector, VentPump, Filter
from .teg import TEG, Turbine


SmLogger   = logging.getLogger('ct.sm')
MACHINERY_TICK = 2 #* u.s # seconds
PIPENET_TICK = 1   #* u.s # seconds

EMITTER_BURST = 4 # or 3 tho?
POWERED_EMITTER_DMG = 142.22 #* hp

EMITTER_POWER_CONS = 100e3 #* u.W # W
EMITTER_FULL_TIME = 6.4    #* u.s # sec

class SmModel:
    he_i = None
    teg_pwr = 0#*u.J

    core_room: GasMixture
    smcore: Supermatter

    cooled : PipeNetwork
    cooling : PipeNetwork
    warm : PipeNetwork
    scalding: PipeNetwork

    filters: list[Filter]

    def __init__(self):
        self.core_room = GasMixture(0, 15*CELL_VOLUME, {})
        self.core_room._gr = 15
        self.smcore =  Supermatter(self.core_room) #
        # pipe networks with pipelines as base
        self.cooled   = PipeNetwork({'cooled':   GasMixture(0,   .070)})
        self.cooling  = PipeNetwork({'cooling':  GasMixture(0, 15.680)})
        self.warm     = PipeNetwork({'warm':     GasMixture(0,  1.890)})
        self.warm.add('nopump',   .2) # empty for waste or whatever
        self.scalding = PipeNetwork({'scalding': GasMixture(0, .315)})
        self.waste = PipeNetwork({'waste': GasMixture(0, .385)})
        # connectors
        self.c_hot =  PipeNetwork({'': mk_canister('N2',engine=True)})
        self.c_cold = PipeNetwork({'': mk_canister('N2',engine=True)})
        self.filters = [Filter(f'F{i}') for i in range(1,3)]
        # pumps
        self.po1 = OutletInjector('PO1', P_max=30e3, target_V=0.700) # Outlet Injector
        self.pi1 = VentPump('PI1', V_out=0.700, ex_P = 100e3) # , V_in=0*u.L, V_out=0.2*u.m**3) # Vent Pump
        self.hp1 = Pump('HP1', P_max=45e3)
        self.pc_hot =  Pump('P[C]hot') # C2-P2 on pic
        self.pc_cold = Pump('P[C]cold') # C5-P6 on pic
        # TEG and HE
        self.t_hot =  Turbine('T(hot)1')
        self.t_cold = Turbine('T(cold)2')
        self.teg = TEG(self.t_hot, self.t_cold)
        self.he = HeatExchanger(191)
        # -- connecting devices --
        self.pc_hot.connect( self.c_hot,     self.warm)
        self.pc_cold.connect(self.c_cold,    self.cooling)
        self.hp1.connect(    self.cooling,   self.cooled)
        self.pi1.connect(    self.core_room, self.scalding)
        self.po1.connect(    self.warm,      self.core_room)
        self.t_hot.connect(  self.scalding,  self.warm)
        self.t_cold.connect( self.cooled,    self.cooling)
        self.filters[0].connect(self.warm, self.warm, {'O2': self.waste})
        self.filters[1].connect(self.warm, self.warm, {'O2': self.waste})

    def energize_core(self, n_pulses: int):
        self.smcore.eer = POWERED_EMITTER_DMG*n_pulses/SM_EER_EMIT

    def pump_cold(self):
        if self.he_i is None: self.he_i = 0
        if self.he_i > 1: 
            self.he.cool(self.cooling.gases['cooling'])
        [x.equalize() for x in [self.c_cold, self.cooling]]
        self.pc_cold.pump()
        self.he_i += 1

    def pump_hot_teg(self):
        pwr_lost = [x.pump() for x in [self.pc_hot, self.pi1, self.po1]]
        self.teg_pwr = self.teg.process()
        [x.equalize() for x in [self.c_hot, self.warm, self.scalding]]

    def full_tick(self):
        self.pump_cold()
        self.hp1.pump()
        self.pump_hot_teg()
        self.smcore.process()
        [v.process() for v in self.filters]
        [x.equalize() for x in [self.c_hot, self.warm, self.scalding, 
                                self.c_cold, self.cooling, self.cooled,
                                self.waste]]
        self.core_room.react()

    @property
    def reftime(self):
        if self.he_i is None: return timedelta(seconds=0)
        dt = timedelta(seconds=self.he_i*2)
        return (dt)

    def _repr_latex_(self) -> str:
        return '\n\n'.join([
            f'Ticks spent: {self.he_i} ; server time: {self.reftime}; ',
            f'Last power: {self.teg_pwr/1e3:.2f} kW',
            self.smcore._repr_latex_(),
            # 'C[hot]: '  + self.c_hot._repr_latex_()[25:],
            # 'C[cold]: ' + self.c_cold._repr_latex_()[25:],
            'scalding: '+ self.scalding._repr_latex_()[25:],
            'warm: ' +    self.warm._repr_latex_()[25:],
            'cooling: ' + self.cooling._repr_latex_()[25:],
            'cooled: ' +  self.cooled._repr_latex_()[25:],
        ])
