"""Supermatter core itself"""

from dataclasses import dataclass
from .gases import GasMixture
from .consts import *

SM_EER_EMIT = 10 #* hp / EER

SM_T_CRIT = 5000 #* u.K
SM_LAMDA = 700   #* (EER**(2/3))
SM_THERM_K = 1e4 #* u.J / EER
SM_OXY_K = 15e3  #* EER / u.mol
SM_USB_K = 1.5e3 #* EER / u.mol
SM_MAX_DMG = 4.5 #* hp
SM_ETA_GAS = 0.25
SM_ETA_N2 = 0.15
SM_ETA_RXN = 1.1
SM_DMG_TEMP = 150   #* u.K / hp
SM_DMG_EER_VAC = 10 #* EER / hp
# SM_REF_HEALTH = 1000 * hp
# SM_MAX_HEALTH = 1000 * hp
SM_REF_EER = 300  # * EER
SM_REF_TEMP = 800 # * u.K
SM_EER_OXY = 450  # * EER
SM_EER_N2 = 200   # * EER

def calc_epr(q: GasMixture):
    if q.V == 0 : return 0 #*atm
    # Vm = q.V/q.Nu
    return CELL_VOLUME**2 / q.V**2 * q.Nu / (23.1) # * u.mol)

@dataclass
class Supermatter:
    room: GasMixture
    eer: float = 0 # Quantity[EER] = 0*EER # `power`
    health = 1000 #* hp # starting
    # --- stats
    def process(self):
        if not self.damage():
            self.react()
            self.decay()
    def damage(self):
        if self.room.p == 0:
            self.health -= max(self.eer-15, 0)/ SM_DMG_EER_VAC # *EER
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
        self.health = min(self.health, 1000) #*hp)
        return self.health < 0
    def react(self):
        f_O = self.room.nus.get('O2', 0) - SM_ETA_N2 * self.room.nus.get('N2',0) #*u.mol
        f_O /= self.room.Nu
        f_O = max(min(f_O, 1), 0)
        P_E = (400 if f_O>0.8 else 250)#*EER
        dP_rxn = (P_E/SM_LAMDA)**3 * (self.room.T / SM_REF_TEMP) * f_O
        self.eer += dP_rxn
        EP = SM_ETA_RXN*self.eer
        ET = (self.room.T - 273.15) #*u.K) * EER / u.K
        nu_Ex = max(EP/SM_USB_K, 0)#*u.mol)
        nu_O2 = max((EP+ET)/SM_OXY_K, 0)#*u.mol)
        self.room.nus['Ex*'] = self.room.nus.get('Ex*', 0) + nu_Ex # *u.mol
        self.room.nus['O2'] = self.room.nus.get('O2',   0) + nu_O2 # *u.mol
        self.room.thermal_E += SM_THERM_K * EP
        if self.room.T > 1e4:# *u.K:
            print('OVERHEATED')
            self.room.T = 1e4 #*u.K
    def decay(self):
        self.eer -= (self.eer / SM_LAMDA)**3
    def _repr_latex_(self) -> str:
        comps = ', '.join([
            f'temp {self.room.T:.2f} K',
            f'pres {self.room.p/1e3:.2f} kPa',
            f'EER {self.eer:.2f}', #.to(EER).round(2)._repr_latex_(),
            f'EPR {calc_epr(self.room):.2f}',#round(2)._repr_latex_(),
            f', '.join([f'{k}: {100*v:.2f}%' for k,v in self.room.X.items()])
        ])
        return "Supermatter Core: " + comps