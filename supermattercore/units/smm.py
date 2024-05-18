"""An attempt to rewrite model to be cantera-less and with proper physical units inlined everywhere"""

import astropy.units as u
from astropy.units import Quantity

# import cantera as ct
import numpy as np
import scipy.constants as C
from typing import NamedTuple, Optional
from dataclasses import dataclass

hp =  u.def_unit('хп')
atm = u.def_unit('atm', 101325*u.Pa)

MACHINERY_TICK = 2 * u.s # seconds

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

IDEAL_GAS_ENTROPY_CONSTANT = 1164 # * (u.mol * u.s / u.kg)**3 / (u.L) # -- looks like shit

