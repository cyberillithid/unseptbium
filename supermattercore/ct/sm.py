import cantera as ct
import numpy as np
import scipy.constants as C
from typing import NamedTuple, Optional
from dataclasses import dataclass
import logging


MACHINERY_TICK = 2 # seconds

POWERED_EMITTER_DMG = 142.22
EMITTER_POWER_CONS = 100e3 # W
EMITTER_FULL_TIME = 6.4 # sec
EMITTER_BURST = 4 # or 3 tho?

def gas_qt_ixr(q: ct.Quantity, mul: float=1) -> dict[str, float]:
    return dict(zip(q.species_names, q.X * mul))

def gas_qt_epr(q: ct.Quantity) -> float:
    if q.moles == 0: return float('nan')
    return q.moles / q.volume**2 * 5**4 / 2.31

class SupermatterCore:
    gas_eff = 0.25 # 1
    decay = 700   # Pa^[2/3]??? Used in radiation loss 
    reaction = 1.1 # `reaction_power_modified`
    k_thermal = 10e3 # J / (MeV/cm3)
    k_exotic = 15e2  
    k_o2 = 15e3
    k_n2_ret = 0.15 # 1; N2 retardation coef.
    k_damage = 4.5e-3 # scaled by 1000
    T_crit = 5000 # K

    damage: float # scaled by 1000 to 0..1 interval; integrity = 1-damage
    power: float # EER, MeV / cm3 -- meaningless!

    def __init__(self, damage = 0, power = 0):
        self.damage = damage
        self.power = power

    def powerup(self, dmg: float):
        self.power += dmg/10.
    
    def tick_damage(self, N: float, T: float) -> bool:
        """N -- moles; T -- temperature"""
        if N == 0:
            self.damage += max((self.power - 15)/10e3, 0)
            return True
        max_dmg = SupermatterCore.k_damage * self.power/300
        techn_dmg = (T - SupermatterCore.T_crit) / 150e3
        self.damage = max(0, self.damage + np.clip(techn_dmg, -SupermatterCore.k_damage, max_dmg))
        return False

    def react(self, o2: float, n2: float, T: float) -> float:
        o2_eff = max(o2 - SupermatterCore.k_n2_ret * n2, 0)
        P_eq = 250 if o2_eff <= 0.8 else 400 # Target equilibrium power?
        fT = T/800 * (P_eq / SupermatterCore.decay)**3
        return o2_eff * fT + self.power

    def tick(self, air: ct.Quantity) -> ct.Quantity:
        if self.damage >= 1:
            air.moles = 0 # exploded
            air.phase.TP = 2.7, 1 # vacuum
            return air
        if self.tick_damage(air.moles, air.T):
            return air
        X = gas_qt_ixr(air)
        P = self.react(X['O2'], X['N2'], air.T)
        
        E = P * SupermatterCore.reaction # ????
        frac_taken = 2.5 * SupermatterCore.gas_eff / air.volume
        new_N = air.moles * frac_taken
        new_X = gas_qt_ixr(air, new_N) # removing a quarter...
        new_Ph = E / SupermatterCore.k_exotic
        new_O2 = max((E + air.T - 273.15) / SupermatterCore.k_o2, 0)
        new_X['Ph'] += new_Ph
        new_X['O2'] += new_O2
        new_N += new_Ph + new_O2

        dQ = E * SupermatterCore.k_thermal

        reag = ct.Solution(air.phase.source)
        reag.TPX = air.T, air.P, new_X
        reag.UV = reag.u + dQ, air.phase.volume_mass
        if reag.T > 10e3: reag.TP = 10e3, reag.P

        reag_q = ct.Quantity(reag, moles=new_N)
        reag_q.constant = 'UV'
        air *= (1 - frac_taken)
        air.constant = 'UV'
        air += reag_q
        
        self.power = P - (P/SupermatterCore.decay)**3
        return air        
    

PUMP_DEFAULT_VOLUME = 0.2 # m^3 == 200 L
SIMPLE_PIPE_VOLUME  = 0.07 # m^3 = 70 L

PIPE_HE_SURFACE = 2 # m^2
SPACE_THC = 1

CANISTER_VOLUME = 1 # m^3 = 1000 L
ENG_CANISTER_PRESSURE = 90*101325 # Pa

AIR_PCRIT = 3771e3 # Pa
AIR_TCRIT = 132.65 # K, 

R_IDEAL_GAS = 8.31
AIR_VMOLAR = R_IDEAL_GAS * AIR_TCRIT / AIR_PCRIT

T_COSMIC_MICROWAVE_BACKGROUND = 3.15
T_MINIMAL = 2.7

class HeatExchanger:
    """ cf.    """
    space_no: int
    plate_no: int
    T_plates: list[float]

    def __init__(self, space: int, plate: int, T_plate: float = 293.15):
        self.space_no = space
        self.plate_no = plate
        self.T_plates = [T_plate]*plate

    def cool_iter(self, q: ct.Quantity):
        """Practically always swappable to cool-lin (even for 3K inital discrepancy <30K, <7 kPa)"""
        s = q.phase # ct.Solution(q.phase.source)
        u1 = s.u
        for j in range(self.space_no):
            lamda = min( AIR_VMOLAR * 1e3 / s.volume_mole , 1 ) 
            dQ1 = 200 * 0.08 - 2 * C.sigma * (s.T - T_COSMIC_MICROWAVE_BACKGROUND)**4
            dQ = (dQ1 * lamda) / q.mass
            s.UV = s.u + dQ, s.v
        u2 = s.u
        for j in range(self.plate_no):
            dT = s.T - self.T_plates[j]
            cp_segm = (0.07 / s.volume_mole * s.cv_mole)
            Q2 = 0.4 * dT * (cp_segm * 10e3 / (cp_segm + 10e3))
            self.T_plates[j] += Q2 / 10e3
            s.UV = s.u - (Q2 / q.mass), s.v
        u3 = s.u
        return ct.Quantity(s, mass = q.mass), (u2-u1), (u3-u2)
    
    def cool_lin(self, q: ct.Quantity):
        s = q.phase
        t_plate = np.sum(self.T_plates) / self.plate_no
        lamda = min( AIR_VMOLAR * 1e3 / q.volume_mole , 1 ) 
        dQ1 = 200 * 0.08 - 2 * C.sigma * (q.T - T_COSMIC_MICROWAVE_BACKGROUND)**4
        dQ = (dQ1 * lamda) / q.mass * self.space_no
        dT = s.T - t_plate
        cp_segm = (0.07 / s.volume_mole * s.cv_mole)
        Q2 = 0.4 * dT * (cp_segm * 10e3 / (cp_segm + 10e3))
        t_plate += Q2 / 10e3
        self.T_plates = [t_plate] * self.plate_no
        dQ2 = - (Q2 / q.mass) * self.plate_no
        s.UV = s.u + dQ + dQ2, s.v
        return ct.Quantity(s, mass = q.mass), dQ, dQ2
    
    def cool_analyt(self, q: ct.Quantity):
        s : ct.Solution = q.phase
        assert self.plate_no == 0
        Td = s.T - T_COSMIC_MICROWAVE_BACKGROUND
        Cve = s.cv_mole * max(q.moles, q.volume / AIR_VMOLAR)
        Tcool = (3*self.space_no*C.sigma*2/Cve + Td**-3)**(-1/3.)
        dTrad = (self.space_no*0.04*2*200)/Cve
        dT = Td - Tcool - dTrad
        dQ = (dT * s.cv_mole * q.moles) / q.mass
        s.UV = s.u - dQ, s.v
        return ct.Quantity(s, mass=q.mass), dT, dQ


IDEAL_GAS_ENTROPY_CONSTANT = 1164    # (mol^3 * s^3) / (kg^3 * L)
@dataclass
class PipeNet: # (NamedTuple):
    """Simulates pipe network; qty == None if there is no gas in it.
    Also simulates segments of it, if V != qty.volume."""
    V: float
    qty: Optional[ct.Quantity]
    _V0: Optional[float] = None
    @property
    def V0(self) -> float:
        return self._V0 if self._V0 else self.V 
    @property
    def T(self) -> float:
        return 0 if self.qty is None else self.qty.phase.T
    @property
    def Nu(self) -> float:
        """kmole"""
        return 0 if self.qty is None else self.qty.moles * (self.V / self.qty.volume)
    @property
    def C(self) -> float:
        if self.qty is None: return 0
        return self.Nu * self.qty.phase.cv_mole
    @property
    def P(self) -> float:
        return 0 if self.qty is None else self.qty.phase.P
    @property
    def S_13(self) -> float:
        """Uses the approximate entropy from xgm-gas-mixture"""
        if self.qty is None: return 150e3
        ret = 0
        T = max(self.T, T_MINIMAL)
        nu = self.Nu * 1e3
        for (i, x) in enumerate(self.qty.phase.X):
            sp = self.qty.phase.species(i)
            scv = sp.thermo.cp(self.T)/1000 - 8.314
            mu = (self.qty.phase.atomic_weights @ np.array([sp.composition.get(x, 0) for x in self.qty.phase.element_names])) * 1e-3
            snu = nu * x
            if snu > 0:
                ss = R_IDEAL_GAS * (np.log( 1 +
                    (mu * scv * T)**(2/3) * IDEAL_GAS_ENTROPY_CONSTANT * self.V * 1e3 / T / snu
                ) + 15)
                ret += snu * ss
        return ret / nu
    def sub(self, moles: float) -> Optional[ct.Quantity]:
        if self.qty is None:
            return None
        q = ct.Quantity(self.qty.phase, moles=moles)
        self.qty.moles -= moles # -- also changes qty.V, so we need to...
        self.V -= q.V
        self.grow(q.V)
        return q

    def segment(self, V: float):
        return PipeNet(V=V, qty=self.qty, _V0=self.V)
    def grow(self, dV: float):
        if self.qty is None:
            self.V += dV
            return self
        oldNu = self.qty.moles
        oldV = self.qty.V
        s = self.qty.phase
        t,p = s.TP
        p /= (oldV + dV)/oldV
        s.TP = t,p
        self.qty = ct.Quantity(s, moles=oldNu)
        self.V += dV
        return self
        
    def __iadd__(self, q: Optional[ct.Quantity]):
        if self.qty is None:
            self.qty = q
            propV = self.V
            self.V = q.V
            self.grow(propV-q.V)
            return self
        if q is None:
            return self
        newV = q.volume
        self.qty.constant = q.constant = 'UV'
        self.qty += q
        return self.grow(-newV)

MIN_KMOLES_PUMP = 1e-5 # 0.01 moles
ATMOS_PUMP_EFFICIENCY = 2.5 # KPD 40%

PumpLogger = logging.getLogger('ct.sm.Pump')

def calc_transfer_moles(source: PipeNet, sink: PipeNet, deltaP: float) -> float:
    """Imitates calc_transfer_moles of original code. 
    Gets segments of whole network as input."""
    if source.T == 0 or source.Nu == 0: 
        return 0 # Nothing to pump
    outV = sink.V0 # full volume
    src_Nu = source.Nu # just the owned source -- `source total moles`
    airT = source.T
    if sink.T > 0 and sink.Nu > 0:
        # approximate mixture temp to re-calc
        eNu = deltaP*outV/sink.T/ct.constants.gas_constant # kMoles
        sinkC = sink.C # total heat cap
        srcC = source.C * eNu / src_Nu
        airT = (sink.T * sinkC + source.C*srcC) / (sinkC + srcC)
    return deltaP*outV/airT/ct.constants.gas_constant # kMoles

def pump_gas(inlet: PipeNet, outlet: PipeNet, transferNu: float = 0, availablePower: float = 0) -> float|None:
    """Imitates pump_gas of original code; modifies nets, returns power draw in W"""
    if inlet.Nu < MIN_KMOLES_PUMP: 
        return None # Nothing to move

    transferNu = min(transferNu, inlet.Nu) if transferNu else inlet.Nu
    # V Calculate specific power
    usedT = outlet.T if outlet.T > 0 else inlet.T
    srcS = inlet.S_13
    sinkS = outlet.S_13
    specific_entropy = sinkS - srcS
    specific_power = -specific_entropy*usedT if specific_entropy<0 else 0
    specific_power /= ATMOS_PUMP_EFFICIENCY # magical coeff.
    PumpLogger.info(f'Source entropy {srcS:.5g} J/mol/K -> Sink entropy {sinkS:.5g} J/mol/K')
    PumpLogger.info(f'Specific entropy change {specific_entropy:.5g} J/mol/K')
    PumpLogger.info(f'Specific power {specific_power:.6g} W/mol')
    transferNu = min(transferNu, availablePower/specific_power/1e3) if specific_power else transferNu
    PumpLogger.info(f'Transferred {transferNu*1e3:.6g} mols')
    return specific_power*transferNu,transferNu

@dataclass
class Pump:
    power_rating = 30e3 # W
    target_P = 15e6 # Pa
    own_V_in = PUMP_DEFAULT_VOLUME # 0.2 m^3
    own_V_out = PUMP_DEFAULT_VOLUME # 0.2 m^3
    idle_power = 150 # W?

    def process(self, inlet: PipeNet, outlet: PipeNet) -> float:
        """Modifies inlet and outlet, returns power draw"""
        deltaP = self.target_P - outlet.P 
        if (deltaP <= 10) or (inlet.T <= 0):
            return 0 # or `inlet.mass > 0`?
        Vin,Vout = self.own_V_in, self.own_V_out
        own_inlet  =  inlet.segment(Vin)  if Vin  else inlet
        own_outlet = outlet.segment(Vout) if Vout else outlet

        deltaNu = calc_transfer_moles(own_inlet, own_outlet, deltaP)
        power_draw, moles = pump_gas(own_inlet, own_outlet, deltaNu, self.power_rating)
        if moles:
            qty = inlet.sub(moles)
            outlet += qty
        return power_draw
    pass

class TurbinePipes(NamedTuple):
    inlet: PipeNet
    outlet: PipeNet
class TEGPipes(NamedTuple):
    hi: TurbinePipes
    lo: TurbinePipes

KAPPA = 0.666 # adiabatic exp - 1 (gamma=1.66 here)

class TurbineCirc:
    kinetic_friction = 5e3 # Pa
    static_friction = 10e3 # Pa
    stored_energy = 0
    eta_kinetic = 0.04 #+kin-to-electric
    volume_ratio = 0.2 # ?
    in_V  = 0.2  # own volumes
    out_V = 0.2
    def extract_air(self, nets: TurbinePipes) -> tuple[TurbinePipes, Optional[ct.Quantity]]:
        # check input vs air1 in DM
        inl, outl = nets
        inP, outP = inl.P, outl.P
        deltaP = inP-outP
        rm = None
        if deltaP>0 and inl.T>0:
            # has things
            deltaNu = (deltaP*inl.V/R_IDEAL_GAS/inl.T)/3
            capUsed = min(1, (deltaP*inl.V/3)/(inP*self.in_V))
            nRT = min(deltaP*inl.V, inP*self.in_V)
            self.stored_energy += nRT/KAPPA * (1-self.volume_ratio**KAPPA) * self.eta_kinetic
            rm = inl.sub(deltaNu)
        return nets, rm
    def getE(self):
        se = self.stored_energy
        self.stored_energy = 0
        return se

class ThermoElectricGen:
    idle_power  = 100 # W
    max_power = 500e3 # W
    eta_thermal = 0.65
    circ_hi: TurbineCirc
    """circ1 -- North [out=E] (or West)"""
    circ_lo: TurbineCirc
    """circ2 -- South [out=W] (or East [])"""

    def process(self, state: TEGPipes) -> TEGPipes:
        hipipes, lopipes = state
        hipipes, hi_proc = self.circ_hi.extract_air(hipipes)
        lopipes, lo_proc = self.circ_lo.extract_air(lopipes)
        prevgen = last_gen
        if hi_proc and lo_proc: #Not None
            hiC,loC = hi_proc.C, lo_proc.C
            dT = abs(lo_proc.T-hi_proc.T)
            if dT>0 and hiC>0 and loC>0:
                dE = dT*loC*hiC/(loC+hiC)
                dQ = dE*(1-self.eta_thermal)
                if lo_proc.T>hi_proc.T:
                    lo_proc.T = lo_proc.T - dE/loC
                    hi_proc.T = hi_proc.T + dQ/hiC
                else:
                    lo_proc.T = lo_proc.T + dQ/loC
                    hi_proc.T = hi_proc.T - dE/hiC
            hiin, hiout = hipipes
            loin, loout = lopipes
            hiout += hi_proc
            loout += lo_proc
        if eff_gen>self.max_power:
            stored_P *=0.5
        storedP += dE*self.eta_thermal + self.circ_hi.getE() + self.circ_lo.getE()
        last_gen = storedP*.4
        storedP -= last_gen
        eff_gen = (last_gen+prevgen)/2
        # gen_power(eff_gen)
        return TEGPipes(TurbinePipes(hiin,hiout),TurbinePipes(loin,loout))
