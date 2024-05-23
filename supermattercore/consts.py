"""Common emulation/mathmodel constants

Includes BYOND-compatible-ish 'Direction' as Python enum. """
from enum import IntEnum
from typing_extensions import Self


ATMOS_TICK = 2
"""Atmosphere module ticks once every 2 seconds"""


ATMOS = '/obj/machinery/atmospherics' 
"""Common atmospheric machinery prefix"""

OBJ_LABELS = {
    ATMOS + '/portables_connector':     'C',
    ATMOS + '/omni/filter':             'F',
    '/obj/machinery/generator':         'G',
    ATMOS + '/unary/heat_exchanger':    'HE',
    '/obj/machinery/meter':             'M',
    ATMOS + '/binary/pump':             'P',
    ATMOS + '/unary/vent_pump':         'PI',
    ATMOS + '/unary/outlet_injector':   'PO',
    ATMOS + '/binary/circulator':       'T',
    ATMOS + '/valve':                   'V',
    ATMOS + '/pipe':                    '_',
}
"""Shorthands for relevant objects"""
OBJ_ARITY = {'C': 1, 'F': 4, 'G': 2, 'HE': 2, 'M': 0, 'P': 2, 'PI': 1, 'PO': 1, 'T': 2, 'V': 2}
"""How many connections there are to consider?"""


VACUUM_FLOOR = '/turf/simulated/floor/reinforced/airless'
AIRLOCK = '/obj/machinery/door/airlock/atmos'
SM_CORE = '/obj/machinery/power/supermatter'

PIPE_VOLUMES = {
    ATMOS + '/pipe/cap': 35,
    ATMOS + '/pipe/simple': 70,
    ATMOS + '/pipe/manifold': 70*1.5,
    ATMOS + '/pipe/manifold4w': 70*2,
    ATMOS + '/pipe/zpipe': 0, # DIRT HACK FOR NOW
}

def match_path(s1: str, s2: str) -> bool:
    """Returns true if s1 is subclass of s2 or vice versa (for BYOND obj paths)"""    
    return all([a==b for a,b in zip(s1.split('/'), s2.split('/'))] )

class ConnectionType(IntEnum):
    REGULAR = 1
    SUPPLY = 2
    SCRUBBER = 4
    HE = 8
    FUEL = 16
    @classmethod
    def all(cls) -> int:
        return ConnectionType.REGULAR | ConnectionType.SUPPLY | ConnectionType.SCRUBBER | ConnectionType.HE | ConnectionType.FUEL
    @classmethod
    def from_path(cls, path: str) -> int:
        endpath = path.split('/')[-1]
        endsarg = CONNECTION_TYPES_ENDSWITH.get(endpath)
        if endsarg is not None:
            return endsarg
        maxp = sorted([(len(k), v) for k,v in CONNECTION_TYPES_CLASSES.items() if match_path(k, path)])
        if len(maxp)>0:
            return maxp[-1][1]
        return ConnectionType.REGULAR

# TODO: autopop that from code!
CONNECTION_TYPES_CLASSES = {
    ATMOS + '/portables_connector':     ConnectionType.REGULAR|ConnectionType.FUEL,
    ATMOS + '/pipe/simple/heat_exchanging/junction':  ConnectionType.REGULAR|ConnectionType.HE|ConnectionType.FUEL,
    ATMOS + '/pipe/simple/heat_exchanging':  ConnectionType.HE,
    ATMOS + '/valve':   ConnectionType.SCRUBBER |  ConnectionType.SUPPLY |  ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/binary/circulator': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/binary/passive_gate': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/binary/pump': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/omni': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/unary/heat_exchanger': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/unary/heater': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/unary/outlet_injector': ConnectionType.REGULAR |  ConnectionType.FUEL,  
    ATMOS + '/unary/tank':   ConnectionType.SCRUBBER |  ConnectionType.SUPPLY |  ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/unary/vent_pump':   ConnectionType.SUPPLY |  ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/unary/vent_scrubber':   ConnectionType.SCRUBBER | ConnectionType.REGULAR,
    ATMOS + '/unary/engine': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/unary/fission_core': ConnectionType.REGULAR |  ConnectionType.FUEL,
    ATMOS + '/binary/stirling': ConnectionType.REGULAR |  ConnectionType.FUEL,
}

CONNECTION_TYPES_ENDSWITH = {
    'universal': ConnectionType.all(),
    'fuel': ConnectionType.FUEL,
    'scrubbers': ConnectionType.SCRUBBER,
    'supply': ConnectionType.SUPPLY,
}

class Direction(IntEnum):
    """BYOND directions enum, as well as a number of helper methods for it."""
    NORTH     =  1
    SOUTH     =  2
    EAST      =  4
    WEST      =  8
    NORTHEAST =  5
    NORTHWEST =  9
    SOUTHEAST =  6
    SOUTHWEST = 10
    UP        = 16 # //reserved for future use
    DOWN      = 32
    @classmethod
    def cardinals(cls) -> list[Self]:
        """List of cardinal directions"""
        return [
            Direction.NORTH, Direction.SOUTH, 
            Direction.EAST, Direction.WEST
        ]
    def __invert__(self) -> Self:
        """Returns reverse direction"""
        ret,arg = 0, int(self)
        if arg & Direction.NORTH: ret |= Direction.SOUTH
        if arg & Direction.SOUTH: ret |= Direction.NORTH
        if arg & Direction.EAST: ret |= Direction.WEST
        if arg & Direction.WEST: ret |= Direction.EAST
        return Direction(ret)
    def __xor__(self, other: Self) -> Self|None:
        if int(self)&other: return None
        return Direction(int(self) ^ int(~other))
    @property
    def dy(self) -> int:
        """per StrongDMM Y"""
        arg = int(self)
        if arg & Direction.NORTH: return 1
        if arg & Direction.SOUTH: return -1
        return 0
    @property
    def dx(self) -> int:
        arg = int(self)
        if arg & Direction.EAST: return 1
        if arg & Direction.WEST: return -1
        return 0
    @property
    def is_composite(self) -> bool:
        return (abs(self.dx) ^ abs(self.dy)) == 0
    def split(self) -> set[Self]:
        """Returns pair of cardinal directions for a composite direction.
        Probably runtimes for a cardinal direction."""
        return set([
            Direction.from_dxdy(self.dx, 0),
            Direction.from_dxdy(0, self.dy),
        ])
    def move1(self, xy: tuple[int, int]) -> tuple[int,int]:
        """Moves from coordinates by 1 into this direction"""
        x,y = xy
        return (x+self.dx, y+self.dy)
    def matches(self, dir0: Self|None):
        if dir0 is None: return True
        if int(self) & ~dir0: return True
        if not dir0.is_composite and self & int(dir0): return True
        return False
    @classmethod
    def from_dxdy(cls, dx: int, dy: int) -> Self:
        """Constructs a direction from X and Y shifts."""
        ret = 0
        if dx==1: ret |= Direction.EAST
        if dx==-1: ret |= Direction.WEST
        if dy==1: ret |= Direction.NORTH
        if dy==-1: ret |= Direction.SOUTH
        return Direction(ret)
    @classmethod
    def from_str(cls, s: str) -> Self:
        if s.upper()=='NORTH'    : return Direction.NORTH       
        if s.upper()=='SOUTH'    : return Direction.SOUTH       
        if s.upper()=='EAST'     : return Direction.EAST       
        if s.upper()=='WEST'     : return Direction.WEST       
        if s.upper()=='NORTHEAST': return Direction.NORTHEAST           
        if s.upper()=='NORTHWEST': return Direction.NORTHWEST           
        if s.upper()=='SOUTHEAST': return Direction.SOUTHEAST           
        if s.upper()=='SOUTHWEST': return Direction.SOUTHWEST           
        if s.upper()=='UP'       : return Direction.UP   
        if s.upper()=='DOWN'     : return Direction.DOWN 
        return None      
    def turn(self, deg: int) -> Self:
        octoroons = ((deg // 45) + 8) % 8
        pass