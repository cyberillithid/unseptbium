"""BYOND DMM analysis and graph parsing utilities"""

from typing import NamedTuple
from enum import Enum, IntEnum
from typing_extensions import Self

class Direction(IntEnum):
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
        return [
            Direction.NORTH, Direction.SOUTH, 
            Direction.EAST, Direction.WEST
        ]
    def __invert__(self) -> Self:
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
        return set([
            Direction.from_dxdy(self.dx, 0),
            Direction.from_dxdy(0, self.dy),
        ])
    def move1(self, xy: tuple[int, int]) -> tuple[int,int]:
        x,y = xy
        return (x+self.dx, y+self.dy)
    def matches(self, dir0: Self|None):
        if dir0 is None: return True
        if int(self) & ~dir0: return True
        if not dir0.is_composite and self & int(dir0): return True
        return False
    @classmethod
    def from_dxdy(cls, dx: int, dy: int) -> Self:
        ret = 0
        if dx==1: ret |= Direction.EAST
        if dx==-1: ret |= Direction.WEST
        if dy==1: ret |= Direction.NORTH
        if dy==-1: ret |= Direction.SOUTH
        return Direction(ret)