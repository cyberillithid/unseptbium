from .ss13_usb import *
from typing import NewType

Coords2 =     NewType('Coords2', tuple[int,int])
PrefabAttrs = NewType('PrefabAttrs', dict[str, ByondConst])
Prefab =      NewType('Prefab', tuple[str, PrefabAttrs])