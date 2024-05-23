"""BYOND DMM analysis and graph parsing utilities"""

from typing import NamedTuple, NewType
from typing_extensions import Self

from dataclasses import dataclass

import pprint
import logging

from ss13_usb import Ss13Map, Coords2, Prefab, PrefabAttrs, ByondConst
from .consts import OBJ_LABELS, OBJ_ARITY, ATMOS, PIPE_VOLUMES, Direction, ConnectionType, match_path

GraphLogger = logging.getLogger('ss13.sm.usb.Graph')

def input_tag(a: PrefabAttrs) -> Direction:
    for k,v in a.items():
        if k.startswith('tag_'):
            if isinstance(v, ByondConst.Float) and int(v) == 1: 
                return Direction.from_str(k[4:]) #direction_names.get(k[4:], Direction.UP)

class StoredObj(NamedTuple):
    """A (non-compact) representation of a parsed object"""
    coords: Coords2
    path: str
    attrs:  PrefabAttrs
    symbol: str
    turf: str|None
    area: str|None
    @classmethod
    def list_from_map_xy(cls, map: Ss13Map, coord: Coords2, symdict: dict[str, str]) -> list[Self]:
        prefabs: list[Prefab] = map[coord]
        turfs = [k for (k,_) in prefabs if match_path(k, '/turf')]
        areas = [k for (k,_) in prefabs if match_path(k, '/area')]
        turf = turfs[0] if len(turfs)==1 else None
        area = areas[0] if len(areas)==1 else None
        return [ StoredObj(coord, path, attrs, sym, turf, area)
            for (path, attrs) in prefabs
            for (k, sym) in symdict.items()
            if match_path(path, k)
        ]
    def __eq__(self, o):
        return self.coords==o.coords and self.path==o.path
    @property
    def dir_attr(self) -> Direction:
        attrdir = self.attrs.get('dir')
        if attrdir is not None:
            return Direction(int(attrdir))
        return Direction(2) # looks default?
    @property
    def conn_type(self) -> int:
        return ConnectionType.from_path(self.path)
    @property
    def conn_dirs(self) -> set[Direction]:
        dirval = self.dir_attr
        arity = OBJ_ARITY.get(self.symbol)
        if self.symbol == '_':
            # pipes are magic
            if 'manifold' in self.path:
                assert not dirval.is_composite
                return set(Direction.cardinals()) - {dirval}
            # else it's binary pipe
            if dirval.is_composite:
                return dirval.split()
            arity = 2 # as normal binary
        if arity == 0: return set()
        if arity == 2: return {dirval, ~dirval}
        if arity == 4: return set(Direction.cardinals())
        return {dirval}
    
# find_amid1 = lambda ssm, xs:  [(StoredObj(k, a, b, OBJ_LABELS[v]),d) for v in OBJ_LABELS for (k,d) in xs for (a, b) in ssm[k] if match_path(a,v)]
find_amid1 = lambda ssm, xds: [(o,*d) for (x,*d) in xds for o in StoredObj.list_from_map_xy(ssm, x, OBJ_LABELS)]

find_amid = lambda ssm, xs:  [k for (k,d) in find_amid1(ssm,[(x, None) for x in xs])]
"""Returns all describable objects amid given xs coords"""


def follow_networks(ssm: Ss13Map, start: list[StoredObj], zeronodes: list[StoredObj]) -> dict[Coords2, list[StoredObj]]:
    """Starting with given object, follows networks until all is consumed."""
    next_objs = [(o, None, ConnectionType.all()) for o in start]
    seen_objs = {}
    for o in zeronodes:
        seen_objs[o.coords] = seen_objs.get(o.coords, []) + [o]
    while len(next_objs) > 0:
        prev_objs = next_objs
        # print(len(prev_objs))
        next_objs = []
        for o, dir, ctype in prev_objs:
            if o.coords in seen_objs: 
                if any(o==k for k in seen_objs[o.coords]):
                    continue
                # print(o, seen_objs[o.coords])
            if dir is not None and len(o.conn_dirs)>0 and dir not in o.conn_dirs: 
                continue # wrong pipe
            if (o.conn_type & ctype) == 0:
                continue # wrong sort of pipe
            coords = [(k.move1(o.coords), ~k, o.conn_type) for k in (o.conn_dirs - {dir})]
            next_objs += [z for z in find_amid1(ssm, coords)]
            # [(StoredObj(k, a, b, OBJ_LABELS[v]), d) for v in OBJS for (k,d) in coords for (a, b) in ssm[k] if match_path(a,v)]
            seen_objs[o.coords] = seen_objs.get(o.coords, []) + [o]
    return seen_objs

def split_to_obj_pipes(seen_objs: dict[Coords2, StoredObj], pipesymbol = '_') -> tuple[dict[str,list[StoredObj]], list[StoredObj]]:
    stored_objs = {}
    known_pipes = []
    # Split objs to pipes and devices
    for k,os in seen_objs.items():
        for o in os:
            if o.symbol == pipesymbol: 
                known_pipes.append(o)
            else:
                stored_objs[o.symbol] = stored_objs.get(o.symbol, []) + [o]
    del seen_objs
    return (stored_objs, known_pipes)

shorten_pipecls = lambda s: '/'.join(s.split('/')[5:])
lengthen_pipecls = lambda s: ATMOS+'/pipe/'+s
simplify_pipecls = lambda s: shorten_pipecls(s).replace('simple','/').replace('visible','')

def gather_pipejoins(pipejoins, remaps, pipekeys):
    """Gathering joinable pipe-segments to pipenets"""
    pipeseqs = []
    for a,b in pipejoins:
        a = remaps.get(a,a) # to account for merges
        b = remaps.get(b,b)
        founds = [s for s in pipeseqs if a in s or b in s]
        rests =  [s for s in pipeseqs if not (a in s or b in s)]
        # print(a,b,pipeseqs)
        if len(founds) == 0:
            pipeseqs.append({a, b})
            continue
        pipeseqs = [*rests, set.union(*founds, {a, b})]
    # and add the unmergeable pipe-segments on their own
    joint = set.union(*pipeseqs)
    pipeseqs += [{p} for p in (set(pipekeys)-joint)]    
    return pipeseqs

SegmentId = NewType('SegmentId', int)
class PipeNetDef(NamedTuple):
    segment_defs: dict[SegmentId, tuple[str, int]]
    """Definitions of identical-path pipeline segments"""
    he_details: dict[SegmentId, tuple[int,list[int]]]
    """Special for heat-exchanging pipes: # of pipetiles in space, coords of ones on turf [for conduction]"""
    @property
    def volume(self) -> float:
        """Returns total volume"""
        ret = 0
        for path, cnt in self.segment_defs.values():
            opts = sorted([(len(k), PIPE_VOLUMES[k]) for k in PIPE_VOLUMES if match_path(k, path)])
            ret += cnt * opts[-1][-1]
        return ret
    @property
    def dominatecls(self) -> str:
        """Returns the path of the longest segment"""
        return sorted([(v,k) for k,v in self.segment_defs.values()])[-1][-1]
    def summarize(self, simplify=simplify_pipecls) -> dict[str, int]:
        ret = {}
        for path, cnt in self.segment_defs.values():
            spath = simplify(path)
            ret[spath] = ret.get(spath, 0) + cnt
        return ret
    def he_summ(self) -> tuple[int, int]: # | None:
        by_d = [(space, len(turfed)) for space,turfed in self.he_details.values()]
        space = sum([space for space,turfed in by_d])
        turfed = sum([turfed for space,turfed in by_d])
        # if space==0 and turfed==0: return None
        return (space, turfed)
    def __str__(self) -> str:
        return f'({self.volume} L, {sum([v for _,v in self.segment_defs.values()])} tiles, {len(self.segment_defs)} segms)'
    def __repr__(self) -> str:
        return self.__str__() + ' : ' + \
            pprint.pformat(self.summarize(), indent=4)
    
@dataclass
class PipeNetsSet:
    segment_map: dict[Coords2, list[tuple[SegmentId, set[Direction], bool]]]
    nets_resolv: dict[SegmentId, str]
    nets_defs: dict[str, PipeNetDef]
    def get_net_in_dir(self, coords: Coords2, dir: Direction|None):
        z = self.segment_map.get(coords)
        if z is None: return None
        for segid, dirs, _ in z:
            if dir is None or ~dir in dirs:
                return self.nets_resolv[segid]
        return None
    def __repr__(self):
        return f'PipeNetsSet(segment_map=<{len(self.segment_map)} coordinates>,\n\t' + \
            f'nets_resolve=<{len(self.nets_resolv)} segments>,\n\tnets_defs={pprint.pformat(self.nets_defs)}\n)'
    def rename(self, old: str, new: str):
        assert old in self.nets_defs, f"{old} net not found"
        self.nets_defs[new] = self.nets_defs[old]
        upds = {k: new for k,v in self.nets_resolv.items() if v==old}
        self.nets_resolv.update(upds)
        del self.nets_defs[old]

def merge_pipes(known_pipes: list[StoredObj], shorten=shorten_pipecls, lengthen=lengthen_pipecls) -> PipeNetsSet:
    seen_pipes: dict[Coords2, list[tuple[SegmentId, set[Direction], bool]]] = {}
    remaps = {}
    pipenets: dict[int, set[Coords2]] = {} 
    pipedefs: dict[int, str] = {} 
    pipejoins: set[tuple[SegmentId,SegmentId]] = set() # pair of pipenet indices
    ix = 1
    # Merge pipes [primary]
    while len(known_pipes) > 0:
        o = known_pipes.pop() # pipe
        k = o.coords
        if k in seen_pipes: GraphLogger.info('?', k)
        # ps = [p for p in set([z for d in o.conn_dirs for z in seen_pipes.get(d.move1(k),[])]) if p]
        ps = [p for p in set([
                z for d in o.conn_dirs for z,r,_ in seen_pipes.get(d.move1(k),[]) if any([d==~d0 for d0 in r])
            ]) if p] # all neighboring
        equ = [p for p in ps if lengthen(pipedefs[p])==o.path] # truly equivalent 
        joint = [*(set(ps) - set(equ))]
        # if 120<=k[0] <=122 and 68<=k[1]<=70:
        #     print(k, ps, equ, joint, o.path)
        if len(equ) == 0: # Truly new pipenet
            pipenets[ix], equ = {k}, [ix]
            pipedefs[ix] = shorten(o.path)
            ix += 1
        if len(equ) == 1: # (incl =0 with new equ) Add pipe-tile to pipenet
            x = equ[0]
            seen_pipes[k] = seen_pipes.get(k,[])+[(x, o.conn_dirs, o.turf=='/turf/space')]
            pipenets[x] |= {k}
        else: #2+ equivalent -- need to merge
            # if 30 in ps: 
            #     print('merge', ps, len(equ))
            sel, *rest = ps
            for p in rest:
                for x in pipenets[p]:
                    seen_pipes[x] = [(sel if v==p else v, dp,q) for v,dp,q in seen_pipes[x]]
                pipenets[sel] |= pipenets[p]
                remaps[p] = sel
                del pipenets[p]
            seen_pipes[k] = seen_pipes.get(k, []) + [(sel, o.conn_dirs, o.turf=='/turf/space')]
            pipenets[sel] |= {k}
        if len(joint) > 0: # Noting the joins, but not merging
            pipejoins |= {tuple(sorted([eq, j])) for j in joint for eq in equ}
    # Gather joins
    pipeseqs = gather_pipejoins(pipejoins,remaps,pipenets.keys())
    # Name segment-seqs [nets]

    namedupes = {}
    nets_resolv = {}
    nets_defs = {}

    he_props = {} # for HE props tuples
    for a,xys in pipenets.items():
        if 'heat_exchanging' not in pipedefs[a]: continue
        is_space = [(q, xy) for xy in xys for v,dp,q in seen_pipes[xy] if v==a]
        tot_spacecnt = sum([q for q,_ in is_space])
        non_spacexs = [xy for q, xy in is_space if not q]
        he_props[a] = (tot_spacecnt, non_spacexs)

    for pipesegments in pipeseqs:
        pipenetdef = PipeNetDef(
            {a: (lengthen(pipedefs[a]), len(pipenets[a])) for a in pipesegments},
            {a: he_props[a] for a in pipesegments if a in he_props}
        )
        candname = pipenetdef.dominatecls.split('/')[-1]
        if candname in namedupes:
            namedupes[candname] += 1
            ix = namedupes[candname]
            candname += f'_{ix}'
        else:
            namedupes[candname] = 1
        nets_defs[candname] = pipenetdef
        nets_resolv.update({a: candname for a in pipesegments})
        
    return PipeNetsSet(seen_pipes, nets_resolv, nets_defs)