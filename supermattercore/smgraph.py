"""Pure Supermatter core [as per ministation] data"""

import networkx as nx
import matplotlib.pyplot as plt 

from .graph import find_amid, input_tag, StoredObj, split_to_obj_pipes, follow_networks, merge_pipes
from .consts import VACUUM_FLOOR, AIRLOCK, SM_CORE

def custom_rename(ssm, stored_objs, pipenetset):
    assert len(stored_objs['PI'])==1
    o = stored_objs['PI'][0]
    pn = pipenetset.get_net_in_dir(o.dir_attr.move1(o.coords), o.dir_attr)
    pipenetset.rename(pn, 'scalding')
    assert len(stored_objs['PO'])==1
    o = stored_objs['PO'][0]
    pn = pipenetset.get_net_in_dir(o.dir_attr.move1(o.coords), o.dir_attr)
    pipenetset.rename(pn, 'warm')
    assert len(stored_objs['T'])==2
    pn,pn1 = [pipenetset.get_net_in_dir(o.dir_attr.move1(o.coords), o.dir_attr) for o in stored_objs['T']]
    o,o1 = stored_objs['T']
    if pn1=='warm': o1,pn1,o,pn = o,pn,o1,pn1
    assert pn=='warm', pn
    del stored_objs['T']
    stored_objs['T(hot)'] = [o]
    stored_objs['T(cold)'] = [o1]
    pipenetset.rename(pn1, 'cooling')
    pn = pipenetset.get_net_in_dir((~o1.dir_attr).move1(o1.coords), ~o1.dir_attr)
    pipenetset.rename(pn, 'cooled')
    gx = tuple([(a+b)//2 for a,b in zip(o.coords, o1.coords)])
    stored_objs['G'] = find_amid(ssm, [gx])

_hedef = lambda herad,hecond: '' if (herad+hecond)==0 else f'(rad {herad}, cond {hecond})'

def gen_graph(stored_objs, pipenetset, reprint = False):
    newnodes = [*[f'{k}({v.volume}L){_hedef(*v.he_summ())}' for k,v in pipenetset.nets_defs.items()],  *[f'{sy}{i+1}' for sy in stored_objs for i in range(len(stored_objs[sy]))]]
    newedges = {('Room', 'SM1'), ('SM1', 'Room'), ('Room', 'PI1'), ('PO1', 'Room')} # , ('T(hot)1', 'G'), ('G', 'T(cold)1')}

    for k in stored_objs:
        for i,o in enumerate(stored_objs[k]):
            syo = f'{k}{i+1}'
            dir_in = input_tag(o.attrs) if k=='F' else ~o.dir_attr
            if reprint: print (f'{syo} @ {o.coords}; d={o.dir_attr}; {dir_in=}')
            conns = [(o.coords,None)] if k == 'M' else [(d.move1(o.coords),d) for d in o.conn_dirs]
            for x,d in conns:
                if z := pipenetset.get_net_in_dir(x,d):
                    hearg = pipenetset.nets_defs[z].he_summ()
                    psn = f'{z}({pipenetset.nets_defs[z].volume}L){_hedef(*hearg)}'
                    if reprint: print(f'\t{d} {"<" if d==dir_in else "-"}- {psn}')
                    newedges |= {(psn,syo) if d==dir_in else (syo, psn)}
                    continue
                o2s = [f'{sy}{jj+1}' for (sy,objs) in stored_objs.items() for (jj,v1) in enumerate(objs) if v1.coords==x]
                if len(o2s):
                    if reprint: print(f'\t{d} {"<" if d==dir_in else "-"}- OBJs: {o2s}')
                    for syo2 in o2s:
                        newedges |= {(syo2,syo) if d==dir_in else (syo, syo2)}
                else:
                    if reprint: print(f'\t{d} -- NOCONN')
    return newnodes, newedges

GREEN_MAINFRAME_FLOOR = "/turf/simulated/floor/greengrid/nitrogen" # dirty hack

def do_all(ssm, do_rename=True):
    core_room = [k for k in [*ssm.find(VACUUM_FLOOR), *ssm.find(GREEN_MAINFRAME_FLOOR)] if AIRLOCK not in dict(ssm[k])]
    smcore = [o for xy in core_room for o in StoredObj.list_from_map_xy(ssm, xy, {SM_CORE: 'SM'})]
    #[StoredObj(k, SM_CORE, dict(ssm[k])[SM_CORE], 'SM') for k in core_room if SM_CORE in dict(ssm[k])]
    assert len(smcore)>0, "SM core not located"
    core_volume = 2.5 * len(core_room)
    # print(f'{core_volume=} m³')
    first = find_amid(ssm, core_room)
    stored_objs, pipes0 = split_to_obj_pipes(follow_networks(ssm, first, smcore))
    pipenetset = merge_pipes(pipes0)
    if do_rename: custom_rename(ssm, stored_objs, pipenetset)
    newnodes,newedges=gen_graph(stored_objs, pipenetset)
    return newnodes, newedges, stored_objs, pipenetset


NODE_SHAPES = {
    'SM1': 'star',
    'PI1': 'triangle',
    'PO1': 'invtriangle',
    '"T(hot)1"': 'trapezium',
    '"T(cold)1"': 'invtrapezium',
    'G1': 'cylinder',
    **{f'T{i}': 'trapezium' for i in range(6)},
    **{f'C{i}': 'doublecircle' for i in range(6)},
    **{f'P{i}': 'pentagon' for i in range(6)},
    **{f'F{i}': 'Msquare' for i in range(5)},
    **{f'HE{i}': 'tripleoctagon' for i in range(5)},
    **{f'M{i}': 'circle' for i in range(5)},
    **{f'V{i}': 'cds' for i in range(5)},
}
NODE_CLRS = {
    'warm': 'cyan',
    'scalding': 'gold',
    'cooling': 'green',
    'cooled': 'chartreuse',
    'T': 'lightgray',
    'heat_exchanging': 'crimson',
}


def plot_graph(label, nodes, edges, ix):
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    Gd = nx.nx_pydot.to_pydot(G)
    for k,v in NODE_SHAPES.items():
        z = Gd.get_node(k)
        if len(z)>0:
            z[0].set_shape(v)
    zl = [k for k in nodes if '(' in k]
    # print(zl)
    for k in zl:
        v = [j for i,j in NODE_CLRS.items() if k.startswith(i)]
        z = Gd.get_node('"'+k+'"')
        # print(k.split('(')[0], v, k, len(z))
        if len(v)>0 and len(z)>0:
            # print(v)
            z[0].set_fillcolor(v[0])
            z[0].set_style('filled')
    Gd.set_label(label)
    Gd.set_colorscheme('X11')
    Gd.write_png(f'./instance/maps/{ix+1}.png')
