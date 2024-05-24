"""Copying code used to create [PR #4050](https://github.com/NebulaSS13/Nebula/pull/4050)
from main Jupyter Notebook to preserve for latter"""

import git
import difflib
import os.path
import re

ROOT_ROOT = "E:\\Code\\opensource\\games\\"
SCAV_ROOT = ROOT_ROOT+ "ScavStation\\"
NEB_ROOT = ROOT_ROOT + "Nebula\\"

repo = git.Repo(NEB_ROOT)

comm1 = repo.commit("b8a70fa94ab4c3fd0904e55f1edd8640646f99e0") # Penny's com
assert len(comm1.parents)==1
comm2 = comm1.parents[0]

possible_edits = {}
lost_files = []
for diff1 in comm1.diff(comm2).iter_change_type('M'):
    a_lns = diff1.a_blob.data_stream.read().decode().split('\n')
    b_lns = diff1.b_blob.data_stream.read().decode().split('\n')
    ops = [x for x in difflib.SequenceMatcher(None,a_lns,b_lns,False).get_opcodes() if x[0]=='replace']
    needs_care = {}
    for _, i1,i2,j1,j2 in ops:
        if j2-j1 != i2-i1:
            print(diff1.a_path, a_lns[i1:i2], b_lns[j1:j2])
            continue
        for i3,s3 in zip(range(i1,i2), b_lns[j1:j2]):
            if 'between' in s3 or 'dd_range' in s3:
                needs_care[i3] = s3
    if len(needs_care) == 0:
        continue # does not need edits
    possible_edits[diff1.a_path] = needs_care
    if not os.path.exists(NEB_ROOT + '/' + diff1.a_path):
        lost_files.append(diff1.a_path)


autoedits = {}
questions = {}
losts = {}
for diff2 in repo.commit().diff(comm1):
    if diff2.b_path not in possible_edits:
        # not a file we care about
        continue
    if diff2.a_blob is None:
        print(diff2.b_path, 'del?', diff2.a_path, diff2.a_blob)
        losts[diff2.b_path] = possible_edits[diff2.b_path]
        continue
    if diff2.b_path in lost_files: 
        print('lost file found')
    a_lns = diff2.a_blob.data_stream.read().decode().split('\n')
    b_lns = diff2.b_blob.data_stream.read().decode().split('\n')
    carelns = possible_edits[diff2.b_path]
    ops = [x for x in difflib.SequenceMatcher(None,a_lns,b_lns,False).get_opcodes()]
    eds = {}
    quests = {}
    for op,i1,i2,j1,j2 in ops:
        for j3 in range(j1,j2):
            if j3 not in carelns: 
                # we don't care about other lines
                continue
            if op == 'equal':
                # line is preserved since comm1
                assert i2-i1==j2-j1
                i3 = i1+(j3-j1)
                eds[i3] = (a_lns[i3], carelns[j3])
            else:
                quests[j3] = ((op,i1,i2,j1,j2), diff2.b_path, j3, carelns[j3])
    if len(eds):
        autoedits[diff2.a_path] = eds
    if len(quests):
        questions[diff2.a_path] = quests


# generate autoedits
autopatches = {}
for fn, linedic in autoedits.items():
    pats = {}
    for i, (nline, oline) in linedic.items():
        dops = [x[8:] for x in re.findall(r'dd_range\(.+?,.+?,.+?\)', oline)]
        bops = [x[7:] for x in re.findall(r'between\(.+?,.+?,.+?\)', oline)]
        nops = [x[5:] for x in re.findall(r'clamp\(.+?,.+?,.+?\)', nline)]
        if len(bops) == 0 and len(dops) == 0: 
            continue
        assert bops == nops or dops == nops, (oline, bops, dops)
        if len(dops) == 0:
            assert 'between' in oline
            nl = re.sub(r'clamp\(\s*(.+?)\s*,\s*(.+?)\s*,\s*(.+?)\s*\)', r'clamp(\2, \1, \3)', nline)
        else:
            assert 'dd_range' in oline, oline
            nl = re.sub(r'clamp\(\s*(.+?)\s*,\s*(.+?)\s*,\s*(.+?)\s*\)', r'clamp(\3, \1, \2)', nline)
        pats[i] = (nl,nline,oline)
    autopatches[fn] = pats

# apply autopatches
for fn,repdic in autopatches.items():
    with open(NEB_ROOT+'/'+fn, 'r') as f:
        ls = f.readlines()
    for i, (neu,old,oldest) in repdic.items():
        assert ls[i] == old+'\n', (ls[i], old)
        ls[i] = neu+'\n'
    with open(NEB_ROOT+'/'+fn, 'w') as f:
        f.writelines(ls)