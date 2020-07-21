#!/usr/bin/env python3

import os, sys
import parmed as pmd
import numpy as np

################# Arguments #################
trimer_ref_pdb = sys.argv[1]
monomer_pdb = sys.argv[2]
monomer_ref_pdb = sys.argv[3]
ppn = int(sys.argv[4])
temp = 310
if len(sys.argv) == 6:
    temp = float(sys.argv[5])
dist_restaint_res_list = [[25, 25], [32, 32], [150, 150], [157, 157], [25, 157], [32, 150]]
fit_sel_1 = ':1-24,34-83,106-121,137-149,158-168,183-213@CA'
fit_sel_3 = ':1-24,34-83,106-121,137-149,158-168,183-213,214-237,247-296,319-334,350-362,371-381,396-426,427-450,460-509,532-547,563-575,584-594,609-639@CA'
interface_sel_1 = ':25-33,150-157@CA'
interface_sel_3 = ':25-33,150-157,238-246,363-370,451-459,576-583@CA'

################# Functions #################
def clean_pdb(pdb, out_dir):
    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP'];
    name = pdb.split('/')[-1].split('.pdb')[0]
    struct = pmd.load_file(pdb)
    sel_idx = np.zeros(len(struct.atoms))
    for idx, res in enumerate(struct.residues):
        res.number = idx+1
        if res.name in AA_name_list:
            for atm in res.atoms:
                sel_idx[atm.idx] = 1
    struct[sel_idx].save(out_dir+'/'+name+'_clean.pdb', overwrite=True)
    return name+'_clean.pdb'

def get_trimerization_results(out):
    interface_area = None
    ACE = None
    f = open(out, 'r')
    for line in f:
        line = line.strip()
        if line.startswith('1 | '):
            words = line.split(' | ')
            interface_area = float(words[3])
            ACE = float(words[7])
            break
    f.close()
    return (interface_area, ACE)

def test_trimer(pdb):
    test_resid = [140, 401]
    test_dist = 20
    tar_pdb = pmd.load_file(pdb)
    test_atmid = []
    for resid in test_resid:
        res = tar_pdb.residues[resid]
        for atm in res.atoms:
            if atm.name == 'CA':
                test_atmid.append(atm.idx)
    dist = np.sum((tar_pdb.coordinates[test_atmid[0]] - tar_pdb.coordinates[test_atmid[1]])**2)**0.5
    if dist > test_dist:
        new_tar_pdb = tar_pdb[':1-213']
        new_tar_pdb = new_tar_pdb + tar_pdb[':427-639']
        new_tar_pdb = new_tar_pdb + tar_pdb[':214-426']
        new_tar_pdb.save(pdb, overwrite=True)

def Amber_TGT_mini_Vacuum(pdb, ref_pdb, steps, fit_sel, interface_sel, frc):
    # Run Amber targeted minimization
    global ppn
    name = pdb.split('/')[-1].split('.pdb')[0] + '_vacuum'
    f = open('leap.in', 'w')
    f.write('''source leaprc.protein.ff14SB
a = loadpdb '''+pdb+'''
saveamberparm a '''+name+'''.prmtop '''+name+'''.inpcrd
quit''')
    f.close()
    os.system('tleap -s -f leap.in > /dev/null')
    
    ref_pdb = clean_pdb(ref_pdb, './')
    f = open('leap.in', 'w')
    f.write('''source leaprc.protein.ff14SB
a = loadpdb '''+ref_pdb+'''
saveamberparm a ref.prmtop ref.inpcrd
quit''')
    f.close()
    os.system('tleap -s -f leap.in > /dev/null')
    
    f = open('min1.in', 'w')
    f.write('''Targeted minimization in Vacuum
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps)+''',
  ncyc   = '''+str(int(steps/2))+''',
  ntb    = 0,
  cut    = 999,
  dielc  = 78.5,
  itgtmd = 1,
  tgtmdfrc = '''+str(frc)+'''
  tgtfitmask = \''''+fit_sel+'''\'
  tgtrmsmask = \''''+interface_sel+'''\'
 /
''')
    f.close()
    os.system('mpirun -np '+str(ppn)+' sander.MPI -O -i min1.in -p '+name+'.prmtop -c '+name+'.inpcrd -r '+name+'_min1.ncrst -o '+name+'_min1.out -ref ref.inpcrd -inf mdinfo.log')

def trans_crd(top, crd, out, mask):
    struct = pmd.load_file(top)
    cor = pmd.load_file(crd)
    struct.coordinates = cor.coordinates
    struct[mask].save(out, overwrite=True)
    
def Amber_equilibrium_water(pdb, steps_1, steps_2, steps_3):
    global ppn, temp
    
    name = pdb.split('/')[-1].split('.pdb')[0] + '_water'
    
    f = open('leap.in', 'w')
    f.write('''source leaprc.protein.ff14SB
source leaprc.water.tip3p
a = loadpdb '''+pdb+'''
solvatebox a TIP3PBOX 10.0
addions a Na+ 0.0
saveamberparm a '''+name+'''.prmtop '''+name+'''.inpcrd
quit''')
    f.close()
    os.system('tleap -s -f leap.in > /dev/null')
    
    crd = name+'.inpcrd'
    
    f = open('min2.in', 'w')
    f.write('''minimization for trimer in TIP3P water box
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps_1)+''',
  ncyc   = '''+str(int(steps_1/2))+''',
  cut=10.0, 
  ntb=1,
  ntc=2,
  ntf=2,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=10.0,
 /
''')
    f.close()
    os.system('mpirun -np '+str(ppn)+' pmemd.MPI -O -i min2.in -p '+name+'.prmtop -c '+crd+' -r '+name+'_min2.ncrst -ref '+crd+' -o '+name+'_min2.out -inf mdinfo.log')
    crd = name+'_min2.ncrst'
    
    f = open('heat.in', 'w')
    f.write('''heating for trimer in TIP3P water box
 &cntrl
  imin=0, ntx=1, irest=0,
  nstlim='''+str(steps_2)+''', dt=0.002,
  ntf=2, ntc=2,
  tempi=0.0, temp0='''+str(temp)+''',
  ntpr=100, ntwx=500,
  cut=10.0, ntb=1, ntc=2, ntf=2,
  ntp=0, ntt=3,
  gamma_ln=2.0,
  ntr=1,
  restraint_wt=10,
  restraintmask='@CA',
  nmropt=1,
  ig=-1,
 /
&wt type='TEMP0', istep1=0, istep2='''+str(int(steps_2/10*9))+''', value1=0.0, value2='''+str(temp)+''' /
&wt type='TEMP0', istep1='''+str(int(steps_2/10*9)+1)+''', istep2='''+str(steps_2)+''', value1='''+str(temp)+''', value2='''+str(temp)+''' /
&wt type='END' /
''')
    f.close()
    os.system('mpirun -np '+str(ppn)+' pmemd.MPI -O -i heat.in -p '+name+'.prmtop -c '+crd+
              ' -r '+name+'_heat.ncrst -o '+name+'_heat.out -x '+name+'_heat.nc -ref '+
              crd+' -inf mdinfo.log')
    crd = name+'_heat.ncrst'
    
    f = open('equil.in', 'w')
    f.write('''equilibration for trimer in TIP3P water box
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim='''+str(steps_3)+''', dt=0.002,
  ntf=2, ntc=2,
  temp0='''+str(temp)+''',
  ntpr=1000, ntwx=5000,
  cut=10.0, ntb=2, ntc=2, ntf=2,
  ntp=1, ntt=3,
  gamma_ln=2.0,
  ig=-1,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=5.0,
 /
''')
    f.close()
    os.system('mpirun -np '+str(ppn)+' pmemd.MPI -O -i equil.in -p '+name+
              '.prmtop -c '+crd+' -r '+name+'_equil.ncrst -o '+name+'_equil.out -x '+
              name+'_equil.nc -ref '+crd+' -inf mdinfo.log')
    crd = name+'_equil.ncrst'
    
    f = open('min3.in', 'w')
    f.write('''minimization for trimer in TIP3P water box
 &cntrl
  imin   = 1,
  maxcyc = 500,
  ncyc   = 250,
  cut=10.0, 
  ntb=1, 
  ntr=1,
  restraintmask='@CA',
  restraint_wt=10.0,
 /
''')
    f.close()
    os.system('mpirun -np '+str(ppn)+' pmemd.MPI -O -i min3.in -p '+name+'.prmtop -c '+crd+' -r '+name+'_min3.ncrst -o '+name+'_min3.out -ref '+crd+' -inf mdinfo.log')

################# MAIN ####################
os.system('mkdir trimerization')

dist_restraint_dist = [[0,0] for i in range(len(dist_restaint_res_list))]
struct = pmd.load_file(trimer_ref_pdb)
res_start_number = struct[0].residue.number
idx_list = [[0,0] for i in range(len(dist_restaint_res_list))]
for i in range(len(dist_restaint_res_list)):
    for res in struct.residues:
        if res.chain == 'A' and (res.number - res_start_number + 1) == dist_restaint_res_list[i][0]:
            for atm in res.atoms:
                if atm.name == 'CA':
                    idx_list[i][0] = atm.idx
                    break
        if res.chain == 'B' and (res.number - res_start_number + 1) == dist_restaint_res_list[i][1]:
            for atm in res.atoms:
                if atm.name == 'CA':
                    idx_list[i][1] = atm.idx
                    break

for i in range(len(idx_list)):
    atm_1 = struct[idx_list[i][0]]
    atm_2 = struct[idx_list[i][1]]
    dist = (atm_1.xx - atm_2.xx)**2 + (atm_1.xy - atm_2.xy)**2 + (atm_1.xz - atm_2.xz)**2
    dist = dist**0.5
    if dist_restaint_res_list[i][0] == dist_restaint_res_list[i][1]:
        dist_restraint_dist[i] = [dist-0.5, 999]
    else:
        dist_restraint_dist[i] = [0, dist+1]

idx_list = [[0,0] for i in range(len(dist_restaint_res_list))]
struct = pmd.load_file(monomer_pdb)
for j in range(len(dist_restaint_res_list)):
    for res in struct.residues:
        if res.number == dist_restaint_res_list[j][0]:
            for atm in res.atoms:
                if atm.name == 'CA':
                    idx_list[j][0] = atm.idx+1
                    break
        if res.number == dist_restaint_res_list[j][1]:
            for atm in res.atoms:
                if atm.name == 'CA':
                    idx_list[j][1] = atm.idx+1
                    break
f = open('trimerization/dist_restraint.dat', 'w')
for j in range(len(idx_list)):
    f.write('%d %d %.4f %.4f\n'%(idx_list[j][0], idx_list[j][1], 
                                 dist_restraint_dist[j][0], dist_restraint_dist[j][1]))
f.close()
        

os.chdir('trimerization')
#refine monomer interface
Amber_TGT_mini_Vacuum('../'+monomer_pdb, '../'+monomer_ref_pdb, 2000, fit_sel_1, interface_sel_1, 100)
name = monomer_pdb.split('/')[-1].split('.pdb')[0]
trans_crd(name+'_vacuum.prmtop', name+'_vacuum_min1.ncrst', 'mono.pdb', '!@/H')

# predict trimer structure
os.system('buildParams.pl 3 mono.pdb')
f = open('params.txt', 'r')
fo = open('params_1.txt', 'w')
for line in f:
    if line.strip().startswith('#distanceConstraintsFile file_name'):
        fo.write('distanceConstraintsFile dist_restraint.dat\n')
    else:
        fo.write(line)
f.close()
fo.close()
os.system('rm -f params.txt')
os.system('mv params_1.txt params.txt')
        
os.system('symm_dock.Linux params.txt build_trimer.out > /dev/null 2>&1')
(int_area, ACE) = get_trimerization_results('build_trimer.out')
if int_area != None and ACE != None:
    os.system('transOutput.pl build_trimer.out 1 1')
    os.system('mv result.1.pdb trimer.pdb')
    test_trimer('trimer.pdb')
    Amber_TGT_mini_Vacuum('trimer.pdb', '../'+trimer_ref_pdb, 2000, fit_sel_3, interface_sel_3, 10)
    trans_crd('trimer_vacuum.prmtop', 'trimer_vacuum_min1.ncrst', 'trimer_refv.pdb', ':*')
    Amber_equilibrium_water('trimer_refv.pdb', 2000, 10000, 250000)
    trans_crd('trimer_refv_water.prmtop', 'trimer_refv_water_min3.ncrst', '../trimer.pdb', '(!@/H)&(!:WAT,Na+)')
else:
    print('Error: cannot form trimer')
