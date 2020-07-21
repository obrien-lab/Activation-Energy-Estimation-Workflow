#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 18:04:12 2019

@author: Yang Jiang @ PSU
"""

import os, time, math, traceback, io, sys
import parmed as pmd
import mdtraj as mdt
import numpy as np

################# Arguments #################
if len(sys.argv) != 2:
    print('Error: Wrong number of arguments')
    print('Usage: python auto_backmap_docking.py cntrl_file')
    sys.exit()
cntrl_file = sys.argv[1]

################# Functions #################
def parse_input(cntrl_file):
    ppn = 1
    aa_ref_pdb = ''
    cg_psf_file = ''
    cg_traj_file = ''
    frame_sel = -1
    temp = 310.0
    if_trimerize = False
    trimer_ref_pdb = ''
    docking_cntrl = ''
    temp_seq = ''
    
    f = open(cntrl_file, 'r')
    for line in f.readlines():
        if line.startswith('ppn'):
            ppn = int(line.strip().split('=')[-1].strip())
        elif line.startswith('aa_ref_pdb'):
            aa_ref_pdb = line.strip().split('=')[-1].strip()
        elif line.startswith('cg_psf_file'):
            cg_psf_file = line.strip().split('=')[-1].strip()
        elif line.startswith('cg_traj_file'):
            cg_traj_file = line.strip().split('=')[-1].strip()
        elif line.startswith('frame_sel'):
            frame_sel = int(line.strip().split('=')[-1].strip())
        elif line.startswith('temperature'):
            temp = float(line.strip().split('=')[-1].strip())
        elif line.startswith('if_trimerize'):
            if_trimerize = int(line.strip().split('=')[-1].strip())
            if if_trimerize == 0:
                if_trimerize = False
            else:
                if_trimerize = True
        elif line.startswith('trimer_ref_pdb'):
            trimer_ref_pdb = line.strip().split('=')[-1].strip()
        elif line.startswith('temp_seq'):
            temp_seq = line.strip().split('=')[-1].strip()
        elif line.strip() != '':
            docking_cntrl += line
            
    tag_error = False
    if aa_ref_pdb == '':
        print('Error: No aa_ref_pdb assigned')
        tag_error = True
    if cg_psf_file == '':
        print('Error: No cg_psf_file assigned')
        tag_error = True
    if cg_traj_file == '':
        print('Error: No cg_traj_file assigned')
        tag_error = True
    if if_trimerize:
        if trimer_ref_pdb == '':
            print('Error: No trimer_ref_pdb assigned')
            tag_error = True
    if temp_seq == '':
        print('Error: No temp_seq assigned')
        tag_error = True
    if docking_cntrl == '':
        print('Error: No docking_cntrl assigned')
        tag_error = True
        
    if tag_error:
        sys.exit()
    else:
        return(ppn, aa_ref_pdb, cg_psf_file, cg_traj_file, frame_sel, temp, if_trimerize, trimer_ref_pdb, temp_seq, docking_cntrl)

def check_log_file():
    if not os.path.exists('auto_backmap_docking.log'):
        return 1
    elif not os.path.getsize('auto_backmap_docking.log'):
        return 1
    
    f = open('auto_backmap_docking.log')
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.find('Error:') != -1:
            print('LOG file contains errors. Please delete "Error" lines to restart or check your system to debug.')
            sys.exit()
    if lines[-1].startswith('Done. Elapsed time'):
        print('All done.')
        sys.exit()
    stage = 0
    if lines[-1].startswith('Backmapping AA model'):
        stage = 1
    elif lines[-1].startswith('Refining rebuilt AA model'):
        stage = 2
    elif lines[-1].startswith('Building trimer model'):
        stage = 3
    elif lines[-1].startswith('Docking reactants'):
        stage = 4
    
    f = open('auto_backmap_docking.log', 'w')
    for i in range(len(lines)-1):
        f.write(lines[i])
    f.close()
    return stage

def align_pdb(pdb, ref):
    tar_sel = []
    ref_sel = []
    pmd_pdb = pmd.load_file(pdb)
    for atm in pmd_pdb.atoms:
        if atm.name == 'CA':
            tar_sel.append(atm.idx)
    
    pmd_ref = pmd.load_file(ref)
    for atm in pmd_ref.atoms:
        if atm.name == 'CA':
            ref_sel.append(atm.idx)
            
    traj_pdb = mdt.load_pdb(pdb, standard_names=False)
    traj_ref = mdt.load_pdb(ref, standard_names=False)
    traj_pdb.superpose(traj_ref, atom_indices=tar_sel, ref_atom_indices=ref_sel)
    
    traj_pdb.save_pdb(pdb, force_overwrite=True)

def Amber_minimization_Vacuum(pdb, steps):
    global ppn
    name = pdb.split('/')[-1].split('.pdb')[0] + '_vacuum'
    f = open('leap.in', 'w')
    f.write('''source leaprc.protein.ff14SB
a = loadpdb '''+pdb+'''
saveamberparm a '''+name+'''.prmtop '''+name+'''.inpcrd
quit''')
    f.close()
    os.system('tleap -s -f leap.in > /dev/null')
    
    f = open('min1.in', 'w')
    f.write('''minimization in Vacuum
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps)+''',
  ncyc   = '''+str(int(steps/2))+''',
  ntb    = 0,
  cut    = 999,
  dielc  = 78.5,
  ntr = 1, 
  restraintmask='@CA',
  restraint_wt=30.0,
 /
''')
    f.close()
    os.system('mpirun -np '+str(ppn)+' sander.MPI -O -i min1.in -p '+name+'.prmtop -c '+name+'.inpcrd -r '+name+'_min1.ncrst -o '+name+'_min1.out -ref '+name+'.inpcrd -inf mdinfo.log')

def trans_crd(top, crd, out, mask):
    struct = pmd.load_file(top)
    cor = pmd.load_file(crd)
    struct.coordinates = cor.coordinates
    struct[mask].save(out, overwrite=True)
    
def Amber_minimization_water(pdb, steps_1, steps_2):
    global ppn
    
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
    
    if steps_1 > 0:
        f = open('min3.in', 'w')
        f.write('''minimization for trimer in TIP3P water box with restraints
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps_1)+''',
  ncyc   = '''+str(int(steps_1/2))+''',
  cut=10.0, 
  ntb=1, 
  ntr=1,
  restraintmask='@CA',
  restraint_wt=100.0,
 /
''')
        f.close()
        os.system('mpirun -np '+str(ppn)+' pmemd.MPI -O -i min3.in -p '+name+'.prmtop -c '+crd+' -r '+name+'_min3.ncrst -o '+name+'_min3.out -ref '+crd+' -inf mdinfo.log')
        crd = name+'_min3.ncrst'
    
    if steps_2 > 0:
        f = open('min4.in', 'w')
        f.write('''minimization for trimer in TIP3P water box with restraints
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps_2)+''',
  ncyc   = '''+str(int(steps_2/2))+''',
  cut=10.0, 
  ntb=1,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=10.0,
 /
''')
        f.close()
        os.system('mpirun -np '+str(ppn)+' pmemd.MPI -O -i min4.in -p '+name+'.prmtop -c '+crd+' -r '+name+'_min_final.ncrst -o '+name+'_min_final.out -ref '+crd+' -inf mdinfo.log')

################# MAIN ####################
start_time = time.time()
(ppn, aa_ref_pdb, cg_psf_file, cg_traj_file, frame_sel, temp, if_trimerize, 
 trimer_ref_pdb, temp_seq, docking_cntrl) = parse_input(cntrl_file)

stage = check_log_file()
if stage == 0:
    print('Error: Cannot find checkpoint. Please check the last line of your LOG file')
    sys.exit()

if stage == 1:
    f_log = open('auto_backmap_docking.log', 'w')
    f_log.write('Backmapping AA model...\n')
    f_log.close()
    if os.path.exists('rebuild_cg'):
        os.system('rm -rf rebuild_cg')
    traj = mdt.load(cg_traj_file, top=cg_psf_file)
    traj = traj[frame_sel].center_coordinates()
    traj.save('cg.pdb', force_overwrite=True)
    os.system('backmap.py -i '+aa_ref_pdb+' -c cg.pdb > backmapping.log 2>&1')
    if not os.path.exists('rebuild_cg/cg_mini_PD2_min_pulchra.pdb'):
        f_log = open('auto_backmap_docking.log', 'a')
        f_log.write('Error: Failed to backmap AA model from cg.pdb\n')
        f_log.close()
        sys.exit()

if stage <= 2:
    pdb = pmd.load_file('rebuild_cg/cg_mini_PD2_min_pulchra.pdb')
    # Modify protonated states according to the res name in temp_seq #
    f = open(temp_seq)
    seq = f.readlines()
    f.close()
    seq = [s.strip() for s in seq]
    for res in pdb.residues:
        res.name = seq[res.idx]
    pdb.save('cg_rebuilt.pdb', overwrite=True)
    
    f_log = open('auto_backmap_docking.log', 'a')
    f_log.write('Refining rebuilt AA model...\n')
    f_log.close()
    if os.path.exists('refine'):
        os.system('rm -rf refine')
    os.system('mkdir refine')
    os.chdir('refine')
    Amber_minimization_Vacuum('../cg_rebuilt.pdb', 1000)
    trans_crd('cg_rebuilt_vacuum.prmtop', 'cg_rebuilt_vacuum_min1.ncrst', 'cg_rebuilt_refv.pdb', ':*')
    Amber_minimization_water('cg_rebuilt_refv.pdb', 2000, 2000)
    trans_crd('cg_rebuilt_refv_water.prmtop', 'cg_rebuilt_refv_water_min_final.ncrst', '../protein.pdb', '(!@/H)&(!:WAT,Na+)')
    os.chdir('../')
    

protein_pdb_file = 'protein.pdb'
if if_trimerize:
    if stage <= 3:
        f_log = open('auto_backmap_docking.log', 'a')
        f_log.write('Building trimer model...\n')
        f_log.close()
        if os.path.exists('trimerization'):
            os.system('rm -rf trimerization')
        os.system('trimerization_sd.py '+trimer_ref_pdb+' protein.pdb '+aa_ref_pdb+' '+str(ppn)+' '+('%.2f'%temp))
        if not os.path.exists('trimer.pdb'):
            f_log = open('auto_backmap_docking.log', 'a')
            f_log.write('Error: Fail to build trimer model from protein.pdb\n')
            f_log.close()
            sys.exit()
        else:
            protein_pdb_file = 'trimer.pdb'
        align_pdb(protein_pdb_file, trimer_ref_pdb)
    else:
        protein_pdb_file = 'trimer.pdb'
elif stage <= 3:    
    align_pdb(protein_pdb_file, aa_ref_pdb)

if stage <= 4:
    f_log = open('auto_backmap_docking.log', 'a')
    f_log.write('Docking reactants to %s...\n'%protein_pdb_file)
    f_log.close()
    f = open('docking.cntrl', 'w')
    f.write('protein_pdb_file = %s\n'%(protein_pdb_file))
    f.write(docking_cntrl+'\n')
    f.write('ppn = %d\n'%ppn)
    f.write('temperature = %.2f\n'%temp)
    f.close()
    os.system('docking_mr_sd_noP.py docking.cntrl')
    
    if not os.path.exists('reactant_complex.ncrst'):
        f_log = open('auto_backmap_docking.log', 'a')
        f_log.write('Error: Failed to predict reactants\' docking pose\n')
        f_log.close()
        sys.exit()
    #elif not os.path.exists('product_complex.ncrst'):
        #f_log = open('auto_backmap_docking.log', 'a')
        #f_log.write('Error: Failed to predict products\' docking pose\n')
        #f_log.close()
        #sys.exit()
    
    end_time = time.time()
    m, s = divmod(end_time-start_time, 60)
    h, m = divmod(m, 60)
    f_log = open('auto_backmap_docking.log', 'a')
    f_log.write('Done. Elapsed time: %d:%d:%d\n'%(h, m, s))
    f_log.close()
