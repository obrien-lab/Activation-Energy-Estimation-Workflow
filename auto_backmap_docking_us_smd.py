#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 18:04:12 2019

@author: Yang Jiang @ PSU
"""

import os, time, math, traceback, io, sys
import parmed as pmd
import numpy as np

################# Arguments #################
if len(sys.argv) != 2:
    print('Error: Wrong number of arguments')
    print('Usage: python auto_backmap_docking_us_smd.py cntrl_file')
    sys.exit()
cntrl_file = sys.argv[1]

################# Functions #################
def parse_input(cntrl_file):
    tpn = 20
    aa_ref_pdb = ''
    cg_psf_file = ''
    cg_traj_file = ''
    frame_sel = -1
    temp = 310
    if_trimerize = False
    trimer_ref_pdb = ''
    adt_home = ''
    ligand_file = []
    box_center = []
    box_dimension = []
    refine_reactant_restraints = ''
    QM_reactant_restraints = ''
    temp_seq = ''
    qm_mask = ''
    qm_charge = None
    if_us = False
    us_reaction_coordinate = ''
    us_rc_list = ''
    us_min_step = 500
    us_sampling_step = 10000
    us_force_constant = 250
    wham_exec = 'wham'
    wham_rc_interval = 0.05
    us_extra_restraint = ''
    us_position_restraint = ''
    if_smd = False
    smd_ligand_list = []
    num_smd_replica = 5
    num_smd_stage = 4
    pull_length = 50
    pull_force_constant = 100
    pull_lig_mask = []
    pull_prot_mask = []
    smd_steps_per_stage = 200000
    
    f = open(cntrl_file, 'r')
    for line in f.readlines():
        if line.startswith('#'):
            continue
        elif line.startswith('tpn'):
            tpn = int(line.strip().split('=')[-1].strip())
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
        elif line.startswith('adt_home'):
            adt_home = line.strip().split('=')[-1].strip()
        elif line.startswith('ligand_file'):
            ligand_file.append(line.strip().split('=')[-1].strip())
        elif line.startswith('box_center'):
            box_center.append(line.strip().split('=')[-1].strip())
        elif line.startswith('box_dimension'):
            box_dimension.append(line.strip().split('=')[-1].strip())
        elif line.startswith('refine_reactant_restraints'):
            refine_reactant_restraints = '='.join(line.strip().split('=')[1:])
        elif line.startswith('QM_reactant_restraints'):
            QM_reactant_restraints = line.strip().split('=')[-1].strip()
        elif line.startswith('qm_mask'):
            qm_mask = '='.join(line.strip().split('=')[1:])
        elif line.startswith('qm_charge'):
            qm_charge = float(line.strip().split('=')[-1].strip())
        elif line.startswith('temp_seq'):
            temp_seq = line.strip().split('=')[-1].strip()
            
        elif line.startswith('if_us'):
            if_us = int(line.strip().split('=')[-1].strip())
            if if_us == 0:
                if_us = False
            else:
                if_us = True
        elif line.startswith('us_reaction_coordinate'):
            us_reaction_coordinate = '='.join(line.strip().split('=')[1:])
        elif line.startswith('us_rc_list'):
            us_rc_list = line.strip().split('=')[-1].strip()
        elif line.startswith('us_min_step'):
            us_min_step = int(line.strip().split('=')[-1].strip())
        elif line.startswith('us_sampling_step'):
            us_sampling_step = int(line.strip().split('=')[-1].strip())
        elif line.startswith('us_force_constant'):
            us_force_constant = float(line.strip().split('=')[-1].strip())
        elif line.startswith('cpptraj_exec'):
            cpptraj_exec = line.strip().split('=')[-1].strip()
        elif line.startswith('wham_exec'):
            wham_exec = line.strip().split('=')[-1].strip()
        elif line.startswith('wham_rc_interval'):
            wham_rc_interval = float(line.strip().split('=')[-1].strip())
        elif line.startswith('extra_restraint'):
            us_extra_restraint = '='.join(line.strip().split('=')[1:])
        elif line.startswith('position_restraint'):
            us_position_restraint = '='.join(line.strip().split('=')[1:])
        
        elif line.startswith('if_smd'):
            if_smd = int(line.strip().split('=')[-1].strip())
            if if_smd == 0:
                if_smd = False
            else:
                if_smd = True
        elif line.startswith('num_smd_replica'):
            num_smd_replica = int(line.strip().split('=')[-1].strip())
        elif line.startswith('num_smd_stage'):
            num_smd_stage = int(line.strip().split('=')[-1].strip())
        elif line.startswith('pull_length'):
            pull_length = float(line.strip().split('=')[-1].strip())
        elif line.startswith('pull_force_constant'):
            pull_force_constant = float(line.strip().split('=')[-1].strip())
        elif line.startswith('smd_ligand_idx'):
            words = line.strip().split('=')[-1].strip().split()
            words = [int(w) for w in words]
            smd_ligand_list.append(words)
        elif line.startswith('pull_lig_mask'):
            pull_lig_mask.append(line.strip().split('=')[-1].strip())
        elif line.startswith('pull_prot_mask'):
            pull_prot_mask.append(line.strip().split('=')[-1].strip())
        elif line.startswith('smd_steps_per_stage'):
            smd_steps_per_stage = int(line.strip().split('=')[-1].strip())
            
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
    if adt_home == '':
        print('Error: No adt_home assigned')
        tag_error = True
    if ligand_file == []:
        print('Error: No ligand_file assigned')
        tag_error = True
    if box_center == []:
        print('Error: No box_center assigned')
        tag_error = True
    if box_dimension == []:
        print('Error: No box_dimension assigned')
        tag_error = True
    if refine_reactant_restraints == '':
        print('Error: No refine_reactant_restraints assigned')
        tag_error = True
    if QM_reactant_restraints == '':
        print('Error: No QM_reactant_restraints assigned')
        tag_error = True
    if qm_mask == '':
        print('Error: No qm_mask assigned')
        tag_error = True
    if qm_charge == None:
        print('Error: No qm_charge assigned')
        tag_error = True
    if us_reaction_coordinate == '' and if_us == True:
        print('Error: No us_reaction_coordinate assigned')
        tag_error = True
    if us_rc_list == '' and if_us == True:
        print('Error: No us_rc_list assigned')
        tag_error = True
    if if_smd == True:
        if smd_ligand_list == []:
            print('Error: No smd_ligand_idx assigned')
            tag_error = True
        if pull_lig_mask == []:
            print('Error: No pull_lig_mask assigned')
            tag_error = True
        if pull_prot_mask == []:
            print('Error: No pull_prot_mask assigned')
            tag_error = True
        if len(smd_ligand_list) != len(pull_lig_mask) or len(pull_lig_mask) != len(pull_prot_mask):
            print('Error: Number of smd_ligand_idx, pull_lig_mask and pull_prot_mask must be the same (%d, %d, %d)'%(
                    len(smd_ligand_list), len(pull_lig_mask), len(pull_prot_mask)))
            tag_error = True
        
    if tag_error:
        sys.exit()
    else:
        return(tpn, aa_ref_pdb, cg_psf_file, cg_traj_file, frame_sel, temp,
               if_trimerize, trimer_ref_pdb, adt_home, ligand_file, box_center,
               box_dimension, refine_reactant_restraints, QM_reactant_restraints, 
               temp_seq, qm_mask, qm_charge, if_us, us_reaction_coordinate, 
               us_rc_list, us_min_step, us_sampling_step, us_force_constant, 
               cpptraj_exec, wham_exec, wham_rc_interval, us_extra_restraint, 
               us_position_restraint, if_smd, num_smd_replica, 
               num_smd_stage, pull_length, pull_force_constant, smd_ligand_list, 
               pull_lig_mask, pull_prot_mask, smd_steps_per_stage)

def check_log_file(): 
    if not os.path.exists('workflow.log'):
        return 1
    elif not os.path.getsize('workflow.log'):
        return 1
    
    f = open('workflow.log')
    lines = f.readlines()
    f.close()
    stage = 0
    line_idx = 0
    for idx, line in enumerate(lines):
        if line.find('Error:') != -1:
            print('LOG file contains errors. Please delete "Error" lines to restart or check your system to debug.')
            sys.exit()
        elif line.startswith('Auto backmapping and docking'):
            stage += 1
            line_idx = idx
        elif line.startswith('Running US'):
            stage += 1
            line_idx = idx
        elif line.startswith('Running ASMD for ligand '):
            stage += 1
            line_idx = idx
    if lines[-1].startswith('Done.'):
        print('All done.')
        sys.exit()
    
    f = open('workflow.log', 'w')
    for i in range(line_idx):
        f.write(lines[i])
    f.close()
    return stage

################# MAIN ####################
start_time = time.time()
(tpn, aa_ref_pdb, cg_psf_file, cg_traj_file, frame_sel, temp,
 if_trimerize, trimer_ref_pdb, adt_home, ligand_file, box_center,
 box_dimension, refine_reactant_restraints, QM_reactant_restraints, 
 temp_seq, qm_mask, qm_charge, if_us, us_reaction_coordinate, 
 us_rc_list, us_min_step, us_sampling_step, us_force_constant, 
 cpptraj_exec, wham_exec, wham_rc_interval, us_extra_restraint, 
 us_position_restraint, if_smd, num_smd_replica, 
 num_smd_stage, pull_length, pull_force_constant, smd_ligand_list, 
 pull_lig_mask, pull_prot_mask, smd_steps_per_stage) = parse_input(cntrl_file)

struct = pmd.load_file(cg_psf_file)
chain_len = len(struct.residues)
if if_trimerize:
    chain_len = 3*chain_len

stage = check_log_file()
if stage == 0:
    print('Error: Cannot find checkpoint. Please check the last line of your LOG file')
    sys.exit()

if stage == 1:
    f_log = open('workflow.log', 'w')
    f_log.write('Auto backmapping and docking\n')
    f_log.close()
    os.system('mkdir auto_backmap_docking')
    os.chdir('auto_backmap_docking')
    f = open('auto_backmap_docking.cntrl', 'w')
    f.write('ppn = %d\n'%tpn)
    f.write('aa_ref_pdb = ../%s\n'%aa_ref_pdb)
    f.write('cg_psf_file = ../%s\n'%cg_psf_file)
    f.write('cg_traj_file = ../%s\n'%cg_traj_file)
    f.write('frame_sel = %d\n'%frame_sel)
    f.write('temperature = %.2f\n'%temp)
    if if_trimerize:
        f.write('if_trimerize = 1\n')
        f.write('trimer_ref_pdb = ../%s\n'%trimer_ref_pdb)
    f.write('temp_seq = ../%s\n'%temp_seq)
    f.write('adt_home = %s\n'%adt_home)
    for i, r in enumerate(ligand_file):
        f.write('ligand_file_%d = ../%s\n'%(i+1, r))
    for i, r in enumerate(box_center):
        f.write('box_center_%d = %s\n'%(i+1, r))
    for i, r in enumerate(box_dimension):
        f.write('box_dimension_%d = %s\n'%(i+1, r))
    f.write('refine_reactant_restraints = %s\n'%refine_reactant_restraints)
    f.write('QM_reactant_restraints = %s\n'%QM_reactant_restraints)
    f.write('temp_seq_file = ../%s\n'%temp_seq)
    f.write('qm_mask = %s\n'%qm_mask)
    f.write('qm_charge = %d\n'%qm_charge)
    f.close()
    os.system('auto_backmap_docking_sd.py auto_backmap_docking.cntrl')
    os.chdir('../')

    if (not os.path.exists('auto_backmap_docking/complex.prmtop') 
        or not os.path.exists('auto_backmap_docking/reactant_complex.ncrst')):
        f_log = open('workflow.log', 'a')
        f_log.write('Error: cannot predict docking poses\n')
        f_log.close()
        sys.exit()

if stage <= 2 and if_us == True:    
    f_log = open('workflow.log', 'a')
    f_log.write('Running US\n')
    f_log.close()
    os.system('mkdir US')
    os.chdir('US')
    f = open('us.cntrl', 'w')
    f.write('''top = ../auto_backmap_docking/complex.prmtop
reactant_cor = ../auto_backmap_docking/reactant_complex.ncrst
''')
    f.write('np = %d\n'%tpn)
    f.write('qm_mask = %s\n'%qm_mask)
    f.write('qm_charge = %d\n'%qm_charge)
    f.write('reaction_coordinate = %s\n'%us_reaction_coordinate)
    f.write('rc_list = %s\n'%us_rc_list)
    f.write('temp = %.2f\n'%temp)
    f.write('min_step = %d\n'%us_min_step)
    f.write('sampling_step = %d\n'%us_sampling_step)
    f.write('force_constant = %.4f\n'%us_force_constant)
    f.write('cpptraj_exec = %s\n'%cpptraj_exec)
    f.write('wham_exec = %s\n'%wham_exec)
    f.write('wham_rc_interval = %.4f\n'%wham_rc_interval)
    if us_extra_restraint != '':
        f.write('extra_restraint = %s\n'%us_extra_restraint)
    if us_position_restraint != '':
        f.write('position_restraint = %s\n'%us_position_restraint)
    f.close()
    os.system('run_ess_us.py us.cntrl')
    os.chdir('../')

    if not os.path.exists('US/pmf.dat') or os.path.getsize('US/pmf.dat') <= 0:
        f_log = open('workflow.log', 'a')
        f_log.write('Error: unfinished US calculations\n')
        f_log.close()
        sys.exit()
    else:
        f = open('US/pmf.dat')
        lines = f.readlines()
        f.close()
        G_list = [float(l.strip().split()[1]) for l in lines if not 'inf' in l]
        G_list = np.array(G_list)
        G_max = np.max(G_list)
        idx = np.argwhere(G_list==G_max)
        idx = idx[0][0]
        G_min = np.min(G_list[:idx+1])
        G_max_err = float(lines[idx].strip().split()[-1])
        idx = np.argwhere(G_list==G_min)
        idx = idx[0][0]
        G_min_err = float(lines[idx].strip().split()[-1])
        G_err = (G_max_err**2 + G_min_err**2)**0.5
        
        fo = open('dE_dG.dat', 'w')
        fo.write('%6s %12.4f +/- %.4f\n'%('dGa', G_max-G_min, G_err))
        fo.close()

if if_smd == True:
    for i, smd_ligand_idx in enumerate(smd_ligand_list):
        if (if_us == True and stage <= 3 + i) or (if_us == False and stage <= 2 + i):
            f_log = open('workflow.log', 'a')
            f_log.write('Running ASMD for ligand %d\n'%(i+1))
            f_log.close()
            os.system('mkdir ASMD_%d'%(i+1))
            os.chdir('ASMD_%d'%(i+1))
            os.system('mkdir setup')
            f = open('asmd.cntrl', 'w')
            f.write('''top = ../auto_backmap_docking/complex.prmtop
cor = ../auto_backmap_docking/reactant_complex.ncrst
''')
            f.write('n_replica = %d\n'%num_smd_replica)
            f.write('n_stage = %d\n'%num_smd_stage)
            f.write('ppn = %d\n'%int(tpn/num_smd_replica))
            sel_syntax = ''
            for j, idx in enumerate(smd_ligand_idx):
                lig_file_name = ligand_file[idx].split()[0].strip()
                f.write('lig_file_%d = ../%s\n'%(j, lig_file_name))
                lfn = lig_file_name.split('/')[-1].split('.mol2')[0].strip()
                if lfn.endswith('_dummy'):
                    lfn = lfn.split('_dummy')[0]
                sel_syntax += lfn
                if j != len(smd_ligand_idx)-1:
                    sel_syntax += ','
            f.write('new_sys_mask = :%s | :1-%d\n'%(sel_syntax, chain_len))
            f.write('force_constant = %f\n'%pull_force_constant)
            f.write('dist_restraints = %s; %s\n'%(pull_lig_mask[i], pull_prot_mask[i]))
            f.write('temp = %f\n'%temp)
            f.write('pull_length = %f\n'%pull_length)
            f.write('steps_per_stage = %d\n'%smd_steps_per_stage)
            f.write('cpptraj_exec = %s\n'%cpptraj_exec)
            f.close()
            os.system('run_smd.py asmd.cntrl')
            os.chdir('../')
            
            if not os.path.exists('ASMD_%d/PMF.dat'%(i+1)) or os.path.getsize('ASMD_%d/PMF.dat'%(i+1)) <= 0:
                f_log = open('workflow.log', 'a')
                f_log.write('Error: unfinished asmd calculations\n')
                f_log.close()
                sys.exit()
            else:
                os.system('cp ASMD_%d/PMF.dat PMF_%d.dat'%(i+1, i+1))
                f = open('PMF_%d.dat'%(i+1))
                lines = f.readlines()
                f.close()
                ene_list = [float(l.strip().split()[-1]) for l in lines]
                min_ene = min(ene_list)
                last_ene = sum(ene_list[-5:])/5
                dG = min_ene - last_ene
                fo = open('dE_dG.dat', 'a')
                fo.write('%6s %12.4f\n'%('dG_%d'%(i+1), dG))
                fo.close()
    
end_time = time.time()
m, s = divmod(end_time-start_time, 60)
h, m = divmod(m, 60)
f_log = open('workflow.log', 'a')
f_log.write('Done. Elapsed time: %d:%d:%d\n'%(h, m, s))
f_log.close()
