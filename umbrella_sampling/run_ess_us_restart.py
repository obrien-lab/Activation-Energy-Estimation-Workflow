#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:55:37 2019

@author: yuj179
"""
import sys, os
import parmed as pmd
import mdtraj as mdt

################# Arguments #################
if len(sys.argv) != 3:
    print('Error: Wrong number of arguments')
    print('Usage: python run_ess_us_restart.py cntrl_file restart_idx')
    sys.exit()
cntrl_file = sys.argv[1]
restart_idx = int(sys.argv[2])

################# Functions #################
def parse_input(cntrl_file):
    top = ''
    reactant_cor = ''
    np = 20
    qm_mask = ''
    qm_charge = 0
    reaction_coordinate = []
    extra_restraint = []
    position_restraint = []
    rc_list = []
    temp = 310.0
    min_step = 500
    sampling_step = 10000
    force_constant = 250
    wham_exec = 'wham'
    wham_rc_interval = 0.05
    if_wham = 1
    only_wham = 0
    
    f = open(cntrl_file, 'r')
    for line in f.readlines():
        if line.startswith('top'):
            top = line.strip().split('=')[-1].strip()
        elif line.startswith('reactant_cor'):
            reactant_cor = line.strip().split('=')[-1].strip()
        elif line.startswith('np'):
            np = int(line.strip().split('=')[-1].strip())
        elif line.startswith('qm_mask'):
            qm_mask = '='.join(line.strip().split('=')[1:])
        elif line.startswith('qm_charge'):
            qm_charge = int(line.strip().split('=')[-1].strip())
        elif line.startswith('reaction_coordinate'):
            words = '='.join(line.strip().split('=')[1:])
            words = words.strip().split(';')
            for w in words:
                if w != '':
                    w1 = w.strip().split()
                    w1 = [ww.strip() for ww in w1]
                    reaction_coordinate.append(w1)
        elif line.startswith('extra_restraint'):
            words = '='.join(line.strip().split('=')[1:])
            words = words.strip().split(';')
            for w in words:
                if w.find('>') != -1:
                    w = w.strip().split('>')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    extra_restraint.append([w1[0], w1[1], 'r1=%.4f, r2=%.4f, r3=%.4f, r4=%.4f, rk2=%.4f, rk3=%.4f'%(
                            float(w2[0])-0.5, float(w2[0]), float(w2[0]), float(w2[0]), float(w2[1]), 0)])
                elif w.find('<') != -1:
                    w = w.strip().split('<')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    extra_restraint.append([w1[0], w1[1], 'r1=%.4f, r2=%.4f, r3=%.4f, r4=%.4f, rk2=%.4f, rk3=%.4f'%(
                            float(w2[0]), float(w2[0]), float(w2[0]), float(w2[0])+0.5, 0, float(w2[1]))])
                elif w.find('=') != -1:
                    w = w.strip().split('=')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    extra_restraint.append([w1[0], w1[1], 'r1=%.4f, r2=%.4f, r3=%.4f, r4=%.4f, rk2=%.4f, rk3=%.4f'%(
                            float(w2[0])-0.5, float(w2[0])-0.1, float(w2[0])+0.1, float(w2[0])+0.5, float(w2[1]), float(w2[1]))])
        elif line.startswith('position_restraint'):
            words = '='.join(line.strip().split('=')[1:])
            position_restraint = words.split()
        elif line.startswith('rc_list'):
            rc_list = line.strip().split('=')[-1].strip().split()
            rc_list = [float(r) for r in rc_list]
        elif line.startswith('temp'):
            temp = float(line.strip().split('=')[-1].strip())
        elif line.startswith('min_step'):
            min_step = int(line.strip().split('=')[-1].strip())
        elif line.startswith('sampling_step'):
            sampling_step = int(line.strip().split('=')[-1].strip())
        elif line.startswith('force_constant'):
            force_constant = float(line.strip().split('=')[-1].strip())
        elif line.startswith('wham_exec'):
            wham_exec = line.strip().split('=')[-1].strip()
        elif line.startswith('wham_rc_interval'):
            wham_rc_interval = float(line.strip().split('=')[-1].strip())
        elif line.startswith('if_wham'):
            if_wham = int(line.strip().split('=')[-1].strip())
        elif line.startswith('only_wham'):
            only_wham = int(line.strip().split('=')[-1].strip())
            
    tag_error = False
    if top == '':
        print('Error: No top assigned')
        tag_error = True
    if reactant_cor == '':
        print('Error: No reactant_cor assigned')
        tag_error = True
    if qm_mask == '':
        print('Error: No qm_mask assigned')
        tag_error = True
    if reaction_coordinate == []:
        print('Error: No reaction_coordinate assigned')
        tag_error = True
    if rc_list == []:
        print('Error: No rc_list assigned')
        tag_error = True
    else:
        if sorted(rc_list) != rc_list:
            print('Error: rc_list must be in ascending order')
            tag_error = True
        #if len(rc_list) <= 2:
            #print('Error: rc_list must have more than two elements')
            #tag_error = True
        
    if tag_error:
        sys.exit()
    else:
        return(top, reactant_cor, np, qm_mask, qm_charge, reaction_coordinate, 
               rc_list, temp, min_step, sampling_step, force_constant, only_wham, if_wham, 
               wham_exec, wham_rc_interval, extra_restraint, position_restraint)

def check_log_file():
    global log_file_name
    if not os.path.exists(log_file_name):
        return 1
    elif not os.path.getsize(log_file_name):
        return 1
    
    f = open(log_file_name)
    lines = f.readlines()
    f.close()
    stage = 0
    line_idx = 0
    for idx, line in enumerate(lines):
        if line.find('Error:') != -1:
            print('LOG file contains errors. Please delete "Error" lines to restart or check your system to debug.')
            sys.exit()
        elif line.startswith('us_'):
            stage += 1
            line_idx = idx
        elif line.startswith('wham'):
            stage += 1
            line_idx = idx
    if lines[-1].startswith('Done.'):
        print('All done.')
        sys.exit()
    
    f = open(log_file_name, 'w')
    for i in range(line_idx):
        f.write(lines[i])
    f.close()
    return stage

def us_jobs(idx, rc_restraints, extra_restraints, pos_restraints):
    global top, np, temp, qm_mask, qm_charge, sampling_step, restart_idx, log_file_name
    
    f_log = open('../'+log_file_name, 'a')
    f_log.write('us_%d\n'%(idx+1))
    f_log.close()
    
    if restart_idx == 1:
        previous_name_postfix = ''
    else:
        previous_name_postfix = '_r%d'%(restart_idx-1)
    name_postfix = '_r%d'%(restart_idx)
    
    crd = 'us_%d%s.ncrst'%(idx+1, previous_name_postfix)
    
    f = open('RST.dist', 'w')
    f.write('#\n# restraint on RC\n')
    f.write(' &rst\n  %s\n &end\n'%rc_restraints)
    if len(extra_restraints) > 0:
        for er in extra_restraints:
            f.write('#\n# extra restraints\n')
            f.write(' &rst\n  %s\n &end\n'%er)
    f.close()
    
    f = open('us_r.in', 'w')
    f.write('''Umbrella Sampling in QM/MM
 &cntrl
  imin=0, ntx=5, irest=1
  nstlim='''+str(int(sampling_step))+''', dt=0.001,
  ntf=2, ntc=2,
  temp0='''+str(temp)+''',
  ntpr=50, ntwx=100, ntwr=100,
  cut=10.0, ntb=2, ntc=2, ntf=2,
  ntp=1, ntt=3,
  gamma_ln=2.0,
  ig=-1,''')
    if len(pos_restraints) > 0:
        f.write('''
  ntr=1,
  restraintmask=\''''+ pos_restraints[0] +'''\',
  restraint_wt='''+ pos_restraints[1] +',')
    f.write('''
  ifqnt=1,
  nmropt=1,
 /
 &qmmm
  qmmask="'''+qm_mask+'''"
  qmcharge='''+str(qm_charge)+''',
  qm_theory='DFTB3',
  dftb_telec=100.0,
  qmshake=0,
  qm_ewald=1, qm_pme=1,
 /
&wt type='DUMPFREQ', istep1=10,/
&wt type='END' /
DISANG=RST.dist
DUMPAVE=dist_'''+str(idx+1)+name_postfix+'''.dat
''')
    
    f.close()
    if pos_restraints == []:
        os.system('mpirun -np '+str(int(np))+' sander.MPI -O -i us_r.in -p ../'+top+' -c '+crd+' -x us_'+str(idx+1)+name_postfix+'.nc -r us_'+str(idx+1)+name_postfix+'.ncrst -o us_'+str(idx+1)+name_postfix+'.out -inf mdinfo.log')
    else:
        os.system('mpirun -np '+str(int(np))+' sander.MPI -O -i us_r.in -p ../'+top+' -c '+crd+' -x us_'+str(idx+1)+name_postfix+'.nc -r us_'+str(idx+1)+name_postfix+'.ncrst -o us_'+str(idx+1)+name_postfix+'.out -ref '+crd+' -inf mdinfo.log')

def get_PMF(rc_list, force_constant, skip_steps):
    global wham_exec, temp, wham_rc_interval, restart_idx, log_file_name
    tol = 1e-4
    
    f_log = open(log_file_name, 'a')
    f_log.write('wham\n')
    f_log.close()
    
    os.system('mkdir wham')
    os.chdir('wham')
    f = open('metadata.dat', 'w')
    for i in range(len(rc_list)):
        ff = open('../us/dist_%d.dat'%(i+1), 'r')
        lines = ff.readlines()[1:]
        ff.close()
        for ii in range(restart_idx):
            ff = open('../us/dist_%d_r%d.dat'%(i+1, ii+1), 'r')
            lines += ff.readlines()[1:]
            ff.close()
        ff = open('%d.dat'%(i+1), 'w')
        for li, line in enumerate(lines):
            frame = 10*li
            dist = float(line.strip().split()[1])
            if frame >= skip_steps:
                ff.write('%-8d%-10.3f\n'%(frame, dist))
        ff.close()
        f.write('%d.dat %.4f %.4f\n'%(i+1, rc_list[i], 2*force_constant))
    f.close()
    
    nbins = int((rc_list[-1] - rc_list[0] + 0.5) / wham_rc_interval)
    os.system('%s %.4f %.4f %d %f %.2f 0 metadata.dat pmf.dat 100 1'%(wham_exec, rc_list[0]-0.25, rc_list[-1]+0.25, nbins, tol, temp))
    os.chdir('../')
        
################# MAIN ####################
log_file_name = 'us_r%d.log'%(restart_idx)
(top, reactant_cor, np, qm_mask, qm_charge, reaction_coordinate, 
 rc_list, temp, min_step, sampling_step, force_constant, only_wham, if_wham, wham_exec,
 wham_rc_interval, extra_restraint, position_restraint) = parse_input(cntrl_file)

n_windows = len(rc_list)
if n_windows < 5:
    print('Warning: n_windows is less than 5. You may have a sampling issue.')

if only_wham == 0:
    stage = check_log_file()
    if stage == 0:
        print('Error: Cannot find checkpoint. Please check the last line of your LOG file')
        sys.exit()
else:
    stage = 9999

struct = pmd.load_file(top)
atm_idx = []
w_list = []
for d in reaction_coordinate:
    for i in range(2):
        sel = d[i].split(':')[-1].split('@')
        for res in struct.residues:
            if res.name == sel[0] or str(res.number+1) == sel[0]:
                for atm in res.atoms:
                    if atm.name == sel[1]:
                        atm_idx.append(atm.idx+1)
                        break
                break
    w_list.append(d[2])
restraints_head = 'iat='
for ai in atm_idx:
    restraints_head += str(ai)+','
restraints_head += ' rstwt='
for w in w_list:
    restraints_head += str(w)+','

extra_restraints = []    
for d in extra_restraint:
    atm_idx = []
    for i in range(2):
        sel = d[i].split(':')[-1].split('@')
        for res in struct.residues:
            if res.name == sel[0] or str(res.number+1) == sel[0]:
                for atm in res.atoms:
                    if atm.name == sel[1]:
                        atm_idx.append(atm.idx+1)
                        break
                break
    res_str = 'iat='
    for ai in atm_idx:
        res_str += str(ai)+','
    res_str += ' '+d[2]
    extra_restraints.append(res_str)

traj = mdt.load('ess/ess_1.ncrst', top=top)
for i in range(1,n_windows):
    traj += mdt.load('ess/ess_%d.ncrst'%(i+1), top=top)
CA_idx = traj.topology.select("name CA")
traj = traj.superpose(traj, atom_indices=CA_idx)
traj.save('RC.nc', force_overwrite=True)

for i in range(n_windows):
    if stage <= i+1:
        rc = rc_list[i]
        restraints = '%s r1=%.4f, r2=%.4f, r3=%.4f, r4=%.4f, rk2=%.4f, rk3=%.4f'%(
                     restraints_head, rc-999, rc, rc, rc+999, force_constant, force_constant)
        os.chdir('us')
        us_jobs(i, restraints, extra_restraints, position_restraint)
        os.chdir('../')        

if (stage <= n_windows+1 and if_wham == 1) or only_wham == 1:
    get_PMF(rc_list, force_constant, int((restart_idx+1)*sampling_step/2))
    f = open('wham/pmf.dat')
    ff = open('pmf.dat', 'w')
    tag = 0
    pmf_list = []
    min_pmf = 1e6
    for line in f:
        if line.startswith('#Coor'):
            continue
        elif line.startswith('#Window'):
            break
        elif 'nan' in line:
            continue
        
        words = line.strip().split()
        data = [float(w) for w in words]
        pmf_list.append(data)
        if min_pmf > data[1]:
            min_pmf = data[1]
    
    for data in pmf_list:
        ff.write('%12.6f %12.6f %12.6f\n'%(data[0], data[1]-min_pmf, 1.91*data[2]))
            
    f.close()
    ff.close()
    os.system('rm -rf wham')

f_log = open(log_file_name, 'a')
f_log.write('Done.\n')
f_log.close()
