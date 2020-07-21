#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:55:37 2019

@author: yuj179
"""
import sys, os, numpy
import mdtraj as mdt

################# Arguments #################
if len(sys.argv) != 2:
    print('Error: Wrong number of arguments')
    print('Usage: python run_smd.py cntrl_file')
    sys.exit()
cntrl_file = sys.argv[1]

################# Functions #################
def parse_input(cntrl_file):
    top = ''
    cor = ''
    n_replica = 10
    n_stage = 4
    ppn = 1
    lig_file = []
    new_sys_mask = ''
    force_constant = 100.0
    dist_restraints = []
    temp = 310
    pull_length = 20
    steps_per_stage = 125000
    cpptraj_exec = 'cpptraj'
    
    f = open(cntrl_file, 'r')
    for line in f.readlines():
        if line.startswith('top'):
            top = line.strip().split('=')[-1].strip()
        elif line.startswith('cor'):
            cor = line.strip().split('=')[-1].strip()
        elif line.startswith('lig_file'):
            lig_file.append(line.strip().split('=')[-1].strip())
        elif line.startswith('n_replica'):
            n_replica = int(line.strip().split('=')[-1].strip())
        elif line.startswith('n_stage'):
            n_stage = int(line.strip().split('=')[-1].strip())
        elif line.startswith('ppn'):
            ppn = int(line.strip().split('=')[-1].strip())
        elif line.startswith('new_sys_mask'):
            new_sys_mask = '='.join(line.strip().split('=')[1:])
            new_sys_mask = ''.join(new_sys_mask.strip().split())
        elif line.startswith('force_constant'):
            force_constant = float(line.strip().split('=')[-1].strip())
        elif line.startswith('dist_restraints'):
            words = line.strip().split('=')[-1].strip().split(';')
            for w in words:
                if w != '':
                    dist_restraints.append(w.strip())
        elif line.startswith('temp'):
            temp = float(line.strip().split('=')[-1].strip())
        elif line.startswith('pull_length'):
            pull_length = float(line.strip().split('=')[-1].strip())
        elif line.startswith('steps_per_stage'):
            steps_per_stage = int(line.strip().split('=')[-1].strip())
        elif line.startswith('cpptraj_exec'):
            cpptraj_exec = line.strip().split('=')[-1].strip()
            
    tag_error = False
    if top == '':
        print('Error: No top assigned')
        tag_error = True
    if cor == '':
        print('Error: No cor file assigned')
        tag_error = True
    if lig_file == []:
        print('Error: No lig_file assigned')
        tag_error = True
    if new_sys_mask == '':
        print('Error: No new_sys_mask assigned')
        tag_error = True
    if len(dist_restraints) == 0:
        print('Error: No dist_restraints assigned')
        tag_error = True
        
    if tag_error:
        sys.exit()
    else:
        return(top, cor, lig_file, n_replica, n_stage, ppn, new_sys_mask, force_constant, 
               dist_restraints, temp, pull_length, steps_per_stage, cpptraj_exec)

def check_log_file():
    if not os.path.exists('smd.log'):
        return 1
    elif not os.path.getsize('smd.log'):
        return 1
    
    f = open('smd.log')
    lines = f.readlines()
    f.close()
    stage = 0
    line_idx = 0
    for idx, line in enumerate(lines):
        if line.find('Error:') != -1:
            print('LOG file contains errors. Please delete "Error" lines to restart or check your system to debug.')
            sys.exit()
        elif line.startswith('min system'):
            stage += 1
            line_idx = idx
        elif line.startswith('heat system'):
            stage += 1
            line_idx = idx
        elif line.startswith('equil_1 system'):
            stage += 1
            line_idx = idx
        elif line.startswith('equil_2 system'):
            stage += 1
            line_idx = idx
        elif line.startswith('smd for stage'):
            stage += 1
            line_idx = idx
    if lines[-1].startswith('Done.'):
        print('All done.')
        sys.exit()
    
    f = open('smd.log', 'w')
    for i in range(line_idx):
        f.write(lines[i])
    f.close()
    return stage

def smd_jobs(top, start_crd, start_len, end_len, np, stage):
    global n_replica, force_constant 
    global dist_restraints, temp, steps_per_stage
    global sel_1, sel_2

    f_log = open('../smd.log', 'a')
    f_log.write('smd for stage %d\n'%(stage+1))
    f_log.close()
    
    # smd
    for i in range(n_replica):
        f = open('asmd_%03d.in'%(i+1), 'w')
        f.write('''ASMD simulation
 &cntrl
   imin = 0, nstlim = '''+str(int(steps_per_stage))+''', dt = 0.002,
   ntx = 1, temp0 = '''+str(temp)+''',
   ntt = 3, gamma_ln=2.0,
   ntc = 2, ntf = 2, ntb =1,
   ntwx =  1000, ntwr = 1000, ntpr = 1000,
   cut = 10.0, ig=-1,
   irest = 0, jar=1, 
 /
 &wt type='DUMPFREQ', istep1=1000 /
 &wt type='END'   /
DISANG=dist.RST.dat.'''+('%03d'%(stage+1))+'''
DUMPAVE=asmd_'''+('%03d'%(i+1))+'''.work.dat.'''+('%03d'%(stage+1))+'''
LISTIN=POUT
LISTOUT=POUT
''')
        f.close()
    
    f = open('dist.RST.dat.%03d'%(stage+1), 'w')
    f.write(''' &rst
        iat=-1,-1,
        r2='''+str(start_len)+''',
        r2a='''+str(end_len)+''',
        rk2='''+str(force_constant)+''',
        igr1='''+sel_1+''',
        igr2='''+sel_2+''',
 &end
 ''')
    f.close()

    f = open('groupfile', 'w')
    for i in range(n_replica):
        f.write('-O -p ../%s -c ../%s -i asmd_%03d.in -x asmd_%03d.nc -o asmd_%03d.out -inf asmd_%03d.info -r asmd_%03d.ncrst\n'%(
                top, start_crd, i+1, i+1, i+1, i+1, i+1))
    f.close()
    
    os.system('mpirun -np '+str(np)+' pmemd.MPI -ng '+str(n_replica)+' -groupfile groupfile')

def find_closest_replica(stage, temp):
    file_list = os.popen('ls stage_%d/asmd_*.work.dat.*'%(stage+1)).readlines()
    file_list = [f.strip() for f in file_list]
    work_list = []
    dist_list = []
    for i, file_name in enumerate(file_list):
        dist = []
        work = []
        f = open(file_name)
        for line in f:
            dist.append(float(line.strip().split()[0]))
            work.append(float(line.strip().split()[-1]))
        f.close()
        
        out_file_name = 'stage_%d/asmd_%03d.out'%(stage+1, i+1)
        f = open(out_file_name)
        lines = f.read()
        f.close()
        if lines.find('A V E R A G E S   O V E R') == -1: 
            f = open('smd.log', 'a')
            f.write('Error: unfinished simulation for %s\n'%file_name)
            f.close()
            sys.exit()
        else:
            dist_list = dist
            work_list.append(work)
    # calculate Jarzynski average
    work_list = numpy.array(work_list)
    beta = 1 / (temp * 1.98722e-3)
    JarAvg_list = -numpy.log(numpy.sum(numpy.exp(-beta*work_list), axis=0)/len(file_list))/beta
    
    dd = numpy.inf
    idx = 0
    for i in range(len(work_list)):
        d = abs(work_list[i,-1] - JarAvg_list[-1])
        if d < dd:
            dd = d
            idx = i
    
    f = open('JarAvg.dat.%d'%(stage+1), 'w')
    for i in range(len(JarAvg_list)):
        f.write('%.4f %.4f\n'%(dist_list[i], JarAvg_list[i]))
    f.close()
    return idx+1
        
    
        
################# MAIN ####################
(top, cor, lig_file, n_replica, n_stage, ppn, new_sys_mask, force_constant, 
 dist_restraints, temp, pull_length, steps_per_stage, cpptraj_exec) = parse_input(cntrl_file)

stage = check_log_file()

if stage == 0:
    print('Error: Cannot find checkpoint. Please check the last line of your LOG file')
    sys.exit()

np = int(ppn*n_replica)

if stage <= 1:
    f = open('cpptraj.in', 'w')
    f.write('trajin %s\n'%cor)
    f.write('strip "!(%s)"\n'%new_sys_mask)
    f.write('trajout setup/new.pdb\n')
    f.close()
    os.system('%s -p %s -i cpptraj.in > /dev/null 2>&1'%(cpptraj_exec, top))
    
    for lf in lig_file:
        lig_name = lf.split('/')[-1].split('.mol2')[0]
        if lig_name.endswith('_dummy'):
            lfp = lf.split('.mol2')[0]
            os.system('cp '+lfp+'.prepi setup/'+lig_name+'.prepi')
            os.system('cp '+lfp+'.frcmod setup/'+lig_name+'.frcmod')
        else:
            os.system('antechamber -i '+lf+' -fi mol2 -o setup/'+lig_name+'.prepi -fo prepi -pf y -dr n > /dev/null')
            os.system('parmchk2 -i setup/'+lig_name+'.prepi -f prepi -o setup/'+lig_name+'.frcmod > /dev/null')
    f = open('leap.in', 'w')
    f.write('''source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2\n''')
    for lf in lig_file:
        lig_name = lf.split('/')[-1].split('.mol2')[0]
        f.write('loadamberprep setup/'+lig_name+'.prepi\n')
        f.write('loadamberparams setup/'+lig_name+'.frcmod\n')
    f.write('''a = loadpdb setup/new.pdb
solvatebox a TIP3PBOX 10.0
addions a Na+ 0.0
saveamberparm a setup/new.prmtop setup/new.inpcrd
quit\n''')
    f.close()
    os.system('tleap -s -f leap.in > /dev/null')

# get restraints selection
amber_top = mdt.load_prmtop('setup/new.prmtop')
sel_1 = amber_top.select(dist_restraints[0])
sel_1 = [str(s+1) for s in sel_1]
sel_1 = ','.join(sel_1)
sel_2 = amber_top.select(dist_restraints[1])
sel_2 = [str(s+1) for s in sel_2]
sel_2 = ','.join(sel_2)

# Calculate initial distance
f = open('cpptraj.in', 'w')
f.write('trajin setup/new.inpcrd\n') 
f.write('distance @%s @%s out dist.dat\n'%(sel_1, sel_2))
f.close()
os.system('%s -p setup/new.prmtop -i cpptraj.in > /dev/null 2>&1'%cpptraj_exec)
f = open('dist.dat')
lines = f.readlines()
f.close()
start_dist = float(lines[-1].strip().split()[-1])

#min
if stage <= 1:
    f_log = open('smd.log', 'w')
    f_log.write('min system\n')
    f_log.close()
    
    os.system('mkdir prepare')
    os.chdir('prepare')
    
    f = open('RST.dist', 'w')
    f.write(' &rst\n')
    f.write('  iat=-1,-1,\n  r1=%.4f, r2=%.4f, r3=%.4f, r4=%.4f, rk2=0, rk3=10,\n  igr1=%s, igr2=%s,\n'%(
            start_dist, start_dist, start_dist, start_dist+5, sel_1, sel_2))
    f.write(' &end\n')
    f.close()
    
    f = open('min.in', 'w')
    f.write('''minimization
 &cntrl
  imin   = 1,
  maxcyc = 2000,
  ncyc   = 1000,
  cut=10.0, 
  ntb=1,
  nmropt=1,
 /
 &wt type='END' /
DISANG=RST.dist
''')
    f.close()
    os.system('mpirun -np %d pmemd.MPI -O -i min.in -p ../setup/new.prmtop -c ../setup/new.inpcrd -r min.ncrst -o min.out -inf mdinfo.log'%np)
    os.chdir('../')

# heat
if stage <= 2:
    f_log = open('smd.log', 'a')
    f_log.write('heat system\n')
    f_log.close()
    
    os.chdir('prepare')
    f = open('heat.in', 'w')
    f.write('''heating in TIP3P water box
 &cntrl
  imin=0, ntx=1, irest=0,
  nstlim=50000, dt=0.002,
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
&wt type='TEMP0', istep1=0, istep2=40000, value1=0.0, value2='''+str(temp)+''' /
&wt type='TEMP0', istep1=40000, istep2=50000, value1='''+str(temp)+''', value2='''+str(temp)+''' /
&wt type='END' /
DISANG=RST.dist
''')
    f.close()
    os.system('mpirun -np %d pmemd.MPI -O -i heat.in -p ../setup/new.prmtop -c min.ncrst -r heat.ncrst -o heat.out -x heat.nc -ref min.ncrst -inf mdinfo.log'%np)
    os.chdir('../')
    
# equil_1
if stage <= 3:
    f_log = open('smd.log', 'a')
    f_log.write('equil_1 system\n')
    f_log.close()
    
    os.chdir('prepare')
    f = open('equil_1.in', 'w')
    f.write('''equilibration in TIP3P water box
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim=250000, dt=0.002,
  ntf=2, ntc=2,
  temp0='''+str(temp)+''',
  ntpr=1000, ntwx=5000,
  cut=10.0, ntb=2, ntc=2, ntf=2,
  ntp=1, ntt=3,
  gamma_ln=2.0,
  ig=-1,
  nmropt=1,
 /
 &wt type='END' /
DISANG=RST.dist
''')
    f.close()
    os.system('mpirun -np %d pmemd.MPI -O -i equil_1.in -p ../setup/new.prmtop -c heat.ncrst -r equil_1.ncrst -o equil_1.out -x equil_1.nc -inf mdinfo.log'%np)
    os.chdir('../')

# equil_2
if stage <= 4:
    f_log = open('smd.log', 'a')
    f_log.write('equil_2 system\n')
    f_log.close()
    
    os.chdir('prepare')
    f = open('equil_2.in', 'w')
    f.write('''equilibration in TIP3P water box
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim=250000, dt=0.002,
  ntf=2, ntc=2,
  temp0='''+str(temp)+''',
  ntpr=1000, ntwx=5000,
  cut=10.0, ntb=1, ntc=2, ntf=2,
  ntp=0, ntt=3,
  gamma_ln=2.0,
  ig=-1,
  nmropt=1,
 /
 &wt type='END' /
DISANG=RST.dist
''')
    f.close()
    os.system('mpirun -np %d pmemd.MPI -O -i equil_2.in -p ../setup/new.prmtop -c equil_1.ncrst -r equil_2.ncrst -o equil_2.out -x equil_2.nc -inf mdinfo.log'%np)
    os.chdir('../')

start_crd = 'prepare/equil_2.ncrst'    
for i in range(n_stage):
    if stage <= i+5:
        os.system('mkdir stage_%d'%(i+1))
        os.chdir('stage_%d'%(i+1))
        start_len = 9*start_dist/10 + pull_length / n_stage * i
        end_len = 9*start_dist/10 + pull_length / n_stage * (i+1)
        smd_jobs('setup/new.prmtop', start_crd, start_len, end_len, np, i)
        os.chdir('../')
    min_idx = find_closest_replica(i, temp)
    start_crd = 'stage_%d/asmd_%03d.ncrst'%(i+1, min_idx)
    f_log = open('smd.log', 'a')
    f_log.write('closest: %d\n'%min_idx)
    f_log.close()
        
fo = open('PMF.dat', 'w')
pmf_init = 0
for i in range(n_stage):
    f = open('JarAvg.dat.%d'%(i+1))
    lines = f.readlines()
    f.close()
    for line in lines:
        rc = float(line.strip().split()[0])
        pmf = float(line.strip().split()[1]) + pmf_init
        fo.write('%10.4f %10.4f\n'%(rc, pmf))
    pmf_init = pmf
fo.close()

f_log = open('smd.log', 'a')
f_log.write('Done.\n')
f_log.close()
