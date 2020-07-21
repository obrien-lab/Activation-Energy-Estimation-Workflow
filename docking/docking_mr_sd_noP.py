#!/usr/bin/env python3

import os, sys
import parmed as pmd
import mdtraj as mdt

################# Arguments #################
if len(sys.argv) != 2:
    print('Error: Wrong number of arguments')
    print('Usage: python docking.py cntrl_file')
    sys.exit()
cntrl_file = sys.argv[1]

################# Functions #################
def parse_input(cntrl_file):
    adt_home = ''
    protein_pdb_file = ''
    ligand_file_list = []
    box_center_list = []
    box_dimension_list = []
    temp_seq = ''
    ppn = 1
    refine_reactant_restraints = []
    QM_reactant_restraints = []
    qm_mask = ''
    qm_charge = 0
    temp = 310.0
    
    f = open(cntrl_file, 'r')
    for line in f.readlines():
        if line.startswith('adt_home'):
            adt_home = line.strip().split('=')[-1].strip()
        elif line.startswith('protein_pdb_file'):
            protein_pdb_file = line.strip().split('=')[-1].strip()
        elif line.startswith('ligand_file'):
            ligand_file_list.append(line.strip().split('=')[-1].strip().split())
        elif line.startswith('box_center'):
            words = line.strip().split('=')[-1].strip().split()
            box_center_list.append([float(c) for c in words])
        elif line.startswith('box_dimension'):
            words = line.strip().split('=')[-1].strip().split()
            box_dimension_list.append([float(c) for c in words])
        elif line.startswith('temp_seq_file'):
            temp_seq = line.strip().split('=')[-1].strip()
        elif line.startswith('ppn'):
            ppn = int(line.strip().split('=')[-1].strip())
        elif line.startswith('refine_reactant_restraints'):
            words = '='.join(line.strip().split('=')[1:])
            words = words.strip().split(';')
            for w in words:
                if w.find('>') != -1:
                    w = w.strip().split('>')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    refine_reactant_restraints.append([w1[0], w1[1], 'r1=%.2f, r2=%.2f, r3=%.2f, r4=%.2f, rk2=%.2f, rk3=%.2f'%(
                            float(w2[0])-0.5, float(w2[0]), float(w2[0]), float(w2[0]), float(w2[1]), 0)])
                elif w.find('<') != -1:
                    w = w.strip().split('<')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    refine_reactant_restraints.append([w1[0], w1[1], 'r1=%.2f, r2=%.2f, r3=%.2f, r4=%.2f, rk2=%.2f, rk3=%.2f'%(
                            float(w2[0]), float(w2[0]), float(w2[0]), float(w2[0])+0.5, 0, float(w2[1]))])
                elif w.find('=') != -1:
                    w = w.strip().split('=')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    refine_reactant_restraints.append([w1[0], w1[1], 'r1=%.2f, r2=%.2f, r3=%.2f, r4=%.2f, rk2=%.2f, rk3=%.2f'%(
                            float(w2[0])-0.5, float(w2[0])-0.1, float(w2[0])+0.1, float(w2[0])+0.5, float(w2[1]), float(w2[1]))])
        elif line.startswith('QM_reactant_restraints'):
            words = line.strip().split('=')[-1].strip().split(';')
            for w in words:
                if w.find('>') != -1:
                    w = w.strip().split('>')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    QM_reactant_restraints.append([w1[0], w1[1], float(w2[0])-0.5, float(w2[0]), float(w2[0]), float(w2[0]), float(w2[1]), 0])
                elif w.find('<') != -1:
                    w = w.strip().split('<')
                    w1 = w[0].strip().split()
                    w2 = w[1].strip().split()
                    QM_reactant_restraints.append([w1[0], w1[1], float(w2[0]), float(w2[0]), float(w2[0]), float(w2[0])+0.5, 0, float(w2[1])])
        elif line.startswith('qm_mask'):
            qm_mask = line.strip().split('=')[-1].strip()
        elif line.startswith('qm_charge'):
            qm_charge = int(line.strip().split('=')[-1].strip())
        elif line.startswith('temperature'):
            temp = float(line.strip().split('=')[-1].strip())
            
    tag_error = False
    if adt_home == '':
        print('Error: No adt_home assigned')
        tag_error = True
    if protein_pdb_file == '':
        print('Error: No protein_pdb_file assigned')
        tag_error = True
    if ligand_file_list == []:
        print('Error: No ligand_file assigned')
        tag_error = True
    if box_center_list == []:
        print('Error: No box_center assigned')
        tag_error = True
    if box_dimension_list == []:
        print('Error: No box_dimension assigned')
        tag_error = True
    if len(ligand_file_list) != len(box_center_list):
        print('Error: The numbers of ligand_file (%d) and box_center (%d) mismatched'%(len(ligand_file_list), len(box_center_list)))
        tag_error = True
    if len(ligand_file_list) != len(box_dimension_list):
        print('Error: The numbers of ligand_file (%d) and box_dimension (%d) mismatched'%(len(ligand_file_list), len(box_dimension_list)))
        tag_error = True
    if len(box_center_list) != len(box_dimension_list):
        print('Error: The numbers of box_center (%d) and box_dimension (%d) mismatched'%(len(box_center_list), len(box_dimension_list)))
        tag_error = True
    if qm_mask == '':
        print('Error: No qm_mask assigned')
        tag_error = True
        
    if tag_error:
        sys.exit()
    else:
        return(adt_home, protein_pdb_file, ligand_file_list, box_center_list, box_dimension_list, 
               temp_seq, ppn, refine_reactant_restraints, QM_reactant_restraints, qm_mask, qm_charge, temp)

def check_log_file():
    if not os.path.exists('docking.log'):
        return 1
    elif not os.path.getsize('docking.log'):
        return 1
    
    f = open('docking.log')
    lines = f.readlines()
    f.close()
    stage = 0
    line_idx = 0
    for idx, line in enumerate(lines):
        if line.find('Error:') != -1:
            print('LOG file contains errors. Please delete "Error" lines to restart or check your system to debug.')
            sys.exit()
        elif line.startswith('Running docking for reactant'):
            stage += 1
            line_idx = idx
        elif line.startswith('Running MD docking refinement of reactants'):
            stage += 1
            line_idx = idx
        elif line.startswith('Running QMMM docking refinement of reactants'):
            stage += 1
            line_idx = idx
    if lines[-1].startswith('Done.'):
        print('All done.')
        sys.exit()
    
    f = open('docking.log', 'w')
    for i in range(line_idx):
        f.write(lines[i])
    f.close()
    return stage
    
def refine_docking_MD(complex_pdb, ligand_file_list, restraint_pairs, steps_1, steps_2, steps_3, steps_4):
    global ppn, temp
    f = open('leap.in', 'w')
    f.write('''source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2\n''')
    lig_name_list = []
    for lig in ligand_file_list:
        lig_name = lig.split('/')[-1].split('.mol2')[0]
        if lig_name == 'WAT':
            continue
        if lig_name.endswith('_dummy'):
            ln =  lig_name.split('_dummy')[0]
            if not ln in lig_name_list:
                lig_name_list.append(ln)
        elif not lig_name in lig_name_list:
            lig_name_list.append(lig_name)
        f.write('loadamberprep '+lig_name+'.prepi\n')
        f.write('loadamberparams '+lig_name+'.frcmod\n')
    f.write('''a = loadpdb '''+complex_pdb+'''
savepdb a complex.pdb
solvatebox a TIP3PBOX 10.0
addions a Na+ 0.0
saveamberparm a complex.prmtop complex.inpcrd
quit\n''')
    f.close()
    os.system('tleap -s -f leap.in > /dev/null')
    
    if restraint_pairs != []:
        struct = pmd.load_file('complex.prmtop')
        f = open('RST.dist', 'w')
        for rest in restraint_pairs:
            atm_idx = []
            for i in range(2):
                sel = rest[i].split(':')[-1].split('@')
                for res in struct.residues:
                    if res.name == sel[0] or str(res.number+1) == sel[0]:
                        for atm in res.atoms:
                            if atm.name == sel[1]:
                                atm_idx.append(atm.idx+1)
                                break
                        break
            f.write('#\n# '+rest[0]+' '+rest[1]+'\n')
            f.write(''' &rst
  iat= '''+str(atm_idx[0])+''', '''+str(atm_idx[1])+''', '''+rest[2]+''',
 &end\n''')
        f.close()
        
    f = open('min1.in', 'w')
    f.write('''minimization in MM with restraints
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps_1)+''',
  ncyc   = '''+str(int(steps_1/2))+''',
  cut=10.0, 
  ntb=1, 
  ntc=1, 
  ntf=1,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=100.0,
 /
''')
    f.close()
    os.system('mpirun -np '+str(int(ppn))+' pmemd.MPI -O -i min1.in -p complex.prmtop -c complex.inpcrd -r min1.ncrst -o min1.out -ref complex.inpcrd -inf mdinfo.log')
    
    f = open('min2.in', 'w')
    f.write('''minimization in MM without restraints
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps_2)+''',
  ncyc   = '''+str(int(steps_2/2))+''',
  cut=10.0, 
  ntb=1, 
  ntc=1, 
  ntf=1,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=5.0,
 /
''')
    f.close()
    os.system('mpirun -np '+str(int(ppn))+' pmemd.MPI -O -i min2.in -p complex.prmtop -c min1.ncrst -r min2.ncrst -o min2.out -ref min1.ncrst -inf mdinfo.log')
    
    res_lig_str = '(@CA) | (:'
    for ln in lig_name_list:
        res_lig_str += ln+','
    res_lig_str = res_lig_str[:-1] + ')'
    f = open('heat.in', 'w')
    f.write('''heating in TIP3P water box
 &cntrl
  imin=0, ntx=1, irest=0,
  nstlim='''+str(int(steps_3))+''', dt=0.002,
  ntf=2, ntc=2,
  tempi=0.0, temp0='''+str(temp)+''',
  ntpr=100, ntwx=500,
  cut=10.0, ntb=1, ntc=2, ntf=2,
  ntp=0, ntt=3,
  gamma_ln=2.0,
  ntr=1,
  restraint_wt=3.0,
  restraintmask=\''''+res_lig_str+'''\',
  nmropt=1,
  ig=-1,
 /
&wt type='TEMP0', istep1=0, istep2='''+str(int(steps_3*9/10))+''', value1=0.0, value2='''+str(temp)+''' /
&wt type='TEMP0', istep1='''+str(int(steps_3*9/10+1))+''', istep2='''+str(int(steps_3))+''', value1='''+str(temp)+''', value2='''+str(temp)+''' /
&wt type='END' /
''')
    f.close()
    os.system('mpirun -np '+str(int(ppn))+' pmemd.MPI -O -i heat.in -p complex.prmtop -c min2.ncrst -r heat.ncrst -o heat.out -x heat.nc -ref min2.ncrst -inf mdinfo.log')
    
    f = open('equil.in', 'w')
    f.write('''equilibration in TIP3P water box with increasing restraints
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim='''+str(int(steps_4))+''', dt=0.002,
  ntf=2, ntc=2,
  temp0='''+str(temp)+''',
  ntpr=1000, ntwx=5000,
  cut=10.0, ntb=2, ntc=2, ntf=2,
  ntp=1, ntt=3,
  gamma_ln=2.0,
  ig=-1,
  ntr=1,
  restraintmask='@CA',
  restraint_wt=1.0,''')
    if restraint_pairs != []:
        f.write('''
  nmropt=1,
 /
 &wt type='REST', istep1=0, istep2='''+str(int(steps_4/2))+''', value1=0.1, value2=1.0,/
 &wt type='REST', istep1='''+str(int(steps_4/2)+1)+''', istep2='''+str(steps_4)+''', value1=1.0, value2=1.0,/
 &wt type='END' /
DISANG=RST.dist
''')
    else:
        f.write('''
 /
''')
    f.close()
    os.system('mpirun -np '+str(int(ppn))+' pmemd.MPI -O -i equil.in -p complex.prmtop -c heat.ncrst -r equil.ncrst -o equil.out -x equil.nc -ref heat.ncrst -inf mdinfo.log')
    
    f = open('equil_min.in', 'w')
    f.write('''minimization in MM without restraints
 &cntrl
  imin   = 1,
  maxcyc = '''+str(steps_2)+''',
  ncyc   = '''+str(int(steps_2/2))+''',
  cut=10.0, 
  ntb=1, 
  ntc=1, 
  ntf=1,''')
    if restraint_pairs != []:
        f.write('''
  nmropt=1,
 /
 &wt type='END' /
DISANG=RST.dist
''')
    else:
        f.write('''
 /
''')
    f.close()
    os.system('mpirun -np '+str(int(ppn))+' pmemd.MPI -O -i equil_min.in -p complex.prmtop -c equil.ncrst -r equil_min.ncrst -o equil_min.out -inf mdinfo.log')

def refine_docking_QMMM(top, crd, restraint_pairs, qm_mask, qm_charge, steps):
    global ppn
    
    if restraint_pairs != []:
        struct = pmd.load_file(top)
        f = open('RST.dist', 'w')
        for rest in restraint_pairs:
            atm_idx = []
            for i in range(2):
                sel = rest[i].split(':')[-1].split('@')
                for res in struct.residues:
                    if res.name == sel[0] or str(res.number+1) == sel[0]:
                        for atm in res.atoms:
                            if atm.name == sel[1]:
                                atm_idx.append(atm.idx+1)
                                break
                        break
            f.write('#\n# '+rest[0]+' '+rest[1]+'\n')
            f.write(' &rst\n')
            f.write('  iat=%d, %d, r1=%.2f, r2=%.2f, r3=%.2f, r4=%.2f, rk2=%.2f, rk3=%.2f\n'%(atm_idx[0], atm_idx[1], 
                    rest[2], rest[3], rest[4], rest[5], rest[6], rest[7]))
            f.write(' &end\n')
        f.close()
    
    f = open('min3.in', 'w')
    f.write('''minimization in QM/MM with increasing restraints
 &cntrl
  imin   = 1,
  maxcyc = '''+str(int(steps))+''',
  ncyc   = '''+str(int(steps/2))+''',
  cut=10.0, 
  ntb=1, 
  ntc=1, 
  ntf=1,
  ifqnt=1,''')
    if restraint_pairs != []:
        f.write('''
  nmropt=1,''')
    f.write('''
 /
 &qmmm
  qmmask=\''''+qm_mask+'''\'
  qmcharge='''+str(qm_charge)+''',
  qm_theory='DFTB3',
  dftb_telec=100.0,
  qmshake=0,
  qm_ewald=1, qm_pme=1,
 /
''')
    if restraint_pairs != []:
        f.write('''&wt type='REST', istep1=0, istep2='''+str(int(steps/2))+''', value1=0.1, value2=1.0,/
&wt type='REST', istep1='''+str(int(steps/2)+1)+''', istep2='''+str(int(steps))+''', value1=1.0, value2=1.0,/
&wt type='END' /
DISANG=RST.dist
''')
    
    f.close()
    os.system('mpirun -np '+str(int(ppn))+' sander.MPI -O -i min3.in -p '+top+' -c '+crd+' -r min3.ncrst -o min3.out -inf mdinfo.log')

################# MAIN ####################
(adt_home, protein_pdb_file, ligand_file_list, box_center_list, box_dimension_list, 
 temp_seq, ppn, refine_reactant_restraints, QM_reactant_restraints, qm_mask, qm_charge, temp) = parse_input(cntrl_file)

stage = check_log_file()
if stage == 0:
    print('Error: Cannot find checkpoint. Please check the last line of your LOG file')
    sys.exit()

if not os.path.exists('reactant_docking'):
    os.mkdir('reactant_docking')
os.chdir('reactant_docking')
if stage == 1:
    f_log = open('../docking.log', 'w')
    f_log.close()

    pdb = pmd.load_file('../'+protein_pdb_file)
    if temp_seq != '':
        # Modify protonated states according to the res name in temp_seq #
        f = open('../'+temp_seq)
        seq = f.readlines()
        f.close()
        seq = [s.strip() for s in seq]
        for res in pdb.residues:
            res.name = seq[res.idx]
    pdb.save('protein.pdb', overwrite=True)

    f = open('leap.in', 'w')
    f.write('''source leaprc.protein.ff14SB
a = loadpdb protein.pdb
savepdb a protein.pdb
quit\n''')
    f.close()

    os.system('tleap -s -f leap.in > /dev/null')
    pdb = pmd.load_file('protein.pdb')
    # END Modify protonated states according to the res name in temp_seq #

for i in range(len(ligand_file_list)):
    if stage <= i+1:
        f_log = open('../docking.log', 'a')
        f_log.write("Running docking for reactant %d\n"%(i+1))
        f_log.close()
        
        os.system(adt_home+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r protein.pdb -o protein.pdbqt')
        lig_file = ligand_file_list[i][0]
        lig_name = lig_file.split('/')[-1].split('.mol2')[0]
        if not lig_name.endswith('_dummy') and lig_name != 'WAT':
            os.system('antechamber -i ../'+lig_file+' -fi mol2 -o '+lig_name+'.prepi -fo prepi -pf y -dr n > /dev/null')
            os.system('parmchk2 -i '+lig_name+'.prepi -f prepi -o '+lig_name+'.frcmod > /dev/null')
        elif lig_name != 'WAT':
            lig_dir = '/'.join(lig_file.split('/')[:-1])
            os.system('cp ../'+lig_dir+'/'+lig_name+'.prepi ./')
            os.system('cp ../'+lig_dir+'/'+lig_name+'.frcmod ./')
        
        if lig_name.endswith('_dummy'):
            os.system(adt_home+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ../'+lig_file+' -o ligand.pdbqt -C')
        elif len(ligand_file_list[i]) == 1:
            os.system(adt_home+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ../'+lig_file+' -o ligand.pdbqt')
        elif ligand_file_list[i][1] == 'all':
            os.system(adt_home+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ../'+lig_file+' -o ligand.pdbqt -Z')
        else:
            os.system(adt_home+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ../'+lig_file+' -o ligand.pdbqt -I '+ligand_file_list[i][1])
        f = open('config.txt', 'w')
        f.write('''receptor = protein.pdbqt
ligand = ligand.pdbqt

size_x =  '''+('%.2f'%box_dimension_list[i][0])+'''
size_y =  '''+('%.2f'%box_dimension_list[i][1])+'''
size_z =  '''+('%.2f'%box_dimension_list[i][2])+'''
center_x =  '''+('%.3f'%box_center_list[i][0])+'''
center_y =  '''+('%.3f'%box_center_list[i][1])+'''
center_z =  '''+('%.3f'%box_center_list[i][2])+'''

cpu = '''+str(int(ppn))+'\n')
        f.close()
        os.system('vina --config config.txt > docking_%d.log 2>&1'%(i+1))
        
        if not os.path.getsize('ligand_out.pdbqt'):
            f_log = open('../docking.log', 'a')
            f_log.write("Error: No proper binding pose found for %s\n"%lig_file)
            f_log.close()
            sys.exit()
    
        f = open('ligand_out.pdbqt')
        lig_struct = pmd.Structure()
        for line in f.readlines():
            if line.startswith('MODEL 2'):
                break
            if line.startswith('ATOM '):
                atm_name = line[12:16].strip()
                res_name = line[17:20].strip()
                element_name = ''.join(list(filter(str.isalpha, atm_name)))
                if not element_name in ['MC', 'DU']:
                    atom = pmd.topologyobjects.Atom(name=atm_name, atomic_number=pmd.periodic_table.AtomicNum[element_name])
                else:
                    atom = pmd.topologyobjects.Atom(name=atm_name, atomic_number=1)
                lig_struct.add_atom(atom, res_name, 1)
                lig_struct[-1].xx = float(line[30:38].strip())
                lig_struct[-1].xy = float(line[38:46].strip())
                lig_struct[-1].xz = float(line[46:54].strip())
        f.close()
        
        pdb = pdb+lig_struct
        pdb.residues[-1].ter = True
        
        pdb.save('protein.pdb', overwrite=True)
    elif stage <= len(ligand_file_list):
        pdb = pmd.load_file('protein.pdb')
        pdb.residues[-1].ter = True

if stage <= len(ligand_file_list)+1:
    f_log = open('../docking.log', 'a')
    f_log.write('Running MD docking refinement of reactants\n')
    f_log.close()
    lig_list = [ligand[0] for ligand in ligand_file_list]
    refine_docking_MD('protein.pdb', lig_list, refine_reactant_restraints, 2000, 2000, 10000, 200000)

if stage <= len(ligand_file_list)+2:
    f_log = open('../docking.log', 'a')
    f_log.write('Running QMMM docking refinement of reactants\n')
    f_log.close()
    refine_docking_QMMM('complex.prmtop', 'equil_min.ncrst', QM_reactant_restraints, qm_mask, qm_charge, 5000)

    os.system('cp complex.prmtop ../')
    os.system('cp min3.ncrst ../')
    os.chdir('../')
    os.system('mv min3.ncrst reactant_complex.ncrst')
    f_log = open('docking.log', 'a')
    f_log.write("Done.\n")
    f_log.close()
