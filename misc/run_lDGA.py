"""
    Description:
Read through directory tree given by ARG1 and search for lDGA_julia directory. 
Run single core (!) LadderDGA.jl script given in ARG2.
Run ARG3 number of jobs per node.
DIRNAME, ARG4, ARG5, ... are given as additional arguments to the script in the example below res_TEST is given as the second argument

    Example:
python run_lDGA.py /scratch/projects/hhp00048/lDGA/cuprates/CTQMC /scratch/projects/hhp00048/codes/LadderDGA.jl/scripts/ldm_sc.jl 30 res_TEST
"""

import os.path
import sys
import subprocess
import copy


def template_berlin(i,joblist):
    res = """#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --ntasks 96
#SBATCH -p cpu-clx:test
#SBATCH -A hhp00048
#SBATCH -J batch{0}

{1}

wait
""".format(i,"".join(joblist))
    return res
batch_len = int(sys.argv[3])
batch_list = []
batch_list_i = []
i_tmp = 0

for path, directories, files in os.walk(sys.argv[1]):
    if path.endswith("lDGA_julia") and "config.toml" in files:
        cmd = "julia " + sys.argv[2] + " " + str(path) + "/config.toml"
        cmd += " " + str(path) + "/" + sys.argv[4]
        for el in sys.argv[5:]:
            cmd += " " + el
        cmd += " & \n"
        batch_list_i.append(cmd)
        i_tmp += 1
        if i_tmp >= batch_len:
            batch_list.append(copy.deepcopy(batch_list_i))
            batch_list_i.clear()
            i_tmp = 0

batch_list.append(copy.deepcopy(batch_list_i))
original_stdout = sys.stdout


for i,el in enumerate(batch_list):
    print("running job #" + str(i))
    job_f = template_berlin(i,el)

    with open('run_tmp_'+str(i)+'.sh', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print("sbatch <<EOF\n" + job_f + "EOF")
        sys.stdout = original_stdout

    pr = subprocess.run("sbatch <<EOF\n" + job_f + "EOF", shell=True, capture_output=True)
    if not pr.returncode == 0:
        print("Submission of slurm job failed!")
        print(pr.stdout.decode("utf-8"))
        print(pr.stderr.decode("utf-8"))
