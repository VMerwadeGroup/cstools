import os, sys
from glob import glob
import subprocess
import shutil
import pytest

## clean the output folders
print('Clean the output folders.')
test_folders = [f for f in glob('./tests/*') if os.path.isdir(f)]

for f in test_folders:
    output_folder = f'{f}/output'
    temp_folder = f'{f}/temp'
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    if os.path.exists(temp_folder):
        shutil.rmtree(temp_folder)
    if 'URL' in f:
        os.makedirs(f'{output_folder}/lackOb')
        os.mkdir(f'{output_folder}/goodOb')
        os.mkdir(f'{output_folder}/multiCl')
    elif 'Raster' in f:
        continue
    else:
        os.makedirs(f'{output_folder}/line')
        os.mkdir(f'{output_folder}/zigzag')
    
if os.path.exists('./temp'):
    shutil.rmtree(output_folder)

def commonTests(folder, fname):
    zigzag_script = f'{folder}/test_zigzag.py'
    line_script = f'{folder}/test_line.py'
    print('\n====================')
    print(f'Testing {fname} for zigzag survey')
    print('--------------------')
    # subprocess.call(['python3.9', zigzag_script])
    subprocess.call([sys.executable, zigzag_script])
    print('--------------------')
    print(f'Test of {fname} for zigzag survey is done.')
    print('====================\n')

    print('\n====================')
    print(f'Testing {fname} for line survey')
    print('--------------------')
    # subprocess.call(['python3.9', line_script])
    subprocess.call([sys.executable, line_script])
    print('--------------------')
    print(f'Test of {fname} for line survey is done.')
    print('====================\n')

def URLTests(folder, fname):
    for case in ['lackOb', 'goodOb', 'multiCl']:
        print('\n====================')
        print(f'Testing datasets {case} from {fname}')
        print('--------------------')
        subprocess.call([sys.executable, f'{folder}/test_URL_{case}.py'])
        print('--------------------')
        print(f'Test of datasets {case} from {fname} is done.')
        print('====================\n')
        shutil.rmtree('./temp')

## start tests
print('Starting tests...')
for f in test_folders:
    folder_name = os.path.basename(f)
    func_name = folder_name.replace('test_', '')
    if ('RiverReach' in func_name):
        continue
    if (func_name == 'URL'):
        URLTests(f, func_name)
    else:
        commonTests(f, func_name)

rrtest = pytest.main(['./tests/test_RiverReach/test_RiverReach.py'])
