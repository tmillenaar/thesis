import os

os.system('rm outputData/output*.xy')
os.system('gfortran LinDiv.f90')
os.system('./a.out')
os.system('python3 animate.py')
