import numpy as np
import sys

def main():
    infile = sys.argv[1]
    prefix = infile.split('.')[0]
    outfile = prefix + '_tdpot.xyz'

    with open(infile, 'r') as fi:
        with open(outfile, 'w') as fo:
            natom = fi.readline().strip()
            pbc = [float(i) for i in fi.readline().strip().split()]
            fo.write(f'{natom}\n')
            fo.write(f'{pbc[0]:3.5f} {pbc[1]:3.5f} {pbc[2]:3.5f}')
            for l in fi.readlines():
                l = l.strip()
                ndat = [float(i) for i in l.split()]
                fo.write(f'\n{ndat[0]:3.5f}  {ndat[1]:3.5f}  {ndat[2]:2.5f}  {6}  {12}')

if __name__=="__main__":
    main()
