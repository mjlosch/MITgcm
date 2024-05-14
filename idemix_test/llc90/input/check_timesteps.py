import glob
import numpy as np

flist = glob.glob("diags2D.*.meta")

n = len(flist)
tsnumber = []
for f in flist:
    tsnumber.append(int(f.split('.')[1]))

tsnumber.sort()
tsnumber_diff=np.diff(tsnumber)
# print(tsnumber_diff)

ii = np.where(np.logical_or(tsnumber_diff<8760,tsnumber_diff>8784))[0]

print('number of years = %i, %f'%(n,float(n)/62.))
print('difference < 8760 or > 8784 ')
print('============================')
for i in ii:
    print("index %i: %i - %i = %i"%( i, tsnumber[i+1], tsnumber[i],
                                        tsnumber_diff[i]) )
