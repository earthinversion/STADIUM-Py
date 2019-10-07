import matplotlib.pyplot as plt
import numpy as np
from time import process_time
st = process_time()

t0 = 1.45
t1 = 5.4 
t2 = 18.55
vp = 6.3
p = 0.06 #ray parameter
w1=0.75 # Weight Ps to 0.75
w2 = 0.25 # Weight PpPs to 0.25

try:
    if w1+w2 != 1:
        raise ValueError('Weights are not properly defined')
except ValueError as e:
    exit(str(e))

# Measure the difference between theory and data:

numpoints = 1000
hs = np.linspace(20,40,numpoints)
Kappas = np.linspace(1.5, 2.5, numpoints)
H, K = np.meshgrid(hs, Kappas)
depth1 = (t1-t0)/(np.sqrt((K/vp)**2-(p)**2)-np.sqrt((1/vp)**2-(p)**2))
depth2 = (t2-t0)/(np.sqrt((K/vp)**2-(p)**2)+np.sqrt((1/vp)**2-(p)**2))
deltas = np.absolute((w1*depth1 + w2* depth2) - H)

fig, ax = plt.subplots()
cmap = plt.get_cmap('viridis')
delta_lvs = np.linspace(np.amin(deltas),np.amax(deltas),30)
CS = ax.contourf(H, K, deltas, levels=delta_lvs,cmap=cmap)

result = np.where(deltas == np.amin(deltas))
# print(H[result], K[result], deltas[result])
ax.plot(H[result],K[result],'ko')
ax.clabel(CS, inline=1, fontsize=10, fmt='%2.1f', colors='w')
fig.colorbar(CS)
ax.set_title(r'$H$-$\kappa$ grid search')
# plt.contour(ka, h, delta, 10, cmap='RdGy')
ax.set_xlabel('H')
ax.set_ylabel(r'$\kappa$')
fig.tight_layout()
plt.savefig('h-k_outfile.png')
edt = process_time()
print(f'Total process time: {edt-st}s')