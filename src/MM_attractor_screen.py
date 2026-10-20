# -*- coding: utf-8 -*-
"""Attractor-based screen + comparison to the single-pass DAG."""
import numpy as np, csv, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from MM_attractor import analyze
from MM_network import mm_measure, FAULT_LABELS, N_FAULTS

REGIMENS = [
    ('Untreated',[0,0,0,0,0]),('Isa',[1,0,0,0,0]),('Dara',[0,1,0,0,0]),
    ('V (Btz)',[0,0,1,0,0]),('d (Dex)',[0,0,0,1,0]),('R (Len)',[0,0,0,0,1]),
    ('Vd',[0,0,1,1,0]),('Rd',[0,0,0,1,1]),('VR',[0,0,1,0,1]),('DR',[0,1,0,0,1]),
    ('VRd',[0,0,1,1,1]),('DRd',[0,1,0,1,1]),('DVd',[0,1,1,1,0]),
    ('Isa-Rd',[1,0,0,1,1]),('D-VRd',[0,1,1,1,1]),('Isa-VRd',[1,0,1,1,1]),
]
faults=list(range(1,N_FAULTS+1)); fx=[FAULT_LABELS[f] for f in faults]
ATT=np.zeros((len(REGIMENS),len(faults)))     # attractor residual (linear burden)
SP =np.zeros((len(REGIMENS),len(faults)))     # single-pass burden = sqrt(measure)
NA =np.zeros((len(REGIMENS),len(faults)),int)  # n attractors (multistability)
for i,(_,vec) in enumerate(REGIMENS):
    for j,f in enumerate(faults):
        r=analyze([f],vec,n_init=120)
        ATT[i,j]=r['residual']; NA[i,j]=r['n_attractors']
        SP[i,j]=mm_measure([f],vec)**0.5

# save raw comparison
with open('MM_attractor_vs_singlepass.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['regimen','fault','singlepass_burden','attractor_residual','n_attractors'])
    for i,(rn,_) in enumerate(REGIMENS):
        for j,f in enumerate(faults):
            w.writerow([rn,FAULT_LABELS[f],round(SP[i,j],3),round(ATT[i,j],3),NA[i,j]])

def heat(ax,Mx,title,cmap='RdYlGn_r',vmin=0,vmax=None,fmt='{:.2f}'):
    vmax=vmax if vmax is not None else (Mx.max() or 1)
    im=ax.imshow(Mx,cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto')
    ax.set_xticks(range(len(faults))); ax.set_xticklabels(fx,rotation=90,fontsize=6.5)
    ax.set_yticks(range(len(REGIMENS))); ax.set_yticklabels([n for n,_ in REGIMENS],fontsize=8)
    ax.set_title(title,fontsize=11)
    for r in range(Mx.shape[0]):
        for c in range(Mx.shape[1]):
            v=Mx[r,c]
            if abs(v)>1e-9: ax.text(c,r,fmt.format(v).lstrip('0').replace('-0.','-.'),ha='center',va='center',fontsize=5)
    ax.set_xticks(np.arange(-.5,len(faults),1),minor=True); ax.set_yticks(np.arange(-.5,len(REGIMENS),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=0.7); ax.tick_params(which='minor',length=0)
    return im

# Figure 1: attractor residual heatmap
fig,ax=plt.subplots(figsize=(15.5,8.2))
im=heat(ax,ATT,'Attractor-based residual malignant burden (feedback restored; dead cells = 0)\n'
        'green = suppressed, red = uncontrolled  [compare to single-pass heatmap]')
ax.set_xlabel('Fault',fontsize=10)
fig.colorbar(im,ax=ax,fraction=0.025,pad=0.01).set_label('residual burden (0=healthy)',fontsize=9)
fig.tight_layout(); fig.savefig('MM_attractor_heatmap.png',dpi=150,bbox_inches='tight'); plt.close(fig)

# Figure 2: single-pass vs attractor + difference + multistability
fig,axes=plt.subplots(3,1,figsize=(15.5,21))
heat(axes[0],SP,'(A) Single-pass DAG  -  burden = sqrt(measure)')
heat(axes[1],ATT,'(B) Attractor model  -  residual burden')
D=ATT-SP
imd=heat(axes[2],D,'(C) Difference  (attractor - single-pass):  red = attractor MORE malignant, blue = LESS',
         cmap='bwr',vmin=-0.5,vmax=0.5)
# mark multistable cells in panel C with a dot
for r in range(NA.shape[0]):
    for c in range(NA.shape[1]):
        if NA[r,c]>1: axes[2].plot(c,r,'o',mfc='none',mec='black',ms=9,mew=1.2)
axes[2].set_xlabel('Fault  (circles = multistable: >1 attractor)',fontsize=10)
fig.colorbar(imd,ax=axes[2],fraction=0.025,pad=0.01).set_label('attractor - single-pass',fontsize=9)
fig.tight_layout(); fig.savefig('MM_attractor_vs_singlepass.png',dpi=140,bbox_inches='tight'); plt.close(fig)

# summary stats
multi=int((NA>1).sum()); diff=np.abs(D)
print(f"cells total = {NA.size}")
print(f"multistable cells (>1 attractor) = {multi}  ({100*multi/NA.size:.1f}%)")
print(f"mean |attractor - single-pass| = {diff.mean():.3f} ; max = {diff.max():.3f}")
print(f"cells where they differ by >0.1 = {(diff>0.1).sum()}")
# where attractor is MORE malignant (bistable survival the DAG missed)
idx=np.dstack(np.unravel_index(np.argsort(-D,axis=None),D.shape))[0][:6]
print('\nTop cells where attractor shows MORE residual than single-pass:')
for r,c in idx:
    print(f"  {REGIMENS[r][0]:9s} x {FAULT_LABELS[faults[c]]:28s}  SP={SP[r,c]:.2f} ATT={ATT[r,c]:.2f} attr#={NA[r,c]}")
