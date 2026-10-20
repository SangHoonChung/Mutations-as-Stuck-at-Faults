# -*- coding: utf-8 -*-
"""Compact single-pass vs attractor comparison for the paper main text:
(A) per-cell scatter with multistable cells highlighted + y=x line,
(B) histogram of the difference. Replaces the tall 3-panel (kept as supplement)."""
import csv, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAVY='#1F3A5F'; RED='#C0392B'; GREY='#888'
sp=[]; at=[]; na=[]
with open('MM_attractor_vs_singlepass.csv') as fh:
    for r in csv.DictReader(fh):
        sp.append(float(r['singlepass_burden'])); at.append(float(r['attractor_residual'])); na.append(int(r['n_attractors']))
sp=np.array(sp); at=np.array(at); na=np.array(na)
d=at-sp; multi=na>1
N=len(sp)
pct_multi=100*multi.mean()
pct_diff=100*(np.abs(d)>0.1).mean()
mad=np.abs(d).mean()
r=np.corrcoef(sp,at)[0,1]

fig,(axA,axB)=plt.subplots(1,2,figsize=(7.4,3.25),gridspec_kw={'width_ratios':[1.15,1]})

# ---- Panel A: scatter ----
rng=np.random.default_rng(0)
jx=sp+rng.uniform(-0.012,0.012,N)   # small jitter (single-pass values are discrete)
jy=at+rng.uniform(-0.012,0.012,N)
axA.plot([0,1],[0,1],color=GREY,lw=1,ls='--',zorder=1)
axA.scatter(jx[~multi],jy[~multi],s=14,c=NAVY,alpha=0.55,edgecolors='none',label='monostable',zorder=2)
axA.scatter(jx[multi],jy[multi],s=20,facecolors='none',edgecolors=RED,linewidths=1.1,label='multistable',zorder=3)
axA.set_xlim(-0.03,1.0); axA.set_ylim(-0.03,1.0); axA.set_aspect('equal')
axA.set_xlabel('single-pass burden'); axA.set_ylabel('attractor residual')
axA.set_title('(A) per-condition agreement', fontsize=10)
axA.legend(fontsize=7.5, loc='upper left', frameon=False)
axA.text(0.98,0.06,f'$n$={N} cells\n$r$={r:.2f}\nmultistable {pct_multi:.0f}%\n$|\\Delta|>0.1$: {pct_diff:.0f}%',
         transform=axA.transAxes,ha='right',va='bottom',fontsize=7.5,
         bbox=dict(boxstyle='round',fc='white',ec=GREY,alpha=0.9))
axA.grid(alpha=0.2)

# ---- Panel B: histogram of difference ----
axB.hist(d,bins=np.arange(-0.55,0.56,0.05),color=NAVY,alpha=0.8,edgecolor='white')
axB.axvline(0,color=GREY,ls='--',lw=1)
axB.axvline(d.mean(),color=RED,lw=1.5,label=f'mean {d.mean():+.3f}')
axB.set_xlabel(r'difference  (attractor $-$ single-pass)'); axB.set_ylabel('conditions')
axB.set_title('(B) where dynamics change the score', fontsize=10)
axB.text(0.02,0.97,f'mean $|\\Delta|$={mad:.3f}\nred$>$0: attractor more\nmalignant (feedback-\nsustained)',
         transform=axB.transAxes,ha='left',va='top',fontsize=7.5)
axB.legend(fontsize=7.5, loc='upper right', frameon=False)
axB.grid(axis='y',alpha=0.2)

fig.tight_layout()
fig.savefig('MM_compare_compact.png',dpi=200,bbox_inches='tight')
fig.savefig('MM_compare_compact.svg',bbox_inches='tight')
print(f'saved MM_compare_compact.png/.svg | N={N} multistable={pct_multi:.1f}% |diff|>0.1={pct_diff:.1f}% mean|d|={mad:.3f} r={r:.3f}')
