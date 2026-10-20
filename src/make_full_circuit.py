# -*- coding: utf-8 -*-
"""Full MM Boolean-network logic circuit (ANSI gate symbols), laid out by graphviz dot."""
import subprocess, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as mp

PRIMARY = ['IL6','IGF1','RANKL','WNT','HYPOXIA','TNF_in','PTEN','CRBN','Isa','Dara','Btz','Dex','Len']

GATES = {
 'TNF':('BUF',['TNF_in']), 'IGF1R':('BUF',['IGF1']), 'RANK':('BUF',['RANKL']),
 'STAT3':('BUF',['IL6']), 'AKT':('BUF',['PI3K']), 'MTOR':('BUF',['AKT']),
 'ERK':('BUF',['IGF1R']), 'BCAT':('BUF',['WNT']), 'NIK':('BUF',['RANK']),
 'IRF4':('BUF',['IKZF']), 'E2F':('BUF',['CDK46']), 'NOXA':('BUF',['Btz']),
 'nPTEN':('NOT',['PTEN']), 'nBtz':('NOT',['Btz']), 'nDex':('NOT',['Dex']),
 'nSTAT3':('NOT',['STAT3']), 'nAKT':('NOT',['AKT']), 'nBIM':('NOT',['BIM']),
 'nP21':('NOT',['P21']), 'nNOXA':('NOT',['NOXA']), 'nguarded':('NOT',['guarded']),
 'nAPOP':('NOT',['APOP']), 'IKZF':('NOT',['Len_active']), 'P53':('NOT',['MDM2']),
 'antiCD38':('OR',['Isa','Dara']), 'PI3K':('OR',['IGF1R','IL6','nPTEN']),
 'orNFKB':('OR',['TNF','RANK','AKT','NIK']),
 'orMYC':('OR',['ERK','BCAT','STAT3','NFKB','MTOR']),
 'CCND2':('OR',['MYC','ERK','BCAT','NFKB']), 'MDM2':('OR',['AKT','ERK']),
 'P21':('OR',['P53','Len_active']), 'MCL1':('OR',['STAT3','NFKB']),
 'BCLXL':('OR',['NFKB','STAT3']), 'guarded':('OR',['free_MCL1','free_BCLXL']),
 'pro':('OR',['NOXA','BIM','antiCD38']), 'ANTIAPOP':('OR',['MCL1','BCLXL']),
 'VEGF':('OR',['STAT3','HYPOXIA']),
 'Len_active':('AND',['Len','CRBN']), 'NFKB':('AND',['orNFKB','nBtz','nDex']),
 'MYC':('AND',['orMYC','IRF4']), 'CDK46':('AND',['CCND2','nP21']),
 'BIM':('AND',['Dex','nSTAT3','nAKT']), 'free_MCL1':('AND',['MCL1','nNOXA','nBIM']),
 'free_BCLXL':('AND',['BCLXL','nBIM']), 'APOP':('AND',['pro','nguarded']),
 'SURVIVAL':('AND',['ANTIAPOP','nAPOP']),
 'out_CCND2':('BUF',['CCND2']), 'out_MYC':('BUF',['MYC']), 'out_E2F':('BUF',['E2F']),
 'out_NFKB':('BUF',['NFKB']), 'out_ANTIAPOP':('BUF',['ANTIAPOP']), 'out_VEGF':('BUF',['VEGF']),
 'out_CDK46':('BUF',['CDK46']), 'out_SURVIVAL':('BUF',['SURVIVAL']),
}
OUTLABEL = {'out_CCND2':'X1 CCND2','out_MYC':'X2 MYC','out_E2F':'X3 E2F',
            'out_NFKB':'X4 NF-kB','out_ANTIAPOP':'X5 ANTIAPOP','out_VEGF':'X6 VEGF',
            'out_CDK46':'X7 CDK4/6','out_SURVIVAL':'X8 SURVIVAL'}

# ---- build DOT ----
lines = ['digraph G {','rankdir=LR; splines=ortho; nodesep=0.28; ranksep=0.75;',
         'node [shape=box, fixedsize=true, width=0.58, height=0.40];']
for p in PRIMARY: lines.append(f'"{p}";')
for g in GATES:   lines.append(f'"{g}";')
for g,(t,ins) in GATES.items():
    for s in ins: lines.append(f'"{s}" -> "{g}";')
lines.append('}')
open('circuit.dot','w').write('\n'.join(lines))

plain = subprocess.run(['dot','-Tplain','circuit.dot'],capture_output=True,text=True).stdout

# ---- parse plain ----
pos={}; edges=[]; W=H=0
for ln in plain.splitlines():
    t=ln.split()
    if t[0]=='graph': W=float(t[2]); H=float(t[3])
    elif t[0]=='node': pos[t[1]]=(float(t[2]),float(t[3]),float(t[4]),float(t[5]))
    elif t[0]=='edge':
        n=int(t[3]); xy=list(map(float,t[4:4+2*n]))
        edges.append(list(zip(xy[0::2],xy[1::2])))

# ---- draw ----
fig,ax=plt.subplots(figsize=(min(W,40),min(H,32)))
ax.set_xlim(-0.5,W+1.6); ax.set_ylim(-0.4,H+0.4); ax.axis('off'); ax.set_aspect('equal')

for pts in edges:
    ax.plot([p[0] for p in pts],[p[1] for p in pts],color='#C0392B',lw=0.7,zorder=1)

def and_poly(x0,y0,w,h):
    r=h/2; cy=y0+h/2; pts=[(x0,y0),(x0,y0+h),(x0+w-r,y0+h)]
    for a in [math.radians(90-10*i) for i in range(19)]:
        pts.append((x0+w-r+r*math.cos(a), cy+r*math.sin(a)))
    pts.append((x0+w-r,y0)); return pts
def or_path(x0,y0,w,h):
    cy=y0+h/2
    verts=[(x0,y0+h),(x0+0.6*w,y0+h),(x0+w,cy),(x0+0.6*w,y0),(x0,y0),(x0+0.22*w,cy),(x0,y0+h)]
    codes=[Path.MOVETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,Path.CURVE3,Path.CURVE3,Path.CURVE3]
    return Path(verts,codes)

def draw_gate(cx,cy,bw,bh,kind):
    x0,y0=cx-bw/2,cy-bh/2
    if kind in ('BUF','NOT'):
        tip = bw-(0.14*bh if kind=='NOT' else 0)
        ax.add_patch(mp.Polygon([(x0,y0),(x0,y0+bh),(x0+tip,cy)],closed=True,
                     facecolor='white',edgecolor='black',lw=1.0,zorder=3))
        if kind=='NOT':
            ax.add_patch(mp.Circle((x0+tip+0.07*bh,cy),0.07*bh,facecolor='white',edgecolor='black',lw=1.0,zorder=3))
    elif kind=='AND':
        ax.add_patch(mp.Polygon(and_poly(x0,y0,bw,bh),closed=True,
                     facecolor='white',edgecolor='black',lw=1.0,zorder=3))
    elif kind=='OR':
        ax.add_patch(mp.PathPatch(or_path(x0,y0,bw,bh),facecolor='white',edgecolor='black',lw=1.0,zorder=3))

for name,(x,y,w,h) in pos.items():
    if name in PRIMARY:
        ax.add_patch(mp.FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.01',
                     facecolor='#DCEAF7',edgecolor='#2C5777',lw=0.9,zorder=3))
        ax.text(x,y,name,ha='center',va='center',fontsize=4.6,zorder=4)
    else:
        kind=GATES[name][0]
        draw_gate(x,y,w*0.92,h*0.86,kind)
        ax.text(x,y-h*0.62,name,ha='center',va='center',fontsize=3.7,color='#222',zorder=4)
        if name in OUTLABEL:
            ax.annotate('',xy=(x+w/2+0.45,y),xytext=(x+w/2,y),
                        arrowprops=dict(arrowstyle='->',color='black',lw=0.9),zorder=3)
            ax.text(x+w/2+0.5,y,OUTLABEL[name],ha='left',va='center',fontsize=6,fontweight='bold',zorder=4)

ax.set_title('Multiple Myeloma Boolean network as a logic-gate circuit  '
             '(BUFFER / NOT / AND / OR ; inputs left, 8 outputs right)',fontsize=12,pad=10)
fig.tight_layout()
fig.savefig('MM_full_circuit.png',dpi=170,bbox_inches='tight')
fig.savefig('MM_full_circuit.svg',bbox_inches='tight')
print('saved MM_full_circuit.png/.svg  canvas(in)=',round(W,1),round(H,1),' gates=',len(GATES))
