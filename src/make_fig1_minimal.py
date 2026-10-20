# -*- coding: utf-8 -*-
"""Minimalist Figure-1 overview: 4-stage narrative as clean concept icons;
detailed circuit/state-graph deferred to Fig. 5 (attractor) and Fig. 6 (gate core)."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Arc

NAVY='#1F3A5F'; HEAD='#22345B'; PANEL='#EEF1F5'; PEDGE='#C9D2DD'
RED='#C0392B'; GREEN='#2E7D46'; BLUE='#2E7CB8'; GREY='#8A94A0'; DARK='#222'

fig,ax=plt.subplots(figsize=(15,6.3)); ax.set_xlim(0,16); ax.set_ylim(0,7)
ax.set_aspect('equal'); ax.axis('off')

L,R,GAP=0.3,0.3,0.36; PW=(16-L-R-3*GAP)/4; PTOP,PBOT=6.05,1.05
def px(i): return L+i*(PW+GAP)

def panel(i,title,subtitle,footer):
    x0=px(i)
    ax.add_patch(FancyBboxPatch((x0,PBOT),PW,PTOP-PBOT,boxstyle='round,pad=0.02,rounding_size=0.12',
                 fc=PANEL,ec=PEDGE,lw=1.2,zorder=1))
    hh=0.62
    ax.add_patch(FancyBboxPatch((x0+0.12,PTOP-hh-0.02),PW-0.24,hh,boxstyle='round,pad=0.02,rounding_size=0.10',
                 fc=HEAD,ec='none',zorder=2))
    ax.text(x0+PW/2,PTOP-hh/2-0.02,title,ha='center',va='center',color='white',fontsize=15,fontweight='bold',zorder=3)
    ax.text(x0+PW/2,PTOP-hh-0.30,subtitle,ha='center',va='center',color=NAVY,fontsize=10.5,style='italic',zorder=3)
    for k,ln in enumerate(footer):
        ax.text(x0+PW/2,PBOT+0.52-0.28*k,ln,ha='center',va='center',color=GREY,fontsize=9,style='italic',zorder=3)
    return x0

def chip(cx,cy,w,h,txt,fc,tc='white',fs=10.5,ec='none',lw=1.1):
    ax.add_patch(FancyBboxPatch((cx-w/2,cy-h/2),w,h,boxstyle='round,pad=0.02,rounding_size=0.09',
                 fc=fc,ec=ec,lw=lw,zorder=4))
    ax.text(cx,cy,txt,ha='center',va='center',color=tc,fontsize=fs,fontweight='bold',zorder=5)

def arrow(x1,y1,x2,y2,color=DARK,lw=2.0,ms=11,rad=0.0):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=ms,color=color,lw=lw,
                 shrinkA=0,shrinkB=0,connectionstyle=f'arc3,rad={rad}',zorder=4))

def tri(i):
    xg=px(i)+PW+GAP/2; y=(PTOP+PBOT)/2
    ax.add_patch(FancyArrowPatch((xg-0.12,y),(xg+0.12,y),arrowstyle='-|>',mutation_scale=13,color=NAVY,lw=2.6,zorder=6))

# ---- Panel 1: Genotype ----
x0=panel(0,'1 · Genotype','stuck-at faults',['mutation pins a node','to 0 / 1  (stuck-at)'])
for cy,tag,node in [(4.45,'sa1','MYC'),(3.45,'sa0','TP53')]:
    chip(x0+0.95,cy,0.74,0.5,tag,RED,fs=11)
    ax.plot([x0+1.32,x0+1.75],[cy,cy],color=RED,lw=2.2,zorder=3)
    ax.add_patch(Circle((x0+2.18,cy),0.42,fc='white',ec=RED,lw=2.2,zorder=4))
    ax.text(x0+2.18,cy,node,ha='center',va='center',color=RED,fontsize=9.5,fontweight='bold',zorder=5)
ax.add_patch(Circle((x0+2.18,2.45),0.42,fc='white',ec=NAVY,lw=2.2,zorder=4))
ax.text(x0+2.18,2.45,'NF-κB',ha='center',va='center',color=NAVY,fontsize=8.5,fontweight='bold',zorder=5)
tri(0)

# ---- Panel 2: Therapy ----
x0=panel(1,'2 · Therapy','masking inputs',['drug masks faults;','combo = test pattern'])
for d,cy in zip(['Btz','Dex','Len','anti-CD38'],[4.55,3.75,2.95,2.15]):
    chip(x0+1.30,cy,1.9,0.55,d,BLUE,fs=11)
    arrow(x0+2.30,cy,x0+3.02,cy,color=BLUE,lw=2.4,ms=12)
tri(1)

# ---- Panel 3: Single-pass (minimal) ----
x0=panel(2,'3 · Single-pass','combinational circuit',['feed-forward eval =','pure Boolean function'])
my=3.55
for dy in (0.5,0,-0.5):
    arrow(x0+0.40,my+dy,x0+0.95,my+dy*0.35,color=GREY,lw=1.8,ms=9)
chip(x0+1.62,my,1.28,1.05,'∧ ∨ ¬',fc='white',tc=NAVY,fs=16,ec=NAVY,lw=1.6)
ax.text(x0+1.62,my-0.68,'combinational logic',ha='center',va='center',color=NAVY,fontsize=8,zorder=5)
arrow(x0+2.28,my,x0+2.66,my,color=DARK,lw=2.2,ms=12)
chip(x0+3.06,my,0.86,0.86,'burden\n0–1',fc=NAVY,tc='white',fs=9.5)
tri(2)

# ---- Panel 4: Attractor (minimal) ----
x0=panel(3,'4 · Attractor','sequential circuit / FSM',['registers + feedback = FSM;','attractors = cell fates'])
my=3.55
chip(x0+1.28,my,1.35,1.05,'∧ ∨ ¬',fc='white',tc=NAVY,fs=15,ec=NAVY,lw=1.6)
ax.text(x0+1.28,my-0.68,'next-state logic',ha='center',va='center',color=NAVY,fontsize=8,zorder=5)
ax.add_patch(FancyArrowPatch((x0+1.75,my+0.60),(x0+0.82,my+0.60),arrowstyle='-|>',mutation_scale=12,
             color=RED,lw=2.1,connectionstyle='arc3,rad=0.75',zorder=5))
ax.text(x0+1.28,my+1.20,'feedback + clock',ha='center',va='center',color=RED,fontsize=8.5,
        fontweight='bold',style='italic',zorder=6)
arrow(x0+1.96,my,x0+2.44,my,color=DARK,lw=2.2,ms=12)
fx=x0+3.02
ax.text(fx,my+1.20,'cell fates',ha='center',va='center',color=NAVY,fontsize=9,fontweight='bold',zorder=5)
ax.add_patch(Circle((fx,my+0.45),0.30,fc='white',ec=RED,lw=2.2,zorder=4))
ax.text(fx,my+0.45,'mal',ha='center',va='center',color=RED,fontsize=8.5,fontweight='bold',zorder=5)
ax.add_patch(Arc((fx+0.02,my+0.80),0.36,0.36,theta1=-40,theta2=210,color=RED,lw=1.5,zorder=4))
ax.add_patch(Circle((fx,my-0.5),0.30,fc='white',ec=GREEN,lw=2.2,zorder=4))
ax.text(fx,my-0.5,'apo',ha='center',va='center',color=GREEN,fontsize=8.5,fontweight='bold',zorder=5)
ax.add_patch(Arc((fx+0.02,my-0.85),0.36,0.36,theta1=150,theta2=400,color=GREEN,lw=1.5,zorder=4))
ax.text(x0+PW/2,PBOT+1.02,'(gate-level detail in supplement)',ha='center',va='center',color=GREY,fontsize=8,style='italic',zorder=5)

# ---- title + bottom strip ----
ax.text(8,6.72,'Mutations as stuck-at faults: from combinational logic to a finite-state machine',
        ha='center',va='center',color=NAVY,fontsize=16.5,fontweight='bold')
eq=['mutation = stuck-at fault','drug = input','combination = test pattern','cell fate = attractor']
xs=[px(i)+PW/2 for i in range(4)]
for xx,t in zip(xs,eq):
    ax.text(xx,0.5,t,ha='center',va='center',color=DARK,fontsize=12)

fig.savefig('MM_figure1_overview.png',dpi=400,bbox_inches='tight',facecolor='white')
print('saved')
