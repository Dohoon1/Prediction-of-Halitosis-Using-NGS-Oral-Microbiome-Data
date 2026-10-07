from common import *
import json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
plt.rcParams.update({'font.family':'Liberation Sans','font.size':9,'axes.linewidth':0.8,'savefig.dpi':400,'axes.spines.top':False,'axes.spines.right':False})
R=json.load(open(R_('R.json'))); models=json.load(open(R_('models.json')))
DA=pd.read_csv(R_('DA_full.csv')); T3=pd.read_csv(R_('T3_cv.csv')); Q=pd.read_csv(R_('qpcr.csv')); FI=pd.read_csv(R_('FI.csv')); NULL=pd.read_csv(R_('null_pg_partners.csv'))
C0,C1='#7F7F7F','#C0392B'
def L(ax,s,x=-0.12,y=1.04): ax.text(x,y,s,transform=ax.transAxes,fontsize=13,fontweight='bold',va='bottom',ha='left')
def pfmt(p):
    if p<0.001: 
        e=int(np.floor(np.log10(p))); m=p/10**e; return f'{m:.1f}×10$^{{{e}}}$'
    return f'{p:.3f}'
def short(t):
    t=t.replace('[','').replace(']',''); g,s=t.split()[0],' '.join(t.split()[1:]); return f'{g[0]}. {s}'
def itl(t): # mathtext italic
    return '$\\it{'+t.replace(' ','\\ ')+'}$'
# ---------------- Figure 1 workflow
fig,ax=plt.subplots(figsize=(9.0,4.2)); ax.axis('off'); ax.set_xlim(0,100); ax.set_ylim(0,50)
def box(x,y,w,h,txt,fc='#EAF2F8',fs=6.6,bold=False):
    ax.add_patch(mpatches.FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.3,rounding_size=1.2',fc=fc,ec='#34495E',lw=0.7))
    ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=fs,fontweight='bold' if bold else 'normal')
def arr(x1,y1,x2,y2): ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='-|>',lw=0.8,color='#34495E'))
box(1,21,13,8,'104 participants\n(72 Non-Halitosis,\n32 Halitosis)')
box(18,21,15,8,'Unstimulated saliva\n16S rRNA V3–V4 NGS\n+ targeted qPCR')
box(36.5,20,15,10,'Sequencing QC\n2 failed samples\nexcluded\nn = 102 (72 / 30)',fc='#FDEBD0')
arr(14.3,25,17.7,25); arr(33.3,25,36.7,25)
box(56,40,19,8,'Community ecology\nα-diversity; Bray–Curtis\nPCoA; PERMANOVA')
box(56,28.5,19,8,'Differential abundance\nWilcoxon rank-sum, prevalence\n≥ 10%, BH-FDR (202 taxa)')
box(56,17,19,8,'Co-occurrence groups\nwithin-group correlation\n(|r| ≥ 0.6)')
box(56,5.5,19,8,'qPCR validation\ntotal bacterial load;\n16 target taxa')
for yy in [44,32.5,21,9.5]: arr(51.3,25,55.7,yy)
box(80,17,19,19,'Logistic regression\nindividual taxa vs\nco-occurrence groups\n\nrepeated stratified\n5-fold CV (×50)\nclass-weighted;\nAUROC, balanced accuracy;\npermutation test',fc='#E8F8F5')
arr(75.3,21,79.7,24); arr(75.3,32.5,79.7,30)
box(80,3,19,10,'Sensitivity analyses\nperiodontal status;\n$\\it{P.\\ gingivalis}$ + random\npartner taxa',fc='#E8F8F5')
arr(89.5,16.7,89.5,13.3)
fig.savefig(F_('Figure1.png'),bbox_inches='tight'); plt.close()
# ---------------- Figure 2
sig=DA[DA.q_primary<0.05].sort_values('q_primary')
fig=plt.figure(figsize=(7.2,6.6)); gs=fig.add_gridspec(2,1,height_ratios=[1.15,1],hspace=0.55)
ax=fig.add_subplot(gs[0]); L(ax,'A',-0.08)
pos=np.arange(len(sig))
for i,t in enumerate(sig.taxon):
    for g,c,dx in [(0,C0,-0.18),(1,C1,0.18)]:
        v=np.log10(SP.loc[y==g,t].values+1)
        bp=ax.boxplot(v,positions=[i+dx],widths=0.3,patch_artist=True,showfliers=False,medianprops=dict(color='k',lw=0.9))
        for b in bp['boxes']: b.set(facecolor=c,alpha=0.35,edgecolor=c)
        ax.scatter(i+dx+np.random.default_rng(i).uniform(-0.08,0.08,len(v)),v,s=4,color=c,alpha=0.8,lw=0)
    q=sig.q_primary.iloc[i]; ax.text(i,4.05,'q = '+pfmt(q),ha='center',fontsize=6.3,rotation=0)
ax.set_xticks(pos); ax.set_xticklabels([itl(short(t)) for t in sig.taxon],rotation=35,ha='right',fontsize=7.5)
ax.set_ylabel('log$_{10}$(abundance + 1)'); ax.set_ylim(-0.1,4.4)
ax.legend(handles=[mpatches.Patch(color=C0,alpha=0.6,label=f'Non-Halitosis (n = {R["n"][0]})'),mpatches.Patch(color=C1,alpha=0.6,label=f'Halitosis (n = {R["n"][1]})')],frameon=False,fontsize=7,loc='upper right',bbox_to_anchor=(1,1.22),ncol=2)
ax=fig.add_subplot(gs[1]); L(ax,'B',-0.08)
T=DA[DA.q_primary.notna()].copy()
T['lfc']=np.log2((T.mean1+1)/(T.mean0+1)); T['nl']=-np.log10(T.q_primary)
ns=T.q_primary>=0.05
ax.scatter(T.lfc[ns],T.nl[ns],s=8,color='#B0B0B0',lw=0,label='q ≥ 0.05')
ax.scatter(T.lfc[~ns],T.nl[~ns],s=16,color=C1,lw=0,label='q < 0.05')
ax.axhline(-np.log10(0.05),ls='--',lw=0.7,color='k'); ax.axvline(0,lw=0.5,color='k')
from adjustText import adjust_text
tx=[ax.text(r.lfc,r.nl,itl(short(r.taxon)),fontsize=6.3) for _,r in T[~ns].iterrows()]
adjust_text(tx,ax=ax,x=T.lfc.values,y=T.nl.values,arrowprops=dict(arrowstyle='-',lw=0.4,color='grey'),expand=(1.3,1.6))
ax.set_xlabel('log$_{2}$ fold change in mean abundance (Halitosis / Non-Halitosis)'); ax.set_ylabel('−log$_{10}$ FDR q-value')
ax.legend(frameon=False,fontsize=7,loc='upper left')
fig.savefig(F_('Figure2.png'),bbox_inches='tight'); plt.close()
# ---------------- Figure 3 / Supp Fig 3 heatmaps
def heat(meth,fname,letters=('A','B')):
    fig,axs=plt.subplots(1,2,figsize=(8.6,4.3),gridspec_kw=dict(wspace=0.55))
    for k,g in enumerate([0,1]):
        Cm=pd.read_csv(R_(f'corr_{meth}_{g}.csv'),index_col=0); Pm=pd.read_csv(R_(f'corrp_{meth}_{g}.csv'),index_col=0)
        ax=axs[k]; im=ax.imshow(Cm.values,cmap='RdBu_r',vmin=-1,vmax=1)
        n=len(Cm)
        for i in range(n):
            for j in range(n):
                if i==j or np.isnan(Cm.values[i,j]): continue
                p=Pm.values[i,j]; s='***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
                ax.text(j,i,f'{Cm.values[i,j]:.2f}\n{s}' if s else f'{Cm.values[i,j]:.2f}',ha='center',va='center',fontsize=3.6,color='w' if abs(Cm.values[i,j])>0.6 else 'k')
        lab=[itl(c) for c in Cm.columns]
        ax.set_xticks(range(n)); ax.set_xticklabels(lab,rotation=90,fontsize=5.5); ax.set_yticks(range(n)); ax.set_yticklabels(lab,fontsize=5.5)
        for s in ax.spines.values(): s.set_visible(True)
        L(ax,letters[k],-0.32,1.01)
    cb=fig.colorbar(im,ax=axs,shrink=0.6,pad=0.02); cb.set_label(('Pearson' if meth=='pearson' else 'Spearman')+' correlation coefficient',fontsize=7); cb.ax.tick_params(labelsize=6)
    fig.savefig(fname,bbox_inches='tight'); plt.close()
heat('pearson',F_('Figure4.png')); heat('spearman',F_('SupplementaryFigure3.png'))
# ---------------- Figure 4 networks
order=['NonHal_Group_1','NonHal_Group_2','Hal_Group_1','Hal_Group_2','Hal_Group_3']
fig,axs=plt.subplots(2,3,figsize=(8.4,5.8)); axs=axs.ravel(); axs[5].axis('off')
cols={'P. gingivalis':'#F1948A','P. intermedia':'#C39BD3','T. forsythia':'#F5B7B1','P. nigrescens':'#E0B98A','P. endodontalis':'#82E0AA','T. denticola':'#CCD1D1','E. corrodens':'#85C1E9','F. nucleatum':'#F8C471'}
for k,nm in enumerate(order):
    g=0 if nm.startswith('Non') else 1
    Cm=pd.read_csv(R_(f'corr_pearson_{g}.csv'),index_col=0); Pm=pd.read_csv(R_(f'corrp_pearson_{g}.csv'),index_col=0)
    F=models[nm]; ax=axs[k]; ax.axis('off'); ax.set_xlim(-1.75,1.75); ax.set_ylim(-1.6,1.6); ax.set_aspect('equal')
    ang=np.pi/2+np.arange(len(F))*2*np.pi/len(F) if len(F)==3 else np.array([np.pi*0.75,-np.pi*0.25])
    xy={f:(np.cos(a),np.sin(a)) for f,a in zip(F,ang)}
    for i,a in enumerate(F):
        for b in F[i+1:]:
            r=Cm.loc[a,b]; p=Pm.loc[a,b]
            lw=0.6+4*max(abs(r),0); st='-' if abs(r)>=0.6 else ':'
            ax.plot([xy[a][0],xy[b][0]],[xy[a][1],xy[b][1]],st,color='#566573' if r>0 else '#2E86C1',lw=lw if abs(r)>=0.6 else 0.8,zorder=1)
            mx,my=(xy[a][0]+xy[b][0])/2,(xy[a][1]+xy[b][1])/2
            ox,oy=(mx*0.35,my*0.35) if len(F)==3 else (0.25,0.25)
            ax.text(mx+ox,my+oy,f'r = {r:.2f}\np = {pfmt(p)}',fontsize=6,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=0.5,alpha=0.9),zorder=3)
    for f in F:
        ax.add_patch(plt.Circle(xy[f],0.42,fc=cols[f],ec='#424949',lw=0.7,zorder=2)); ax.text(*xy[f],itl(f.split()[0])+'\n'+itl(f.split()[1]),ha='center',va='center',fontsize=7,zorder=4)
    ax.text(-1.75,1.6,'ABCDE'[k],fontsize=13,fontweight='bold',va='top')
fig.savefig(F_('Figure5.png'),bbox_inches='tight'); plt.close()
# ---------------- Figure 5
roc=np.load(R_('roc.npy'),allow_pickle=True).item(); rsd=np.load(R_('roc_sd.npy'),allow_pickle=True).item(); grid=np.linspace(0,1,101)
T3i=T3.set_index('model')
fig=plt.figure(figsize=(8.4,7.2)); gs=fig.add_gridspec(2,2,hspace=0.42,wspace=0.62)
def rocpanel(ax,names,colors,let):
    for nm,c in zip(names,colors):
        m=roc[nm]; s=rsd[nm]; r=T3i.loc[nm]
        ax.plot(grid,m,color=c,lw=1.2,label=f'{itl(nm.replace("Single_","")) if nm.startswith("Single") else nm}\nAUROC {r.auc:.3f} (perm. p {"< 0.005" if r.perm_p<0.005 else "= "+format(r.perm_p,".3f")})')
        ax.fill_between(grid,np.clip(m-s,0,1),np.clip(m+s,0,1),color=c,alpha=0.12,lw=0)
    ax.plot([0,1],[0,1],'--',color='grey',lw=0.7); ax.set_xlabel('1 − specificity'); ax.set_ylabel('Sensitivity'); ax.set_xlim(0,1); ax.set_ylim(0,1.02)
    ax.legend(fontsize=5.6,frameon=False,loc='lower right'); L(ax,let)
rocpanel(fig.add_subplot(gs[0,0]),['NonHal_Group_1','Hal_Group_2','Single_P. gingivalis'],['#1F618D','#C0392B','#D68910'],'A')
rocpanel(fig.add_subplot(gs[0,1]),['Hal_Group_1','Hal_Group_3','NonHal_Group_2'],['#17A589','#7D3C98','#839192'],'B')
ax=fig.add_subplot(gs[1,0]); L(ax,'C',-0.42)
T=T3.sort_values('auc'); yy=np.arange(len(T))
ax.errorbar(T.auc,yy,xerr=[T.auc-T.auc_lo,T.auc_hi-T.auc],fmt='o',ms=3.5,color='#1F618D',ecolor='#7FB3D5',elinewidth=1,capsize=1.5)
ax.axvline(0.5,ls='--',color='grey',lw=0.7)
for i,(_,r) in enumerate(T.iterrows()): ax.text(1.03,i,('p < 0.005' if r.perm_p<0.005 else f'p = {r.perm_p:.3f}'),fontsize=5.6,va='center',transform=ax.get_yaxis_transform())
ax.set_yticks(yy); ax.set_yticklabels([ (itl(m.replace('Single_','')) if m.startswith('Single') else m) for m in T.model],fontsize=6.2)
ax.set_xlabel('Cross-validated AUROC (2.5th–97.5th percentile)'); ax.set_xlim(0,1.0)
ax=fig.add_subplot(gs[1,1]); L(ax,'D')
ax.hist(NULL.auc,bins=20,color='#BFC9CA',edgecolor='white')
pgA=T3i.loc['Single_P. gingivalis','auc']; ax.axvline(pgA,color='#D68910',lw=1.1,ls='--')
for nm,c in [('NonHal_Group_1','#1F618D'),('Hal_Group_2','#C0392B')]:
    a=T3i.loc[nm,'auc']; pe=(np.sum(NULL.auc>=a)+1)/(len(NULL)+1)
    ax.axvline(a,color=c,lw=1.3); ax.text(a,ax.get_ylim()[1]*(0.97 if c=='#1F618D' else 0.80),f'{nm} \nempirical p = {pe:.2f} ',color=c,fontsize=5.8,va='top',ha='right')
ax.text(pgA,ax.get_ylim()[1]*0.55,'$\\it{P.\\ gingivalis}$ alone ',color='#D68910',fontsize=5.6,ha='right')
ax.set_xlabel('AUROC, $\\it{P.\\ gingivalis}$ + two partner taxa'); ax.set_ylabel(f'Number of partner pairs (n = {len(NULL)})')
fig.savefig(F_('Figure6.png'),bbox_inches='tight'); plt.close()
# ---------------- Figure 6 feature importance (OR per SD)
fig,axs=plt.subplots(3,2,figsize=(7.6,6.6),gridspec_kw=dict(wspace=0.75,hspace=0.6)); axs=axs.ravel(); axs[5].axis('off')
for k,nm in enumerate(order):
    d=FI[FI.model==nm].sort_values('coef'); ax=axs[k]
    yy=np.arange(len(d)); c=[C1 if v>0 else '#2E86C1' for v in d.coef]
    ax.errorbar(d.OR,yy,xerr=[d.OR-d.ORlo,d.ORhi-d.OR],fmt='none',ecolor='#7F8C8D',elinewidth=1,capsize=2)
    ax.scatter(d.OR,yy,c=c,s=22,zorder=3); ax.axvline(1,ls='--',color='grey',lw=0.7); ax.set_xscale('log')
    ax.set_yticks(yy); ax.set_yticklabels([itl(f) for f in d.feature],fontsize=7)
    for i,(_,r) in enumerate(d.iterrows()): ax.text(1.02,i,'p = '+pfmt(r.p),transform=ax.get_yaxis_transform(),fontsize=6,va='center')
    ax.set_xlabel('Odds ratio per 1-SD increase (95% CI)',fontsize=7)
    from matplotlib.ticker import LogLocator,FuncFormatter,NullFormatter
    lo,hi=d.ORlo.min(),d.ORhi.max(); cand=[0.01,0.1,1,10,100,1000] if hi/lo>50 else [0.25,0.5,1,2,3]; ticks=[t for t in cand if lo*0.8<=t<=hi*1.25]
    ax.set_xticks(ticks); ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_: f'{v:g}')); ax.xaxis.set_minor_formatter(NullFormatter()); ax.tick_params(labelsize=7); ax.set_ylim(-0.6,len(d)-0.4); L(ax,'ABCDE'[k],-0.3)
fig.savefig(F_('Figure7.png'),bbox_inches='tight'); plt.close()
# ---------------- Figure 7 qPCR
fig=plt.figure(figsize=(7.2,6.0)); gs=fig.add_gridspec(2,3,hspace=0.45,wspace=0.45)
def bx(ax,v,lab,p,let):
    for g,c in [(0,C0),(1,C1)]:
        vv=np.log10(v[y==g]+1); b=ax.boxplot(vv,positions=[g],widths=0.5,patch_artist=True,showfliers=False,medianprops=dict(color='k'))
        for bb in b['boxes']: bb.set(facecolor=c,alpha=0.35,edgecolor=c)
        ax.scatter(g+np.random.default_rng(g).uniform(-0.12,0.12,len(vv)),vv,s=5,color=c,lw=0)
    ax.set_xticks([0,1]); ax.set_xticklabels(['Non-\nHalitosis','Halitosis'],fontsize=7); ax.set_ylabel(lab,fontsize=7)
    top=ax.get_ylim()[1]; ax.plot([0,0,1,1],[top*1.01,top*1.03,top*1.03,top*1.01],color='k',lw=0.7); ax.text(0.5,top*1.04,'p = '+pfmt(p),ha='center',fontsize=7); L(ax,let,-0.35)
Qi=Q.set_index('target')
bx(fig.add_subplot(gs[0,0]),TB.values.astype(float),'log$_{10}$ total bacteria (copies/mL)',Qi.loc['Total bacteria','p'],'A')
bx(fig.add_subplot(gs[0,1]),PCR['Porphyromonas gingivalis'].values.astype(float),'log$_{10}$ $\\it{P.\\ gingivalis}$ (copies/mL, qPCR)',Qi.loc['Porphyromonas gingivalis','p'],'B')
ax=fig.add_subplot(gs[0,2]); L(ax,'C',-0.35)
a=np.log10(PCR['Porphyromonas gingivalis'].values.astype(float)+1); b=np.log10(SP['Porphyromonas gingivalis'].values+1)
ax.scatter(a[y==0],b[y==0],s=8,color=C0,lw=0,label='Non-Halitosis'); ax.scatter(a[y==1],b[y==1],s=8,color=C1,lw=0,label='Halitosis')
r=Qi.loc['Porphyromonas gingivalis']; ax.text(0.03,0.97,f'Spearman ρ = {r.rho:.2f}\np = {pfmt(r.p_rho)}',transform=ax.transAxes,va='top',fontsize=6.5)
ax.set_xlabel('qPCR, log$_{10}$(copies/mL + 1)',fontsize=7); ax.set_ylabel('NGS, log$_{10}$(abundance + 1)',fontsize=7); ax.legend(fontsize=6,frameon=False,loc='upper left',bbox_to_anchor=(0,0.85))
ax=fig.add_subplot(gs[1,:]); L(ax,'D',-0.06)
QQ=Q.copy(); QQ['lab']=[t if t=='Total bacteria' else short(t.replace('Fusobacterium alocis','Filifactor alocis')) for t in QQ.target]
QQ['dir']=np.sign(QQ.med1-QQ.med0); QQ=QQ.sort_values('p')
v=-np.log10(QQ.q); ax.bar(range(len(QQ)),v,color=[C1 if q<0.05 else '#BFC9CA' for q in QQ.q])
ax.axhline(-np.log10(0.05),ls='--',lw=0.7,color='k')
ax.set_xticks(range(len(QQ))); ax.set_xticklabels([l if l=='Total bacteria' else itl(l) for l in QQ.lab],rotation=40,ha='right',fontsize=6.5)
ax.set_ylabel('−log$_{10}$ FDR q-value',fontsize=7)
fig.savefig(F_('Figure3.png'),bbox_inches='tight'); plt.close()
# ---------------- Supp Fig 1 composition
fig,axs=plt.subplots(1,3,figsize=(7.4,4.0))
n0,n1=R["n"]
ordr=np.r_[np.where(y==0)[0],np.where(y==1)[0]]
for k,lvl in enumerate(['Phylum','Class','Order']):
    G=pd.read_csv(R_(f'comp_{lvl}.csv')); top=G.mean().sort_values(ascending=False).index[:8]
    H=G[top].copy(); H['Others']=100-H.sum(1); ax=axs[k]; bottom=np.zeros(len(y)); cm=plt.get_cmap('tab20')
    handles=[]
    for i,c in enumerate(H.columns):
        p=stats.mannwhitneyu(H[c][y==0],H[c][y==1]).pvalue
        ax.bar(range(len(y)),H[c].values[ordr],bottom=bottom,width=1.0,color=cm(i),lw=0); bottom+=H[c].values[ordr]
        handles.append(mpatches.Patch(color=cm(i),label=f'{c} (p = {pfmt(p)})'))
    ax.axvline(n0-0.5,color='k',lw=1); ax.set_xlim(-0.5,len(y)-0.5); ax.set_ylim(0,100)
    ax.set_xticks([n0/2,n0+n1/2]); ax.set_xticklabels(['Non-Halitosis','Halitosis'],fontsize=7); ax.set_ylabel('Relative abundance (%)' if k==0 else '')
    ax.legend(handles=handles[::-1],fontsize=4.9,frameon=False,loc='upper center',bbox_to_anchor=(0.5,-0.09),ncol=1)
    L(ax,'ABC'[k],-0.2)
n0,n1=R['n']
fig.savefig(F_('SupplementaryFigure1.png'),bbox_inches='tight'); plt.close()
# ---------------- Supp Fig 2 PCoA
co=np.load(R_('pcoa.npy')); pc=R['pcoa']
fig,axs=plt.subplots(1,2,figsize=(7.2,3.4))
ax=axs[0]
for g,c,l in [(0,C0,'Non-Halitosis'),(1,C1,'Halitosis')]: ax.scatter(co[y==g,0],co[y==g,1],s=12,color=c,alpha=0.8,lw=0,label=l)
ax.text(0.02,0.98,f'PERMANOVA R$^{{2}}$ = {pc["R2"]:.3f}\np = {pc["p"]:.3f}',transform=ax.transAxes,va='top',fontsize=7)
ax=axs[1]
for s,c,l in [('healthy','#27AE60','Healthy'),('gingivitis','#F39C12','Gingivitis'),('periodontitis','#8E44AD','Periodontitis')]:
    m=pstat==s; axs[1].scatter(co[m,0],co[m,1],s=12,color=c,alpha=0.8,lw=0,label=l)
ax.text(0.02,0.98,f'PERMANOVA R$^{{2}}$ = {pc["R2p"]:.3f}\np = {pc["pp"]:.3f}',transform=ax.transAxes,va='top',fontsize=7)
for k,ax in enumerate(axs):
    ax.set_xlabel(f'PCo1 ({pc["var"][0]*100:.1f}%)'); ax.set_ylabel(f'PCo2 ({pc["var"][1]*100:.1f}%)'); ax.legend(fontsize=6.5,frameon=False,loc='lower right'); L(ax,'AB'[k],-0.2)
fig.tight_layout(); fig.savefig(F_('SupplementaryFigure2.png'),bbox_inches='tight'); plt.close()
print('done')

# TIFF copies for submission
from PIL import Image
import glob
for f in glob.glob(F_('*.png')):
    Image.open(f).convert('RGB').save(f[:-4]+'.tif',compression='tiff_lzw',dpi=(400,400))
