from common import *
import json, itertools
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, recall_score, f1_score, average_precision_score, roc_curve
import statsmodels.api as sm
R={}
n0,n1=(y==0).sum(),(y==1).sum(); R['n']=[int(n0),int(n1)]
# ---------- Table 1 recomputable rows
def mw(a,b): return stats.mannwhitneyu(a,b,alternative='two-sided').pvalue
T1={}
ct=pd.crosstab(pstat,y); T1['perio_ct']=ct.to_dict(); T1['perio_p']=stats.chi2_contingency(ct)[1]
ct2=pd.crosstab(perio,y); T1['perio_bin_p']=stats.fisher_exact(ct2.values)[1]
T1['teeth']=[teeth[y==0].mean(),teeth[y==0].std(ddof=1),teeth[y==1].mean(),teeth[y==1].std(ddof=1),mw(teeth[y==0],teeth[y==1])]
pg=SP['Porphyromonas gingivalis'].values
ctp=pd.crosstab(pg>0,y); T1['pg_det']=[int(((pg>0)&(y==0)).sum()),int(((pg>0)&(y==1)).sum()),stats.fisher_exact(ctp.values)[1]]
tb=TB.values.astype(float)
T1['tb']=[np.median(tb[y==0]),*np.percentile(tb[y==0],[25,75]),np.median(tb[y==1]),*np.percentile(tb[y==1],[25,75]),mw(tb[y==0],tb[y==1])]
# diversity, rarefied
rng=np.random.default_rng(2026); D0=22000; its=10
def rare(v):
    idx=np.repeat(np.arange(len(v)),v.astype(int)); s=rng.choice(idx,D0,replace=False); return np.bincount(s,minlength=len(v))
div=np.zeros((len(y),4))
for i in range(len(y)):
    acc=np.zeros(4)
    for _ in range(its):
        r=rare(A[i]); S=(r>0).sum(); p=r[r>0]/D0; H=-(p*np.log(p)).sum(); f1=(r==1).sum(); f2=(r==2).sum()
        ch=S+(f1**2/(2*f2) if f2>0 else f1*(f1-1)/2); acc+=[S,H,H/np.log(S),ch]
    div[i]=acc/its
T1['div']={nm:[div[y==0,k].mean(),div[y==0,k].std(ddof=1),div[y==1,k].mean(),div[y==1,k].std(ddof=1),mw(div[y==0,k],div[y==1,k])] for k,nm in enumerate(['Richness','Shannon','Evenness','Chao1'])}
T1['age']=[age[y==0].mean(),age[y==0].std(ddof=1),age[y==1].mean(),age[y==1].std(ddof=1),stats.ttest_ind(age[y==0],age[y==1],equal_var=False).pvalue]
T1['sex']=[int(male[y==0].sum()),int((1-male[y==0]).sum()),int(male[y==1].sum()),int((1-male[y==1]).sum()),stats.chi2_contingency(pd.crosstab(male,y))[1] if False else stats.fisher_exact(pd.crosstab(male,y).values)[1]]
T1['smk']=[int(smoking[y==0].sum()),int(smoking[y==1].sum()),stats.fisher_exact(pd.crosstab(smoking,y).values)[1]]
def welch(v): return [v[y==0].mean(),v[y==0].std(ddof=1),v[y==1].mean(),v[y==1].std(ddof=1),stats.ttest_ind(v[y==0],v[y==1],equal_var=False).pvalue]
T1['clin']={'Age':welch(age),**{k:welch(v) for k,v in VSC.items()}}
T1['clin_cat']={'Sex':stats.chi2_contingency(pd.crosstab(male,y))[1],'Smoking':stats.chi2_contingency(pd.crosstab(smoking,y))[1],'Oral hygiene':stats.chi2_contingency(pd.crosstab(hyg_mod,y))[1]}
T1['hyg']=[int((hyg_mod[y==0]==0).sum()),int((hyg_mod[y==1]==0).sum())]
R['T1']=T1
# ---------- Differential abundance
prev=(SP>0).mean(); keep=prev[prev>=0.10].index
rows=[]
for t in SP.columns:
    a=SP.loc[y==0,t].values; b=SP.loc[y==1,t].values
    rows.append(dict(taxon=t,prev=prev[t],prev0=(a>0).mean(),prev1=(b>0).mean(),med0=np.median(a),q10=np.percentile(a,25),q30=np.percentile(a,75),
        med1=np.median(b),q11=np.percentile(b,25),q31=np.percentile(b,75),mean0=a.mean(),sd0=a.std(ddof=1),mean1=b.mean(),sd1=b.std(ddof=1),
        p_w=mw(a,b),p_t=stats.ttest_ind(a,b,equal_var=False).pvalue if (a.std()+b.std())>0 else 1.0,
        direction='Higher in Halitosis' if np.mean(stats.rankdata(np.r_[a,b])[len(a):])>np.mean(stats.rankdata(np.r_[a,b])[:len(a)]) else 'Lower in Halitosis'))
DA=pd.DataFrame(rows); DA['p_t']=DA.p_t.fillna(1)
DA['q_primary']=np.nan; m=DA.taxon.isin(keep)
DA.loc[m,'q_primary']=multipletests(DA.loc[m,'p_w'],method='fdr_bh')[1]
DA['q_all409']=multipletests(DA.p_w,method='fdr_bh')[1]
DA['q_welch']=multipletests(DA.p_t,method='fdr_bh')[1]
DA['p_w_ord']=DA.p_w
DA=DA.sort_values('p_w').drop(columns='p_w_ord')
DA.to_csv(R_('DA_full.csv'),index=False)
R['DA']=dict(n_all=int(len(DA)),n_prev=int(m.sum()),sig_primary=int((DA.q_primary<0.05).sum()),sig_all=int((DA.q_all409<0.05).sum()),sig_welch=int((DA.q_welch<0.05).sum()),
             sig_list=DA[DA.q_primary<0.05].taxon.tolist(),sig_all_list=DA[DA.q_all409<0.05].taxon.tolist())
# zero / skew stats
Z=(SP[keep]==0).mean(); R['dist']=dict(med_zero_all=float((SP==0).mean().median()),frac_gt50=float(((SP==0).mean()>0.5).mean()),
     med_skew=float(SP[keep].apply(lambda c: stats.skew(c)).median()),
     shapiro_frac=float(np.mean([stats.shapiro(SP[c])[1]<0.05 for c in keep])))
# ---------- Correlations
X=pd.DataFrame({k:(SP[v] if v in SP else np.zeros(len(y))) for k,v in panel.items()})
X=X.loc[:,X.sum()>0]
for g in [0,1]:
    for meth in ['pearson','spearman']:
        Xi=X[y==g]; cols=Xi.columns; C=np.full((len(cols),len(cols)),np.nan); Pm=C.copy()
        for i,a in enumerate(cols):
            for j,b in enumerate(cols):
                if Xi[a].std()==0 or Xi[b].std()==0: continue
                r,p=(stats.pearsonr if meth=='pearson' else stats.spearmanr)(Xi[a],Xi[b]); C[i,j]=r; Pm[i,j]=p
        pd.DataFrame(C,index=cols,columns=cols).to_csv(R_(f'corr_{meth}_{g}.csv')); pd.DataFrame(Pm,index=cols,columns=cols).to_csv(R_(f'corrp_{meth}_{g}.csv'))
# ---------- CV models
Xl=np.log1p(X)
models={'Hal_Group_1':['E. corrodens','P. intermedia','P. nigrescens'],'Hal_Group_2':['F. nucleatum','P. gingivalis','T. forsythia'],
 'Hal_Group_3':['P. endodontalis','T. denticola'],'NonHal_Group_1':['P. gingivalis','P. intermedia','T. forsythia'],
 'NonHal_Group_2':['P. endodontalis','P. nigrescens','T. denticola']}
for s in ['E. corrodens','F. nucleatum','P. endodontalis','P. gingivalis','P. intermedia','T. denticola','P. nigrescens','T. forsythia']: models['Single_'+s]=[s]
json.dump(models,open(R_('models.json'),'w'))
cv=RepeatedStratifiedKFold(n_splits=5,n_repeats=50,random_state=42); splits=list(cv.split(Xl,y))
grid=np.linspace(0,1,101)
def run(F,extra=None,spl=splits,yy=y):
    M=Xl[F].values if extra is None else np.c_[Xl[F].values,extra] if F else extra
    out=[];tprs=[]
    for tr,te in spl:
        mdl=make_pipeline(StandardScaler(),LogisticRegression(class_weight='balanced',max_iter=5000))
        mdl.fit(M[tr],yy[tr]); p=mdl.predict_proba(M[te])[:,1]; pr=(p>=0.5).astype(int)
        out.append([roc_auc_score(yy[te],p),balanced_accuracy_score(yy[te],pr),recall_score(yy[te],pr),recall_score(yy[te],pr,pos_label=0),f1_score(yy[te],pr,zero_division=0),average_precision_score(yy[te],p)])
        fpr,tpr,_=roc_curve(yy[te],p); t=np.interp(grid,fpr,tpr); t[0]=0; tprs.append(t)
    return np.array(out),np.array(tprs)
res={};ROC={}
for nm,F in models.items():
    o,t=run(F); res[nm]=o; ROC[nm]=t
def perm_p(F,extra=None,B=200):
    obs=res_mean=None
    cvp=RepeatedStratifiedKFold(n_splits=5,n_repeats=2,random_state=7); sp=list(cvp.split(Xl,y))
    obs=run(F,extra,sp)[0][:,0].mean(); r=np.random.default_rng(11); cnt=0
    for b in range(B):
        yp=r.permutation(y); cnt+= run(F,extra,sp,yp)[0][:,0].mean()>=obs
    return (cnt+1)/(B+1)
PP={nm:perm_p(F) for nm,F in models.items()}
T3=[]
for nm,F in models.items():
    o=res[nm]; T3.append(dict(model=nm,features=', '.join(F),k=len(F),auc=o[:,0].mean(),auc_lo=np.percentile(o[:,0],2.5),auc_hi=np.percentile(o[:,0],97.5),
        ba=o[:,1].mean(),ba_lo=np.percentile(o[:,1],2.5),ba_hi=np.percentile(o[:,1],97.5),sens=o[:,2].mean(),spec=o[:,3].mean(),f1=o[:,4].mean(),prauc=o[:,5].mean(),perm_p=PP[nm]))
T3=pd.DataFrame(T3); T3.to_csv(R_('T3_cv.csv'),index=False)
np.save(R_('roc.npy'),{k:v.mean(0) for k,v in ROC.items()},allow_pickle=True); np.save(R_('roc_sd.npy'),{k:v.std(0) for k,v in ROC.items()},allow_pickle=True)
np.save(R_('fold_auc.npy'),{k:v[:,0] for k,v in res.items()},allow_pickle=True)
def nb(a,b):
    d=a-b; k=len(d); t=d.mean()/np.sqrt((1/k+1/4)*d.var(ddof=1)); return d.mean(),2*stats.t.sf(abs(t),k-1)
pairs=[('NonHal_Group_1','Hal_Group_1'),('NonHal_Group_1','Hal_Group_3'),('NonHal_Group_1','NonHal_Group_2'),('NonHal_Group_2','Hal_Group_3'),('NonHal_Group_1','Hal_Group_2'),
 ('Hal_Group_2','Hal_Group_1'),('Hal_Group_2','Hal_Group_3'),('Hal_Group_2','NonHal_Group_2'),('NonHal_Group_2','Hal_Group_1'),('Hal_Group_3','Hal_Group_1'),
 ('NonHal_Group_1','Single_P. gingivalis'),('Hal_Group_2','Single_P. gingivalis')]
T4=[]
for a,b in pairs:
    d,p=nb(res[a][:,0],res[b][:,0]); T4.append(dict(m1=a,m2=b,auc1=res[a][:,0].mean(),auc2=res[b][:,0].mean(),diff=d,p=p))
T4=pd.DataFrame(T4); T4['p_holm']=multipletests(T4.p,method='holm')[1]; T4.to_csv(R_('T4_pairs.csv'),index=False)
# ---------- Pg + random partners null
others=[c for c in X.columns if c!='P. gingivalis']
nullr=[]
for a,b in itertools.combinations(others,2):
    o,_=run(['P. gingivalis',a,b]); nullr.append(dict(partners=f'{a}+{b}',auc=o[:,0].mean(),ba=o[:,1].mean()))
NULL=pd.DataFrame(nullr); NULL.to_csv(R_('null_pg_partners.csv'),index=False)
# ---------- periodontitis-adjusted
PA=[]
o,_=run([],extra=perio.reshape(-1,1).astype(float)); PA.append(dict(model='Periodontitis only',features='Periodontitis',auc=o[:,0].mean(),auc_lo=np.percentile(o[:,0],2.5),auc_hi=np.percentile(o[:,0],97.5),ba=o[:,1].mean(),sens=o[:,2].mean(),spec=o[:,3].mean(),f1=o[:,4].mean(),prauc=o[:,5].mean()))
o,_=run([],extra=np.c_[perio,teeth].astype(float)); PA.append(dict(model='Periodontitis + number of teeth',features='Periodontitis, number of teeth',auc=o[:,0].mean(),auc_lo=np.percentile(o[:,0],2.5),auc_hi=np.percentile(o[:,0],97.5),ba=o[:,1].mean(),sens=o[:,2].mean(),spec=o[:,3].mean(),f1=o[:,4].mean(),prauc=o[:,5].mean()))
for nm,F in models.items():
    o,_=run(F,extra=np.c_[perio,teeth].astype(float))
    PA.append(dict(model=nm,features=', '.join(F)+', periodontitis, number of teeth',auc=o[:,0].mean(),auc_lo=np.percentile(o[:,0],2.5),auc_hi=np.percentile(o[:,0],97.5),ba=o[:,1].mean(),sens=o[:,2].mean(),spec=o[:,3].mean(),f1=o[:,4].mean(),prauc=o[:,5].mean(),d_auc=o[:,0].mean()-res[nm][:,0].mean()))
pd.DataFrame(PA).to_csv(R_('S_perio_models.csv'),index=False)
# ---------- Table 5: demographic-augmented models (age, sex, smoking, oral hygiene)
DEM=np.c_[age,male,smoking,hyg_mod]
T5=[]
o,_=run([],extra=DEM); T5.append(dict(model='Demographics only',features='Age, sex, smoking, oral hygiene',auc=o[:,0].mean(),auc_lo=np.percentile(o[:,0],2.5),auc_hi=np.percentile(o[:,0],97.5),ba=o[:,1].mean(),sens=o[:,2].mean(),spec=o[:,3].mean(),f1=o[:,4].mean(),prauc=o[:,5].mean(),d_auc=np.nan))
for nm,F in models.items():
    o,_=run(F,extra=DEM)
    T5.append(dict(model=nm,features=', '.join(F)+', age, sex, smoking, oral hygiene',auc=o[:,0].mean(),auc_lo=np.percentile(o[:,0],2.5),auc_hi=np.percentile(o[:,0],97.5),ba=o[:,1].mean(),sens=o[:,2].mean(),spec=o[:,3].mean(),f1=o[:,4].mean(),prauc=o[:,5].mean(),d_auc=o[:,0].mean()-res[nm][:,0].mean()))
pd.DataFrame(T5).to_csv(R_('T5_demographic_models.csv'),index=False)
# adjusted association Pg ~ halitosis + periodontitis (+ demographics)
zpg=(Xl['P. gingivalis']-Xl['P. gingivalis'].mean())/Xl['P. gingivalis'].std()
Xa=sm.add_constant(pd.DataFrame({'Pg':zpg,'Perio':perio,'Teeth':(teeth-teeth.mean())/teeth.std()}))
fa=sm.Logit(y,sm.add_constant(pd.DataFrame({'Pg':zpg,'Perio':perio,'Teeth':(teeth-teeth.mean())/teeth.std(),'Age':(age-age.mean())/age.std(),'Male':male,'Smoking':smoking}))).fit(disp=0)
R['adj_full']=dict(or_pg=float(np.exp(fa.params['Pg'])),ci=np.exp(fa.conf_int().loc['Pg']).tolist(),p=float(fa.pvalues['Pg']))
f=sm.Logit(y,Xa).fit(disp=0); R['adj']=dict(or_pg=float(np.exp(f.params['Pg'])),ci=np.exp(f.conf_int().loc['Pg']).tolist(),p=float(f.pvalues['Pg']),or_perio=float(np.exp(f.params['Perio'])),p_perio=float(f.pvalues['Perio']))
strat={}
for s in ['healthy','gingivitis','periodontitis']:
    mm=pstat==s; strat[s]=[int(((pg>0)&mm&(y==1)).sum()),int((mm&(y==1)).sum()),int(((pg>0)&mm&(y==0)).sum()),int((mm&(y==0)).sum()),mw(pg[mm&(y==0)],pg[mm&(y==1)])]
R['strat']=strat
# ---------- Table 6: linear regression for log-transformed VSC concentrations (repeated 5-fold CV x 50)
# VSC concentrations are highly right-skewed (many zeros, few values > 1,000 ppb); the outcome is log(x + 1)-transformed.
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import r2_score, mean_absolute_error
rk=list(RepeatedKFold(n_splits=5,n_repeats=50,random_state=42).split(Xl))
T6=[]
for tgt,v in VSC.items():
    lv=np.log1p(v)
    for nm,F in models.items():
        r2=[];mae=[]
        for tr,te in rk:
            mdl=make_pipeline(StandardScaler(),LinearRegression()).fit(Xl[F].values[tr],lv[tr]); pr=mdl.predict(Xl[F].values[te])
            r2.append(r2_score(lv[te],pr)); mae.append(mean_absolute_error(v[te],np.expm1(pr)))
        T6.append(dict(model=nm,target=tgt,r2=np.mean(r2),r2_lo=np.percentile(r2,2.5),r2_hi=np.percentile(r2,97.5),mae_ppb=np.mean(mae)))
pd.DataFrame(T6).to_csv(R_('T6_vsc_regression.csv'),index=False)
R['vsc_cor']={k:list(stats.spearmanr(Xl['P. gingivalis'],v)) for k,v in VSC.items()}
R['vsc_cor_withinH']={k:list(stats.spearmanr(Xl['P. gingivalis'][y==1],v[y==1])) for k,v in VSC.items()}
# ---------- feature importance (full data, standardized log1p)
FI=[]
for nm,F in models.items():
    if nm.startswith('Single') and nm!='Single_P. gingivalis': continue
    Z=(Xl[F]-Xl[F].mean())/Xl[F].std()
    try:
        f=sm.Logit(y,sm.add_constant(Z)).fit(disp=0,maxiter=200); ci=f.conf_int()
        for c in F: FI.append(dict(model=nm,feature=c,coef=f.params[c],lo=ci.loc[c,0],hi=ci.loc[c,1],OR=np.exp(f.params[c]),ORlo=np.exp(ci.loc[c,0]),ORhi=np.exp(ci.loc[c,1]),p=f.pvalues[c]))
    except Exception as e: print('FI fail',nm,e)
pd.DataFrame(FI).to_csv(R_('FI.csv'),index=False)
# ---------- qPCR
Q=[]
for c in PCR.columns:
    v=pd.to_numeric(PCR[c],errors='coerce').values.astype(float)
    a,b=v[y==0],v[y==1]
    sh=[k for k,vv in panel.items() if vv==c or (c=='Fusobacterium alocis' and k=='F. alocis')]
    rho=pr=np.nan
    if sh and sh[0] in X: rho,pr=stats.spearmanr(v,X[sh[0]])
    Q.append(dict(target=c,det0=(a>0).mean(),det1=(b>0).mean(),med0=np.median(a),q10=np.percentile(a,25),q30=np.percentile(a,75),med1=np.median(b),q11=np.percentile(b,25),q31=np.percentile(b,75),p=mw(a,b),rho=rho,p_rho=pr))
Q.append(dict(target='Total bacteria',det0=1,det1=1,med0=np.median(tb[y==0]),q10=np.percentile(tb[y==0],25),q30=np.percentile(tb[y==0],75),med1=np.median(tb[y==1]),q11=np.percentile(tb[y==1],25),q31=np.percentile(tb[y==1],75),p=mw(tb[y==0],tb[y==1]),rho=np.nan,p_rho=np.nan))
Q=pd.DataFrame(Q); Q['q']=multipletests(Q.p,method='fdr_bh')[1]; Q.to_csv(R_('qpcr.csv'),index=False)
# ---------- PCoA / PERMANOVA (ASV-level relative abundance)
RA=A/A.sum(1,keepdims=True)
from scipy.spatial.distance import pdist,squareform
Dm=squareform(pdist(RA,'braycurtis')); n=len(y)
Jc=np.eye(n)-1/n; Bm=-0.5*Jc@(Dm**2)@Jc; ev,evec=np.linalg.eigh(Bm); o=np.argsort(ev)[::-1]; ev,evec=ev[o],evec[:,o]
coords=evec[:,:2]*np.sqrt(ev[:2]); varexp=ev[:2]/ev[ev>0].sum()
def pseudoF(lab):
    D2=Dm**2; SST=D2[np.triu_indices(n,1)].sum()/n; SSW=0
    for g in np.unique(lab):
        idx=np.where(lab==g)[0]; SSW+=D2[np.ix_(idx,idx)][np.triu_indices(len(idx),1)].sum()/len(idx)
    k=len(np.unique(lab)); return ((SST-SSW)/(k-1))/(SSW/(n-k)),(SST-SSW)/SST
Fo,R2=pseudoF(y); r=np.random.default_rng(3); cnt=sum(pseudoF(r.permutation(y))[0]>=Fo for _ in range(999))
Fp,R2p=pseudoF(pstat); cntp=sum(pseudoF(r.permutation(pstat))[0]>=Fp for _ in range(999))
R['pcoa']=dict(var=varexp.tolist(),F=Fo,R2=R2,p=(cnt+1)/1000,Fp=Fp,R2p=R2p,pp=(cntp+1)/1000)
np.save(R_('pcoa.npy'),coords)
# composition
comp={}
for lvl in ['Phylum','Class','Order']:
    G=pd.DataFrame(RA,columns=taxo[lvl].values).T.groupby(level=0).sum().T*100; comp[lvl]=G
    G.to_csv(R_(f'comp_{lvl}.csv'),index=False)
R['comp_top']={lvl:comp[lvl].mean().sort_values(ascending=False).head(4).round(1).to_dict() for lvl in comp}
np.save(R_('div.npy'),div)
json.dump(R,open(R_('R.json'),'w'),indent=1,default=float)
print('analysis complete; outputs written to results/')
