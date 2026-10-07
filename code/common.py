"""Shared data loading for the revised analysis (n = 102).
Inputs : ../data/PCR_NGS_Data.xlsx (qPCR + NGS ASV counts), ../data/sample_metadata.csv (clinical and sample metadata)."""
import os, re, warnings, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
DATA=os.path.join(ROOT,'data'); RES=os.path.join(ROOT,'results'); FIG=os.path.join(ROOT,'figures')
os.makedirs(RES,exist_ok=True); os.makedirs(FIG,exist_ok=True)
raw=pd.read_excel(os.path.join(DATA,'PCR_NGS_Data.xlsx'),header=None)
hdr=raw.iloc[6].tolist(); tax=raw.iloc[0:6]
ngs_idx=[i for i,h in enumerate(hdr) if isinstance(h,str) and h.startswith('__') and h!='__']
D=raw.iloc[7:].reset_index(drop=True)
ID_all=D[0].astype(str).str.strip(); y_all=D[1].astype(int).values
A_all=D[ngs_idx].apply(pd.to_numeric,errors='coerce').fillna(0).values          # ASV counts
species=[hdr[i][2:].strip() for i in ngs_idx]
taxo=pd.DataFrame({'Phylum':tax.iloc[2,ngs_idx].values,'Class':tax.iloc[3,ngs_idx].values,'Order':tax.iloc[4,ngs_idx].values,'Family':tax.iloc[5,ngs_idx].values,'Species':species})
for c in ['Phylum','Class','Order','Family']: taxo[c]=taxo[c].astype(str).str.replace('__','').str.strip()
PCR_all=D[list(range(2,18))].apply(pd.to_numeric,errors='coerce'); PCR_all.columns=hdr[2:18]
TB_all=pd.to_numeric(D[18],errors='coerce')
reads_all=A_all.sum(1)
# ---- sample QC: exclude failed sequencing (< 1,000 reads): 2S-16 (0 reads), 2S-20 (1 read)
ok=reads_all>=1000
A=A_all[ok]; ID=ID_all[ok].reset_index(drop=True); y=y_all[ok]
PCR=PCR_all[ok].reset_index(drop=True); TB=TB_all[ok].reset_index(drop=True); reads=reads_all[ok]
SP=pd.DataFrame(A,columns=species).T.groupby(level=0).sum().T          # species-level (putative) abundance
SP=SP.loc[:,SP.sum()>0]
META=pd.read_csv(os.path.join(DATA,'sample_metadata.csv')).set_index('participant_id').loc[ID]
assert (META.halitosis.values==y).all()
pstat=META.periodontal_status.values; perio=(pstat=='periodontitis').astype(int); teeth=META.n_teeth.values.astype(float)
age=META.age.values.astype(float); male=(META.sex.values=='M').astype(float); smoking=META.smoking.values.astype(float)
hyg_mod=(META.oral_hygiene.values=='Moderate').astype(float)
VSC={'H2S (ppb)':META.h2s_ppb.values.astype(float),'CH3SH (ppb)':META.ch3sh_ppb.values.astype(float),'VSCs (ppb)':META.vscs_ppb.values.astype(float)}
panel={'A. actinomycetemcomitans':'Aggregatibacter actinomycetemcomitans','P. intermedia':'Prevotella intermedia','P. nigrescens':'Prevotella nigrescens',
 'L. casei':'Lactobacillus casei','F. nucleatum':'Fusobacterium nucleatum','S. mutans':'Streptococcus mutans','S. sobrinus':'Streptococcus sobrinus',
 'T. denticola':'Treponema denticola','P. gingivalis':'Porphyromonas gingivalis','T. forsythia':'Tannerella forsythia','E. corrodens':'Eikenella corrodens',
 'P. micra':'Parvimonas micra','C. rectus':'Campylobacter rectus','E. nodatum':'Eubacterium nodatum','F. alocis':'Filifactor alocis',   # qPCR panel lists F. alocis under its former name Fusobacterium alocis
 'P. endodontalis':'Porphyromonas endodontalis'}
def R_(f): return os.path.join(RES,f)
def F_(f): return os.path.join(FIG,f)
