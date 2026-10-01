#!/usr/bin/env python3
"""Presentation only: plot frozen summary tables; no model fitting or inference."""
import os
from pathlib import Path
import os
os.environ['MPLCONFIGDIR']=str(Path(os.environ['AGE_FIGURE_TEMP'])/'mpl')
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R=Path(__file__).resolve().parents[1]
E=R.parent/'09_EXTERNAL_VALIDATION'
TMP=Path(os.environ['AGE_FIGURE_TEMP']);TMP.mkdir(exist_ok=True)
INPUTS=[*sorted((R/'tables').glob('*.tsv')),*sorted((E/'tables').glob('*.tsv'))]
before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUTS if p.name!='figure_2_external_source.tsv'}
p=pd.read_csv(R/'tables/patient_phenotype_WGS_chemistry.tsv',sep='\t')
stat=pd.read_csv(R/'tables/retained_statistics.tsv',sep='\t').set_index('result')
null=pd.read_csv(R/'tables/five_CpG_matched_null.tsv',sep='\t')
local=pd.read_csv(R/'tables/local_discordance_statistics.tsv',sep='\t')
normal=pd.read_csv(E/'tables/normal_brain_validation.tsv',sep='\t').set_index('CpGmarker')
external=pd.read_csv(E/'tables/GBM_validation.tsv',sep='\t').set_index(['cohort','CpGmarker'])
summary=pd.read_csv(E/'tables/five_locus_summary.tsv',sep='\t').set_index('CpGmarker')
LOCI=['cg07158339','cg22947000','cg01511567','cg19761273','cg13460409']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
 'axes.titlesize':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,
 'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
blue='#236987';orange='#BC7038';teal='#2F8075';gray='#7A8188'
colors={s:(blue if s.startswith('GBM') else orange) for s in p.sample_id}
def title(ax,letter,text):ax.set_title(letter+'  '+text,loc='left',fontweight='bold',pad=12)
def save(fig,name):
    fig.savefig(R/'figures'/(name+'.pdf'),bbox_inches='tight',metadata={'Title':name,'Author':'','CreationDate':None,'ModDate':None})
    fig.savefig(TMP/(name+'.png'),dpi=170,bbox_inches='tight');plt.close(fig)

# Figure 1: retain phenotype panels and uncertainty; correct clipped paediatric x.
fig,(a,b)=plt.subplots(1,2,figsize=(10.3,6.5),gridspec_kw={'width_ratios':[1,1.08]})
fig.subplots_adjust(wspace=.45,bottom=.14,top=.87)
for _,r in p.iterrows():
    a.errorbar(r.chronological_age,r.Age_BS,yerr=[[r.Age_BS-r.Age_BS_CI_low],[r.Age_BS_CI_high-r.Age_BS]],fmt='o',ms=4.5,color=colors[r.sample_id],alpha=.8,lw=.8)
a.plot([0,175],[0,175],ls='--',color=gray,lw=1,label='Equality')
a.set(xlim=(10,90),ylim=(0,175),xlabel='Chronological age (years)',ylabel='Projected DNAm age (years)')
for label,col in [('Glioblastoma',blue),('Meningioma',orange)]:a.scatter([],[],color=col,s=24,label=label)
a.legend(frameon=False,loc='upper left',fontsize=8);title(a,'a','Apparent age acceleration')
pp=p.sort_values('DNAm_age_acceleration',ascending=False)
for i,(_,r) in enumerate(pp.iterrows()):
    b.errorbar(r.DNAm_age_acceleration,i,xerr=[[r.DNAm_age_acceleration-r.DeltaAge_CI_low],[r.DeltaAge_CI_high-r.DNAm_age_acceleration]],fmt='o',ms=4,color=colors[r.sample_id],lw=1)
b.set(yticks=np.arange(20),yticklabels=pp.sample_id.tolist(),xlabel='DNAm age minus chronological age (years)',ylim=(19.8,-.8))
b.axvline(0,color=gray,ls='--',lw=1);b.axvline(stat.loc['mean_age_offset','estimate'],color=teal,ls=':',lw=1.2)
title(b,'b','Patient-level offsets')
fig.suptitle('Mean apparent acceleration +32.59 years (95% CI 21.94–43.81)',fontsize=12,fontweight='bold',y=.98)
fig.text(.07,.035,'19/20 positive estimates; 16/20 conditional sampling intervals above zero.\nIntervals quantify molecule sampling, not cross-platform calibration error.',fontsize=9)
assert a.get_xlim()[0]<p.chronological_age.min()
save(fig,'FIGURE_1_phenotype')

# Figure 2: no new models/CI calculations; all forests consume saved estimates.
fig=plt.figure(figsize=(8.2,8.7))
gs=fig.add_gridspec(2,3,height_ratios=[1,1.65],left=.105,right=.975,top=.88,bottom=.29,wspace=.65,hspace=.70)
a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[0,2])
d=fig.add_subplot(gs[1,0]);e=fig.add_subplot(gs[1,1:])
for _,r in p[p.primary_purity.notna()].iterrows():
    a.scatter(100*r.primary_purity,r.DNAm_age_acceleration,s=22,facecolor='white' if r.purity_at_search_boundary else colors[r.sample_id],edgecolor=colors[r.sample_id],lw=1)
a.set(xlabel='Qualified WGS purity (%)',ylabel='Apparent acceleration (years)')
title(a,'a','Genomic context')
a.text(.03,.97,'n = 13; ρ = 0.622\n95% CI 0.157–0.896',transform=a.transAxes,va='top',fontsize=7)
a.margins(y=.30)
for _,r in p.iterrows():b.plot([0,1],[r.DNAm_age_acceleration,r.remaining_DeltaAge],color=colors[r.sample_id],alpha=.4,lw=.65)
b.scatter([0,1],[stat.loc['mean_age_offset','estimate'],stat.loc['mean_offset_after_top5_reference_replacement','estimate']],s=32,color=teal,zorder=5)
b.set(xticks=[0,1],xticklabels=['Full\nclock','Five loci\nreplaced'],xlim=(-.25,1.25),ylabel='Apparent acceleration (years)')
b.axhline(0,color=gray,lw=.7);title(b,'b','Five-locus influence')
b.text(.5,1.01,'Mean +32.59 → +17.72',transform=b.transAxes,ha='center',fontsize=7)
c.hist(null.mean_replacement_effect_years.to_numpy(),bins=38,color='#A8BDC5',edgecolor='white',lw=.25)
c.axvline(stat.loc['five_CpG_mean_reference_replacement_effect','estimate'],color=orange,lw=1.7)
c.set(xlabel='Replacement effect (years)',ylabel='Matched five-locus sets',xlim=(-1,17))
title(c,'c','Conditional matched null')
c.text(.98,.96,'Selected: 14.87\nMedian null: 4.21\nTail: 1/10,001',transform=c.transAxes,ha='right',va='top',fontsize=6.8,bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1))
categories={
'cg07158339':('Method-sensitive\nTCGA unavailable','#fff4dd'),
'cg22947000':('Age-concordant\nexaggerated','#e6f2ed'),
'cg01511567':('Opposes normal age\ntest-dependent','#efeaf5'),
'cg19761273':('Opposes normal age\ninconsistent','#efeaf5'),
'cg13460409':('Age-concordant\nexaggerated','#e6f2ed')}
cc={'GSE60274':'#167d8d','GSE195684':'#7950a3','TCGA_GBM':'#c67829'}
offset={'GSE60274':-.20,'GSE195684':0,'TCGA_GBM':.20}
source=[]
for j,l in enumerate(LOCI):
    r=normal.loc[l];x=r.robust_beta_slope_per_decade*100;lo=r.slope_CI_low*100;hi=r.slope_CI_high*100
    d.errorbar(x,j,xerr=[[x-lo],[hi-x]],fmt='o',color='#333333',ms=4,capsize=2,lw=1)
    source.append(dict(panel='d',cohort='GSE74193',CpGmarker=l,n=r.n_analysed,estimate=x,CI_low=lo,CI_high=hi,unit='beta percentage points per decade',BH_q=r.age_association_BH_q,source_table='../09_EXTERNAL_VALIDATION/tables/normal_brain_validation.tsv',source_field='robust_beta_slope_per_decade;slope_CI_low;slope_CI_high',formal_classification=summary.loc[l,'classification'],display_annotation='normal age direction'))
    label,bg=categories[l];e.axhspan(j-.44,j+.44,color=bg,zorder=-2)
    for cohort in cc:
        r=external.loc[(cohort,l)];x=r.normal_trajectory_delta*100;lo=r.normal_trajectory_low*100;hi=r.normal_trajectory_high*100
        if np.isfinite(x):e.errorbar(x,j+offset[cohort],xerr=[[x-lo],[hi-x]],fmt='o',ms=3.7,color=cc[cohort],capsize=2,lw=1)
        source.append(dict(panel='e',cohort=cohort,CpGmarker=l,n=r.normal_trajectory_n,estimate=x,CI_low=lo,CI_high=hi,unit='beta percentage points beyond normal age trajectory',BH_q=r.BH_q,source_table='../09_EXTERNAL_VALIDATION/tables/GBM_validation.tsv',source_field='normal_trajectory_delta;normal_trajectory_low;normal_trajectory_high',formal_classification=summary.loc[l,'classification'],display_annotation=label.replace('\n','; ')))
    e.text(33,j,label,fontsize=6.5,va='center',ha='left',linespacing=1.35)
d.axvline(0,color=gray,ls='--',lw=.75);e.axvline(0,color=gray,ls='--',lw=.75)
d.set(yticks=range(5),yticklabels=LOCI,ylim=(4.5,-.5),xlim=(-1.2,1.9),xlabel='Normal age slope\n(percentage points/decade)')
e.set(yticks=range(5),yticklabels=[],ylim=(4.5,-.5),xlim=(-25,62),xticks=[-20,0,20],xlabel='Excess beyond normal age prediction\n(percentage points)')
title(d,'d','Normal adult DLPFC');title(e,'e','Independent GBM cohorts: direction and evidence')
d.text(.5,1.015,'n = 226; all five q ≤ 0.00345',transform=d.transAxes,ha='center',fontsize=6.7)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=cc[k],label=v) for k,v in [('GSE60274','GSE60274 (n=63)'),('GSE195684','Nordic (n=116)'),('TCGA_GBM','TCGA (n=138)')]],frameon=False,ncol=3,fontsize=7,loc='upper center',bbox_to_anchor=(.54,.19),borderaxespad=0)
fig.suptitle('Clock distortion and external validation\nAge-concordant excess coexists with opposing, heterogeneous shifts',fontsize=10,fontweight='bold',y=.98,linespacing=1.5)
fig.text(.105,.115,'a: Adjusted purity CI crosses zero. b–c: Same-cohort locus selection; the matched tail is not selection-adjusted.\nd–e: Saved 95% bootstrap CIs. Cross-study/platform/tissue uncertainty is not included.',fontsize=6.9,linespacing=1.45)
fig.text(.105,.06,'Strongest reproducible support: cg22947000 and cg13460409. Other loci retain inconclusive formal replication.\nOpposing trends do not establish tumour-cell specificity. Arrays do not resolve 5mC/5hmC or molecule architecture.',fontsize=6.9,linespacing=1.45)
src=pd.DataFrame(source)
src['q_endpoint']=np.where(src.panel.eq('d'),'normal-age Spearman association','primary tumour-control permutation; NOT a test of the plotted trajectory residual')
src.to_csv(R/'tables/figure_2_external_source.tsv',sep='\t',index=False,float_format='%.12g')
save(fig,'FIGURE_2_clock_distortion')

# Figure 3: replace weak component-partition panel with measured M/H signal shares.
fig=plt.figure(figsize=(10.2,6.8))
gs=fig.add_gridspec(2,3,width_ratios=[1.05,1.05,1.35],height_ratios=[1,.95],left=.09,right=.97,top=.82,bottom=.20,wspace=.55,hspace=.36)
a=fig.add_subplot(gs[:,0]);b=fig.add_subplot(gs[:,1]);c=fig.add_subplot(gs[0,2]);note=fig.add_subplot(gs[1,2]);note.axis('off')
pp=p.sort_values('sample_id');y=np.arange(len(pp));m=pp.absolute_measured_5mC_signal_fraction.to_numpy()*100
a.barh(y,m,height=.64,color=blue,label='5mC');a.barh(y,100-m,left=m,height=.64,color='#cf9a5a',label='5hmC')
a.set(yticks=y,yticklabels=pp.sample_id,ylim=(19.8,-.8),xlim=(0,100),xlabel='Absolute measured clock signal (%)')
fig.legend(handles=[Line2D([],[],color=blue,lw=5,label='5mC'),Line2D([],[],color='#cf9a5a',lw=5,label='5hmC')],frameon=False,ncol=2,fontsize=8,loc='upper left',bbox_to_anchor=(.09,.12));title(a,'a','Modification identity')
for i,(_,r) in enumerate(pp.iterrows()):
    x=r.hydroxymethylation_clock_effect
    b.errorbar(x,i,xerr=[[x-r.hydroxymethylation_clock_effect_CI_low],[r.hydroxymethylation_clock_effect_CI_high-x]],fmt='o',color=teal,ms=3.8,lw=.85)
b.axvline(0,color=gray,ls='--',lw=.8);b.axvline(stat.loc['mean_5hmC_clock_effect','estimate'],color=orange,ls=':',lw=1)
b.set(yticks=y,yticklabels=[],ylim=(19.8,-.8),xlabel='Age(M+H) − age(M), years');title(b,'b','Signed 5hmC effect')
for i,l in enumerate(['cg19761273','cg01511567']):
    r=local[(local.scenario=='primary')&(local.metric=='PDR')&(local.locus==l)].iloc[0]
    c.errorbar(r.effect*100,i,xerr=[[(r.effect-r.CI_low)*100],[(r.CI_high-r.effect)*100]],fmt='o',color=teal,ms=5,capsize=3,lw=1.3)
    c.text(8.05,i,f'n={int(r.n_patients)}',va='center',ha='right',fontsize=7)
c.set(yticks=[0,1],yticklabels=['cg19761273','cg01511567'],ylim=(1.6,-.65),xlim=(-.3,8.3),xlabel='PDR excess (percentage points)')
c.axvline(0,color=gray,ls='--',lw=.75);title(c,'c','Local beta-controlled discordance')
note.text(0,1,'Local four-CpG discordance\nFive-locus q = 0.041 at both regions.\nJoint ten-endpoint q = 0.082.\nSignificance is not universal under\npatient omission or dependence stress.\n\nNo robust binary-entropy excess.\nThe supported signal is focal, not\na general tumour-disorder phenotype.',fontsize=8,va='top',linespacing=1.55)
fig.suptitle('Predominantly 5mC-based clock inputs and focal local methylation discordance',fontsize=12,fontweight='bold',y=.97)
fig.text(.09,.884,'Median 5mC share: 91.4%. Mean signed 5hmC effect: +1.80 years (95% CI −0.31–3.96).',fontsize=9)
fig.text(.09,.035,'a: Shares exclude imputation and do not represent fractions of age. b: Conditional molecule-sampling intervals.\nc: 95% patient-bootstrap CIs; complete four-CpG windows with ≥20 molecules; matched beta, coverage, density, span and context.',fontsize=8)
save(fig,'FIGURE_3_molecular_basis')
for path,h in before.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h,path
print('Three main PDFs updated from frozen values; original canonical and external result tables unchanged.')
