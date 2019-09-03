# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import codecs
from collections import Counter
from collections import OrderedDict


## part 1
cnep_df = pd.read_excel('./output/cntr-ep差异分析结果_filtbygene.xlsx', encoding='utf-8')
cnep_hwe_df = cnep_df[(cnep_df['HWE_p.value'] > 0.05) & (cnep_df['allelic_chisq_p.value'] < 0.05)]
cnep_genelist = cnep_hwe_df['Gene.refGene'].to_list()
cnep_genelist_dict = dict(Counter(cnep_genelist))
cnep_genegt5 = {}
for k in sorted(cnep_genelist_dict,key=cnep_genelist_dict.__getitem__, reverse=True):
    if cnep_genelist_dict[k] < 5:
        break
    cnep_genegt5[k] = cnep_genelist_dict[k]

cnidi_df = pd.read_excel('./output/cntr-idi差异分析结果_filtbygene.xlsx', encoding='utf-8')
cnidi_hwe_df = cnidi_df[(cnidi_df['HWE_p.value'] > 0.05) & (cnidi_df['allelic_chisq_p.value'] < 0.05)]
cnidi_genelist = cnidi_hwe_df['Gene.refGene'].to_list()
cnidi_genelist_dict = dict(Counter(cnidi_genelist))
cnidi_genegt5 = {}
for k in sorted(cnidi_genelist_dict,key=cnidi_genelist_dict.__getitem__, reverse=True):
    if cnidi_genelist_dict[k] < 5:
        break
    cnidi_genegt5[k] = cnidi_genelist_dict[k]


cnsym_df = pd.read_excel('./output/cntr-sym差异分析结果_filtbygene.xlsx', encoding='utf-8')
cnsym_hwe_df = cnsym_df[(cnsym_df['HWE_p.value'] > 0.05) & (cnsym_df['allelic_chisq_p.value'] < 0.05)]
cnsym_genelist = cnsym_hwe_df['Gene.refGene'].to_list()
cnsym_genelist_dict = dict(Counter(cnsym_genelist))
cnsym_genegt5 = {}
for k in sorted(cnsym_genelist_dict,key=cnsym_genelist_dict.__getitem__, reverse=True):
    if cnsym_genelist_dict[k] < 5:
        break
    cnsym_genegt5[k] = cnsym_genelist_dict[k]

for gn in cnep_genelist_dict:
    #     print(gn)
    if gn not in cnidi_genelist_dict:
        cnidi_genelist_dict[gn] = 0
    if gn not in cnsym_genelist_dict:
        cnsym_genelist_dict[gn] = 0

for gn in cnidi_genelist_dict:
    if gn not in cnep_genelist_dict:
        cnep_genelist_dict[gn] = 0
    if gn not in cnsym_genelist_dict:
        cnsym_genelist_dict[gn] = 0

for gn in cnsym_genelist_dict:
    if gn not in cnep_genelist_dict:
        cnep_genelist_dict[gn] = 0
        
    if gn not in cnidi_genelist_dict:
        cnidi_genelist_dict[gn] = 0

cnidi_genelist_ordereddict = OrderedDict()
cnsym_genelist_ordereddict = OrderedDict()
for gn in cnep_genelist_dict:
    cnidi_genelist_ordereddict[gn] = cnidi_genelist_dict[gn]
    cnsym_genelist_ordereddict[gn] = cnsym_genelist_dict[gn]

## part 3
## 突变类型，同义/非同义/移码等
def ext_vartype(epdf):
    epvdf = dict(Counter(epdf['ExonicFunc.refGene']))
    epvdf.update({'other': epvdf.pop('.')})
    return epvdf

cnep_vartype = ext_vartype(cnep_hwe_df)
cnidi_vartype = ext_vartype(cnidi_hwe_df)
cnsym_vartype = ext_vartype(cnsym_hwe_df)
for t in cnep_vartype:
    if t not in cnidi_vartype:
        cnidi_vartype[t] = 0
    if t not in cnsym_vartype:
        cnsym_vartype[t] = 0
        
cnidi_vartype_ordered = OrderedDict()
cnsym_vartype_ordered = OrderedDict()
for t in cnep_vartype:
    cnidi_vartype_ordered[t] = cnidi_vartype[t]
    cnsym_vartype_ordered[t] = cnsym_vartype[t]

## 突变区域，外显子/内含子/UTR等
def ext_region(epdf):
    epdf_region = dict(Counter(epdf['Func.refGene']))
    return epdf_region

cnep_region = ext_region(cnep_hwe_df)
cnidi_region = ext_region(cnidi_hwe_df)
cnsym_region = ext_region(cnsym_hwe_df)

cnidi_region_ordered = OrderedDict()
cnsym_region_ordered = OrderedDict()
for r in cnep_region:
    cnidi_region_ordered[r] = cnidi_region[r]
    cnsym_region_ordered[r] = cnsym_region[r]

## 千人基因组频率
def ext_1000gdict(ep_hwe_df):
    ep_1000gdict = {'.':0, '0 ~ 0.1':0, '0.1 ~ 0.3':0, '0.3 ~ 0.5':0, '0.5 ~ 1':0}
    ep_1000g_list = ep_hwe_df['X1000g2015aug_all'].tolist()
    for emb in ep_1000g_list:
        if emb == '.':
            ep_1000gdict['.'] += 1
        elif 0 < float(emb) < 0.1:
            ep_1000gdict['0 ~ 0.1'] += 1
        elif 0.1 <= float(emb) < 0.3:
            ep_1000gdict['0.1 ~ 0.3'] += 1
        elif 0.3 <= float(emb) < 0.5:
            ep_1000gdict['0.3 ~ 0.5'] += 1
        elif 0.5 <= float(emb) < 1:
            ep_1000gdict['0.5 ~ 1'] += 1
            
    return ep_1000gdict

cnep_1000g = ext_1000gdict(cnep_hwe_df)
cnidi_1000g = ext_1000gdict(cnidi_hwe_df)
cnsym_1000g = ext_1000gdict(cnsym_hwe_df)

def ext_1000g(epdf):
    ep_1000g_list = epdf['X1000g2015aug_all'].tolist()
    for i in range(len(ep_1000g_list)):
        if ep_1000g_list[i] == '.':
            ep_1000g_list[i] = 0
        else:
            ep_1000g_list[i] = float(ep_1000g_list[i])

    ep_1000g_sortedlist = sorted(ep_1000g_list)
    return ep_1000g_sortedlist

cnep_1000g_sortedlist = ext_1000g(cnep_hwe_df)
cnidi_1000g_sortedlist = ext_1000g(cnidi_hwe_df)
cnsym_1000g_sortedlist = ext_1000g(cnsym_hwe_df)

##---------------------------matplotlib.plot 1--------------------------------##
plt.close('all')

def autolabel(rects, ax):
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
            xy=(rect.get_x() +  rect.get_width()/2, height),
            xytext=(0,3),
            textcoords='offset points',
            ha='center', va='bottom', fontsize=10)

fig = plt.figure(figsize=(27,35))

ax1 = plt.subplot(411)
x1 = list(cnep_genegt5.keys())
y1 = list(cnep_genegt5.values())
rects1 = ax1.bar(x1, y1, color='#9DB63B')
plt.xticks(rotation=45, fontsize=14)
plt.xlabel('Gene', fontdict={'family':'Times New Roman','size':16})
plt.ylabel('Count', fontdict={'family':'Times New Roman','size':16})
plt.title('CN-EP-GeneCount', fontdict={'family':'Times New Roman','size':22})
autolabel(rects1, ax1)

ax2 = plt.subplot(412)
x2 = list(cnidi_genegt5.keys())
y2 = list(cnidi_genegt5.values())
rects2 = ax2.bar(x2, y2, color='#D92521')
plt.xticks(rotation=45, fontsize=14)
plt.xlabel('Gene', fontdict={'family':'Times New Roman','size':16})
plt.ylabel('Count', fontdict={'family':'Times New Roman','size':16})
plt.title('CN-IDI-GeneCount', fontdict={'family':'Times New Roman', 'size':22})
autolabel(rects2, ax2)

ax3 = plt.subplot(413)
x3 = list(cnsym_genegt5.keys())
y3 = list(cnsym_genegt5.values())
rects3 = ax3.bar(x3, y3, color='#4F81BD')
plt.xticks(rotation=45, fontsize=14)
plt.xlabel('Gene', fontdict={'family':'Times New Roman','size':16})
plt.ylabel('Count', fontdict={'family':'Times New Roman','size':16})
plt.title('CN-SYM-GeneCount', fontdict={'family':'Times New Roman', 'size':22})
autolabel(rects3, ax3)
        
plt.subplots_adjust(hspace=0.3)
# plt.show()
plt.savefig('gene3grp1.png', dpi=400, bbox_inches='tight')
##----------------------------------------------------------------------------##

##---------------------------matplotlib.plot 2--------------------------------## 
plt.close('all')
labels = list(cnep_genelist_dict.keys())
width = 0.24

fig = plt.figure(figsize=(22,14))
ax = plt.subplot(211)
x1 = np.arange(len(labels))
y1 = list(cnep_genelist_dict.values())
y2 = list(cnidi_genelist_ordereddict.values())
y3 = list(cnsym_genelist_ordereddict.values())

rects1 = ax.bar(x1-width, y1, width, color='#9DB63B',label='CN-EP')
rects2 = ax.bar(x1, y2, width, color='#D92521', label='CN-IDI')
rects3 = ax.bar(x1+width, y3, width, color='#4F81BD', label='CN-SYM')

plt.xticks(x1, labels, family='Times New Roman', rotation=90, fontsize=5)
ax.set_ylabel('Count', fontdict={'family':'Times New Roman','size':14})
ax.legend(loc='upper right', prop={'family':'Times New Roman', 'size':12})
plt.grid(axis='y', linewidth=1, linestyle='--')

ax2 = plt.subplot(212)
rects11 = ax2.plot(x1, y1, '-', color='#9DB63B', label='CN-EP')
rects22 = ax2.plot(x1, y2, '-.', color='#D92521', label='CN-IDI')
rects33 = ax2.plot(x1, y3, '--', color='#4F81BD', label='CN-SYM')
plt.xticks(x1, labels, family='Times New Roman', rotation=90, fontsize=5)
ax2.set_ylabel('Count', fontdict={'family':'Times New Roman','size':14})
ax2.legend(loc='upper right', prop={'family':'Times New Roman', 'size':12})
plt.grid(axis='y', linewidth=1, linestyle='--')

plt.suptitle('Count by GeneGroups', family='Times New Roman', fontsize=20, y=0.94)
# plt.show()
plt.savefig('gene3grp2.png', dpi=400, bbox_inches='tight')
##----------------------------------------------------------------------------##

##---------------------------matplotlib.plot 3--------------------------------##
plt.close('all')
fig, axes = plt.subplots(2,3,figsize=(40,20))
def autolabel(rects, ax):
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
            xy=(rect.get_x() +  rect.get_width()/2, height),
            xytext=(0,3),
            textcoords='offset points',
            ha='center', va='bottom', fontsize=10)
## axes[0,0]
labels1 = list(cnep_vartype.keys())
width = 0.24 

x1 = np.arange(len(labels1))
y1 = list(cnep_vartype.values())
y2 = list(cnidi_vartype_ordered.values())
y3 = list(cnsym_vartype_ordered.values())

rects1 = axes[0,0].bar(x1-width, y1, width, color='#9DB63B', label='CN-EP')
rects2 = axes[0,0].bar(x1, y2, width, color='#D92521', label='CN-IDI')
rects3 = axes[0,0].bar(x1+width, y3, width, color='#4F81BD', label='CN-SYM')

axes[0,0].set_xticks(x1)
axes[0,0].set_xticklabels(labels1, fontdict={'rotation':15}, fontsize=15, family='Times New Roman')
axes[0,0].set_yticklabels(np.arange(0,450,50), fontdict={'family':'Times New Roman','size':15})
axes[0,0].set_xlabel('Variant Type', fontdict={'family':'Times New Roman','size':18})
axes[0,0].set_ylabel('Count', fontdict={'family':'Times New Roman','size':18})

autolabel(rects1, axes[0,0])
autolabel(rects2, axes[0,0])
autolabel(rects3, axes[0,0])
axes[0,0].legend(loc='upper left', prop={'family':'Times New Roman','size':14})

## axes[0,1]
labels2 = list(cnep_region.keys())

x2 = np.arange(len(labels2))
z1 = list(cnep_region.values())
z2 = list(cnidi_region_ordered.values())
z3 = list(cnsym_region_ordered.values())

rects4 = axes[0,1].bar(x2-width, z1, width, color='#9DB63B', label='CN-EP')
rects5 = axes[0,1].bar(x2, z2, width, color='#D92521', label='CN-IDI')
rects6 = axes[0,1].bar(x2+width, z3, width, color='#4F81BD', label='CN-SYM')
axes[0,1].set_xticks(x2)
axes[0,1].set_xticklabels(labels2, fontdict={'rotation':15}, fontsize=15, family='Times New Roman')
axes[0,1].set_yticklabels(np.arange(0,450,50), fontdict={'family':'Times New Roman','size':15})
axes[0,1].set_xlabel('Transport Region', fontdict={'family':'Times New Roman','size':18}, labelpad=30)
axes[0,1].set_ylabel('Count', fontdict={'family':'Times New Roman','size':18})

autolabel(rects4, axes[0,1])
autolabel(rects5, axes[0,1])
autolabel(rects6, axes[0,1])
axes[0,1].legend(loc='upper right', prop={'family':'Times New Roman','size':14})

## axes[0,2]
labels3 = list(cnep_1000g.keys())
width = 0.24

x3 = np.arange(len(labels3))
ft1 = cnep_1000g.values()
ft2 = cnidi_1000g.values()
ft3 = cnsym_1000g.values()

rects7 = axes[0,2].bar(x3-width, ft1, width, color='#9DB63B', label='CN-EP')
rects8 = axes[0,2].bar(x3, ft2, width, color='#D92521', label='CN-IDI')
rects9 = axes[0,2].bar(x3+width, ft3, width, color='#4F81BD', label='CN-SYM')

axes[0,2].set_xticks(x3)
axes[0,2].set_xticklabels(labels3, fontdict={'rotation':15}, fontsize=15, family='Times New Roman')
axes[0,2].set_yticklabels(np.arange(0,350,50), fontdict={'family':'Times New Roman','size':15})
axes[0,2].set_xlabel('1000G FreqGroup', fontdict={'family':'Times New Roman','size':18}, labelpad=30) 
axes[0,2].set_ylabel('Count', fontdict={'family':'Times New Roman','size':18})

autolabel(rects7, axes[0,2])
autolabel(rects8, axes[0,2])
autolabel(rects9, axes[0,2])
axes[0,2].legend(loc='upper right', prop={'family':'Times New Roman','size':14})


## axes[1,0]
outer_colors = [(0,0,0.5), (0.6,0.42,0), (0,0.25,0), (0.5,0,0), (0.5,0,0.5)]
middle_colors = [(0,0,0.7), (0.8,0.52,0), (0,0.35,0), (0.7,0,0), (0.7,0,0.7)]
inner_colors = [(0,0,0.9), (1,0.65,0), (0,0.45,0), (0.9,0,0), (0.9,0,0.9)]
outer_colors2 = [(0,0,0.5), (0.6,0.42,0), (0,0.25,0), (0.5,0,0), (0.5,0,0.5), (0.5,0.5,0), (145/255,24/255,24/255)]
middle_colors2 = [(0,0,0.7), (0.8,0.52,0), (0,0.35,0), (0.7,0,0), (0.7,0,0.7), (0.7,0.7,0), (165/255,42/255,42/255)]
inner_colors2 = [(0,0,0.9), (1,0.65,0), (0,0.45,0), (0.9,0,0), (0.9,0,0.9), (0.9,0.9,0), (185/255,60/255,60/255)]

size = 0.25
x11 = list(cnep_vartype.values())
x12 = list(cnidi_vartype_ordered.values())
x13 = list(cnsym_vartype_ordered.values())
l11 = list(cnep_vartype.keys())
l12 = list(cnidi_vartype_ordered.keys())
l13 = list(cnsym_vartype_ordered.keys())

axes[1,0].pie(x11, labels=l11, radius=1, labeldistance=None, colors=outer_colors, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,0].pie(x12, labels=l12, radius=0.75, labeldistance=None, colors=middle_colors, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,0].pie(x13, labels=l13, radius=0.5, labeldistance=None, colors= inner_colors, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,0].legend(labels=l13, loc=(0.8,0.8), prop={'family':'Times New Roman', 'size':14})

## axes[1,1]
x21 = list(cnep_region.values())
x22 = list(cnidi_region_ordered.values())
x23 = list(cnsym_region_ordered.values())
l21 = list(cnep_region.keys())
l22 = list(cnidi_region_ordered.keys())
l23 = list(cnsym_region_ordered.keys())

axes[1,1].pie(x21, labels=l21, radius=1, labeldistance=None, colors=outer_colors2, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,1].pie(x22, labels=l22, radius=0.75, labeldistance=None, colors=middle_colors2, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,1].pie(x23, labels=l23, radius=0.5, labeldistance=None, colors=inner_colors2, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,1].legend(labels=l23, loc=(0.8,0.8), prop={'family':'Times New Roman', 'size':14})

## axes[1,2]
x31 = list(cnep_1000g.values())
x32 = list(cnidi_1000g.values())
x33 = list(cnsym_1000g.values())
l31 = list(cnep_1000g.keys())
l32 = list(cnidi_1000g.keys())
l33 = list(cnsym_1000g.keys())

axes[1,2].pie(x31, labels=l31, radius=1, labeldistance=None, colors=outer_colors, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,2].pie(x32, labels=l32, radius=0.75, labeldistance=None, colors=middle_colors, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,2].pie(x33, labels=l33, radius=0.5, labeldistance=None, colors= inner_colors, wedgeprops=dict(width=size, edgecolor='w'))
axes[1,2].legend(labels=l33, loc=(0.8,0.8), prop={'family':'Times New Roman', 'size':14})

plt.suptitle('Other Characteristic Statistics', fontsize=30, fontdict={'family':'Times New Roman'}, y=0.92)
# plt.show()
plt.savefig('gene3grp4.png', dpi=400, bbox_inches='tight')
##----------------------------------------------------------------------------## 