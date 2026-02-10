import os,shutil
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import special
import sympy
import decimal
import mpmath
import statistics
import array
import scipy
from scipy.optimize import fmin_l_bfgs_b
from scipy.stats import norm
from scipy.stats import exponnorm
from mpmath import mp

def Generator(number):

    """
    Generate all 3-combinations (unordered triples) from {1, 2, ..., number}.

    Returns
    -------
    list[list[int]]
        A list of triples [i, j, k] with 1 <= i < j < k <= number.
        The output size is C(number, 3).
    """

    Threecomb=list()
    for i in range(number):
        for j in range(number):
            for k in range(number):
                if (i<j<k):
                    Threecomb.append([i+1,j+1,k+1])
                else:
                    continue
    return Threecomb


def Re_sorted(w_count):

    # Sort the inner dictionary of each w_count entry by value in descending order.
    nw=dict()
    for i in w_count:
        nw[i] = dict()
        for j in range(len(w_count[i])):
            nw[i][str(sorted(w_count[i].items(),key=lambda d:d[1],reverse=True)[j][0])]=sorted(w_count[i].items(),key=lambda d:d[1],reverse=True)[j][1]


    return nw



# User-defined: provide the path to the input *.withoutweight_num.tre file.

path1='E:/Stefan/2025.5.29/TriMouNet/Tnet_script/'
f1=open(path1+'gamma_test_6.704_0.14916467780429596_M_quartet_adjust_deleteoutlier_modify_withoutweight_num.tre','r+',encoding='utf-8')


# User-defined: specify the list of taxa to be analyzed and the outgroup taxon here.
w_count=dict()
temp=0
taxa={'1':'A','2':'B','3':'C','4':'D','5':'E','6':'F','7':'G','8':'H','9':'O'}

# There are C(8,3)=56 three-taxon subsets; initialize w_count for each subset with the three corresponding triplet keys.
for i in range(56):
    w_count[str(i+1)]=dict()
    w_count[str(i + 1)][str(Generator(8)[temp][0]) +'_'+ str(Generator(8)[temp][1])+'_' + '|' + str(Generator(8)[temp][2])+'_'] = 0
    w_count[str(i + 1)][str(Generator(8)[temp][0]) +'_'+ str(Generator(8)[temp][2])+'_' + '|' + str(Generator(8)[temp][1])+'_']  = 0
    w_count[str(i + 1)][str(Generator(8)[temp][1]) +'_'+ str(Generator(8)[temp][2])+'_' + '|' + str(Generator(8)[temp][0])+'_'] = 0
    temp+=1



lines=f1.readlines()
temp=0
wtemp=1

# Compute w_count by scanning the tree file: for each ingroup triple, count how many trees support each of the three rooted triplet topologies.


for line in lines:
    if str(Generator(8)[temp][0]) in line and str(Generator(8)[temp][1]) in line and str(
            Generator(8)[temp][2]) in line:
        if '('+'('+ str(Generator(8)[temp][0]) + ',' + str(Generator(8)[temp][1]) + '),' + str(
                Generator(8)[temp][2]) in line or '('+'('+ str(Generator(8)[temp][1]) + ',' + str(Generator(8)[temp][0]) + '),' + str(
                Generator(8)[temp][2])  in line or '(' + str(
                Generator(8)[temp][2]) + ',(' + str(Generator(8)[temp][0]) + ',' + str(
                Generator(8)[temp][1]) + '))' in line or '(' + str(Generator(8)[temp][2]) + ',(' + str(
                Generator(8)[temp][1]) + ',' + str(Generator(8)[temp][0]) + '))' in line:
            w_count[str(wtemp)][str(Generator(8)[temp][0]) + '_' + str(Generator(8)[temp][1]) + '_' + '|' + str(
                Generator(8)[temp][2]) + '_'] += 1

        elif '('+'('+ str(Generator(8)[temp][0]) + ',' + str(Generator(8)[temp][2]) + '),' + str(
                Generator(8)[temp][1]) in line or '('+'('+ str(Generator(8)[temp][2]) + ',' + str(Generator(8)[temp][0]) + '),' + str(
                Generator(8)[temp][1])  in line or '(' + str(
                Generator(8)[temp][1]) + ',(' + str(Generator(8)[temp][0]) + ',' + str(
                Generator(8)[temp][2]) + '))' in line or '(' + str(Generator(8)[temp][1]) + ',(' + str(
                Generator(8)[temp][2]) + ',' + str(Generator(8)[temp][0]) + '))' in line:
            w_count[str(wtemp)][str(Generator(8)[temp][0]) + '_' + str(Generator(8)[temp][2]) + '_' + '|' + str(
                Generator(8)[temp][1]) + '_'] += 1

        elif '('+'('+ str(Generator(8)[temp][1]) + ',' + str(Generator(8)[temp][2]) + '),' + str(
                Generator(8)[temp][0]) in line or '('+'('+ str(Generator(8)[temp][2]) + ',' + str(Generator(8)[temp][1]) + '),' + str(
                Generator(8)[temp][0])  in line or '(' + str(
                Generator(8)[temp][0]) + ',(' + str(Generator(8)[temp][1]) + ',' + str(
                Generator(8)[temp][2]) + '))' in line or '(' + str(Generator(8)[temp][0]) + ',(' + str(
                Generator(8)[temp][2]) + ',' + str(Generator(8)[temp][1]) + '))' in line:
            w_count[str(wtemp)][str(Generator(8)[temp][1]) + '_' + str(Generator(8)[temp][2]) + '_' + '|' + str(
                Generator(8)[temp][0]) + '_'] += 1

    else:
        wtemp += 1
        temp += 1
        if str(Generator(8)[temp][0]) in line and str(Generator(8)[temp][1]) in line and str(
                Generator(8)[temp][2]) in line:
            if '(' + '(' + str(Generator(8)[temp][0]) + ',' + str(Generator(8)[temp][1]) + '),' + str(
                    Generator(8)[temp][2]) in line or '(' + '(' + str(Generator(8)[temp][1]) + ',' + str(
                Generator(8)[temp][0]) + '),' + str(
                Generator(8)[temp][2]) in line or '(' + str(
                Generator(8)[temp][2]) + ',(' + str(Generator(8)[temp][0]) + ',' + str(
                Generator(8)[temp][1]) + '))' in line or '(' + str(Generator(8)[temp][2]) + ',(' + str(
                Generator(8)[temp][1]) + ',' + str(Generator(8)[temp][0]) + '))' in line:
                w_count[str(wtemp)][str(Generator(8)[temp][0]) + '_' + str(Generator(8)[temp][1]) + '_' + '|' + str(
                    Generator(8)[temp][2]) + '_'] += 1

            elif '(' + '(' + str(Generator(8)[temp][0]) + ',' + str(Generator(8)[temp][2]) + '),' + str(
                    Generator(8)[temp][1]) in line or '(' + '(' + str(Generator(8)[temp][2]) + ',' + str(
                Generator(8)[temp][0]) + '),' + str(
                Generator(8)[temp][1]) in line or '(' + str(
                Generator(8)[temp][1]) + ',(' + str(Generator(8)[temp][0]) + ',' + str(
                Generator(8)[temp][2]) + '))' in line or '(' + str(Generator(8)[temp][1]) + ',(' + str(
                Generator(8)[temp][2]) + ',' + str(Generator(8)[temp][0]) + '))' in line:
                w_count[str(wtemp)][str(Generator(8)[temp][0]) + '_' + str(Generator(8)[temp][2]) + '_' + '|' + str(
                    Generator(8)[temp][1]) + '_'] += 1

            elif '(' + '(' + str(Generator(8)[temp][1]) + ',' + str(Generator(8)[temp][2]) + '),' + str(
                    Generator(8)[temp][0]) in line or '(' + '(' + str(Generator(8)[temp][2]) + ',' + str(
                Generator(8)[temp][1]) + '),' + str(
                Generator(8)[temp][0]) in line or '(' + str(
                Generator(8)[temp][0]) + ',(' + str(Generator(8)[temp][1]) + ',' + str(
                Generator(8)[temp][2]) + '))' in line or '(' + str(Generator(8)[temp][0]) + ',(' + str(
                Generator(8)[temp][2]) + ',' + str(Generator(8)[temp][1]) + '))' in line:
                w_count[str(wtemp)][str(Generator(8)[temp][1]) + '_' + str(Generator(8)[temp][2]) + '_' + '|' + str(
                    Generator(8)[temp][0]) + '_'] += 1

print("The count of split are:",w_count)







w_count=Re_sorted(w_count)
print("The count of split after sorting are:",w_count)



second=list()
third=list()
first=list()

Un=list()
sigma_k_index=list()
sigma_kc_index=list()

bio_p=list()
bio_p_first=list()


for key, sub_dict in w_count.items():
    # 按值排序，得到排名第二和第三的键值对
    sorted_items = sorted(sub_dict.items(), key=lambda x: x[1], reverse=True)
    second.append(sorted_items[1][1])  # 排名第二
    third.append(sorted_items[2][1])
    first.append(sorted_items[0][1])




# Assign indices to sigma_kc and sigma_k based on binomial-test p-values (Bonferroni-corrected),
# and record the corresponding p-values (pf) for the 2nd-vs-3rd and 1st-vs-2nd counts.

for i in range(len(second)):
    if second[i]==0:
        continue
    else:
        bio_p.append(scipy.stats.binomtest(second[i], second[i]+third[i], p=0.5, alternative='two-sided').pvalue)
        if scipy.stats.binomtest(second[i], second[i]+third[i], p=0.5, alternative='two-sided').pvalue<0.05/56:
            sigma_kc_index.append(str(i+1))

for i in range(len(first)):
    if first[i]==0:
        continue
    else:
        bio_p_first.append(scipy.stats.binomtest(first[i], first[i]+second[i], p=0.5, alternative='two-sided').pvalue)


for i in range(56):
    Un.append(str(i+1))

sigma_k_index = list(set(Un) - set(sigma_kc_index))
bio_p_second_value=dict()
bio_p_first_value=dict()



for i in range(len(bio_p)):
    bio_p_second_value[str(i+1)]=bio_p[i]
    bio_p_first_value[str(i+1)]=bio_p_first[i]

print('This is bio_p_second_value\n',bio_p_second_value)
print('This is bio_p_first_value\n',bio_p_first_value)

sigma_k_index = sorted(sigma_k_index, key=int)
sigma_kc_index = sorted(sigma_kc_index, key=int)

print('This is sigma_k_index\n',sigma_k_index)
print("This is sigma_kc_index\n",sigma_kc_index)




f2=open(path1+'gamma_test_6.704_0.14916467780429596_M_quartet_adjust_deleteoutlier_modify.tre','r+')
lines=f2.readlines()

w_index=dict()
w_length=dict()
for i in w_count:
    w_index[i]=list()
    for j in w_count[i]:
        w_index[i].append(j)
    w_length[i]=0
    for j in w_count[i]:
        w_length[i]+=w_count[i][j]

w_page=dict()
w_page['1'] = 0
for i in range(56):
    w_page[str(i+2)]=0
    for j in range(i+1):
        w_page[str(i+2)]+=w_length[str(j+1)]

Un=list()


second=list()
third=list()
first=list()



DA_k=dict()
DB_k=dict()
DC_k=dict()
DO_k=dict()
Di_k=dict()
DXY_k=dict()
DXZ_k=dict()
DYZ_k=dict()
DXO_k=dict()
DYO_k=dict()
DZO_k=dict()

DXY_O_k=dict()
DXZ_O_k=dict()
DYZ_O_k=dict()
DZ_ZO_k=dict()

DADB_CACB_k=dict()
DC_CACB_k=dict()
DADB_OAOB_k=dict()
DO_OAOB_k=dict()
DC_CO_k=dict()

DB_kc_first=dict()
Di_kc_first=dict()
DO_kc_first=dict()
DB_OB_kc_first=dict()
DB_OB_kc_avg_first=dict()

DB_kc_second=dict()
Di_kc_second=dict()
DO_kc_second=dict()
DB_OB_kc_second=dict()
DB_OB_kc_avg_second=dict()

DB_kc_third=dict()
Di_kc_third=dict()
DO_kc_third=dict()
DB_OB_kc_third=dict()
DB_OB_kc_avg_third=dict()

DB_OB_kc=dict()





dab=list()
dac=list()
dbc=list()

# For each index in sigma_k, parse branch lengths from the corresponding quartet trees and compute distance-based ratios
# (e.g., DXY_O_k as a cherry edge ratio (V1 in the TriMouNet) and DXZ_O_k (V2 in the TriMouNet) as a pending ratio).
for i in sigma_k_index:
    DA_k[str(i)] = list()
    DB_k[str(i)] = list()
    DC_k[str(i)] = list()
    DO_k[str(i)] = list()
    Di_k[str(i)] = list()
    DXY_k[str(i)] = list()
    DXZ_k[str(i)] = list()
    DYZ_k[str(i)] = list()
    DXO_k[str(i)] = list()
    DYO_k[str(i)] = list()
    DZO_k[str(i)] = list()
    for k in range(w_page[i], w_page[str(int(i) + 1)]):
        for j in range(len(lines[k].strip(';\n').split(':'))):
            if taxa[str(w_index[str(i)][0].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[j] and '(' in lines[k].strip(';\n').split(':')[j - 1] and '(' in lines[k].strip(';\n').split(':')[j + 1] or taxa[str(w_index[str(i)][0].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[j] and ')' in lines[k].strip(';\n').split(':')[j - 1] and ')' in lines[k].strip(';\n').split(':')[j + 1]:
                for m in lines[k].strip(';\n').split(','):
                    if 'O' in m:
                        DO_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[i][0].split('_')[2].split('|')[1])] in m:
                        DC_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][0].split('_')[0])] in m:
                        DA_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][0].split('_')[1])] in m:
                        DB_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if '):' in m:
                        Di_k[str(i)].append(float(m.split(':')[2].strip(')').strip(')')))
                DXY_k[str(i)].append(DA_k[str(i)][-1] + DB_k[str(i)][-1])
                DXZ_k[str(i)].append(DA_k[str(i)][-1] + Di_k[str(i)][-1] + DC_k[str(i)][-1])
                DYZ_k[str(i)].append(DB_k[str(i)][-1] + Di_k[str(i)][-1] + DC_k[str(i)][-1])
                DXO_k[str(i)].append(DA_k[str(i)][-1] + Di_k[str(i)][-1] + DO_k[str(i)][-1])
                DYO_k[str(i)].append(DB_k[str(i)][-1] + Di_k[str(i)][-1] + DO_k[str(i)][-1])
                DZO_k[str(i)].append(DC_k[str(i)][-1] + DO_k[str(i)][-1])


            elif taxa[str(w_index[str(i)][1].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[
                j] and '(' in lines[k].strip(';\n').split(':')[j - 1] and '(' in lines[k].strip(';\n').split(':')[
                j + 1] or taxa[str(w_index[str(i)][1].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[
                j] and ')' in lines[k].strip(';\n').split(':')[j - 1] and ')' in lines[k].strip(';\n').split(':')[
                j + 1]:
                for m in lines[k].strip(';\n').split(','):
                    if 'O' in m:
                        DO_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[i][1].split('_')[2].split('|')[1])] in m:
                        DB_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][1].split('_')[0])] in m:
                        DA_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][1].split('_')[1])] in m:
                        DC_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if '):' in m:
                        Di_k[str(i)].append(float(m.split(':')[2].strip(')').strip(')')))
                DXY_k[str(i)].append(DA_k[str(i)][-1] + Di_k[str(i)][-1] + DB_k[str(i)][-1])
                DXZ_k[str(i)].append(DA_k[str(i)][-1] + DC_k[str(i)][-1])
                DYZ_k[str(i)].append(DB_k[str(i)][-1] + Di_k[str(i)][-1] + DC_k[str(i)][-1])
                DXO_k[str(i)].append(DA_k[str(i)][-1] + Di_k[str(i)][-1] + DO_k[str(i)][-1])
                DYO_k[str(i)].append(DB_k[str(i)][-1] + DO_k[str(i)][-1])
                DZO_k[str(i)].append(DC_k[str(i)][-1] + Di_k[str(i)][-1] + DO_k[str(i)][-1])


            elif taxa[str(w_index[str(i)][2].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[
                j] and '(' in lines[k].strip(';\n').split(':')[j - 1] and '(' in lines[k].strip(';\n').split(':')[
                j + 1] or taxa[str(w_index[str(i)][2].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[
                j] and ')' in lines[k].strip(';\n').split(':')[j - 1] and ')' in lines[k].strip(';\n').split(':')[
                j + 1]:
                for m in lines[k].strip(';\n').split(','):
                    if 'O' in m:
                        DO_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[i][2].split('_')[2].split('|')[1])] in m:
                        DA_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][2].split('_')[0])] in m:
                        DB_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][2].split('_')[1])] in m:
                        DC_k[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if '):' in m:
                        Di_k[str(i)].append(float(m.split(':')[2].strip(')').strip(')')))
                DXY_k[str(i)].append(DA_k[str(i)][-1] + Di_k[str(i)][-1] + DB_k[str(i)][-1])
                DXZ_k[str(i)].append(DA_k[str(i)][-1] + Di_k[str(i)][-1] + DC_k[str(i)][-1])
                DYZ_k[str(i)].append(DB_k[str(i)][-1] + DC_k[str(i)][-1])
                DXO_k[str(i)].append(DA_k[str(i)][-1] + DO_k[str(i)][-1])
                DYO_k[str(i)].append(DB_k[str(i)][-1] + Di_k[str(i)][-1] + DO_k[str(i)][-1])
                DZO_k[str(i)].append(DC_k[str(i)][-1] + Di_k[str(i)][-1] + DO_k[str(i)][-1])
            else:
                continue
for i in sigma_k_index:
    DXY_O_k[i] = list()
    DXZ_O_k[i] = list()
    DYZ_O_k[i] = list()
    DZ_ZO_k[i] = list()
    for j in range(w_page[str(int(i) + 1)] - w_page[i]):
        DXY_O_k[i].append(float(DXY_k[i][j]) / (float((DXO_k[i][j]))  + float(DYO_k[i][j])+float(DXY_k[i][j])))
        #DXZ_O_k[i].append(0.5*(float(DXZ_k[i][j])+float(DZO_k[i][j])-float(DXO_k[i][j])) / (float(DXO_k[i][j])+float(DZO_k[i][j])+float(DXZ_k[i][j]))  +  0.5*(float(DYZ_k[i][j])+float(DZO_k[i][j])-float(DYO_k[i][j])) / (float(DYO_k[i][j])+float(DZO_k[i][j])+float(DYZ_k[i][j])))
        DXZ_O_k[i].append(0.5 * float(DXZ_k[i][j])/ (float(DXO_k[i][j]) + float(DZO_k[i][j]) + float(DXZ_k[i][j])) + 0.5 * float(DYZ_k[i][j]) / (float(DYO_k[i][j]) + float(DZO_k[i][j]) + float(DYZ_k[i][j])))

        DYZ_O_k[i].append(float(DYZ_k[i][j]) / float((DYO_k[i][j])))
        DZ_ZO_k[i].append((float(DC_k[i][j])+0.5*float(Di_k[i][j]))/float(DZO_k[i][j]))

# For each sigma_kc index, collect branch lengths supporting the first split (w_index[i][0]),
# then compute the common edge ratio DB_OB_kc (V3 in the TriMouNet) relative to the outgroup O and average it across trees.

for i in sigma_kc_index:
    DB_kc_first[str(i)]=list()
    DO_kc_first[str(i)]=list()
    Di_kc_first[str(i)]=list()
    for k in range(w_page[i],w_page[str(int(i)+1)]):
        for j in range(len(lines[k].strip(';\n').split(':'))):
            if taxa[str(w_index[str(i)][0].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[j] and '(' in lines[k].strip(';\n').split(':')[j-1] and '(' in  lines[k].strip(';\n').split(':')[j+1] or taxa[str(w_index[str(i)][0].split('_')[2].split('|')[1])]  in lines[k].strip(';\n').split(':')[j] and ')' in lines[k].strip(';\n').split(':')[j-1] and ')' in  lines[k].strip(';\n').split(':')[j+1]:
                for m in lines[k].strip(';\n').split(','):
                    if 'O' in m:
                        DO_kc_first[str(i)].append(float(m.split(':')[1].strip(')').strip(')').strip('(').strip('(')))
                    if taxa[str(w_index[str(i)][2].split('_')[2].split('|')[1])] in m:
                        DB_kc_first[str(i)].append(float(m.split(':')[1].strip(')').strip(')').strip('(').strip('(')))
                    if '):' in m:
                        Di_kc_first[str(i)].append(float(m.split(':')[2].strip(')').strip(')').strip('(').strip('(')))
            else:
                continue
for i in sigma_kc_index:
    DB_OB_kc_first[str(i)]=list()
    for j in range(len(DO_kc_first[i])):
        DB_OB_kc_first[str(i)].append(float(DB_kc_first[str(i)][j]/(DB_kc_first[str(i)][j]+Di_kc_first[str(i)][j]+DO_kc_first[str(i)][j])))
    DB_OB_kc_avg_first[i] = sum(DB_OB_kc_first[i]) / len(DB_OB_kc_first[i])

# For each sigma_kc index, collect branch lengths supporting the second split (w_index[i][1]),
# then compute the common edge ratio DB_OB_kc (V3 in the TriMouNet) relative to the outgroup O and average it across trees.
for i in sigma_kc_index:
    DB_kc_second[str(i)]=list()
    DO_kc_second[str(i)]=list()
    Di_kc_second[str(i)]=list()
    for k in range(w_page[i],w_page[str(int(i)+1)]):
        for j in range(len(lines[k].strip(';\n').split(':'))):
            if taxa[str(w_index[str(i)][1].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[j] and '(' in lines[k].strip(';\n').split(':')[j-1] and '(' in  lines[k].strip(';\n').split(':')[j+1] or taxa[str(w_index[str(i)][1].split('_')[2].split('|')[1])]  in lines[k].strip(';\n').split(':')[j] and ')' in lines[k].strip(';\n').split(':')[j-1] and ')' in  lines[k].strip(';\n').split(':')[j+1]:
                for m in lines[k].strip(';\n').split(','):
                    if 'O' in m:
                        DO_kc_second[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][2].split('_')[2].split('|')[1])] in m:
                        DB_kc_second[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if '):' in m:
                        Di_kc_second[str(i)].append(float(m.split(':')[2].strip(')').strip(')')))
            else:
                continue
for i in sigma_kc_index:
    DB_OB_kc_second[str(i)]=list()
    for j in range(len(DO_kc_second[i])):
        DB_OB_kc_second[str(i)].append(float(DB_kc_second[str(i)][j]/(DB_kc_second[str(i)][j]+Di_kc_second[str(i)][j]+DO_kc_second[str(i)][j])))
    DB_OB_kc_avg_second[i] = sum(DB_OB_kc_second[i]) / len(DB_OB_kc_second[i])

# For each sigma_kc index, collect branch lengths supporting the third split (w_index[i][2]),
# then compute the common edge ratio DB_OB_kc (V3 in the TriMouNet) relative to the outgroup O and average it across trees.
for i in sigma_kc_index:
    DB_kc_third[str(i)]=list()
    DO_kc_third[str(i)]=list()
    Di_kc_third[str(i)]=list()
    for k in range(w_page[i],w_page[str(int(i)+1)]):
        for j in range(len(lines[k].strip(';\n').split(':'))):
            if taxa[str(w_index[str(i)][2].split('_')[2].split('|')[1])] in lines[k].strip(';\n').split(':')[j] and '(' in lines[k].strip(';\n').split(':')[j-1] and '(' in  lines[k].strip(';\n').split(':')[j+1] or taxa[str(w_index[str(i)][2].split('_')[2].split('|')[1])]  in lines[k].strip(';\n').split(':')[j] and ')' in lines[k].strip(';\n').split(':')[j-1] and ')' in  lines[k].strip(';\n').split(':')[j+1]:
                for m in lines[k].strip(';\n').split(','):
                    if 'O' in m:
                        DO_kc_third[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if taxa[str(w_index[str(i)][2].split('_')[2].split('|')[1])] in m:
                        DB_kc_third[str(i)].append(float(m.split(':')[1].strip(')').strip(')')))
                    if '):' in m:
                        Di_kc_third[str(i)].append(float(m.split(':')[2].strip(')').strip(')')))
            else:
                continue
for i in sigma_kc_index:
    DB_OB_kc_third[str(i)]=list()
    DB_OB_kc[str(i)] = list()
    for j in range(len(DO_kc_third[i])):
        DB_OB_kc_third[str(i)].append(float(DB_kc_third[str(i)][j]/(DB_kc_third[str(i)][j]+Di_kc_third[str(i)][j]+DO_kc_third[str(i)][j])))
    DB_OB_kc_avg_third[i] = sum(DB_OB_kc_third[i]) / len(DB_OB_kc_third[i])

    DB_OB_kc[str(i)] = DB_OB_kc_first[str(i)] + DB_OB_kc_second[str(i)]



dab = np.array(dab)
dac = np.array(dac)
dbc = np.array(dbc)

log_value_pending_kc_m=dict()
pams_pending_kc_m=dict()
log_value_pending_kc_1=dict()
pams_pending_kc_1=dict()
log_value_pending_kc_2=dict()
pams_pending_kc_2=dict()


# ---------------------------------------------------------------------
# LV3 objective computation for triples in Sigma_kC (sigma_kc_index)
#
# For each triple t = {x, y, z} in Sigma_kC, we evaluate two competing models:
#   LV3,1: fit two separate EMG (exGaussian) distributions to two groups of V3 values
#          (e.g., corresponding to the two candidate sister-pair edges), and sum their NLLs.
#   LV3,2: fit a single EMG distribution to the combined V3 values (merged group) and take its NLL.
#
# Later (not shown here), we assign a trinet type by comparing the two objective values:
#   if |LV3,1 - LV3,2| <= kappa: the common edge is consistent with a unimodal distribution -> S1(x, y; z)
#   else: evidence for a mixture / bimodality -> S2(x; z; y)
#
# In the code below:
#   - log_value_pending_kc_m[m] stores the NLL for the combined-group fit  (used for LV3,2)
#   - log_value_pending_kc_1[m] stores the NLL for the first subgroup fit
#   - log_value_pending_kc_2[m] stores the NLL for the second subgroup fit
#   - LV3,1 is obtained by log_value_pending_kc_1[m] + log_value_pending_kc_2[m]
#
# Each EMG fit is performed by L-BFGS-B on the negative log-likelihood, using a two-stage strategy:
#   Stage 1: loose tolerances for a quick fit
#   Stage 2: if the projected gradient is still large, restart with stricter tolerances
# ---------------------------------------------------------------------

for m in sigma_kc_index:

    log_value_pending_kc_m[m] = float(0)
    pams_pending_kc_m[m] = list()

    data = DB_OB_kc[m]

    # -----------------------------------
    # 定义负对数似然函数
    def likelihood(params):
        u1, sigma, lm = params
        Lh = 0.0
        for x in data:
            pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
            if pdf_val > 0:
                Lh -= math.log(pdf_val)
            else:
                Lh += 1000
        return Lh

    # -----------------------------------
    # 数值梯度函数
    def likelihood_grad(params, eps=1e-8):
        grad = np.zeros(3)
        for i in range(3):
            p1 = params.copy()
            p2 = params.copy()
            p1[i] += eps
            p2[i] -= eps
            f1 = likelihood(p1)
            f2 = likelihood(p2)
            grad[i] = (f1 - f2) / (2 * eps)
        return grad

    # -----------------------------------
    # 初始参数和边界
    initial_params = np.array([0.5*np.mean(data), 0.5*np.var(data), 20])
    bounds = [(0.001, 1.2*np.mean(data)), (0.001,1), (0.001,1000)]

    # -----------------------------------
    # Stage 1 优化（宽松设置）
    x1, f1, info1 = fmin_l_bfgs_b(
        func=likelihood,
        x0=initial_params,
        fprime=likelihood_grad,
        bounds=bounds,
        factr=1e7,        # 宽松收敛
        pgtol=1e-4,
        maxfun=5000,
        disp=1
    )

    proj_g1 = np.linalg.norm(likelihood_grad(x1))
    print(f"[Stage1] f = {f1}, projected_grad = {proj_g1}, iters = {info1['nit']}")

    # -----------------------------------
    # 如果梯度太大，触发 Stage2 严格优化
    if proj_g1 > 1.0:
        print(" -> Triggering stage2 restart (stricter settings)...")
        x2, f2, info2 = fmin_l_bfgs_b(
            func=likelihood,
            x0=x1,
            fprime=likelihood_grad,
            bounds=bounds,
            factr=1e5,        # stricter
            pgtol=1e-8,
            maxfun=10000,
            disp=1
        )
        optimized_params = x2
    else:
        optimized_params = x1

    u1, sigma, lm = optimized_params

    # -----------------------------------
    # 打印参数
    print('u1=' + str(format(u1, '.3f')) + ',' +
          'sigma=' + str(format(sigma, '.3f')) + ',' +
          'lm=' + str(format(lm, '.3f')))

    pams_pending_kc_m[m].extend([u1, sigma, lm])

    # 计算最终的 log-likelihood
    Lh = 0.0
    for x in data:
        pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
        if pdf_val > 0:
            Lh -= math.log(pdf_val)
        else:
            Lh += 1000

    print('The log_likelihood in this loop is:', Lh)
    log_value_pending_kc_m[m] = format(float(Lh), '.3f')
for m in sigma_kc_index:

    log_value_pending_kc_1[m] = float(0)
    pams_pending_kc_1[m] = list()

    data = DB_OB_kc_first[m]

    # -----------------------------------
    # 定义负对数似然函数
    def likelihood(params):
        u1, sigma, lm = params
        Lh = 0.0
        for x in data:
            pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
            if pdf_val > 0:
                Lh -= math.log(pdf_val)
            else:
                Lh += 1000
        return Lh

    # -----------------------------------
    # 数值梯度函数
    def likelihood_grad(params, eps=1e-8):
        grad = np.zeros(3)
        for i in range(3):
            p1 = params.copy()
            p2 = params.copy()
            p1[i] += eps
            p2[i] -= eps
            f1 = likelihood(p1)
            f2 = likelihood(p2)
            grad[i] = (f1 - f2) / (2 * eps)
        return grad

    # -----------------------------------
    # 初始参数和边界
    initial_params = np.array([0.5*np.mean(data), 0.5*np.var(data), 20])
    bounds = [(0.001, 1.2*np.mean(data)), (0.001,1), (0.001,1000)]

    # -----------------------------------
    # Stage 1 优化（宽松设置）
    x1, f1, info1 = fmin_l_bfgs_b(
        func=likelihood,
        x0=initial_params,
        fprime=likelihood_grad,
        bounds=bounds,
        factr=1e7,        # 宽松收敛
        pgtol=1e-4,
        maxfun=5000,
        disp=1
    )

    proj_g1 = np.linalg.norm(likelihood_grad(x1))
    print(f"[Stage1] f = {f1}, projected_grad = {proj_g1}, iters = {info1['nit']}")

    # -----------------------------------
    # 如果梯度太大，触发 Stage2 严格优化
    if proj_g1 > 1.0:
        print(" -> Triggering stage2 restart (stricter settings)...")
        x2, f2, info2 = fmin_l_bfgs_b(
            func=likelihood,
            x0=x1,
            fprime=likelihood_grad,
            bounds=bounds,
            factr=1e5,        # stricter
            pgtol=1e-8,
            maxfun=10000,
            disp=1
        )
        optimized_params = x2
    else:
        optimized_params = x1

    u1, sigma, lm = optimized_params

    # -----------------------------------
    # 打印参数
    print('u1=' + str(format(u1, '.3f')) + ',' +
          'sigma=' + str(format(sigma, '.3f')) + ',' +
          'lm=' + str(format(lm, '.3f')))

    pams_pending_kc_1[m].extend([u1, sigma, lm])

    # 计算最终的 log-likelihood
    Lh = 0.0
    for x in data:
        pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
        if pdf_val > 0:
            Lh -= math.log(pdf_val)
        else:
            Lh += 1000

    print('The log_likelihood in this loop is:', Lh)
    log_value_pending_kc_1[m] = format(float(Lh), '.3f')
for m in sigma_kc_index:

    log_value_pending_kc_2[m] = float(0)
    pams_pending_kc_2[m] = list()

    data = DB_OB_kc_second[m]

    # -----------------------------------
    # 定义负对数似然函数
    def likelihood(params):
        u1, sigma, lm = params
        Lh = 0.0
        for x in data:
            pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
            if pdf_val > 0:
                Lh -= math.log(pdf_val)
            else:
                Lh += 1000
        return Lh

    # -----------------------------------
    # 数值梯度函数
    def likelihood_grad(params, eps=1e-8):
        grad = np.zeros(3)
        for i in range(3):
            p1 = params.copy()
            p2 = params.copy()
            p1[i] += eps
            p2[i] -= eps
            f1 = likelihood(p1)
            f2 = likelihood(p2)
            grad[i] = (f1 - f2) / (2 * eps)
        return grad

    # -----------------------------------
    # 初始参数和边界
    initial_params = np.array([0.5*np.mean(data), 0.5*np.var(data), 20])
    bounds = [(0.001, 1.2*np.mean(data)), (0.001,1), (0.001,1000)]

    # -----------------------------------
    # Stage 1 优化（宽松设置）
    x1, f1, info1 = fmin_l_bfgs_b(
        func=likelihood,
        x0=initial_params,
        fprime=likelihood_grad,
        bounds=bounds,
        factr=1e7,        # 宽松收敛
        pgtol=1e-4,
        maxfun=5000,
        disp=1
    )

    proj_g1 = np.linalg.norm(likelihood_grad(x1))
    print(f"[Stage1] f = {f1}, projected_grad = {proj_g1}, iters = {info1['nit']}")

    # -----------------------------------
    # 如果梯度太大，触发 Stage2 严格优化
    if proj_g1 > 1.0:
        print(" -> Triggering stage2 restart (stricter settings)...")
        x2, f2, info2 = fmin_l_bfgs_b(
            func=likelihood,
            x0=x1,
            fprime=likelihood_grad,
            bounds=bounds,
            factr=1e5,        # stricter
            pgtol=1e-8,
            maxfun=10000,
            disp=1
        )
        optimized_params = x2
    else:
        optimized_params = x1

    u1, sigma, lm = optimized_params

    # -----------------------------------
    # 打印参数
    print('u1=' + str(format(u1, '.3f')) + ',' +
          'sigma=' + str(format(sigma, '.3f')) + ',' +
          'lm=' + str(format(lm, '.3f')))

    pams_pending_kc_2[m].extend([u1, sigma, lm])

    # 计算最终的 log-likelihood
    Lh = 0.0
    for x in data:
        pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
        if pdf_val > 0:
            Lh -= math.log(pdf_val)
        else:
            Lh += 1000

    print('The log_likelihood in this loop is:', Lh)
    log_value_pending_kc_2[m] = format(float(Lh), '.3f')


# ---------------------------------------------------------------------
# Objective-function evaluation on sigmak (sigma_k_index): LV1 and LV2
#
# For each triple t = {x, y, z} in sigmak, we evaluate whether the relative-distance signals
# derived from quartet trees are better explained by:
#   (i) a unimodal EMG/exGaussian distribution (single mean μ), or
#   (ii) a bimodal / mixture model (two means μ1 and μ2) represented as a two-component EMG mixture.
#
# We compute two sets of objective values (negative log-likelihoods):
#
#   V1 := DXY_O_k[m]   (cherry-related ratio; captures the relative distance between the sister pair x,y)
#       LV1,1(m) = Σ_i  -log f(V1_i | μ, σ, λ)                           [single-component EMG]
#       LV1,2(m) = Σ_i  -log ( w f(V1_i | μ1, σ1, λ) + (1-w) f(V1_i | μ2, σ2, λ) )  [two-component mixture]
#
#   V2 := DXZ_O_k[m]   (pending/attachment-related ratio; captures distances between (x or y) and z)
#       LV2,1(m) = Σ_i  -log f(V2_i | μ, σ, λ)                           [single-component EMG]
#       LV2,2(m) = Σ_i  -log ( w f(V2_i | μ1, σ1, λ) + (1-w) f(V2_i | μ2, σ2, λ) )  [two-component mixture]
#
# Model selection is then based on the gaps:
#   Δ1(m) = |LV1,1(m) - LV1,2(m)|   and   Δ2(m) = |LV2,1(m) - LV2,2(m)|
# using a threshold κ:
#   - If Δ1(m) ≤ κ, V1 is effectively unimodal (supports T1 or N2); otherwise V1 is bimodal (supports N3 or N4).
#   - If Δ2(m) ≤ κ, V2 is effectively unimodal (supports T1 or N3); otherwise V2 is bimodal (supports N2 or N4).
# Combining the two decisions (Δ1 and Δ2) uniquely determines the trinet type for t.
#
# Implementation notes:
#   - The mixture EMG fits for V1 and V2 are performed via constrained L-BFGS-B, with multiple restarts
#     (adjusting initial y and w) to avoid degenerate solutions and poor local optima.
#   - The single-component EMG fits ("*_taxa" blocks) provide the unimodal baselines LV1,1 and LV2,1.
#   - The outputs stored in log_value_* dictionaries are NLL objective values used later for type assignment.
# ---------------------------------------------------------------------


log_value_cherry = dict()
pams_cherry = dict()
for m in sigma_k_index:
    log_value_cherry[m] = float(0)
    pams_cherry[m] = list()

    Sm = np.mean(DXY_O_k[m])
    Sv = np.var(DXY_O_k[m])


    def likelihood(params, Lh=float(0)):
        x, y, w = params

        lm = 1 / x
        y1 = y * (Sv - x ** 2) / (Sv + Sm ** 2 - 2 * Sm * x)
        yw = 1 - y1
        w1 = 1 - w * ((Sm - x) ** 2 / ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1)))
        w2 = 1 - w1
        z1 = w1
        z2 = 1 - z1

        sigma1 = math.sqrt(z1 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / w1)
        sigma2 = math.sqrt(z2 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / (1 - w1))

        if (((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)) < 0:
            u1 = Sm - x
            u2 = Sm - x
        else:
            u1 = Sm - x - w2 * (math.sqrt(
                ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))  ## if y1=0 sqrt(Sv-x^2)
            u2 = Sm - x + w1 * (math.sqrt(((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))

        for j in range(len(DXY_O_k[m])):
            # print('test')
            # print(Lh)
            # print(x1, y1, z1, w1, x2, y2, z2)
            # print(exponnorm.pdf(DXY_O_k[j], 1 / math.sqrt(z2/(1-z2)), loc=(x2-math.sqrt((1-z2)*y2)),scale=math.sqrt(y2*z2)))
            # print((1-w1)* exponnorm.pdf(DXY_O_k[j], 1 / math.sqrt(z2/(1-z2)), loc=(x2-math.sqrt((1-z2)*y2)),scale=math.sqrt(y2*z2)))

            term1 = w1 * exponnorm.pdf(DXY_O_k[m][j], 1 / (sigma1 * lm), loc=u1, scale=sigma1)
            term2 = (1 - w1) * exponnorm.pdf(DXY_O_k[m][j], 1 / (sigma2 * lm), loc=u2, scale=sigma2)
            if term1 + term2 > 0:
                Lh -= (math.log(term1 + term2))
            else:
                Lh += 1000
        print(Lh)
        print(u1, u2, sigma1, sigma2, lm)

        return Lh


    def bfgs_with_bounds():
        initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.5, 0.5])
        print(f"\n[Stage1] Initial params = {initial_params}")
        bounds = [
            (0.01 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.999 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv))),
            (0.00001, 0.999), (0.00001, 0.999)]
        result1 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds, options={
            'disp': True,
            'iprint': 1,  # 显示迭代细节
            'ftol': 1e-10,
            'gtol': 1e-8,
            'maxiter': 3000,
            'maxfun': 5000
        })

        optimized = result1
        grad_norm = np.linalg.norm(result1.jac)
        print(f"[Stage0] f = {result1.fun:.6f}, |proj_g| = {grad_norm:.6f}")
        if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
            print(" we need a lower initial value of y and w in the 1 loop.")
            initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.3, 0.3])
            result2 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds, options={
                'disp': True,
                'iprint': 1,  # 显示迭代细节
                'ftol': 1e-10,
                'gtol': 1e-8,
                'maxiter': 3000,
                'maxfun': 5000
            })
            optimized = result2
            grad_norm = np.linalg.norm(result2.jac)
            print(f"[Stage1] f = {result2.fun:.6f}, |proj_g| = {grad_norm:.6f}")
            print(grad_norm)
            if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
                print(" we need a lower initial value of y and w in the 2 loop.")
                print("The optimized params2 is:", optimized.x)
                initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.1, 0.1])
                result3 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds,
                                   options={'disp': True, 'iprint': 1, 'ftol': 1e-10, 'gtol': 1e-8, 'maxiter': 3000,
                                            'maxfun': 5000})
                optimized = result3
                grad_norm = np.linalg.norm(result3.jac)
                print(f"[Stage2] f = {result3.fun:.6f}, |proj_g| = {grad_norm:.6f}")
                print(grad_norm)
                if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
                    print(" we need a lower initial value of y and w in the 3 loop.")
                    initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.05, 0.05])
                    result4 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds,
                                       options={'disp': True, 'iprint': 1, 'ftol': 1e-10, 'gtol': 1e-8, 'maxiter': 3000,
                                                'maxfun': 5000})
                    optimized = result4
                    grad_norm = np.linalg.norm(result4.jac)
                    print(f"[Stage3] f = {result4.fun:.6f}, |proj_g| = {grad_norm:.6f}")
                    print(grad_norm)
                    if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
                        print(" we need a lower initial value of y and w in the 4 loop.")
                        initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.75, 0.75])
                        result5 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds,
                                           options={'disp': True, 'iprint': 1, 'ftol': 1e-10, 'gtol': 1e-8,
                                                    'maxiter': 3000, 'maxfun': 5000})
                        optimized = result5
                        grad_norm = np.linalg.norm(result5.jac)
                        print(f"[Stage4] f = {result5.fun:.6f}, |proj_g| = {grad_norm:.6f}")
                        print(grad_norm)
                        print("The optimized params5* is:", optimized.x)
                    else:
                        optimized = result4
                        print("The optimized params4* is:", optimized.x)
                else:
                    optimized = result3
                    print("The optimized params3* is:", optimized.x)
            else:
                optimized = result2
                print("The optimized params2* is:", optimized.x)
        else:
            optimized = result1
            print("The optimized params1* is:", optimized.x)

        return optimized.x


    optimized_params = []
    optimized_params = bfgs_with_bounds()

    x = optimized_params[0]
    y = optimized_params[1]
    w = optimized_params[2]

    y1 = y * (Sv - x ** 2) / (Sv + Sm ** 2 - 2 * Sm * x)
    y2 = 1 - y1

    lm = 1 / x

    w1 = 1 - w * ((Sm - x) ** 2 / ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1)))
    w2 = 1 - w1
    z1 = w1
    z2 = 1 - z1
    print(x, y, z1, w)

    sigma1 = math.sqrt(z1 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / w1)
    sigma2 = math.sqrt(z2 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / (1 - w1))

    u1 = Sm - x - w2 * (math.sqrt(((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))

    u2 = Sm - x + w1 * (math.sqrt(((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))
    Lh = float(0)

    print('x=' + str(format(x, '.3f')) + 'y=' + str(format(y, '.3f')) + 'z1=' + str(format(z1, '.3f')) + 'w=' + str(
        format(w, '.3f')) + 'w1=' + str(format(w1, '.3f')) + ','
          + 'u1=' + str(format(u1, '.3f')) + ',' + 'u2=' + str(format(u2, '.3f')) + ',' + 'sigma1=' + str(
        format(sigma1, '.3f')) + ',' + 'sigma2=' + str(format(sigma2, '.3f')) + ',''lm=' + str(format(lm, '.3f')))

    pams_cherry[m].append(x)
    pams_cherry[m].append(y)
    pams_cherry[m].append(z1)
    pams_cherry[m].append(w)
    pams_cherry[m].append(str("||||||"))
    pams_cherry[m].append(w1)
    pams_cherry[m].append(u1)
    pams_cherry[m].append(u2)
    pams_cherry[m].append(sigma1)
    pams_cherry[m].append(sigma2)
    pams_cherry[m].append(lm)

    x_new = 1 / lm
    y1_new = (w1 * sigma1 ** 2 + w2 * sigma2 ** 2) / (Sv + Sm ** 2 - 2 * Sm * x)
    y_new = (w1 * sigma1 ** 2 + w2 * sigma2 ** 2) / (Sv - x ** 2)
    z1_new = w1 * sigma1 ** 2 / (w1 * sigma1 ** 2 + w2 * sigma2 ** 2)
    w_new = (1 - w1) / ((Sm - x) ** 2 / ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1_new)))

    a = 0
    b = 0

    for j in range(len(DXY_O_k[m])):

        a = w1 * exponnorm.pdf(DXY_O_k[m][j], 1 / (sigma1 * lm), loc=u1, scale=sigma1)
        b = (1 - w1) * exponnorm.pdf(DXY_O_k[m][j], 1 / (sigma2 * lm), loc=u2, scale=sigma2)

        # print(a,b)

        if a + b > 0:
            Lh -= (math.log(a + b))
        else:
            print('no')
            Lh += 1000

    print('The log_likelihood in this loop is:', Lh)
    log_value_cherry[m] = format(float(Lh), '.3f')
mp.dps = 500  # 设置精度为 50 位
log_value_cherry_taxa=dict()
pams_cherry_taxa=dict()
for m in sigma_k_index:

    log_value_cherry_taxa[m] = float(0)
    pams_cherry_taxa[m] = list()

    data = DXY_O_k[m]
    data = np.array(data)


    # -----------------------------
    # 定义负对数似然函数
    def likelihood(params):
        u1, sigma, lm = params
        Lh = 0.0
        for x in data:
            pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
            if pdf_val > 0:
                Lh -= math.log(pdf_val)
            else:
                Lh += 1000  # 罚项，避免溢出或负概率
        return Lh


    # -----------------------------
    # 数值梯度函数（中央差分）
    def likelihood_grad(params, eps=1e-6):
        grad = np.zeros(3)
        for i in range(3):
            p1, p2 = params.copy(), params.copy()
            p1[i] += eps
            p2[i] -= eps
            grad[i] = (likelihood(p1) - likelihood(p2)) / (2 * eps)
        return grad


    # -----------------------------
    # 初始参数和边界
    initial_params = np.array([
        0.4 * np.mean(data),
        0.5 * np.var(data),
        20
    ])
    bounds = [
        (0.02 * np.mean(data), 1.2 * np.mean(data)),
        (0.001, 0.999),
        (0.002, 300)
    ]

    print(f"\n===== Start optimizing dataset {m} =====")

    # -----------------------------
    # Stage 1 优化（宽松）
    x1, f1, info1 = fmin_l_bfgs_b(
        func=likelihood,
        x0=initial_params,
        fprime=likelihood_grad,
        bounds=bounds,
        factr=1e7,  # 宽松收敛
        pgtol=1e-4,
        maxfun=2000,
        disp=True
    )

    proj_g1 = np.linalg.norm(likelihood_grad(x1))
    print(f"[Stage 1] f = {f1:.3f}, grad_norm = {proj_g1:.3e}, iters = {info1['nit']}")

    # -----------------------------
    # 若梯度仍大，进入 Stage 2 严格优化
    if proj_g1 > 1.0:
        print(" → Stage 1 未充分收敛，进入 Stage 2 严格优化 ...")
        x2, f2, info2 = fmin_l_bfgs_b(
            func=likelihood,
            x0=x1,
            fprime=likelihood_grad,
            bounds=bounds,
            factr=1e5,  # 严格收敛
            pgtol=1e-8,
            maxfun=5000,
            disp=True
        )
        optimized_params = x2
        final_Lh = f2
    else:
        optimized_params = x1
        final_Lh = f1

    # -----------------------------
    # 输出最终参数
    u1, sigma, lm = optimized_params
    print(f"[Final Params] u1 = {u1:.5f}, sigma = {sigma:.5f}, lm = {lm:.5f}")
    print(f"[Final LogLikelihood] {final_Lh:.3f}")

    # 存储结果
    pams_cherry_taxa[m].extend([u1, sigma, lm])
    log_value_cherry_taxa[m] = format(float(final_Lh), '.3f')
print("\n===== Optimization finished for all datasets =====")
log_value_pending=dict()
pams_pending=dict()
for m in sigma_k_index:
    log_value_pending[m]=float(0)
    pams_pending[m]=list()

    Sm = np.mean(DXZ_O_k[m])
    Sv = np.var(DXZ_O_k[m])

    def likelihood(params, Lh=float(0)):
        x, y, w = params

        lm = 1 / x
        y1 = y * (Sv - x ** 2) / (Sv + Sm ** 2 - 2 * Sm * x)
        yw = 1 - y1
        w1 = 1 - w * ((Sm - x) ** 2 / ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1)))
        w2 = 1 - w1
        z1 = w1
        z2 = 1 - z1

        sigma1 = math.sqrt(z1 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / w1)
        sigma2 = math.sqrt(z2 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / (1 - w1))

        if (((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)) < 0:
            u1 = Sm - x
            u2 = Sm - x
        else:
            u1 = Sm - x - w2 * (math.sqrt(
                ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))  ## if y1=0 sqrt(Sv-x^2)
            u2 = Sm - x + w1 * (math.sqrt(((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))

        for j in range(len(DXZ_O_k[m])):
            # print('test')
            # print(Lh)
            # print(x1, y1, z1, w1, x2, y2, z2)
            # print(exponnorm.pdf(DXZ_O_k[j], 1 / math.sqrt(z2/(1-z2)), loc=(x2-math.sqrt((1-z2)*y2)),scale=math.sqrt(y2*z2)))
            # print((1-w1)* exponnorm.pdf(DXZ_O_k[j], 1 / math.sqrt(z2/(1-z2)), loc=(x2-math.sqrt((1-z2)*y2)),scale=math.sqrt(y2*z2)))

            term1 = w1 * exponnorm.pdf(DXZ_O_k[m][j], 1 / (sigma1 * lm), loc=u1, scale=sigma1)
            term2 = (1 - w1) * exponnorm.pdf(DXZ_O_k[m][j], 1 / (sigma2 * lm), loc=u2, scale=sigma2)
            if term1 + term2 > 0:
                Lh -= (math.log(term1 + term2))
            else:
                Lh += 1000
        print(Lh)
        print(u1, u2, sigma1, sigma2, lm)

        return Lh


    def bfgs_with_bounds():
        initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.5, 0.5])
        print(f"\n[Stage1] Initial params = {initial_params}")
        bounds = [
            (0.01 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.999 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv))),
            (0.00001, 0.999), (0.00001, 0.999)]
        result1 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds, options={
            'disp': True,
            'iprint': 1,  # 显示迭代细节
            'ftol': 1e-10,
            'gtol': 1e-8,
            'maxiter': 3000,
            'maxfun': 5000
        })

        optimized = result1
        grad_norm = np.linalg.norm(result1.jac)
        print(f"[Stage0] f = {result1.fun:.6f}, |proj_g| = {grad_norm:.6f}")
        if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
            print(" we need a lower initial value of y and w in the 1 loop.")
            initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.3, 0.3])
            result2 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds, options={
                'disp': True,
                'iprint': 1,  # 显示迭代细节
                'ftol': 1e-10,
                'gtol': 1e-8,
                'maxiter': 3000,
                'maxfun': 5000
            })
            optimized = result2
            grad_norm = np.linalg.norm(result2.jac)
            print(f"[Stage1] f = {result2.fun:.6f}, |proj_g| = {grad_norm:.6f}")
            print(grad_norm)
            if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
                print(" we need a lower initial value of y and w in the 2 loop.")
                print("The optimized params2 is:", optimized.x)
                initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.1, 0.1])
                result3 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds,
                                   options={'disp': True, 'iprint': 1, 'ftol': 1e-10, 'gtol': 1e-8, 'maxiter': 3000,
                                            'maxfun': 5000})
                optimized = result3
                grad_norm = np.linalg.norm(result3.jac)
                print(f"[Stage2] f = {result3.fun:.6f}, |proj_g| = {grad_norm:.6f}")
                print(grad_norm)
                if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
                    print(" we need a lower initial value of y and w in the 3 loop.")
                    initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.05, 0.05])
                    result4 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds,
                                       options={'disp': True, 'iprint': 1, 'ftol': 1e-10, 'gtol': 1e-8, 'maxiter': 3000,
                                                'maxfun': 5000})
                    optimized = result4
                    grad_norm = np.linalg.norm(result4.jac)
                    print(f"[Stage3] f = {result4.fun:.6f}, |proj_g| = {grad_norm:.6f}")
                    print(grad_norm)
                    if optimized.x[2] < 0.002 or optimized.x[2] > 0.998 or grad_norm > 1:
                        print(" we need a lower initial value of y and w in the 4 loop.")
                        initial_params = np.array([0.3 * min(Sv / (2 * Sm) + Sm / 2, math.sqrt(Sv)), 0.75, 0.75])
                        result5 = minimize(likelihood, initial_params, method='L-BFGS-B', bounds=bounds,
                                           options={'disp': True, 'iprint': 1, 'ftol': 1e-10, 'gtol': 1e-8,
                                                    'maxiter': 3000, 'maxfun': 5000})
                        optimized = result5
                        grad_norm = np.linalg.norm(result5.jac)
                        print(f"[Stage4] f = {result5.fun:.6f}, |proj_g| = {grad_norm:.6f}")
                        print(grad_norm)
                        print("The optimized params5* is:", optimized.x)
                    else:
                        optimized = result4
                        print("The optimized params4* is:", optimized.x)
                else:
                    optimized = result3
                    print("The optimized params3* is:", optimized.x)
            else:
                optimized = result2
                print("The optimized params2* is:", optimized.x)
        else:
            optimized = result1
            print("The optimized params1* is:", optimized.x)

        return optimized.x


    optimized_params = []
    optimized_params = bfgs_with_bounds()

    x = optimized_params[0]
    y = optimized_params[1]
    w = optimized_params[2]

    y1 = y * (Sv - x ** 2) / (Sv + Sm ** 2 - 2 * Sm * x)
    y2 = 1 - y1

    lm = 1 / x

    w1 = 1 - w * ((Sm - x) ** 2 / ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1)))
    w2 = 1 - w1
    z1 = w1
    z2 = 1 - z1
    print(x, y, z1, w)

    sigma1 = math.sqrt(z1 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / w1)
    sigma2 = math.sqrt(z2 * y1 * (Sv + Sm ** 2 - 2 * Sm * x) / (1 - w1))

    u1 = Sm - x - w2 * (math.sqrt(((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))

    u2 = Sm - x + w1 * (math.sqrt(((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1) - (Sm - x) ** 2) / (w1 * w2)))
    Lh = float(0)

    print('x=' + str(format(x, '.3f')) + 'y=' + str(format(y, '.3f')) + 'z1=' + str(format(z1, '.3f')) + 'w=' + str(format(w, '.3f')) + 'w1=' + str(format(w1, '.3f')) + ','
          + 'u1=' + str(format(u1, '.3f')) + ',' + 'u2=' + str(format(u2, '.3f')) + ',' + 'sigma1=' + str(format(sigma1, '.3f')) + ',' + 'sigma2=' + str(format(sigma2, '.3f')) + ',''lm=' + str(format(lm, '.3f')))

    pams_pending[m].append(x)
    pams_pending[m].append(y)
    pams_pending[m].append(z1)
    pams_pending[m].append(w)
    pams_pending[m].append(str("||||||"))
    pams_pending[m].append(w1)
    pams_pending[m].append(u1)
    pams_pending[m].append(u2)
    pams_pending[m].append(sigma1)
    pams_pending[m].append(sigma2)
    pams_pending[m].append(lm)

    x_new = 1 / lm
    y1_new = (w1 * sigma1 ** 2 + w2 * sigma2 ** 2) / (Sv + Sm ** 2 - 2 * Sm * x)
    y_new = (w1 * sigma1 ** 2 + w2 * sigma2 ** 2) / (Sv - x ** 2)
    z1_new = w1 * sigma1 ** 2 / (w1 * sigma1 ** 2 + w2 * sigma2 ** 2)
    w_new = (1 - w1) / ((Sm - x) ** 2 / ((Sv + Sm ** 2 - 2 * Sm * x) * (1 - y1_new)))

    a = 0
    b = 0

    for j in range(len(DXZ_O_k[m])):

        a = w1 * exponnorm.pdf(DXZ_O_k[m][j], 1 / (sigma1 * lm), loc=u1, scale=sigma1)
        b = (1 - w1) * exponnorm.pdf(DXZ_O_k[m][j], 1 / (sigma2 * lm), loc=u2, scale=sigma2)

        # print(a,b)

        if a + b > 0:
            Lh -= (math.log(a + b))
        else:
            print('no')
            Lh += 1000

    print('The log_likelihood in this loop is:', Lh)
    log_value_pending[m] = format(float(Lh), '.3f')
log_value_pending_taxa=dict()
pams_pending_taxa=dict()
for m in sigma_k_index:


    log_value_pending_taxa[m] = float(0)
    pams_pending_taxa[m] = list()

    data = DXZ_O_k[m]
    data = np.array(data)

    # -----------------------------
    # 定义负对数似然函数
    def likelihood(params):
        u1, sigma, lm = params
        Lh = 0.0
        for x in data:
            pdf_val = exponnorm.pdf(x, 1 / (sigma * lm), loc=u1, scale=sigma)
            if pdf_val > 0:
                Lh -= math.log(pdf_val)
            else:
                Lh += 1000  # 罚项，避免溢出或负概率
        return Lh

    # -----------------------------
    # 数值梯度函数（中央差分）
    def likelihood_grad(params, eps=1e-6):
        grad = np.zeros(3)
        for i in range(3):
            p1, p2 = params.copy(), params.copy()
            p1[i] += eps
            p2[i] -= eps
            grad[i] = (likelihood(p1) - likelihood(p2)) / (2 * eps)
        return grad

    # -----------------------------
    # 初始参数和边界
    initial_params = np.array([
        0.4 * np.mean(data),
        0.5 * np.var(data),
        20
    ])
    bounds = [
        (0.02 * np.mean(data), 1.2 * np.mean(data)),
        (0.001, 0.999),
        (0.002, 1000)
    ]

    print(f"\n===== Start optimizing dataset {m} =====")

    # -----------------------------
    # Stage 1 优化（宽松）
    x1, f1, info1 = fmin_l_bfgs_b(
        func=likelihood,
        x0=initial_params,
        fprime=likelihood_grad,
        bounds=bounds,
        factr=1e7,        # 宽松收敛
        pgtol=1e-4,
        maxfun=2000,
        disp=True
    )

    proj_g1 = np.linalg.norm(likelihood_grad(x1))
    print(f"[Stage 1] f = {f1:.3f}, grad_norm = {proj_g1:.3e}, iters = {info1['nit']}")

    # -----------------------------
    # 若梯度仍大，进入 Stage 2 严格优化
    if proj_g1 > 1.0:
        print(" → Stage 1 未充分收敛，进入 Stage 2 严格优化 ...")
        x2, f2, info2 = fmin_l_bfgs_b(
            func=likelihood,
            x0=x1,
            fprime=likelihood_grad,
            bounds=bounds,
            factr=1e5,        # 严格收敛
            pgtol=1e-8,
            maxfun=5000,
            disp=True
        )
        optimized_params = x2
        final_Lh = f2
    else:
        optimized_params = x1
        final_Lh = f1

    # -----------------------------
    # 输出最终参数
    u1, sigma, lm = optimized_params
    print(f"[Final Params] u1 = {u1:.5f}, sigma = {sigma:.5f}, lm = {lm:.5f}")
    print(f"[Final LogLikelihood] {final_Lh:.3f}")

    # 存储结果
    pams_pending_taxa[m].extend([u1, sigma, lm])
    log_value_pending_taxa[m] = format(float(final_Lh), '.3f')
print("\n===== Optimization finished for all datasets =====")


# ---------------------------------------------------------------------
# Final step: write the inferred trinets to a .tnets/.tnet file.
# For each triple index, determine its trinet type (T1/N2/N3/N4 for Σ_k, and S1/S2 for Σ_k^C)
# based on the objective-function gaps, then output the trinet label together with diagnostics
# (gap values, fitted parameters, and binomial-test p-values).
# The last line reports the average number of loci per trinet.
# ---------------------------------------------------------------------


f1=open(path1+'gamma_test_6.704_0.149_M.tnets','w+')

gap_value_pending_taxa=dict()
gap_value_cherry_taxa=dict()
gap_value_kc_taxa=dict()


for i in w_count:
    w_index[i]=list()
    for j in w_count[i]:
        w_index[i].append(j)

for i in sigma_k_index:
    gap_value_pending_taxa[i]=0
    gap_value_pending_taxa[i]=float(log_value_pending_taxa[i])-float(log_value_pending[i])

    gap_value_cherry_taxa[i]=0
    gap_value_cherry_taxa[i]=float(log_value_cherry_taxa[i])-float(log_value_cherry[i])

for i in log_value_pending_kc_m:
    gap_value_kc_taxa[i]=0
    gap_value_kc_taxa[i]=(float(log_value_pending_kc_m[i])-(float(log_value_pending_kc_1[i])+float(log_value_pending_kc_2[i])))




for i in w_index:
    if i in sigma_k_index:
        if  gap_value_pending_taxa[i]<=32 and gap_value_cherry_taxa[i]<=32:
            f1.write(taxa[str(w_index[i][0].split('_')[0])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[1])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[2].split('|')[1])]+' ')
            f1.write('T1 ')
            f1.write(str(gap_value_pending_taxa[i]))
            f1.write(' '+str(gap_value_cherry_taxa[i]))
            f1.write(' ' + str(pams_cherry[i][5]))
            f1.write(' ' + str(pams_cherry[i][6]))
            f1.write(' ' + str(pams_cherry[i][7]))
            f1.write(' ' + str(pams_pending[i][5]))
            f1.write(' ' + str(pams_pending[i][6]))
            f1.write(' ' + str(pams_pending[i][7]))
            f1.write(' ' + str(bio_p_second_value[i]))
            f1.write(' ' + str(bio_p_first_value[i]))
            f1.write('\n')
        elif gap_value_pending_taxa[i]>32 and gap_value_cherry_taxa[i]<=32:
            f1.write(taxa[str(w_index[i][0].split('_')[0])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[1])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[2].split('|')[1])]+' ')
            f1.write('N2 ')
            f1.write(str(gap_value_pending_taxa[i]))
            f1.write(' '+str(gap_value_cherry_taxa[i]))
            f1.write(' ' + str(pams_cherry[i][5]))
            f1.write(' ' + str(pams_cherry[i][6]))
            f1.write(' ' + str(pams_cherry[i][7]))
            f1.write(' ' + str(pams_pending[i][5]))
            f1.write(' ' + str(pams_pending[i][6]))
            f1.write(' ' + str(pams_pending[i][7]))
            f1.write(' ' + str(bio_p_second_value[i]))
            f1.write(' ' + str(bio_p_first_value[i]))
            f1.write('\n')
        elif gap_value_pending_taxa[i]<=32 and gap_value_cherry_taxa[i]>32:
            f1.write(taxa[str(w_index[i][0].split('_')[0])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[1])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[2].split('|')[1])]+' ')
            f1.write('N3 ')
            f1.write(str(gap_value_pending_taxa[i]))
            f1.write(' '+str(gap_value_cherry_taxa[i]))
            f1.write(' ' + str(pams_cherry[i][5]))
            f1.write(' ' + str(pams_cherry[i][6]))
            f1.write(' ' + str(pams_cherry[i][7]))
            f1.write(' ' + str(pams_pending[i][5]))
            f1.write(' ' + str(pams_pending[i][6]))
            f1.write(' ' + str(pams_pending[i][7]))
            f1.write(' ' + str(bio_p_second_value[i]))
            f1.write(' ' + str(bio_p_first_value[i]))
            f1.write('\n')
        elif gap_value_pending_taxa[i]>32 and gap_value_cherry_taxa[i]>32:
            f1.write(taxa[str(w_index[i][0].split('_')[0])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[1])]+' ')
            f1.write(taxa[str(w_index[i][0].split('_')[2].split('|')[1])]+' ')
            f1.write('N4 ')
            f1.write(str(gap_value_pending_taxa[i]))
            f1.write(' '+str(gap_value_cherry_taxa[i]))
            f1.write(' ' + str(pams_cherry[i][5]))
            f1.write(' ' + str(pams_cherry[i][6]))
            f1.write(' ' + str(pams_cherry[i][7]))
            f1.write(' ' + str(pams_pending[i][5]))
            f1.write(' ' + str(pams_pending[i][6]))
            f1.write(' ' + str(pams_pending[i][7]))
            f1.write(' ' + str(bio_p_second_value[i]))
            f1.write(' ' + str(bio_p_first_value[i]))
            f1.write('\n')
    elif i in sigma_kc_index:
        if  gap_value_kc_taxa[i]<32:
            f1.write(taxa[str(w_index[i][2].split('_')[0])]+' ')
            f1.write(taxa[str(w_index[i][2].split('_')[1])]+' ')
            f1.write(taxa[str(w_index[i][2].split('_')[2].split('|')[1])]+' ')
            f1.write('S1 ')
            f1.write(str(gap_value_kc_taxa[i]))
            f1.write(' ' + str(bio_p_second_value[i]))
            f1.write(' ' + str(bio_p_first_value[i]))
            f1.write(' ' + str((list(w_count[i].values())[0] - list(w_count[i].values())[2]) / (
                        (list(w_count[i].values())[0] - list(w_count[i].values())[2]) + (
                            list(w_count[i].values())[1] - list(w_count[i].values())[2]))))
            f1.write(' ' + str((list(w_count[i].values())[1] - list(w_count[i].values())[2]) / (
                        (list(w_count[i].values())[0] - list(w_count[i].values())[2]) + (
                            list(w_count[i].values())[1] - list(w_count[i].values())[2]))))

            f1.write('\n')
        elif gap_value_kc_taxa[i]>=32:
            if DB_OB_kc_avg_first[i]>=DB_OB_kc_avg_second[i]:
                x = w_index[i][0].split('_')
                z=w_index[i][2].split('_')[2].split('|')[1]
                x.remove(str(z))
                f1.write(taxa[x[0]]+' ')
                f1.write(taxa[str(w_index[i][2].split('_')[2].split('|')[1])]+' ')
                f1.write(taxa[x[1].split('|')[1]]+' ')
                f1.write('S2 ')
                f1.write(str(gap_value_kc_taxa[i]))
                f1.write(' ' + str(bio_p_second_value[i]))
                f1.write(' ' + str(bio_p_first_value[i]))
                f1.write(' ' + str((list(w_count[i].values())[0] - list(w_count[i].values())[2]) / (
                            (list(w_count[i].values())[0] - list(w_count[i].values())[2]) + (
                                list(w_count[i].values())[1] - list(w_count[i].values())[2]))))
                f1.write(' ' + str((list(w_count[i].values())[1] - list(w_count[i].values())[2]) / (
                            (list(w_count[i].values())[0] - list(w_count[i].values())[2]) + (
                                list(w_count[i].values())[1] - list(w_count[i].values())[2]))))

                f1.write('\n')
            else:
                y = w_index[i][1].split('_')
                z = w_index[i][2].split('_')[2].split('|')[1]
                y.remove(str(z))
                f1.write(taxa[y[0]]+' ')
                f1.write(taxa[str(w_index[i][2].split('_')[2].split('|')[1])]+' ')
                f1.write(taxa[y[1].split('|')[1]]+' ')
                f1.write('S2 ')
                f1.write(str(gap_value_kc_taxa[i]))
                f1.write(' ' + str(bio_p_second_value[i]))
                f1.write(' ' + str(bio_p_first_value[i]))
                f1.write(' ' + str((list(w_count[i].values())[0] - list(w_count[i].values())[2]) / (
                            (list(w_count[i].values())[0] - list(w_count[i].values())[2]) + (
                                list(w_count[i].values())[1] - list(w_count[i].values())[2]))))
                f1.write(' ' + str((list(w_count[i].values())[1] - list(w_count[i].values())[2]) / (
                            (list(w_count[i].values())[0] - list(w_count[i].values())[2]) + (
                                list(w_count[i].values())[1] - list(w_count[i].values())[2]))))

                f1.write('\n')
    else:
        continue

# Append summary information to the output file (includes total number of quartet trees).

f2=open(path1+'gamma_test_6.704_0.14916467780429596_M_quartet_adjust_deleteoutlier_modify.tre','r+',encoding='utf-8')
lines=f2.readlines()
loci_num=len(lines)
print('The number of loci is',loci_num/len(w_count))
f1.write("N U M Num "+str(loci_num))

f1.close()

