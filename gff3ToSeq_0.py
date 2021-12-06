import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from datetime import datetime as dt
from datetime import timedelta
from scipy.signal import lfilter
import pickle as pickle
import sys
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

# get reference protein sequences for SARS CoV2 (Wuhan-1)
ref_dict = {rec.id : rec.seq for rec in SeqIO.parse("SARS_CoV2_ref_protSeq_NSP.fasta", "fasta")}

# dictionary between protein names and Uniprot
uni2names={'YP_009725318.1':'ORF7b', 'YP_009725295.1':'ORF1a', 'QHI42199.1':'ORF10', 'QHD43423.2':'N', 'QHD43423.1':'N', 'QHD43422.1':'ORF8', 'QHD43421.1':'ORF7a', 'QHD43420.1':'ORF6', 'QHD43419.1':'M', 'QHD43418.1':'E', 'QHD43417.1':'ORF3a', 'QHD43416.1':'S', 'QHD43415.1':'ORF1ab'}

# nsp proteins coordinates within ORF1ab
ORF1ab_ref={'NSP1':[1,180],
            'NSP2':[181,818],
            'NSP3':[819,2763],
            'NSP4':[2764,3263],
            'NSP5':[3264,3569],
            'NSP6':[3570,3859],
            'NSP7':[3860,3942],
            'NSP8':[3943,4140],
            'NSP9':[4141,4253],
            'NSP10':[4254,4392],
            'NSP12':[4393,5324],
            'NSP13':[5325,5925],
            'NSP14':[5926,6452],
            'NSP15':[6453,6798],
            'NSP16':[6799,7096]}

# extract mutations from the annotation field in gff3 - reindex to individual nsp proteins in case of ORF1ab. Return DF.
def MutExtract(block,ORF1ab_ref):
    
    start_l=[]
    end_l=[]
    target_l=[]
    prot_l=[]
    for i in range(len(block)) :
        
        mut=block.iloc[i]['mut']
        if block.iloc[i]['VEP']=='missense_variant':
            target=mut[-1]
            start=end=int(re.findall(r'[A-Za-z]+|\d+', mut.split('>')[0])[0])
        elif block.loc[i,'VEP']=='inframe_deletion':
            target='' if mut.split('>')[1]=='-' else mut.split('>')[1]
            if '-' in mut.split('>')[0]:
                start=int(mut.split('>')[0].split('-')[0])
                end=int(re.findall(r'[A-Za-z]+|\d+', mut.split('>')[0].split('-')[1])[0])
            else:
                start=end=int(re.findall(r'[A-Za-z]+|\d+', mut.split('>')[0])[0])
        elif block.loc[i,'VEP']=='inframe_insertion':
            if '->' in mut:
                target=mut.split('->')[1]
                start=int(mut.split('->')[0].split('-')[0])
                end=int(mut.split('->')[0].split('-')[1])
            else:
                target=mut.split('>')[1]
                start=end=int(re.findall(r'[A-Za-z]+|\d+', mut.split('>')[0])[0])

             
        if block.iloc[i]['prot0']=='ORF1ab':
            for var in ORF1ab_ref.keys():
                if ORF1ab_ref[var][0] <= start <= (ORF1ab_ref[var][1]):
                    prot=var
                    start=start-ORF1ab_ref[var][0]+1
                    end=end-ORF1ab_ref[var][0]+1
                    
        else:
            prot=block.iloc[i]['prot0']
            
        start_l.append(start)
        end_l.append(end)
        target_l.append(target)
        prot_l.append(prot)
        
    return pd.DataFrame(list(zip(prot_l,start_l,end_l,target_l)),
              columns=['prot','start_p','end_p','target'])
            
        
#Modify ref sequence by adding aa-changing mutations in gff3 file 
def SeqMod(block,ref_seq):
    
    block=block.sort_values('start',ascending=False)
    
    mutated_seq=MutableSeq(str(ref_seq))
    for i in range(len(block)) :
        if block.iloc[i]['VEP']=='missense_variant':
            mutated_seq[int(block.iloc[i]['start_p'])-1]=block.iloc[i]['target']
        elif block.iloc[i]['VEP']=='inframe_deletion':
            mutated_seq=mutated_seq[0:block.iloc[i]['start_p']-1] + block.iloc[i]['target'] + mutated_seq[block.iloc[i]['end_p']:] 
        elif block.iloc[i]['VEP']=='inframe_insertion':
            if '->' in block.iloc[i]['mut']:
                mutated_seq=mutated_seq[0:block.iloc[i]['start_p']-1] + block.iloc[i]['target'] + mutated_seq[block.iloc[i]['start_p']-1:]
            else:
                mutated_seq=mutated_seq[0:block.iloc[i]['start_p']-1] + block.iloc[i]['target'] + mutated_seq[block.iloc[i]['start_p']:]
    return mutated_seq       



with open(f"../Covid19/Input_files/jh_daily_ir_08_22_21.pkl",'rb') as file: #Contains only country name for each date
    jh_daily_ir = pickle.load(file)
    
with open(f"../Covid19/Input_files/meta_8_22_21.pkl",'rb') as file:
    meta = pickle.load(file)

date_range = jh_daily_ir.columns[0:100]

all_names = os.listdir('../gff3_cncb_08_23_21')
column_names = ['variant type','start','end', 'info']
missing = []

for date in date_range:
    meta_daily=meta[(meta['Sample Collection Date']==date)]
    
    for identifier in tqdm(meta_daily.index):
    #Searching for correct identifier
    #--------------------------
        #No alternate name is ' '
        file_name = ''
        #Check if accession id in file names, if not check related ids
        if '2019-nCoV_'+identifier+'_variants.gff3' in all_names:
            file_name = '2019-nCoV_'+identifier+'_variants.gff3'
        # checking alternate names
        elif meta_daily.loc[identifier,'Related ID'] != ' ':
            for alt_identifier in meta_daily.loc[identifier,'Related ID'].replace(' ','').split(','):
                if '2019-nCoV_'+alt_identifier+'_variants.gff3' in all_names:
                    file_name = '2019-nCoV_'+alt_identifier+'_variants.gff3'
                    break
            #Added in case alternate names are also not found in gffs
            if file_name == '':
                missing.append(identifier)
                continue
        # If file name has not been updated, then there is no matching identifier, move to next index
        elif file_name == '':
            missing.append(identifier)
            continue
        #--------------------------

        #Filtering files with no variants
        #--------------------------
        with open(f'../gff3_cncb_08_23_21/{file_name}') as text_file:
            lines = text_file.readlines()
            counter = 0
            for l in lines:
                if '#' in l:
                    counter += 1
        #Number of info lines should be less than total, if not then there are no mutations
        #--------------------------

        #List to keep track of which files used already, save and import this in future to avoid redundant search
        #processed_identifiers.append(identifier)

        if counter<len(lines): #and meta.loc[identifier,'Country'] in set(jh_data['Country_Region']): #country needs to be found in jh data for inf/fata rates

            gff = pd.read_csv(f'../gff3_cncb_08_23_21/{file_name}',sep='\t',skiprows=counter,usecols=[1,3,4,8],names=column_names)
            print(file_name)

            info_df = pd.DataFrame(gff['info'].str.split(';').values.tolist(),columns=[0,1,'Ref','Alt','Description']).drop([0,1],axis=1)
            gff = gff.drop(['info'],axis=1)
            gff['Country'] = [meta.loc[identifier,'Country']]*gff.shape[0]
            temp_df = pd.concat([gff,info_df],axis=1)

            # limit to VEP that have AA changes; missense and indels.
            temp_df=temp_df.loc[temp_df['Description'].str.contains('missense_variant|inframe_deletion|inframe_insertion'),].reset_index(drop=True)
            if temp_df.empty:
                continue

            # remove ambiguous cases with multiple values for VEP
            if any(temp_df.Description.str.count(',')!=2):
                temp_df=temp_df.loc[temp_df.Description.str.count(',')==2,].reset_index(drop=True)

            desc_df=pd.DataFrame(temp_df['Description'].str.split(',').values.tolist(),columns=['VEP','p','c'])
            desc_df['VEP'] = desc_df['VEP'].str.replace(r'VEP=','')

            prot_df=pd.DataFrame(desc_df['p'].str.split(':').values.tolist(),columns=['prot0','mut']) #,columns=['protein','p','c',0]).drop([0],axis=1)
            prot_df['mut'] = prot_df['mut'].str.replace(r'p\.','')

            out_df=pd.concat([temp_df,desc_df,prot_df],axis=1,join='inner')
            # conditional replace 
            out_df.loc[out_df['prot0'].isin(uni2names.keys()),'prot0']=out_df['prot0'].map(uni2names)
            # simplified df
            out_df2=out_df.drop(['Ref','Alt','Description','p','c'],axis=1)

            # get detailed table of mutations
            mut_table=pd.concat([out_df2,MutExtract(out_df2,ORF1ab_ref)],axis=1)

            # split df by unique values in 'prot' column
            mut_table_byprot=dict(tuple((mut_table.groupby(['prot']))))

            fasta_list=[]
            for prot in mut_table_byprot.keys():
                seq_out=str(SeqMod(mut_table_byprot[prot],ref_dict[prot]))
                fasta_list.append(SeqRecord(Seq(seq_out),id=prot,description='',name=''))

            SeqIO.write(fasta_list, f"fasta/{identifier}.fasta","fasta")





