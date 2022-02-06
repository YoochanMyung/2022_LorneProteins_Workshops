import sys
import re
import os
import io
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
import time
import csv
import subprocess
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1
from Bio.PDB import *
from Bio import SeqIO
import ast
import numpy as np
from subprocess import Popen, PIPE, check_output
from itertools import permutations
from multiprocessing import Pool

start_time = time.time()
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

def RabinKarp_search(pattern, text):
    d = len(text)
    q = len(text) + len(pattern)
    m = len(pattern)
    n = len(text)
    p = 0
    t = 0
    h = 1
    i = 0
    j = 0

    for i in range(m-1):
        h = (h*d) % q

    # Calculate hash value for pattern and text
    for i in range(m):
        p = (d*p + ord(pattern[i])) % q
        t = (d*t + ord(text[i])) % q

    # Find the match
    for i in range(n-m+1):
        if p == t:
            for j in range(m):
                if text[i+j] != pattern[j]:
                    break

            j += 1
            if j == m:
                # print("Pattern is found at position: " + str(i+1))
                return int(i)

        if i < n-m:
            t = (d*(t-ord(text[i])*h) + ord(text[i+m])) % q

            if t < 0:
                t = t+q

    CDR_dictionary = dict()
    result_list_L_chain = list()
    annotated_seq = runDirectly(fasta, 'chothia')
    buf = StringIO(annotated_seq)
    temp_pd = pd.DataFrame()
    actual_fasta = fasta
    fasta = readFasta(fasta)

    if len(annotated_seq) < 30:
        # print("Given chain doesn't have any ANTIBODY Light sequence.")
        return 'none'
    else:
        try:
            for each in buf.readlines():
                if not each.startswith('#'):
                    stripted = each.strip()
                    trimed = list(filter(None, stripted.split(' ')))

                    if len(trimed) == 3:
                        if trimed[2] != '-':
                            # temp_pd = temp_pd.append([trimed], ignore_index=True)
                            temp_pd = temp_pd.append([[chainID, trimed[1], trimed[2]]], ignore_index=True)
                    if len(trimed) == 4:
                        if trimed[3] != '-':
                            temp_pd = temp_pd.append([[chainID, trimed[1] + trimed[2], trimed[3]]], ignore_index=True)
                            # print(temp_pd)
            temp_pd.index += 1
            temp_pd.columns = ['chain', 'chothia_numbering', 'chothia_fasta']
            chothia_fasta = ''.join(temp_pd['chothia_fasta'].to_list())
            offset = RabinKarp_search(chothia_fasta[0:10],actual_fasta)
            temp_pd.index += offset

            result = pd.merge(temp_pd, fasta, left_index=True, right_index=True, how='outer')
            result['CDR-'] = result.apply(lambda row: numberChothiaCDR_L(row['chothia_numbering']), axis=1)
            result['chothia_numbering'] = result.apply(
               lambda row: addChainLetter('L', row['chothia_numbering']), axis=1)

            CDR_L1 = list()
            CDR_L2 = list()
            CDR_L3 = list()

            for index, row in result.iterrows():
                if row['CDR-'].startswith('CDR'):
                    pdb_numb = fastanPDBMap.query('fasta_numb == {}'.format(row['fasta_numbering']))['pdb_numb'].values[0]
                    result_list_L_chain.extend([pdb_numb])
                    result_list_L_chain.extend([row['fasta_numbering']])

                if row['CDR-'].startswith('CDR-1'):
                    # CDR_L1.append(pdb_numb)
                    CDR_L1.append('{}_{}'.format(chainID,pdb_numb))
                elif row['CDR-'].startswith('CDR-2'):
                    # CDR_L2.append(pdb_numb)
                    CDR_L2.append('{}_{}'.format(chainID,pdb_numb))
                elif row['CDR-'].startswith('CDR-3'):
                    # CDR_L3.append(pdb_numb)
                    CDR_L3.append('{}_{}'.format(chainID,pdb_numb))

            CDR_dictionary['CDR-1'] = CDR_L1
            CDR_dictionary['CDR-2'] = CDR_L2
            CDR_dictionary['CDR-3'] = CDR_L3

            # return result['fasta'].values,result['fasta_numbering'].values,result['chothia_numbering'].values,CDR_dictionary
            return CDR_dictionary

        except (KeyError):
            # print("Given information is wrong")
            return 'error'

def getCDRnumbForHighlight_IMGT(fasta,fastanPDBMap,chainID,pdbID,chainType):
    CDR_dictionary = dict()
    
    annotated_seq = runDirectly(fasta, 'imgt')
    buf = StringIO(annotated_seq)
    temp_pd = pd.DataFrame()
    actual_fasta = fasta
    fasta = readFasta(fasta)

    if len(annotated_seq) < 30:
        return 'none'
    else:
        for each in buf.readlines():
            if not each.startswith('#'):
                stripted = each.strip()
                trimed = list(filter(None, stripted.split(' ')))
                if len(trimed) == 3:
                    if trimed[2] != '-':
                        temp_pd = temp_pd.append({'chainType':trimed[0],'IMGT':trimed[1],'{}_{}'.format(pdbID,chainID):trimed[2]},ignore_index=True)

                if len(trimed) == 4:
                    if trimed[3] != '-':
                        temp_pd = temp_pd.append({'chainType':trimed[0],'IMGT':trimed[1] + trimed[2],'{}_{}'.format(pdbID,chainID):trimed[3]},ignore_index=True)

            if each.startswith('# Domain 2') and not temp_pd.query('chainType == @chainType').empty:
                print("{}_{} has more than single {} chain domains.".format(pdbID,chainID,chainType))
                break

        temp_pd.index += 1
        _fasta = ''.join(temp_pd['{}_{}'.format(pdbID,chainID)].to_list())
        offset = RabinKarp_search(_fasta[0:10],actual_fasta)
        temp_pd.index += offset

        result = pd.merge(temp_pd, fasta, left_index=True, right_index=True, how='outer')
        result['CDR-{}'.format(chainType)] = result.apply(lambda row: numberIMGTCDR(row['IMGT']), axis=1)
        result['imgt_numbering'] = result.apply(lambda row: addChainLetter('{}'.format(chainType), row['IMGT']), axis=1)

        CDR_1 = list()
        CDR_2 = list()
        CDR_3 = list()
        for index, row in result.iterrows():
            if row['CDR-{}'.format(chainType)].startswith('CDR'):
                pdb_numb = fastanPDBMap.query('fasta_numb == {}'.format(row['fasta_numbering']))['pdb_numb'].values[0]

            if row['CDR-{}'.format(chainType)].startswith('CDR-1'):
                CDR_1.append('{}_{}'.format(chainID,pdb_numb))
            elif row['CDR-{}'.format(chainType)].startswith('CDR-2'):
                CDR_2.append('{}_{}'.format(chainID,pdb_numb))
            elif row['CDR-{}'.format(chainType)].startswith('CDR-3'):
                CDR_3.append('{}_{}'.format(chainID,pdb_numb))

        CDR_dictionary['CDR-1'] = CDR_1
        CDR_dictionary['CDR-2'] = CDR_2
        CDR_dictionary['CDR-3'] = CDR_3

        return CDR_dictionary

def getFASTAnPDBnumberMap(input_pdb,chain_id):
    pdb_id = input_pdb.split('/')[-1]
    new_dic = dict()
    nn = list()
    out_list = pd.DataFrame()
    chain_id = str(chain_id)

    for each in new_dic.values():
        nn.append(" ".join(each).split(" "))

    structure = PDBParser(QUIET=True).get_structure('input_pdb', input_pdb)
    tempp_list = list()

    # for model in structure:
    model = structure[0]

    new_number = 1
    for residue in model[chain_id]:
        tempp_list.append([new_number, str(residue)])
        chain_ID = str(model[chain_id]).rsplit('id=')[1][0].strip()
        amino_acid_name = str(residue).split('Residue')[1].strip().split(' ')[0]
        amino_acid_number = str(residue).split('resseq=')[1].split('icode')[0].strip()
        icode_code = str(residue).rsplit('icode=')[1].strip()

        if len(icode_code) != 1:
            icode_code = icode_code[0].strip()
            amino_acid_number = amino_acid_number + icode_code
        out_list = out_list.append({'chain':chain_id,'wild':amino_acid_name,'pdb_numb':amino_acid_number,'fasta_numb':str(new_number)},ignore_index=True)
        new_number += 1

    return out_list

def getFASTAnum(input_pdb, chain_id, amino_acid, res_num):
    pdb_id = input_pdb.split('/')[-1]
    new_dic = {}
    nn = []
    out_list = []

    for each in new_dic.values():
        nn.append(" ".join(each).split(" "))

    structure = PDBParser(QUIET=True).get_structure('input_pdb', input_pdb)
    tempp_list = []

    model = structure[0]
    new_number = 1
    for residue in model[chain_id]:
        tempp_list.append([new_number, str(residue)])
        chain_ID = str(model[chain_id]).rsplit('id=')[1][0].strip()
        amino_acid_name = str(residue).split('<Residue')[1].strip().split(' ')[0]
        amino_acid_number = str(residue).split('resseq=')[1].split('icode')[0].strip()
        icode_code = str(residue).rsplit('icode=')[1].strip()
        if len(icode_code) != 1:
            icode_code = icode_code[0].strip()
            amino_acid_number = amino_acid_number + icode_code
        if res_num == amino_acid_number:
            out_list.extend(
                [pdb_id, chain_id, amino_acid, res_num, "-->", chain_ID, amino_acid_name, new_number])
        new_number += 1

    return out_list[-1]

def numberIMGTCDR(IMGT_number):
    IMGT_number = str(IMGT_number)
    position = int()

    if IMGT_number[-1].isalpha():
        if IMGT_number == 'nan':
            return '-'
        else:
            position = int(IMGT_number[:-1])
    else:
        position = int(IMGT_number)

    if position < 27:
        return "FR1"
    if (27 <= position) & (position <= 38):
        return "CDR-1"
    if (39 <= position) & (position <= 55):
        return "FR2"
    if (56 <= position) & (position <= 65):
        return "CDR-2"
    if (66 <= position) & (position <= 104):
        return "FR3"
    if (105 <= position) & (position <= 117):
        return "CDR-3"
    if 118 <= position:
        return "FR4"
    else:
        return "-"

def addChainLetter(letter, row):
    row = str(row)
    if row == 'nan':
        return '-'
    else:
        return letter.upper() + row

def getAllChains(pdb_file):
    chain_list = []
    p = PDBParser()
    structure = p.get_structure('input_pdb', pdb_file)
    for each in structure[0]:
        chain_list.append(each.get_id())
    return ''.join(chain_list)

def checkAbAg_using_pdbtofasta(pdb_file):
    chain_list = list(getAllChains(pdb_file))
    result_pd = pd.DataFrame()
    H_list = list()
    L_list = list()
    ag_list = list()

    for chain in chain_list:
        fasta = pdbtofasta(pdb_file,chain)
        if fasta:
            annotation_result_temp = runDirectly(fasta, 'chothia')
            buf = StringIO(annotation_result_temp)
            if len(annotation_result_temp) < 30:
                ag_list.extend(chain)
            else:
                try:
                    temp_set = list()
                    for each in buf.readlines():
                        if not each.startswith('#'):
                            stripted = each.strip()
                            trimed = list(filter(None, stripted.split(' ')))
                            temp_set.extend([trimed[0]])
                    temp_set = set(temp_set[0:-1])
                    for each_set in temp_set:
                        if each_set == 'H':
                            H_list.extend(chain)
                        else:
                            L_list.extend(chain)

                except (KeyError):
                    return 'error'
    if len(H_list) > 0:
        if len(H_list) > 0 and len(L_list) == 0:
            if len(ag_list) > 0:
                return ('Nanobody-Ag;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
            else:
                return ('Nanobody;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
        elif H_list == L_list and len(H_list) !=0:
            if len(ag_list) > 0:
                return ('scFv-Ag;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
            else:
                return ('scFv;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
        else:
            if len(ag_list) > 0:
                return ('Fab-Ag;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
            else:
                return ('Fab;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
    elif len(L_list) >0:
        if len(ag_list) > 0:
            return ('Fll-Ag;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
        else:
            return ('Fll;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
    else:
        if len(ag_list) == 1:
            return ('single-chain_Protein;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))
        elif len(ag_list) > 1:
            return ('multi-chain_protein;{};{};{}'.format(','.join(H_list),','.join(L_list),','.join(ag_list)))


def runDirectly(fasta, scheme):
    # only for chothia
    # anarci_run = Popen(["{}/bin/ANARCI".format('/home/local/BHRI/ymyung/anaconda3/envs/mcsm_ab2'), "--sequence", fasta, "--scheme", scheme], stdout=PIPE, stderr=PIPE)
    # anarci_run = Popen(["/usr/local/bin/ANARCI", "--sequence", fasta, "--scheme", scheme], stdout=PIPE, stderr=PIPE) # bio21
    anarci_run = Popen(["/Users/ymyung_m1pro/miniconda3/envs/landscape/bin/ANARCI", "--sequence", fasta, "--scheme", scheme], stdout=PIPE, stderr=PIPE) # m1 mac

    stdout, stderr = anarci_run.communicate();  # print(len(stdout.decode('ascii')))
    return stdout.decode('ascii')

def pdbatomtofasta(pdb_file,selected_chain):

    # for chain in chain_list:
    for record in SeqIO.parse(pdb_file,'pdb-atom'):
        if record.annotations['chain'] == selected_chain:
            fasta = str(record.seq)
            return fasta

def pdbtofasta(pdb_file, selected_chain):
    nucleic_acids=['DA','DC','DG','DT','DI','A','C','G','U','I']
    aa3to1 = {
        'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
        'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
        'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
        'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
        'UNK': 'X'
    }
    chain_pattern = re.compile("^ATOM\s{2,6}\d{1,5}[\sA-Z1-9]{10}([\w])")
    ca_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])")
    nucleic_acid_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s{2}P [\sA]{3}([\s\w]{2})")
    chain_dict = dict()
    chain_list = []

    fp = open(pdb_file, 'r')
    for line in fp.read().splitlines():
        if line.startswith("ENDMDL"):
            break
        if ''.join(chain_pattern.findall(line)) == selected_chain:
            match_list = ca_pattern.findall(line)
            na_match_list = nucleic_acid_pattern.findall(line)
            if na_match_list:
                return 'nucleic_acid'
            elif match_list:
                resn = match_list[0][0]
                chain = match_list[0][1]
                try:
                    if chain in chain_dict:
                        chain_dict[chain] += aa3to1[resn]
                    else:
                        chain_dict[chain] = aa3to1[resn]
                        chain_list.append(chain)
                except:
                    pass
    fp.close()
    result = chain_dict.get(selected_chain[0])
    # print("pdbtofasta:{}".format(result))
    return result

def readFasta(fasta):
    buf = StringIO(fasta)
    new_pd = pd.DataFrame(list(fasta), columns=['fasta'])
    new_pd.index += 1
    new_pd['fasta_numbering'] = new_pd.index

    return new_pd

def getIMGTCDRs_fab(input_pdb):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", PDBConstructionWarning)

        fastanPDBMap = pd.DataFrame()
        result_pd = pd.DataFrame()

        type_of_ab,heavy_chain,light_chain,_ = checkAbAg_using_pdbtofasta(input_pdb).split(';')

        for each_chain in getAllChains(input_pdb):
            fastanPDBMap = fastanPDBMap.append([getFASTAnPDBnumberMap(input_pdb,each_chain)],ignore_index=True, sort=False)
        fastanPDBMap.fasta_numb = fastanPDBMap.fasta_numb.astype(int)
        fastanPDBMap.pdb_numb = fastanPDBMap.pdb_numb.astype(str)
        fastanPDBMap.chain = fastanPDBMap.chain.astype(str)

        result_pd.loc['heavy','ID'] = input_pdb
        result_pd.loc['heavy','type'] = 'Fab'

        CDR_H_dict = getCDRnumbForHighlight_IMGT(pdbtofasta(input_pdb,heavy_chain),fastanPDBMap.query('chain == "{}"'.format(heavy_chain)),heavy_chain,input_pdb,'H')
        if CDR_H_dict != 'none':
            result_pd.loc['heavy','CDR-1'] = ','.join(CDR_H_dict['CDR-1'])
            result_pd.loc['heavy','nCDR-1'] = len(CDR_H_dict['CDR-1'])
            result_pd.loc['heavy','CDR-2'] = ','.join(CDR_H_dict['CDR-2'])
            result_pd.loc['heavy','nCDR-2'] = len(CDR_H_dict['CDR-2'])
            result_pd.loc['heavy','CDR-3'] = ','.join(CDR_H_dict['CDR-3'])
            result_pd.loc['heavy','nCDR-3'] = len(CDR_H_dict['CDR-3'])
            result_pd.loc['heavy','nCDRs'] =  len(CDR_H_dict['CDR-1']) + len(CDR_H_dict['CDR-2']) + len(CDR_H_dict['CDR-3'])
        else:
            result_pd.loc['heavy','heavy_chain'] = 'might_not_H_chain'
            result_pd.loc['heavy','CDR-1'] = 'might_not_H_chain'
            result_pd.loc['heavy','nCDR-1'] = 0
            result_pd.loc['heavy','CDR-2'] = 'might_not_H_chain'
            result_pd.loc['heavy','nCDR-2'] = 0
            result_pd.loc['heavy','CDR-3'] = 'might_not_H_chain'
            result_pd.loc['heavy','nCDR-3'] = 0
            result_pd.loc['heavy','nCDRs'] = 0

        result_pd.loc['light','ID'] = input_pdb
        result_pd.loc['light','type'] = 'Fab'

        CDR_L_dict = getCDRnumbForHighlight_IMGT(pdbtofasta(input_pdb,light_chain),fastanPDBMap.query('chain == "{}"'.format(light_chain)),light_chain,input_pdb,'L')
        if CDR_L_dict != 'none':
            result_pd.loc['light','CDR-1'] = ','.join(CDR_L_dict['CDR-1'])
            result_pd.loc['light','nCDR-1'] = len(CDR_L_dict['CDR-1'])
            result_pd.loc['light','CDR-2'] = ','.join(CDR_L_dict['CDR-2'])
            result_pd.loc['light','nCDR-2'] = len(CDR_L_dict['CDR-2'])
            result_pd.loc['light','CDR-3'] = ','.join(CDR_L_dict['CDR-3'])
            result_pd.loc['light','nCDR-3'] = len(CDR_L_dict['CDR-3'])
            result_pd.loc['light','nCDRs'] =  len(CDR_L_dict['CDR-1']) + len(CDR_L_dict['CDR-2']) + len(CDR_L_dict['CDR-3'])
        else:
            result_pd.loc['light','light_chain'] = 'might_not_L_chain'
            result_pd.loc['light','CDR-1'] = 'might_not_L_chain'
            result_pd.loc['light','nCDR-1'] = 0
            result_pd.loc['light','CDR-2'] = 'might_not_L_chain'
            result_pd.loc['light','nCDR-2'] = 0
            result_pd.loc['light','CDR-3'] = 'might_not_L_chain'
            result_pd.loc['light','nCDR-3'] = 0
            result_pd.loc['light','nCDRs'] = 0

        db_all_pd = pd.DataFrame()
        db_res_pd = pd.DataFrame()
        temp_list = list()

        for line in open(input_pdb,'r').readlines():
            if line.startswith('ATOM'):
                _chain = line.strip()[21:22].strip()
                _resi = line.strip()[22:28].strip()
                _atomi = line.strip()[4:11].strip()
                _pLDDT = line.strip()[61:67].strip()
                temp_list.append({'chain': _chain,'resi':_resi,'atomi':_atomi,'pLDDT':_pLDDT})
        db_all_pd = pd.DataFrame(temp_list)
        db_all_pd.pLDDT = db_all_pd.pLDDT.astype(float)

        for group_id,group in db_all_pd.groupby(by=['chain','resi']):         
            db_res_pd = db_res_pd.append({'chain_resi':'{}_{}'.format(group_id[0],group_id[1]),\
                'pLDDT':group['pLDDT'].mean(),\
                },ignore_index=True)


        CDR_H1_pLDDT = list()
        CDR_H2_pLDDT = list()
        CDR_H3_pLDDT = list()

        if result_pd.loc['heavy','nCDR-1'] > 0:
            query_list = ','.join([each for each in result_pd.loc['heavy','CDR-1'].split(',')])
            for each in getpLDDT(db_res_pd,query_list).values():
                CDR_H1_pLDDT.append(each['avr_pLDDT'])

        if result_pd.loc['heavy','nCDR-2'] > 0:
            query_list = ','.join([each for each in result_pd.loc['heavy','CDR-2'].split(',')])
            for each in getpLDDT(db_res_pd,query_list).values():
                CDR_H2_pLDDT.append(each['avr_pLDDT'])

        if result_pd.loc['heavy','nCDR-3'] > 0:
            query_list = ','.join([each for each in result_pd.loc['heavy','CDR-3'].split(',')])
            for each in getpLDDT(db_res_pd,query_list).values():
                CDR_H3_pLDDT.append(each['avr_pLDDT'])

        CDR_H1_pLDDT = np.array(CDR_H1_pLDDT).astype(float)
        CDR_H2_pLDDT = np.array(CDR_H2_pLDDT).astype(float)
        CDR_H3_pLDDT = np.array(CDR_H3_pLDDT).astype(float)

        result_pd.loc['heavy','CDR-1_min_pLDDT'] = np.min(CDR_H1_pLDDT).round(2)
        result_pd.loc['heavy','CDR-1_avr_pLDDT'] = np.mean(CDR_H1_pLDDT).round(2)
        result_pd.loc['heavy','CDR-1_max_pLDDT'] = np.max(CDR_H1_pLDDT).round(2)
        result_pd.loc['heavy','CDR-2_min_pLDDT'] = np.min(CDR_H2_pLDDT).round(2)
        result_pd.loc['heavy','CDR-2_avr_pLDDT'] = np.mean(CDR_H2_pLDDT).round(2)
        result_pd.loc['heavy','CDR-2_max_pLDDT'] = np.max(CDR_H2_pLDDT).round(2)
        result_pd.loc['heavy','CDR-3_min_pLDDT'] = np.min(CDR_H3_pLDDT).round(2)
        result_pd.loc['heavy','CDR-3_avr_pLDDT'] = np.mean(CDR_H3_pLDDT).round(2)
        result_pd.loc['heavy','CDR-3_max_pLDDT'] = np.max(CDR_H3_pLDDT).round(2)
        result_pd.loc['heavy','CDR-1_std_pLDDT'] = np.std(CDR_H1_pLDDT).round(2)
        result_pd.loc['heavy','CDR-2_std_pLDDT'] = np.std(CDR_H2_pLDDT).round(2)
        result_pd.loc['heavy','CDR-3_std_pLDDT'] = np.std(CDR_H3_pLDDT).round(2)        
        result_pd.loc['heavy','CDR_avr_pLDDT'] = np.mean(np.concatenate((CDR_H1_pLDDT,CDR_H2_pLDDT,CDR_H3_pLDDT),axis=None)).round(2)

        CDR_L1_pLDDT = list()
        CDR_L2_pLDDT = list()
        CDR_L3_pLDDT = list()

        if result_pd.loc['light','nCDR-1'] > 0:
            query_list = ','.join([each for each in result_pd.loc['light','CDR-1'].split(',')])
            for each in getpLDDT(db_res_pd,query_list).values():
                CDR_L1_pLDDT.append(each['avr_pLDDT'])

        if result_pd.loc['light','nCDR-2'] > 0:
            query_list = ','.join([each for each in result_pd.loc['light','CDR-2'].split(',')])
            for each in getpLDDT(db_res_pd,query_list).values():
                CDR_L2_pLDDT.append(each['avr_pLDDT'])

        if result_pd.loc['light','nCDR-3'] > 0:
            query_list = ','.join([each for each in result_pd.loc['light','CDR-3'].split(',')])
            for each in getpLDDT(db_res_pd,query_list).values():
                CDR_L3_pLDDT.append(each['avr_pLDDT'])

        CDR_L1_pLDDT = np.array(CDR_L1_pLDDT).astype(float)
        CDR_L2_pLDDT = np.array(CDR_L2_pLDDT).astype(float)
        CDR_L3_pLDDT = np.array(CDR_L3_pLDDT).astype(float)

        result_pd.loc['light','CDR-1_min_pLDDT'] = np.min(CDR_L1_pLDDT).round(2)
        result_pd.loc['light','CDR-1_avr_pLDDT'] = np.mean(CDR_L1_pLDDT).round(2)
        result_pd.loc['light','CDR-1_max_pLDDT'] = np.max(CDR_L1_pLDDT).round(2)
        result_pd.loc['light','CDR-2_min_pLDDT'] = np.min(CDR_L2_pLDDT).round(2)
        result_pd.loc['light','CDR-2_avr_pLDDT'] = np.mean(CDR_L2_pLDDT).round(2)
        result_pd.loc['light','CDR-2_max_pLDDT'] = np.max(CDR_L2_pLDDT).round(2)
        result_pd.loc['light','CDR-3_min_pLDDT'] = np.min(CDR_L3_pLDDT).round(2)
        result_pd.loc['light','CDR-3_avr_pLDDT'] = np.mean(CDR_L3_pLDDT).round(2)
        result_pd.loc['light','CDR-3_max_pLDDT'] = np.max(CDR_L3_pLDDT).round(2)
        result_pd.loc['light','CDR-1_std_pLDDT'] = np.std(CDR_L1_pLDDT).round(2)
        result_pd.loc['light','CDR-2_std_pLDDT'] = np.std(CDR_L2_pLDDT).round(2)
        result_pd.loc['light','CDR-3_std_pLDDT'] = np.std(CDR_L3_pLDDT).round(2)
        result_pd.loc['light','CDR_avr_pLDDT'] = np.mean(np.concatenate((CDR_L1_pLDDT,CDR_L2_pLDDT,CDR_L3_pLDDT),axis=None)).round(2)

    print(result_pd[['CDR_avr_pLDDT','CDR-1_avr_pLDDT','CDR-2_avr_pLDDT','CDR-3_avr_pLDDT']])
    result_pd.to_csv('CDR_confidence.csv')

def getpLDDT(db_res_pd,query_list):
    output_dict = dict()
    for query in query_list.split(','):
        _temp = db_res_pd.query('chain_resi == @query').to_dict('list')
        output_dict[_temp['chain_resi'][0]] = {'avr_pLDDT':_temp['pLDDT'][0]}
    return output_dict

def print_CDR_pymol(pdb_file):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", PDBConstructionWarning)
        fastanPDBMap = pd.DataFrame()
        type_of_ab,heavy_chain_list,light_chain_list,ag_list = checkAbAg_using_pdbtofasta(pdb_file).split(';')

        for each_chain in getAllChains(pdb_file):
           fastanPDBMap = fastanPDBMap.append([getFASTAnPDBnumberMap(pdb_file,each_chain)],ignore_index=True, sort=False)
        fastanPDBMap.fasta_numb = fastanPDBMap.fasta_numb.astype(int)
        fastanPDBMap.pdb_numb = fastanPDBMap.pdb_numb.astype(str)
        fastanPDBMap.chain = fastanPDBMap.chain.astype(str)
        print("")
        print("PYMOL COMMANDS")
        print("PDB:{}".format(pdb_file))
        print("Type:{}".format(type_of_ab))
        print("Heavy:{}, Light:{}, Ag:{}".format(heavy_chain_list,light_chain_list,ag_list))

        if heavy_chain_list:
            for heavy_chain in heavy_chain_list:
                CDR_H_dict = getCDRnumbForHighlight_IMGT(pdbtofasta(input_pdb,heavy_chain),fastanPDBMap.query('chain == "{}"'.format(heavy_chain)),heavy_chain,input_pdb,'H')
                print("chain {}:".format(heavy_chain))
                for key,values in CDR_H_dict.items():
                    resi = ' or '.join(['resi '+ each.split('_')[1] for each in values])
                    print("sele {}H, chain {} and ({})".format(key,heavy_chain,resi))


        if light_chain_list:
            for light_chain in light_chain_list:
                CDR_L_dict = getCDRnumbForHighlight_IMGT(pdbtofasta(input_pdb,light_chain),fastanPDBMap.query('chain == "{}"'.format(light_chain)),light_chain,input_pdb,'L')
                print("chain {}:".format(light_chain))
                for key,values in CDR_L_dict.items():
                    resi = ' or '.join(['resi '+ each.split('_')[1] for each in values])
                    print("sele {}L, chain {} and ({})".format(key,light_chain,resi))

        print("spectrum b, rainbow_rev, minimum=40, maximum=100")

if __name__ == '__main__':
    input_pdb = sys.argv[1]
    getIMGTCDRs_fab(input_pdb)
    print_CDR_pymol(input_pdb)

