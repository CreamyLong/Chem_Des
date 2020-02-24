
# coding: utf-8

### ==========================参数设定区================================ ###

in_file_name = 'data/SRC.csv'  # 输入文件路径
smi_column = 2  # smiles所在列的位置(从1开始记)
name_column = 1  # name所在列的位置(从1开始记)，没有的话请设为0
header = 1  # 表头行所在行数，没有的话请设为0
des_type = 'rdkit2d'  # 要计算的描述符类型，可以多选(传入list)
sep = ','  # 文件分隔符
out_file_prefix = None  # 输出文件名的前缀，如果定义为None的话前缀为输入文件名

### =================================================================== ###

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.SaltRemover import SaltRemover
import os
import sys
import padelpy
from time import time


def Sec2Time(seconds):  # convert seconds to time
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return ("{:02d}h:{:02d}m:{:02d}s".format(h, m, s))

def CalDes(in_file_name, smi_column=1, name_column=0, header=1, des_type='all', sep=',', out_file_prefix=None):
    """in_file_name: str，输入文件路径
       smi_column: int，smiles所在列的位置(从1开始记)
       name_column: int，name所在列的位置(从1开始记)，没有的话请设为0
       header: int，表头行所在行数，没有的话请设为0
       des_type: str/list，要计算的描述符类型，可以多选(传入list)，可选项有:
                                                                ‘rdkit2d’——计算rdkit2D描述符(199维)
                                                                ‘morgan’——计算Morgan指纹(ECFP4，1024维)
                                                                ‘maccs’——计算MACCS指纹(167维)
                                                                ‘pubchem’——计算PubChem指纹(881维)
                                                                ‘all’——计算所有以上描述符
       sep: str，文件分隔符
       out_file_prefix: str，输出文件名的前缀，如果定义为None的话前缀为输入文件名"""

    # 输入文件处理
    t0 = time()
    init_df = pd.read_csv(in_file_name, header=header - 1, sep=sep)
    in_df = pd.DataFrame(
        init_df.iloc[:, smi_column - 1].values, columns=['Smiles'])

    if name_column:
        in_df['Name'] = init_df.iloc[:, name_column - 1].values
    else:
        in_df['Name'] = in_df.index.astype(str) + '_' + in_df['Smiles']

    parent_dir = os.path.dirname(in_file_name)
    if type(out_file_prefix) == str:
        out_file_prefix = os.path.join(parent_dir, out_file_prefix)
    else:
        out_file_prefix = in_file_name[:in_file_name.rfind('.')]

    if des_type == 'all':
        des_type = ['rdkit2d', 'morgan', 'maccs', 'pubchem']
    elif type(des_type) == str:
        des_type = [des_type]
    elif (type(des_type) == list) or (type(des_type) == tuple):
        pass
    else:
        print("The des_type '{}' is invalid".format(des_type))
        exit(0)

    invalid_id = []
    #  读入输入文件Smiles
    wt = Chem.SmilesWriter('temp.smi', delimiter='\t', includeHeader=False)
    for i, smi in enumerate(in_df['Smiles']):
        m = Chem.MolFromSmiles(smi)
        if m:
            m.SetProp('_Name', str(in_df.loc[i, 'Name']))
            wt.write(m)
        else:
            invalid_id.append(i)
    wt.close()

    print('\nNumber of input molecules: {}\n'.format(len(in_df)))
    print('Number of invalid molecules: {}\n'.format(len(invalid_id)))
    # in_df.drop(labels=invalid_id,axis=0,inplace=True)
    # in_df.reset_index(drop=True,inplace=True)
    # in_df.to_csv('temp.smi',sep='\t',index=False,header=None)
    mols = Chem.SmilesMolSupplier('temp.smi', smilesColumn=0, nameColumn=1, delimiter='\t', titleLine=False)
    t1 = time()
    print('Preprocessing done, time cost: {}'.format(Sec2Time(t1 - t0)))

    # 计算描述符
    if 'rdkit2d' in des_type:
        t2 = time()
        if os.path.exists('temp.csv'):
            os.remove('temp.csv')
        script = './CalRDKitDes/RDKitCalculateMolecularDescriptors.py'
        in_params = 'smilesColumn,1,smilesNameColumn,2,smilesDelimiter,tab,smilesTitleLine,False'
        os.system('''python {} -i temp.smi --infileParams "{}" -o temp.csv'''.format(script, in_params))
        rdkit2d_df = pd.read_csv('temp.csv').drop(labels='Ipc', axis=1).rename(columns={'MolID': 'Name'})
        # rdkit2d_df.insert(1,'Smiles',in_df['Smiles'].values)
        rdkit2d_path = '{}_rdkit2d.csv'.format(out_file_prefix)
        rdkit2d_df.to_csv(rdkit2d_path, index=False)
        print("Calculated \033[36;1mrdkit2D descriptors\033[0m have been saved in '\033[45;1m{}\033[0m'".format(rdkit2d_path))
        print('Time cost: {}'.format(Sec2Time(time() - t2)))

    if 'morgan' in des_type:
        t2 = time()
        mgfps2 = np.array([list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024,useChirality=True)) for mol in mols])
        mg_df2 = pd.DataFrame(mgfps2, columns=['ECFP4_{}'.format(i) for i in range(mgfps2.shape[1])])
        mg_df2 = pd.concat([in_df['Name'], mg_df2], axis=1, sort=False)
        mg2_path = '{}_ecfp4.csv'.format(out_file_prefix)
        mg_df2.to_csv(mg2_path, index=False)
        print("Calculated \033[36;1mMorgan fingerprint\033[0m have been saved in '\033[45;1m{}\033[0m'".format(mg2_path))
        print('Time cost: {}'.format(Sec2Time(time() - t2)))
    # mgfps4 = np.array([list(AllChem.GetMorganFingerprintAsBitVect(mol,4,nBits=1024,useChirality=True)) for mol in mols])
    # mg_df4 = pd.DataFrame(mgfps4,columns=['ECFP8_{}'.format(i) for i in range(mgfps4.shape[1])])

    if 'maccs' in des_type:
        t2 = time()
        maccsfp = np.array([list(MACCSkeys.GenMACCSKeys(mol)) for mol in mols])
        maccs_df = pd.DataFrame(maccsfp, columns=['MACCS_{}'.format(i) for i in range(maccsfp.shape[1])])
        maccs_df = pd.concat([in_df['Name'], maccs_df], axis=1, sort=False)
        maccs_path = '{}_maccs.csv'.format(out_file_prefix)
        maccs_df.to_csv(maccs_path, index=False)
        print("Calculated \033[36;1mMACCS fingerprint\033[0m have been saved in '\033[45;1m{}\033[0m'".format(maccs_path))
        print('Time cost: {}'.format(Sec2Time(time() - t2)))

    if 'pubchem' in des_type:
        t2 = time()
        if os.path.exists('temp.csv'):
            os.remove('temp.csv')
        temp_df = pd.read_csv('temp.smi', sep='\t', header=None)
        if len(temp_df) > 5000:
            loop_num = len(temp_df) // 5000 + 1
            for i in range(loop_num):
                temp_df.iloc[i * 5000:(i + 1) * 5000].to_csv('temp.smi', sep='\t', index=False, header=False)
                pubchem_path = '{}_pubchem_{}.csv'.format(out_file_prefix, i)
                padelpy.padeldescriptor(mol_dir='temp.smi', fingerprints=True, threads=-1, removesalt=False, d_file=pubchem_path, retainorder=True)
                # pubchem_df = pd.read_csv('temp.csv')
                # pubchem_path = '{}_pubchem_{}.csv'.format(out_file_prefix,i)
                # pubchem_df.to_csv(pubchem_path,index=False)
                print("The {}/{}th calculated \033[36;1mPubChem fingerprint\033[0m have been saved in '\033[45;1m{}\033[0m'".format(i, loop_num, pubchem_path))
                # print('Time cost: {}'.format(Sec2Time(time()-t2)))
            pubchem_dfs = []
            for i in range(loop_num):
                pubchem_path = '{}_pubchem_{}.csv'.format(out_file_prefix, i)
                pubchem_dfs.append(pd.read_csv(pubchem_path))
                os.remove(pubchem_path)
            pubchem_df = pd.concat(pubchem_dfs, axis=0, sort=False)
            pubchem_path = '{}_pubchem.csv'.format(out_file_prefix)
            pubchem_df.to_csv(pubchem_path, index=False)
        else:
            padelpy.padeldescriptor(mol_dir='temp.smi', fingerprints=True, threads=-1, removesalt=False, d_file='temp.csv', retainorder=True)
            pubchem_df = pd.read_csv('temp.csv')
            pubchem_path = '{}_pubchem.csv'.format(out_file_prefix)
            pubchem_df.to_csv(pubchem_path, index=False)
        print("Calculated \033[36;1mPubChem fingerprint\033[0m have been saved in '\033[45;1m{}\033[0m'".format(pubchem_path))
        print('Time cost: {}'.format(Sec2Time(time() - t2)))

    # mg_df4.to_csv('{}_ecfp8.csv'.format(out_file_prefix),index=False)

    # des_df1 = rdkit2d_df.merge(mg_df2,left_on='MolID',right_on='Name').drop(labels=['Name'], axis=1)
    # des_df2 = des_df1.merge(puchem_df,left_on='MolID',right_on='Name').drop(labels=['Name'], axis=1)
    # des_df3 = des_df2.merge(maccs_df,left_on='MolID',right_on='Name').drop(labels=['Name'], axis=1)

    if os.path.exists('temp.csv'):
        os.remove('temp.csv')
    # if os.path.exists('temp.smi'):
    #     os.remove('temp.smi')

    # return des_df3

if __name__ == '__main__':
    CalDes(in_file_name, smi_column=smi_column, name_column=name_column, header=header, des_type=des_type, sep=sep, out_file_prefix=out_file_prefix)
