import os
import pandas as pd
import numpy as np
import argparse


"""
1. load expr files for each species and combine them. 
2. load csv of processed set to a df.
2. join the dfs and save the new set to a csv.
"""
parser = argparse.ArgumentParser(description='Homologs windows processing')
parser.add_argument('--data_path', type=str, default='../../../data/datasets/processed_data_SCPECBS3/homologs',
                    help='location of the partition fasta files.')
parser.add_argument('--dataset_name', type=str, default='SCPECBS3',
                    help='location of the partition fasta files.')
parser.add_argument('--expr_dirpath', type=str, default='../../../data/raw_data/SCPECBS/expr/',
                    help='location of the partition fasta files.')

args=parser.parse_args()
expr_paths = [
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/data_files/S_cerevisiae.dat",
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/data_files/S_pombe.dat",
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/expr_files/E_coli.csv",
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/expr_files/B_subtilis.csv"
        ]


def discritize_expr(df, expr_th_perc = [0,0.25,0.5,0.75,0.9]):
    df = df[~df['expr'].isna()]
    print(len(df))
    df.sort_values(by='expr',inplace=True, ascending=False)
    top10 = int(np.ceil(len(df)*(1-expr_th_perc[-1])))
    sec15 = top10 + int(np.ceil((len(df)*0.15)))
    med_high = sec15 + int(np.ceil((len(df)*(expr_th_perc[3]-expr_th_perc[2]))))
    med_low = med_high + int(np.ceil((len(df)*(expr_th_perc[2]-expr_th_perc[1]))))
    print(top10, sec15, med_high, med_low,len(df))
    df['expr'][:top10] = 'expr_top10'
    df['expr'][top10:sec15] = 'expr_pre75_90'
    df['expr'][sec15:med_high] = 'expr_pre50_75'
    df['expr'][med_high:med_low] = 'expr_pre25_50'
    df['expr'][med_low:] = 'expr_low25'
    return df

def combine_expr_data(expr_dirpath):
    """ Assumption:  all expression files in expr_dirpath have the same format and column names.
    """
    expr_paths = [os.path.join(expr_dirpath, p) for p in os.listdir(expr_dirpath)]
    print(expr_paths)
    exdfs=[]
    for i in range(len(expr_paths)):
        exdf = pd.read_csv(expr_paths[i], delimiter=',')
        exdf = exdf[['ids', 'expr']]
        exdf = discritize_expr(exdf)
        exdfs.append(exdf)
    exdf = pd.concat(exdfs)
    return exdf

def merge_exdf_with_csv(exdf, df):
    exdf.rename(columns={"ids":"qseqid"},inplace=True)
    df = df.join(exdf.set_index('qseqid'), on='qseqid', how='left')
    return df

def add_expr2csvs(data_dirpath, exdf):
    csv_paths = [os.path.join(data_dirpath, p) for p in os.listdir(data_dirpath) if p.endswith('.csv') and not p.endswith('.expr.csv')]
    for csvp in csv_paths:
        print(csvp)
        df = pd.read_csv(csvp, delimiter=',')
        df = merge_exdf_with_csv(exdf,df)
        df.to_csv(csvp[:-4]+".expr.csv",index=False)



exdf = combine_expr_data(args.expr_dirpath)
add_expr2csvs(args.data_path, exdf)
