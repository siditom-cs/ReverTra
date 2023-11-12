import argparse
import os
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser(description='get subset of sequences with high expression level for CAI calculation.')
parser.add_argument('--species', type=str, default='S_cerevisiae',
                    help='')
parser.add_argument('--data_path', type=str, default='/home/tomer/COnTRA/data/datasets/processed_data_SCPECBS3',
                    help='')
parser.add_argument('--expr_path', type=str, default='/home/tomer/COnTRA/data/raw_data/SCPECBS/expr/S_cerevisiae.csv',
                    help='')
parser.add_argument('--expr_th_perc', type=float, default='0.05',
                    help='')

args = parser.parse_args()

if __name__ == '__main__':
    training_file_path = os.path.join(os.path.join(args.data_path, args.species),args.species+".0.nt.fasta")
    expr_path = args.expr_path

    
    records_train = list(SeqIO.parse(training_file_path, "fasta"))
    train_ids = [r.id for r in records_train]
    print(records_train[0].id)

    df = pd.read_csv(expr_path)
    print(len(df))
    df = df[~df['expr'].isna()]
    print(len(df))
    df = df[df['ids'].isin(train_ids)]
    print(len(df))
    df.sort_values(by='expr',inplace=True, ascending=False)
    df = df.iloc[:int(len(df)*args.expr_th_perc)]
    high_expr_ids = list(df['ids'])

    high_expr_records = []
    for r in records_train:
        if r.id in high_expr_ids:
            high_expr_records.append(r)
    #print(high_expr_records)
    subset_output = os.path.join(os.path.join(args.data_path, args.species),args.species+".subset.e"+str(args.expr_th_perc)+".0.nt.fasta")
    SeqIO.write(high_expr_records, subset_output,format='fasta') 
