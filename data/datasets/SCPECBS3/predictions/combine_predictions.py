import pandas as pd
import numpy as np

def combine_eval_data(df, path):
    df = pd.read_csv(expression_path)
    df['prot_name'] = df['ids']
    df1 = pd.read_csv(model_path, delimiter=',')
    df1 = df1[df1['inx'].notna()]
    df1['prot_name'] = df1['query_prot_name']
    df1['inx'] = df1['inx']
    print("Before expr: ", model_path, len(df1))
    df_all = df1.join(df.set_index('prot_name'), on='prot_name', how='inner')

def get_df(path, model_num, prediction_type=None,harmonization=False, cai=None):
    df = pd.read_csv(path)
    df = df[~df['inx'].isnull()]
    if harmonization:
        df['cross_entropy_loss']='NaN'
        df['entropy']='NaN'
    df.rename(columns={"qseqid":"prot_name",
        "query_prot_name":"prot_name",
        "query_species":"organism",
        'qstart':'start',
        'qend':'end',
        'prot_len':'len',
        'accuracy':'accuracy',
        'cross_entropy_loss':'cross_entropy_loss',
        'entropy':'entropy',
        'subject_species':'mimick_org',
        'sseqid':'mimick_prot_name',
        'subject_prot_name':'mimick_prot_name',
        'evalue':'eval'},
        inplace=True)
    df['model_num']=model_num
    
    pred_type = ['mask']*len(df)
    if prediction_type == 'mimick':
        pred_type = [prediction_type if df.iloc[i]['prot_name']!=df.iloc[i]['mimick_prot_name'] else 'mask' for i in range(len(df))] 
        df['prediction_type']=pred_type
        df = df[df['prediction_type']==prediction_type]
    else:# prediction_type == 'mask'
        df['prediction_type']=pred_type
    
    df= df[['prot_name','organism','start','end','len','pident','model_num','prediction_type','accuracy','cross_entropy_loss','entropy','mimick_org','mimick_prot_name','eval']]
    if harmonization:
        df = df[df['prot_name']!=df['mimick_prot_name']]
    if cai!=None:
        print("TEST=", len(df))
        cai_filter = (df['organism']==cai)
        df = df[cai_filter]
    print(len(df))
    return df

def add_expr_data(df, expr_paths):
    exdfs=[]
    for i in range(len(expr_paths)):
        exdf = pd.read_csv(expr_paths[i], delimiter=',')
        exdf = exdf[['ids', 'expr']]
        exdfs.append(exdf)
    exdf = pd.concat(exdfs)
    exdf.rename(columns={"ids":"prot_name"},inplace=True)
    df = df.join(exdf.set_index('prot_name'), on='prot_name', how='left')
    return df

mask_paths=["/home/tomer/CodOpTRM/Models/SCPECBS3/FTTR_SCPECBS3_allpairs_ws75_Msk100_SingleSpecies_100ksteps_noisy-fcodons/best_model_40k/mask_model_eval.csv",
        ]

homologs_paths=[
        "./Finetuned_oSiS_10_ExR/best_model/mask_model_eval.csv",
        "./Finetuned2Steps_homologs_10_ExR/best_model/model_eval_mimic.csv",
        "./Finetuned2Steps_homologs_10_ExR/best_model/mask/mask_model_eval.csv",
        "./Finetuned_oSiS_30_ExR/best_model/mask_model_eval.csv",
        "./Finetuned2Steps_homologs_30_ExR/best_model/model_eval_mimic.csv",
        "./Finetuned2Steps_homologs_30_ExR/best_model/mask/mask_model_eval.csv",
        "./Finetuned_oSiS_50_ExR/best_model/mask_model_eval.csv",
        "./Finetuned2Steps_homologs_50_ExR/best_model/model_eval_mimic.csv",
        "./Finetuned2Steps_homologs_50_ExR/best_model/mask/mask_model_eval.csv",
        "./Finetuned_oSiS_75_ExR/best_model/mask_model_eval.csv",
        "./Finetuned2Steps_homologs_75_ExR/best_model/model_eval_mimic.csv",
        "./Finetuned2Steps_homologs_75_ExR/best_model/mask/mask_model_eval.csv",
        "./Finetuned_oSiS_100_ExR/best_model/mask_model_eval.csv",
        "./Finetuned2Steps_homologs_100_ExR/best_model/model_eval_mimic.csv",
        "./Finetuned2Steps_homologs_100_ExR/best_model/mask/mask_model_eval.csv",
        "./Finetuned_oSiS_150_ExR/best_model/mask_model_eval.csv",
        ]

cai_paths=[
        "/home/tomer/COnTRA/CAI_Evals/S_cerevisiae_0.1.csv",
        "/home/tomer/COnTRA/CAI_Evals/S_cerevisiae_all.csv",
        "/home/tomer/COnTRA/CAI_Evals/S_pombe_0.1.csv",
        "/home/tomer/COnTRA/CAI_Evals/S_pombe_all.csv",
        "/home/tomer/COnTRA/CAI_Evals/E_coli_0.1.csv",
        "/home/tomer/COnTRA/CAI_Evals/E_coli_all.csv",
        "/home/tomer/COnTRA/CAI_Evals/B_subtilis_0.1.csv",
        "/home/tomer/COnTRA/CAI_Evals/B_subtilis_all.csv",
        "/home/tomer/COnTRA/models/expr_model/hrm/HRM_eval.csv",
#        "/home/tomer/CodOpTRM/Models/SCPECBS3/SCPECBS3_CAI/SCPECBS3.cai.S_cerevisiae.csv",
#        "/home/tomer/CodOpTRM/Models/SCPECBS3/SCPECBS3_CAI/SCPECBS3.cai.S_pombe.csv",
#        "/home/tomer/CodOpTRM/Models/SCPECBS3/SCPECBS3_CAI/SCPECBS3.cai.E_coli.csv",
#        "/home/tomer/CodOpTRM/Models/SCPECBS3/SCPECBS3_CAI/SCPECBS3.cai.B_subtilis.csv",
#        "/home/tomer/CodOpTRM/Models/SCPECBS3/SCPECBS3_CAI/SCPECBS3.harmonized.wPreds.csv"
        ]
expr_paths = [
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/data_files/S_cerevisiae.dat",
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/data_files/S_pombe.dat",
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/expr_files/E_coli.csv",
        "/home/tomer/CodOpTRM/data/raw_data/SCPECBS/expr_files/B_subtilis.csv"
        ]

df = pd.DataFrame(columns=['prot_name','organism','start','end','len','pident','model_num','prediction_type','accuracy','cross_entropy_loss','entropy','mimick_org','mimick_prot_name','eval','pid'])
df1 = get_df(mask_paths[0],"MBart_MaskFT_noExpr",'mask')
df61 = get_df(homologs_paths[0],"MBart10_MaskFT_wExpr",'mask')
df71 = get_df(homologs_paths[1],"MBart10_MimickFT", 'mimick')
df81 = get_df(homologs_paths[2],"MBart10_MimickFT", 'mask')
df62 = get_df(homologs_paths[3],"MBart30_MaskFT_wExpr",'mask')
df72 = get_df(homologs_paths[4],"MBart30_MimickFT", 'mimick')
df82 = get_df(homologs_paths[5],"MBart30_MimickFT", 'mask')
df63 = get_df(homologs_paths[6],"MBart50_MaskFT_wExpr",'mask')
df73 = get_df(homologs_paths[7],"MBart50_MimickFT", 'mimick')
df83 = get_df(homologs_paths[8],"MBart50_MimickFT", 'mask')
df64 = get_df(homologs_paths[9],"MBart75_MaskFT_wExpr",'mask')
df74 = get_df(homologs_paths[10],"MBart75_MimickFT", 'mimick')
df84 = get_df(homologs_paths[11],"MBart75_MimickFT", 'mask')
df65 = get_df(homologs_paths[12],"MBart100_MaskFT_wExpr",'mask')
df75 = get_df(homologs_paths[13],"MBart100_MimickFT", 'mimick')
df85 = get_df(homologs_paths[14],"MBart100_MimickFT", 'mask')
df66 = get_df(homologs_paths[15],"MBart150_MaskFT_wExpr",'mask')

df21 = get_df(cai_paths[0],"CAI_0.1", cai='S_cerevisiae')
df22 = get_df(cai_paths[1],"CAI_all", cai='S_cerevisiae')

df31 = get_df(cai_paths[2],"CAI_0.1", cai='S_pombe')
df32 = get_df(cai_paths[3],"CAI_all", cai='S_pombe')

df41 = get_df(cai_paths[4],"CAI_0.1", cai='E_coli')
df42 = get_df(cai_paths[5],"CAI_all", cai='E_coli')

df51 = get_df(cai_paths[6],"CAI_0.1", cai='B_subtilis')
df52 = get_df(cai_paths[7],"CAI_all", cai='B_subtilis')
df9 = get_df(cai_paths[-1],"Harmonized", harmonization=True)

#df = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9])
df = pd.concat([df1,
    df21,df22,
    df31,df32,
    df41,df42,
    df51,df52,
    df61,df71,df81,
    df62,df72,df82,
    df63,df73,df83,
    df64,df74,df84,
    df65,df75,df85,
    df66,
    df9,
    ])
print(df.columns)
print(len(df))

df['perplexity'] = df['cross_entropy_loss']
df['perplexity'] = df['perplexity'].apply(lambda x: np.exp(x) if x != 'NaN' else 'NaN')

#df.to_csv('SCPECBS3.combined_eval_data.no_expr.csv', index=False, sep=' ', na_rep='NAN')
#df = pd.read_csv('SCPECBS3.combined_eval_data.csv')
df = add_expr_data(df, expr_paths)
df.to_csv('SCPECBS3.combined_eval_data.csv', index=False, sep=' ', na_rep='NAN')
 



