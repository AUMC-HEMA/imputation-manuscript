# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
import os
import argparse
import pandas as pd
import numpy as np
from fcsy import DataFrame

flow_variable1 = ['FITC-A', 'APC-A', 'BV605-A', 'BV786-A']
flow_variable2 = ['PE-A', 'PE-CF594-A', 'BV711-A', 'PC7-A']

nilsson_variable1 = ['Qdot 655-A', 'PE-Cy7-A', 'DAPI-A']
nilsson_variable2 = ['PE-A', 'QDot 800-A', 'QDot705-A', 'APC-A']

mosmann_variable1 = ['Violet H 450/50-A', 'Blue A 710/50-A', 'Blue B 515/20-A',
                     'Red B 710/50-A']
mosmann_variable2 = ['Green A 780/40-A', 'Violet E 585/42-A', 'Green E 575/25-A',
                     'Red C 660/20-A']

def parse_results(path):
    dfs = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('_exprs.fcs'):
                method = file.split('_')[2]
                if method != 'CytoBackBone':
                    exprs = DataFrame.from_fcs(root+'/'+file)
                    exprs['sample_id'] = root.split('/')[-1]
                    exprs['method'] = method
                    # cyCombine results in NAs for clusters with n < 50.
                    # We annotated these events previously with metacluster label -1
                    # However, we want to annotate these events in the ground-truth data as well.
                    if method == 'cyCombine':
                        gt_data = exprs[exprs['imp_state']==0]
                        imp_data = exprs[exprs['imp_state']==1]
                        imp_data = imp_data.sort_values('original_ID')
                        imp_data['cyCombine_NA'] = np.where(imp_data['fSOM_metacluster'] == -1, True, False)
                        gt_data['cyCombine_NA'] = np.where(imp_data['fSOM_metacluster'] == -1, True, False)
                        exprs = pd.concat([gt_data, imp_data])
                    dfs.append(exprs)
    expression_data = pd.concat(dfs)
    return expression_data


def parse_CytoBackBone(gt_path, results_path, ids, variable1, variable2):
    dfs = []
    for sample_id in ids:
        results_prefix = f"{results_path}{sample_id}/{sample_id}"
        gt_prefix = f"{gt_path}{sample_id}/{sample_id}"
        # Prepare the ground-truth data for backtracking
        gt = DataFrame.from_fcs(gt_prefix + "_gt.fcs")
        # Add fSOM clustering data to ground-truth data
        cbb_exprs =  DataFrame.from_fcs(results_prefix + "_CytoBackBone_exprs.fcs")
        gt_cbb_exprs = cbb_exprs[cbb_exprs['imp_state'] == 0]
        gt['fSOM_cluster'] = list(gt_cbb_exprs['fSOM_cluster'])
        gt['fSOM_metacluster'] = list(gt_cbb_exprs['fSOM_metacluster'])
        # Recreate the split datasets
        ff1 = gt[gt['dataset']==1]
        ff2 = gt[gt['dataset']==2]
        # CytoBackBone drops some events above a certain threshold
        # We remove these beforehand based on their indices
        ff1_excl = pd.read_csv(results_prefix +"_FCS1_exclusions.csv", index_col=0)
        ff2_excl = pd.read_csv(results_prefix +"_FCS2_exclusions.csv", index_col=0)
        ff1['drop_ff1'] = list(ff1_excl['V1'])
        ff2['drop_ff2'] = list(ff2_excl['V1'])
        ff1 = ff1[ff1['drop_ff1']==False]
        ff2 = ff2[ff2['drop_ff2']==False]
        ff1 = ff1.reset_index(drop=True)
        ff2 = ff2.reset_index(drop=True)
        # Prepare the CytoBackBone merged data
        cbb = pd.read_csv(results_prefix +"_CBB_meta.csv", index_col=0)
        cbb['idx_a'] = pd.to_numeric(cbb['idx_a'], downcast='integer')
        cbb['idx_b'] = pd.to_numeric(cbb['idx_b'], downcast='integer')
        cbb['iteration'] = pd.to_numeric(cbb['iteration'], downcast='integer')
        # Convert to 0-indexing
        cbb['idx_a'] = cbb['idx_a'] - 1
        cbb['idx_b'] = cbb['idx_b'] - 1
        # Backtrack the CytoBackBone algorithm by moving over the iterations
        gt_data = []
        for i in cbb['iteration'].unique():
            iteration = cbb[cbb['iteration'] == i]
            temp = iteration[['idx_a', 'idx_b', 'iteration']]
            # Cross-match
            # FF1 contains the filled in value for idx_a, but FF2 contains the "true" value...
            for var in variable1:
                temp[var] = list(ff2.iloc[temp['idx_b']][var])
            for var in variable2:
                temp[var] = list(ff1.iloc[temp['idx_a']][var])
            # Drop the rows from ff1, ff2, and reset index
            ff1 = ff1.drop(iteration['idx_a'], axis=0)
            ff1 = ff1.reset_index(drop=True)
            ff2 = ff2.drop(iteration['idx_b'], axis=0)
            ff2 = ff2.reset_index(drop=True)
            gt_data.append(temp)
        gt_data = pd.concat(gt_data)
        gt_data['imp_state'] = 0
        cbb_exprs = cbb_exprs[cbb_exprs['imp_state']==1]
        assert(len(gt_data) == len(cbb_exprs))
        exprs = pd.concat([gt_data, cbb_exprs])
        exprs['sample_id'] = sample_id
        exprs['method'] = 'CytoBackBone'
        dfs.append(exprs)
    cbb_expression_data = pd.concat(dfs)
    return cbb_expression_data


def load_data(gt_path, results_path):
    flow_data = parse_results(results_path)
    # Parse CytoBackBone data
    flow_cbb = parse_CytoBackBone(gt_path, results_path, 
                              list(flow_data['sample_id'].unique()),
                              flow_variable1, flow_variable2)
    flow_data = pd.concat([flow_data, flow_cbb])
    # Drop some unused columns to free up memory
    flow_data = flow_data.drop(['FSC-A', 'FSC-H', 'FSC-W', 'SSC-A', 'SSC-H', 'SSC-W', 'Time', 
                                'idx_a', 'idx_b', 'iteration'], axis=1)
    return flow_data


def parse_Nilsson_results(path):
    dfs = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('_exprs.fcs'):
                method = file.split('_')[2]
                if method != 'CytoBackBone':
                    exprs = DataFrame.from_fcs(root+file)
                    exprs['sample_id'] = 'Nilsson_rare'
                    exprs['method'] = method
                    # cyCombine results in NAs for clusters with n < 50.
                    # We annotated these events previously with metacluster label -1
                    # However, we want to annotate these events in the ground-truth data as well.
                    if method == 'cyCombine':
                        gt_data = exprs[exprs['imp_state']==0]
                        imp_data = exprs[exprs['imp_state']==1]
                        imp_data = imp_data.sort_values('original_ID')
                        imp_data['cyCombine_NA'] = np.where(imp_data['fSOM_metacluster'] == -1, True, False)
                        gt_data['cyCombine_NA'] = np.where(imp_data['fSOM_metacluster'] == -1, True, False)
                        exprs = pd.concat([gt_data, imp_data])
                    dfs.append(exprs)
    expression_data = pd.concat(dfs)
    return expression_data


def parse_Nilsson_CytoBackBone(path, ids, variable1, variable2):
    dfs = []
    for sample_id in ids:
        results_prefix = path + "Nilsson_rare"
        # Prepare the ground-truth data for backtracking
        gt = DataFrame.from_fcs(results_prefix + "_gt.fcs")
        # Add fSOM clustering data to ground-truth data
        cbb_exprs =  DataFrame.from_fcs(results_prefix + "_CytoBackBone_exprs.fcs")
        gt_cbb_exprs = cbb_exprs[cbb_exprs['imp_state'] == 0]
        gt['fSOM_cluster'] = list(gt_cbb_exprs['fSOM_cluster'])
        gt['fSOM_metacluster'] = list(gt_cbb_exprs['fSOM_metacluster'])
        # Recreate the split datasets
        ff1 = gt[gt['dataset']==1]
        ff2 = gt[gt['dataset']==2]
        # CytoBackBone drops some events above a certain threshold
        # We remove these beforehand based on their indices
        ff1_excl = pd.read_csv(results_prefix +"_FCS1_exclusions.csv", index_col=0)
        ff2_excl = pd.read_csv(results_prefix +"_FCS2_exclusions.csv", index_col=0)
        ff1['drop_ff1'] = list(ff1_excl['V1'])
        ff2['drop_ff2'] = list(ff2_excl['V1'])
        ff1 = ff1[ff1['drop_ff1']==False]
        ff2 = ff2[ff2['drop_ff2']==False]
        ff1 = ff1.reset_index(drop=True)
        ff2 = ff2.reset_index(drop=True)
        # Prepare the CytoBackBone merged data
        cbb = pd.read_csv(results_prefix +"_CBB_meta.csv", index_col=0)
        cbb['idx_a'] = pd.to_numeric(cbb['idx_a'], downcast='integer')
        cbb['idx_b'] = pd.to_numeric(cbb['idx_b'], downcast='integer')
        cbb['iteration'] = pd.to_numeric(cbb['iteration'], downcast='integer')
        # Convert to 0-indexing
        cbb['idx_a'] = cbb['idx_a'] - 1
        cbb['idx_b'] = cbb['idx_b'] - 1
        # Backtrack the CytoBackBone algorithm by moving over the iterations
        gt_data = []
        for i in cbb['iteration'].unique():
            iteration = cbb[cbb['iteration'] == i]
            temp = iteration[['idx_a', 'idx_b', 'iteration']]
            # Cross-match
            # FF1 contains the filled in value for idx_a, but FF2 contains the "true" value...
            for var in variable1:
                temp[var] = list(ff2.iloc[temp['idx_b']][var])
            for var in variable2:
                temp[var] = list(ff1.iloc[temp['idx_a']][var])
            # Drop the rows from ff1, ff2, and reset index
            ff1 = ff1.drop(iteration['idx_a'], axis=0)
            ff1 = ff1.reset_index(drop=True)
            ff2 = ff2.drop(iteration['idx_b'], axis=0)
            ff2 = ff2.reset_index(drop=True)
            gt_data.append(temp)
        gt_data = pd.concat(gt_data)
        gt_data['imp_state'] = 0
        cbb_exprs = cbb_exprs[cbb_exprs['imp_state']==1]
        assert(len(gt_data) == len(cbb_exprs))
        exprs = pd.concat([gt_data, cbb_exprs])
        exprs['sample_id'] = sample_id
        exprs['method'] = 'CytoBackBone'
        dfs.append(exprs)
    cbb_expression_data = pd.concat(dfs)
    return cbb_expression_data


def load_Nilsson_data(gt_path, results_path):
    flow_data = parse_Nilsson_results(results_path)
    # Parse CytoBackBone data
    flow_cbb = parse_Nilsson_CytoBackBone(results_path, 
                                          list(flow_data['sample_id'].unique()),
                                          nilsson_variable1, nilsson_variable2)
    flow_data = pd.concat([flow_data, flow_cbb])
    return flow_data


def parse_Mosmann_results(path):
    dfs = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('_exprs.fcs'):
                method = file.split('_')[2]
                if method != 'CytoBackBone':
                    exprs = DataFrame.from_fcs(root+file)
                    exprs['sample_id'] = 'Mosmann_rare'
                    exprs['method'] = method
                    # cyCombine results in NAs for clusters with n < 50.
                    # We annotated these events previously with metacluster label -1
                    # However, we want to annotate these events in the ground-truth data as well.
                    if method == 'cyCombine':
                        gt_data = exprs[exprs['imp_state']==0]
                        imp_data = exprs[exprs['imp_state']==1]
                        imp_data = imp_data.sort_values('original_ID')
                        imp_data['cyCombine_NA'] = np.where(imp_data['fSOM_metacluster'] == -1, True, False)
                        gt_data['cyCombine_NA'] = np.where(imp_data['fSOM_metacluster'] == -1, True, False)
                        exprs = pd.concat([gt_data, imp_data])
                    dfs.append(exprs)
    expression_data = pd.concat(dfs)
    return expression_data


def parse_Mosmann_CytoBackBone(path, ids, variable1, variable2):
    dfs = []
    for sample_id in ids:
        results_prefix = path + "Mosmann_rare"
        # Prepare the ground-truth data for backtracking
        gt = DataFrame.from_fcs(results_prefix + "_gt.fcs")
        # Add fSOM clustering data to ground-truth data
        cbb_exprs =  DataFrame.from_fcs(results_prefix + "_CytoBackBone_exprs.fcs")
        gt_cbb_exprs = cbb_exprs[cbb_exprs['imp_state'] == 0]
        gt['fSOM_cluster'] = list(gt_cbb_exprs['fSOM_cluster'])
        gt['fSOM_metacluster'] = list(gt_cbb_exprs['fSOM_metacluster'])
        # Recreate the split datasets
        ff1 = gt[gt['dataset']==1]
        ff2 = gt[gt['dataset']==2]
        # CytoBackBone drops some events above a certain threshold
        # We remove these beforehand based on their indices
        ff1_excl = pd.read_csv(results_prefix +"_FCS1_exclusions.csv", index_col=0)
        ff2_excl = pd.read_csv(results_prefix +"_FCS2_exclusions.csv", index_col=0)
        ff1['drop_ff1'] = list(ff1_excl['V1'])
        ff2['drop_ff2'] = list(ff2_excl['V1'])
        ff1 = ff1[ff1['drop_ff1']==False]
        ff2 = ff2[ff2['drop_ff2']==False]
        ff1 = ff1.reset_index(drop=True)
        ff2 = ff2.reset_index(drop=True)
        # Prepare the CytoBackBone merged data
        cbb = pd.read_csv(results_prefix +"_CBB_meta.csv", index_col=0)
        cbb['idx_a'] = pd.to_numeric(cbb['idx_a'], downcast='integer')
        cbb['idx_b'] = pd.to_numeric(cbb['idx_b'], downcast='integer')
        cbb['iteration'] = pd.to_numeric(cbb['iteration'], downcast='integer')
        # Convert to 0-indexing
        cbb['idx_a'] = cbb['idx_a'] - 1
        cbb['idx_b'] = cbb['idx_b'] - 1
        # Backtrack the CytoBackBone algorithm by moving over the iterations
        gt_data = []
        for i in cbb['iteration'].unique():
            iteration = cbb[cbb['iteration'] == i]
            temp = iteration[['idx_a', 'idx_b', 'iteration']]
            # Cross-match
            # FF1 contains the filled in value for idx_a, but FF2 contains the "true" value...
            for var in variable1:
                temp[var] = list(ff2.iloc[temp['idx_b']][var])
            for var in variable2:
                temp[var] = list(ff1.iloc[temp['idx_a']][var])
            # Drop the rows from ff1, ff2, and reset index
            ff1 = ff1.drop(iteration['idx_a'], axis=0)
            ff1 = ff1.reset_index(drop=True)
            ff2 = ff2.drop(iteration['idx_b'], axis=0)
            ff2 = ff2.reset_index(drop=True)
            gt_data.append(temp)
        gt_data = pd.concat(gt_data)
        gt_data['imp_state'] = 0
        cbb_exprs = cbb_exprs[cbb_exprs['imp_state']==1]
        assert(len(gt_data) == len(cbb_exprs))
        exprs = pd.concat([gt_data, cbb_exprs])
        exprs['sample_id'] = sample_id
        exprs['method'] = 'CytoBackBone'
        dfs.append(exprs)
    cbb_expression_data = pd.concat(dfs)
    return cbb_expression_data


def load_Mosmann_data(gt_path, results_path):
    flow_data = parse_Mosmann_results(results_path)
    # Parse CytoBackBone data
    flow_cbb = parse_Mosmann_CytoBackBone(results_path, 
                                          list(flow_data['sample_id'].unique()),
                                          mosmann_variable1, mosmann_variable2)
    flow_data = pd.concat([flow_data, flow_cbb])
    return flow_data