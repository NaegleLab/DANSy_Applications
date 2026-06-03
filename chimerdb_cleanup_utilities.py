import liftover
import datetime
import pandas as pd

def fix_datetime_names(row,fusion_col,partner_col, pos):
    
    fusion_name = row[fusion_col]
    partner = row[partner_col]
    if isinstance(partner, datetime.datetime):
        new_name = fusion_name.split('-')[pos]
    else:
        new_name = partner
    return new_name

def check_for_tophat(df):
    '''
    Checks if there is a tophat fusion and a second source of the breakpoints
    '''

    tophat_presence = 'Tophat-Fusion' in df['Source'].tolist()
    multiple = len(df) > 1
    single_tophat = sum(df['Source'] == 'Tophat-Fusion') == 1
    output = all([tophat_presence, multiple, single_tophat])

    return output

def check_tophat_pos(df):
    '''
    Checks the position of the Tophat fusion and compares to the other ones
    '''

    temp_df = df[df['Genome_Build_Version']=='hg19'].copy()
    temp_df['H_position'] = temp_df['H_position'].astype(int)
    temp_df['T_position'] = temp_df['T_position'].astype(int)
    tophat = temp_df[temp_df['Source'] == 'Tophat-Fusion']
    others = temp_df[temp_df['Source'] != 'Tophat-Fusion']
    tp_hpos = tophat['H_position'].astype(int).values
    tp_tpos = tophat['T_position'].astype(int).values

    if len(tp_hpos) == 1:
        h_check = any(others['H_position'] - tp_hpos == 1)
        t_check = any(others['T_position'] - tp_tpos == 1)

    return h_check and t_check

def compare_hg19_hg38(df, lift, tophat_correction = False):
    hg19_df = df[df['Genome_Build_Version'] == 'hg19']
    hg38_df = df[df['Genome_Build_Version'] == 'hg38']
    
    h_chr = hg38_df['H_chr'].unique()[0]
    t_chr = hg38_df['T_chr'].unique()[0]
    h_pos = hg38_df['H_position'].unique()[0]
    t_pos = hg38_df['T_position'].unique()[0]
    res = False
    idx_res = list(df.index)
    for i in range(len(hg19_df)):
        hchr_19 = hg19_df.iloc[i]['H_chr']
        tchr_19 = hg19_df.iloc[i]['T_chr']
        hpos_19 = hg19_df.iloc[i]['H_position']
        tpos_19 = hg19_df.iloc[i]['T_position']

        if tophat_correction:
            hpos_19 += 1
            tpos_19 += 1

        # Now converting
        h_cnv_pos = lift.query(hchr_19, hpos_19)[0]
        t_cnv_pos = lift.query(tchr_19, tpos_19)[0]

        comp_res = inner_hg_build_comp([h_chr, h_pos, t_chr,t_pos],[h_cnv_pos,t_cnv_pos])
        if comp_res:
            res = True
            overlap_res = [hg19_df.index[i]]
            break
    
    if res:
        idx_res = [x for x in idx_res if x not in overlap_res]

    return res, idx_res

def inner_hg_build_comp(hg38_info:list, cnv_info):

    h_chr,h_pos,t_chr,t_pos = hg38_info
    h_cnv_pos, t_cnv_pos = cnv_info

    # Checking values
    hchr_check = h_cnv_pos[0] == h_chr
    hpos_check = h_cnv_pos[1] == h_pos
    tchr_check = t_cnv_pos[0] == t_chr
    tpos_check = t_cnv_pos[1] == t_pos


    return all([hchr_check,hpos_check,tchr_check, tpos_check])

def check_for_collapsing(df:pd.DataFrame, lift:liftover.lifter):
    '''
    Determines if we can collapse the entire group into a single line or if there was a second "isoform" of the fusion gene also detected.
    '''
    #print(df['Genome_Build_Version'].tolist())
    hg38_check = 'hg38' in df['Genome_Build_Version'].tolist()
    hg19_check = 'hg19' in df['Genome_Build_Version'].tolist()
    idx_2_keep = df.index
    code = 0
    if check_for_tophat(df):
        tophat_collapse = check_tophat_pos(df)
        if tophat_collapse:

            # Now find how many positions were there and if the hg38 was present as well
            notp = df[df['Source'] != 'Tophat-Fusion']
            notp = notp.filter(['H_gene','H_chr','H_position', 'T_gene', 'T_chr', 'T_position','Genome_Build_Version','Cancertype','BarcodeID','Frame'], axis =1).drop_duplicates()
            unique_hpos = notp['H_position'].unique()
            unique_tpos = notp['T_position'].unique()

            if len(unique_hpos) == 1 and len(unique_tpos) == 1:
                
                idx_2_keep = notp.index[0]
                code = 1
                
            elif hg38_check:
                hg38_comp, hg38_index = compare_hg19_hg38(notp, lift)
                if hg38_comp:
                    idx_2_keep = hg38_index
                    if len(hg38_index) > 1:
                        code = 8
                    else:
                        hg38_index = hg38_index[0]
                        code = 2
                elif len(notp) != len(df):
                    idx_2_keep = list(notp.index)
                    code = 7
        elif len(df) == 2 and hg38_check: # Instance where the tophat fusion and hg38 are what are being left behind
            hg38_comp, hg38_index = compare_hg19_hg38(df, lift, tophat_correction = True)
            if hg38_comp:
                
                if len(hg38_index) > 1:
                        code = 8
                else:
                    hg38_index = hg38_index[0]
                idx_2_keep = hg38_index
                code =3
    
    elif hg38_check and hg19_check:
        temp_df = df.filter(['H_gene','H_chr','H_position', 'T_gene', 'T_chr', 'T_position','Genome_Build_Version','Cancertype','BarcodeID','Frame'], axis =1).drop_duplicates()
        hg38_comp, hg38_index = compare_hg19_hg38(temp_df, lift)
        if hg38_comp:
            if len(hg38_index) > 1:
                code = 8
            else:
                hg38_index = hg38_index[0]
            idx_2_keep = hg38_index
            code = 4

    elif len(df) > 1:
        # Check if can collapse the values together
        temp_df = df.filter(['H_gene','H_chr','H_position', 'T_gene', 'T_chr', 'T_position','Genome_Build_Version','Cancertype','BarcodeID','Frame'], axis =1)
        if len(temp_df.drop_duplicates()) != len(df):
            idx_2_keep = list(temp_df.drop_duplicates().index)
            code = 5
    
    elif len(df) == 1:
        code = 6
        idx_2_keep = idx_2_keep[0]

    if isinstance(idx_2_keep, pd.Index):
        idx_2_keep = list(idx_2_keep)
    elif not isinstance(idx_2_keep, list):
        idx_2_keep = [idx_2_keep]
    return idx_2_keep,code