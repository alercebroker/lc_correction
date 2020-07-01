import pandas as pd

COL_DATA_QUALITY = ['aimage', 'aimagerat', 'bimage', 'bimagerat', 'candid', 'chipsf', 'classtar', 'fid', 'fwhm',
                   'mindtoedge', 'nbad', 'nneg', 'objectId', 'scorr', 'seeratio', 'sky', 'sumrat', 'xpos', 'ypos']
COL_SS_ZTF = ['candid', 'objectId', 'ssdistnr', 'ssmagnr', 'ssnamenr']
COL_PS1_ZTF = ['sgscore1', 'szmag1', 'srmag2', 'nmtchps', 'distpsnr3', 'objectidps2', 'szmag2', 'simag3', 'sgmag1',
               'distpsnr1', 'sgmag2', 'sgscore3', 'sgscore2', 'candid', 'distpsnr2', 'sgmag3', 'srmag3', 'simag1',
               'objectidps1', 'simag2', 'objectidps3', 'szmag3', 'srmag1']
COL_CORRECTED = ['objectId', 'candid', 'mjd', 'fid', 'pid', 'diffmaglim', 'isdiffpos', 'nid', 'ra', 'dec', 'magpsf',
                 'sigmapsf', 'magap', 'sigmagap', 'distnr', 'rb', 'magapbig', 'sigmagapbig', 'rfid', 'magpsf_corr',
                 'sigmapsf_corr', 'sigmapsf_corr_ext', 'corrected', 'dubious', 'has_stamp', 'parent_candid', 'step_id_corr']


def get_data_quality(df):
    return df[COL_DATA_QUALITY]


def get_ss_ztf(df):
    return df[COL_SS_ZTF]


def get_ps1_ztf(df):
    def ps1(group):
        first = group.mjd.min()
        data = group[COL_PS1_ZTF][group.mjd == first]
        return pd.Series(data.iloc[0].to_dict())
    return df.groupby(["objectId", "fid"]).apply(ps1).reset_index()


def get_clean_corrected(df, step_name=None):
    df["has_stamp"] = (df["parent_candid"] == 0)
    df['step_id_corr'] = 'corr_bulk_0.0.1' if step_name is None else step_name
    # df.loc[(df["parent_candid"] == 0), "parent_candid"] = np.nan
    return df.reset_index()[COL_CORRECTED]
