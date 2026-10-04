""" analyse.py -- build the observation-log analysis workbook from db/latest.csv.

Run after spider4.py:  python analyse.py [latest.csv] [out.xlsx]
Needs: pandas, xlsxwriter.
"""

import sys
import pandas as pd

HOURS = range(24)


def add_derived(df: pd.DataFrame) -> pd.DataFrame:
    """ Derived columns, same logic as the formula columns in the hand-made workbook. """
    t = pd.to_datetime(df['UTC Start'], errors='coerce')
    sweep = (df['Mode'] == 'correlator') & (df['Sub-mode'] == 'sweep')
    failed = sweep & (df['n_files'] != 385)

    df['QA'] = failed.map({True: 'Failed', False: ''})
    df['Month Bin'] = t.dt.to_period('M').dt.to_timestamp()
    df['QA Status'] = ''
    df.loc[sweep, 'QA Status'] = failed[sweep].map({True: 'Failed', False: 'Success'})
    df['Successful Sweep Flag'] = (df['QA Status'] == 'Success').astype(int)
    local = t + pd.Timedelta(hours=8)
    df['Perth Local Time (AWST)'] = local
    df['Day/Night (AWST)'] = ((local.dt.hour >= 6) & (local.dt.hour < 18)).map({True: 'Day', False: 'Night'})
    df.loc[local.isna(), 'Day/Night (AWST)'] = ''
    df['LST Hour Bin'] = (pd.to_numeric(df['LST start (hr)'], errors='coerce') % 24).floordiv(1)
    df['UTC Hour Bin'] = t.dt.hour
    return df


def volume_summary(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby(['Mode', 'Sub-mode'], sort=False)
    s = pd.DataFrame({
        'Total Size (MB)': g['size_mb'].sum(),
        'Observation Count': g.size(),
        'QA Failed Count': g['QA'].apply(lambda q: (q == 'Failed').sum()),
    }).reset_index()
    s['Share of Total'] = s['Total Size (MB)'] / s['Total Size (MB)'].sum()
    s['Chart Label'] = s['Mode'] + ' (' + s['Sub-mode'] + ')'
    s['QA Fail Rate'] = s['QA Failed Count'] / s['Observation Count']
    s = s.sort_values('Total Size (MB)', ascending=False, ignore_index=True)
    total = {'Mode': 'Total', 'Total Size (MB)': s['Total Size (MB)'].sum(),
             'Observation Count': s['Observation Count'].sum(), 'Share of Total': s['Share of Total'].sum()}
    return pd.concat([s, pd.DataFrame([total])], ignore_index=True)


def monthly_pivot(df: pd.DataFrame) -> pd.DataFrame:
    d = df.assign(**{'Station ID': df['Station ID'].fillna('(blank)')})
    p = d.pivot_table(index='Month Bin', columns='Station ID', values='Successful Sweep Flag',
                      aggfunc='sum', margins=True, margins_name='Grand Total')
    p.index = [i.strftime('%Y-%m') if isinstance(i, pd.Timestamp) else i for i in p.index]
    return p.rename_axis('Month')


def hour_hist(df: pd.DataFrame, col: str) -> pd.DataFrame:
    sw = df[(df['Mode'] == 'correlator') & (df['Sub-mode'] == 'sweep')]
    sw = sw.assign(**{'Station ID': sw['Station ID'].fillna('(blank)')})
    h = pd.crosstab(sw[col], sw['Station ID']).reindex(HOURS, fill_value=0)
    h['Total'] = h.sum(axis=1)
    return h.rename_axis(col.replace(' Bin', ''))


def main(csv='db/latest.csv', out='db/observation_log_analysis.xlsx'):
    df = pd.read_csv(csv, keep_default_na=False, na_values=[''], dtype={'Mode': str, 'Sub-mode': str})
    df['Mode'] = df['Mode'].fillna('')
    df['Sub-mode'] = df['Sub-mode'].fillna('')   # antenna-bandpass sub-mode is ' ' -- preserved
    df = add_derived(df)

    vol, monthly = volume_summary(df), monthly_pivot(df)
    lst, utc = hour_hist(df, 'LST Hour Bin'), hour_hist(df, 'UTC Hour Bin')

    with pd.ExcelWriter(out, engine='xlsxwriter', datetime_format='yyyy-mm-dd hh:mm:ss') as xw:
        wb = xw.book
        df.to_excel(xw, sheet_name='observation_log_cleaned', index=False)
        vol.to_excel(xw, sheet_name='Data Volume Summary', index=False, startrow=2)
        monthly.to_excel(xw, sheet_name='Correlator Sweep Monthly', startrow=2)
        lst.to_excel(xw, sheet_name='Sweep Coverage Histograms', startrow=3)
        utc.to_excel(xw, sheet_name='Sweep Coverage Histograms', startrow=32)

        ws = xw.sheets['Data Volume Summary']
        n = len(vol) - 1  # exclude Total row
        pie = wb.add_chart({'type': 'pie'})
        pie.add_series({'categories': ['Data Volume Summary', 3, 6, 2 + n, 6],
                        'values': ['Data Volume Summary', 3, 2, 2 + n, 2],
                        'data_labels': {'percentage': True}})
        pie.set_title({'name': 'Data Volume by Mode / Sub-mode (MB)'})
        ws.insert_chart('J3', pie)

        wh = xw.sheets['Sweep Coverage Histograms']
        for first, title, cell in ((4, 'Correlator Sweeps by LST Hour', 'V4'),
                                   (33, 'Correlator Sweeps by UTC Hour', 'V33')):
            ncol = lst.shape[1]   # Total is the last column
            ch = wb.add_chart({'type': 'column'})
            ch.add_series({'categories': ['Sweep Coverage Histograms', first, 0, first + 23, 0],
                           'values': ['Sweep Coverage Histograms', first, ncol, first + 23, ncol]})
            ch.set_title({'name': title})
            ch.set_legend({'none': True})
            wh.insert_chart(cell, ch)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
