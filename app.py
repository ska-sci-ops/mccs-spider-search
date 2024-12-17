""" app.py - Simple MCCS single station data search. """

import pandas as pd
import panel as pn
from bokeh.models import HTMLTemplateFormatter
pn.extension('tabulator')

# Read data from file
df = pd.read_csv('db/2024-12-17.csv')

# Convert floats into s string (helps searching)
df['lst_start']    = df['lst_start'].round(1).astype('str')
df['obs_duration'] = df['obs_duration'].round(3).astype('str')

col_names = {
    'lst_start': 'LST start (hr)',
    'obs_duration': 'Duration (s)',
    'obs_id': 'Observation ID',
    'mode': 'Mode',
    'station': 'Station ID',
    'sub_mode': 'Sub-mode',
    'utc_start': 'UTC Start',
    'qa': 'QA',
    'bandwidth': 'Bandwidth (MHz)',
    'observer': 'Observer',
    'reference': 'reference'
    }
df = df.rename(columns=col_names)

# Select out primary columns to show in UI
colsel = ['obs_id', 'station', 'mode', 'sub_mode', 'utc_start', 'lst_start', 'reference', 'observer', 'qa']
colsel = [col_names[c] for c in colsel]
df_view = df[colsel].fillna('-').sort_values('UTC Start', ascending=False)

def show_entry(tab_idx: int=0) -> pd.DataFrame:
    """ Show all information about selected entry.

    Used as a callback via pn.bind()

    Args:
        tab_idx (int): Table index (table.param.selection)

    Returns:
        Full dataframe for given index.
    """
    df_tab = table.selected_dataframe
    obs_col = col_names['obs_id']
    if len(df_tab) >= 1:
        df_sel = df[df[obs_col] == df_tab.iloc[0][obs_col]]
    else:
        df_sel = table.current_view[:1]

    if 'reference' in df_sel.columns:
        href = str(df_sel['reference'].values[0])
        if href.startswith('http'):
            ref_id = href.split('/')[-1]
            ref_str = f"<a href={href}  target='_blank'> {ref_id}</a>"
            df_sel['reference'] = ref_str

    fmts = {
        "reference": HTMLTemplateFormatter(template="<code><%= value %></code>")
    }

    tab_sel = pn.widgets.Tabulator(df_sel.dropna(axis=1), disabled=True, formatters=fmts)
    return tab_sel

tab_config = {
    'theme': 'semantic-ui',
    'show_index': False,
    'selectable': True,
    'disabled': True,
    'header_filters': True,
    'initial_page_size': 12,
    'page_size': 12
}

table = pn.widgets.Tabulator(df_view, **tab_config)

plot = pn.bind(show_entry, tab_idx=table.param.selection)

# #070068 - blueshift navy
# #E70068 - Redshift magenta
tpl = pn.template.FastListTemplate(
     title="MCCS Observations", main=[table, plot], header_background='#FFFFFF', header_color='#070068', logo='assets/logo.png'
)

_ = tpl.servable()
pn.serve(tpl, address='10.151.6.200', port=8765)
