""" app.py - Simple MCCS single station data search."""
import os
import pandas as pd
import panel as pn
from bokeh.models import HTMLTemplateFormatter
pn.extension('tabulator')

def read_db():
    # Read data from file
    df = pd.read_csv('db/latest.csv')

    # Select out primary columns to show in UI
    colsel = ['Observation ID', 'Station ID', 'Mode', 'Sub-mode', 'UTC Start', 'LST start (hr)', 'Observer', 'source_name', 'QA', 'n_files', 'size_mb']
    df_view = df[colsel].fillna('-').sort_values('UTC Start', ascending=False)
    return df, df_view

df, df_view = read_db()

#def update_db(event):
#    print("Updating DB")
#    os.system("rm ./db/latest.csv")
#    os.system("rclone copy acacia:/aa05/mccs-spider-search/db/latest.csv ./db/")
#    df, df_view = read_db()
#    table.value = df_view

def show_entry(tab_idx: int=0) -> pd.DataFrame:
    """ Callback via pn.bind() to show all information about selected entry.

    Args:
        tab_idx (int): Table index (table.param.selection)

    Returns full dataframe for given index (pd.DataFrame)
    """
    df_tab = table.selected_dataframe
    obs_col = 'Observation ID'
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


if __name__ == "__main__":
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

    #button = pn.widgets.Button(name='Refresh DB', button_type='primary', icon='refresh')
    #button.on_click(update_db)

    # SKAO colors: #070068 - blueshift navy, #E70068 - Redshift magenta
    tpl = pn.template.FastListTemplate(
         title="MCCS Observations", main=[table, plot], header_background='#FFFFFF', header_color='#070068', logo='assets/logo.png'
    )

    _ = tpl.servable()
    pn.serve(tpl, address='10.151.6.200', port=5555) # run with python app.py --alow-websocket-origin=10.151.6.200:5555
