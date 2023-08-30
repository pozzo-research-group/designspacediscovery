import streamlit as st
import altair as alt
import pandas as pd

vendables = pd.read_csv('streamlit/DES_redox_components_altairdata.csv')

#vendables['image'] = vendables['CID'].apply(lambda x: f'http://localhost:8501/app/static/altair_images/{x}.png')

chart = alt.Chart(vendables, width=600,
                  height=500).mark_circle(size=100).encode(
                      x='tsne_1',
                      y='tsne_2',
                      tooltip=['cid', 'image'],
                      color='overall_score',
                      href='link')

chart['usermeta'] = {"embedOptions": {'loader': {'target': '_blank'}}}

#st.set_page_config(layout = 'wide')
st.title('Redox component candidates for DES flow battery electrolytes')
st.write(
    'Mouse over chart to explore, click to view pubchem page. Not tested on mobile'
)
st.altair_chart(chart)
