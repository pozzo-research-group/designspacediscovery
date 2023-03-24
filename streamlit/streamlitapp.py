import streamlit as st
import altair as alt
import pandas as pd


vendables = pd.read_csv('streamlit/TMI_demo_altairdata.csv')

chart = alt.Chart(vendables, width = 600, height = 500).mark_circle(size=100).encode(
    x = 'tsne_1',
    y = 'tsne_2',
    tooltip=['CID', 'image'],
    href = 'link')

chart['usermeta'] = {
    "embedOptions": {
        'loader': {'target': '_blank'}
    }
}

#st.set_page_config(layout = 'wide')
st.title('t-SNE Embedding of 2,3,3-Trimethely Indolenine similar candidates')
st.write('Mouse over chart to explore, click to view Pubchem page for molecule. Not recommended for mobile')
st.altair_chart(chart)
