import streamlit as st
import altair as alt
import pandas as pd


vendables = pd.read_csv('purchaseable_with_tsne.csv')

chart = alt.Chart(vendables, width = 600, height = 500).mark_circle(size=100).encode(
    x = 'tsne_1',
    y = 'tsne_2',
    tooltip=['CID', 'image'],
    color='hydro_rank',
    href = 'link')

chart['usermeta'] = {
    "embedOptions": {
        'loader': {'target': '_blank'}
    }
}

#st.set_page_config(layout = 'wide')
st.title('t-SNE Embedding of Methacrylamide candidate monomers')
st.write('Mouse over chart to explore, click to view Pubchem page for molecule. Not recommended for mobile')
st.altair_chart(chart)
