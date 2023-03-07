import streamlit as st
import altair as alt
import pandas as pd


vendables = pd.read_csv('purchaseable_with_tsne.csv')

vendables['image'] = vendables['CID'].apply(lambda x: f'http://localhost:8501/app/static/altair_images/{x}.png')

chart = alt.Chart(vendables, width = 800, height = 800).mark_circle(size=100).encode(
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
st.altair_chart(chart)
