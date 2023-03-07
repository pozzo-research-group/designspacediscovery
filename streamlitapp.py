import streamlit as st
import altair as alt
import pandas as pd


vendables = pd.read_csv('purchaseable_with_tsne.csv')

#vendables['image'] = vendables['CID'].apply(lambda x: f'http://localhost:8501/app/static/altair_images/{x}.png')


"""
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
"""

source = pd.DataFrame.from_records(
    [{'a': 1, 'b': 1, 'image': 'https://altair-viz.github.io/_static/altair-logo-light.png'},
     {'a': 2, 'b': 2, 'image': 'https://avatars.githubusercontent.com/u/11796929?s=200&v=4'}]
)
chart = alt.Chart(source).mark_circle(size=200).encode(
    x='a',
    y='b',
    tooltip=['image']  # Must be a list containing a field called "image"
)

#st.set_page_config(layout = 'wide')
st.title('t-SNE Embedding of Methacrylamide candidate monomers')
st.altair_chart(chart)
