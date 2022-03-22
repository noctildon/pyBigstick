from pyBigstick.nucleus import Nucleus
import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from barChartPlotly import plotly_barcharts_3d


nu = Nucleus('f19')
nu.save_results()

header_container = st.container()
states_container = st.container()
densities_container = st.container()

hide_table_row_index = """
            <style>
            tbody th {display:none}
            .blank {display:none}
            </style>
            """

with header_container:
    st.title('pyBigstick')
    st.markdown("""This streamlit app visualizes the nuclear transitions calculated by [pyBigstick](https://github.com/noctildon/pyBigstick),
        including the energy levels and the density matrices.""")

with states_container:
    st.title('Energy level states')
    st.markdown("""When the scattering happens to a nucleus, the nucleus could be excited to higher state.
        In general, initially the nucleus is in ground state (the state with the lowest energy).
        Then after scattering, it is excited to some higher state with energy higher than ground state.
        We also label the state with a letter n. Ground state has n=1. First excited state has n=2. Second excited has n=3, and so on.""")

    fig = px.bar(nu.states, x='state', y='Ex',
        labels={
            'state': 'State',
            'Ex': 'Excitation energy (MeV)'
    })
    st.plotly_chart(fig, use_container_width=True)

    if st.checkbox('Show all states data'):
        st.text('States')
        st.write(nu.states)


with densities_container:
    st.title('Density matrices')
    st.markdown("""The amp (transition amplitdue) in the last column below is (to some degree) proportional to the probability that
        a nucleon moves from one orbit to another, given the condition that the nucleus jumps from one state to another (say from n=1 to n=2).
        Jt and Tt are the spin and isospin of the transition, respecitvely. They are the attributes of a transition.
        A transition could have multiple values of Jt. Tt can be either 0 or 1. Most of the amp is zero.""")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        statei = st.selectbox('Initial state', nu.states['state'])
    with col2:
        statej = st.selectbox('Final state', nu.states['state'])
    with col3:
        Jt = st.selectbox('Jt', np.unique(nu.densities['Jt']))
    with col4:
        Tt = st.selectbox('Tt', [0,1])


    filter_densities = nu.densities.loc[(nu.densities['statei']==statei) & (nu.densities['statej']==statej) &\
        (nu.densities['Jt']==Jt) & (nu.densities['Tt']==Tt)]

    st.subheader('Non-zero elements')
    st.markdown(hide_table_row_index, unsafe_allow_html=True)
    st.table(filter_densities)

    st.subheader('3D plot of the density matrices')
    st.text('note that the plot only shows non-zero elements.')


    if not filter_densities.empty:
        fig = plotly_barcharts_3d(filter_densities['orba'], filter_densities['orbb'], filter_densities['amp'],
                x_title='orbit a', y_title='orbit b', z_title='amp')
        fig.update_layout(width=700, height=700, yaxis = dict(scaleanchor = 'x'))
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.text('All elements are zero, so the plot is skipped.')


    if st.checkbox('Show all raw densities data'):
        st.text('Density matrices')
        st.write(nu.densities)