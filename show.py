from pyBigstick.nucleus import Nucleus
import streamlit as st
import numpy as np
import plotly.express as px
from barChartPlotly import plotly_barcharts_3d
from PIL import Image


he4_image = Image.open('assets/he4.png')
nucl_image = Image.open('assets/nucl_symbol.png')
table_image = Image.open('assets/table.jpg')
scattering_image = Image.open('assets/scattering.jpeg')
deexcitation_image = Image.open('assets/deexcitation.png')
lvl_image = Image.open('assets/Energy_levels.png')
shells_image = Image.open('assets/shells.png')

bs = '/bigstick/src/bigstick.x'

header_container = st.container()
intro_container = st.container()
bs_container = st.container()
states_container = st.container()
densities_container = st.container()

hide_table_row_index = """
            <style>
            tbody th {display:none}
            .blank {display:none}
            </style>
            """

light_nuclei = ['F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']

with header_container:
    st.title('pyBigstick')
    st.markdown("""This streamlit app visualizes the nuclear transitions calculated by [pyBigstick](https://github.com/noctildon/pyBigstick),
        including the energy levels and the density matrices.""")


with intro_container:
    st.subheader("Basic knowledge about nuclear physics")
    st.markdown('Physicists usually use a symbol + 2 number to represent a unique nucleus')
    st.image(nucl_image, width=500)
    st.markdown('For example, the following is identical to He (Helium) with mass number 4 and atomic number 2, or equivalently 2 protons and 2 neutrons.')
    st.image(he4_image,width=300)
    st.markdown('Atomic number can be determined by element symbol uniquely, so sometimes it is skipped and ignored.')

    st.text('And this is the well-known periodic table')
    st.image(table_image,width=800)

    st.markdown('Experimentalists use neutrinos (an extremely small and light particle) to hit the nucleus. This process is called "scattering".')
    st.image(scattering_image,width=800)
    st.markdown("""Before scattering the nucleus has lowest possbile energy (ground state). After scattering nucleus gain some energy from the neutrinos,
        being called "excited nucleus" or "excited state". Then there is a chance that the excited nucleus would drop back to the ground state by emitting gamma ray.
    """)

    col1, col2 = st.columns(2)
    with col1:
        st.image(deexcitation_image,width=400)
    with col2:
        st.image(lvl_image,width=400)


    st.markdown("""What happen in the scattering is that some of the nucleons get excited to the orbit with high energy.
        The core algorithm of pyBigstick is to iterate all possible combinations and transitions of the nucleons.
        And the density matrices describe how nucleons move among the orbits by talking us a probability-like value.
    """)
    st.image(shells_image,width=700)


with bs_container:
    st.subheader("Control panel")
    st.markdown("""Input the info of the interested nucleus, eg. F19, Be10.
        Not sure which nucleus to pick? check out [this](https://periodictable.com/Isotopes/009.19/index.html).
        (Not all of the nucleus is possible to calculate).""")

    col1_1, col1_2 = st.columns(2)
    with col1_1:
        s1 = st.selectbox('The nucleus to calculate', light_nuclei, index=0)
    with col1_2:
        s2 = st.selectbox('Mass number', range(9,41), index=10)

    col2_1, col2_2 = st.columns(2)
    with col2_1:
        n_states = st.selectbox('Number of states to calculate (more states always need more time)', range(1,7), index=2)
    with col2_2:
        maxiter = st.selectbox('Number of iteration (higher iteration is more accurate on results, but takes longer)', np.arange(50,510,10), index=5)

    s1 = s1.lower()
    nucl_name = f'{s1}{s2}'
    st.text(f'Calculate {nucl_name}...')
    nu = Nucleus(nucl_name, n_states=n_states, maxiter=maxiter)

    if st.button('Clean the previous result and rerun'):
        st.write(f'Successfully clean nucleus {nucl_name}. Running {nucl_name} again...')
        nu.clean()

    if not nu.check():
        nu.script_save()
        nu.prepare()
        nu.run(bs=bs)

    nu.save_results()

with states_container:
    st.subheader('Energy level states')
    st.markdown("""When the scattering happens to a nucleus, the nucleus could be excited to higher state.
        In general, initially the nucleus is in ground state (the state with the lowest energy).
        Then after scattering, it is excited to some higher state with energy higher than ground state.
        We also label the state with n. Ground state has n=1. First excited state has n=2. Second excited has n=3, and so on.""")

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
    st.subheader('Density matrices')
    st.markdown("""The amp (transition amplitdue) in the last column below is (to some degree) proportional to the probability that
        a nucleon moves from one orbit to another, given the condition that the nucleus jumps from one state to another (say from n=1 to n=2).
        Jt and Tt are the spin and isospin of the transition, respectively. They are the attributes of a transition.
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
    st.text('The plot only shows non-zero elements.')


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