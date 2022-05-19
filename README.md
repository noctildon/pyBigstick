# pyBigstick

A python script that generates the input for [BIGSTICK](https://github.com/cwjsdsu/BigstickPublick), which uses shell model (theory & algorithm) and fortran (implementation) to calculate the nuclear structure.

See main.py for the example code, and show.py for data visualization (using streamlit).
Note that [BIGSTICK](https://github.com/cwjsdsu/BigstickPublick) only supports unix-like system, so does pyBigstick
But one can execute it using docker.


```python
# Run pyBigstick calculating the nucleus
python3 main.py
```

```python
# Launch streamlit app for data visualization
streamlit run show.py
```

![project-demo](./assets/project-demo.gif)


# Todos
- [ ] use docker pack BIGSTICK fortran code along with streamlit
- [ ] Write test
- [ ] Write documentation
- [ ] Setup database
- [ ] Deploy the streamlit app
- [x] Write some intros
- [x] Add gif demo
- [x] Add clean button
- [x] Set the pre-selected nuclei
- [x] Visualize the results (.res, .dres) by streamlit