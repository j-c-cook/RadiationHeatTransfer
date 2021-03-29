# Jack C. Cook
# Sunday, March 28, 2021

"""
viewgr.py
A Python module that dispalys a VSfile geometry.


Walton, G. N.: "Calculation of obstructed view factors by adaptive integration",
Technical Report NISTIR–6925, National Institute of Standards and Technology,
Gaithersburg, MD, 2002.

https://github.com/jasondegraw/View3D
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import plotly.graph_objects as go


def read_vs3_file(path_to_vs3, delimiter='\t'):

    file = open(path_to_vs3, 'r+')
    data = file.readlines()
    file.close()

    vertices = {}
    shapes = {}

    def process_line(chars):
        current_line_data = []
        for i in range(len(chars)):
            try:
                value = float(chars[i].replace('\n', ''))
                current_line_data.append(value)
            except:
                pass
        return current_line_data

    for line in data:
        characters = line.split(delimiter)
        # collect all of the vertices
        if characters[0] == 'V':
            current_line_data = process_line(characters)
            vertice_number = current_line_data[0]
            x = current_line_data[1]
            y = current_line_data[2]
            z = current_line_data[3]
            vertices[int(vertice_number)] = (x, y, z)
        # collect the shapes
        if characters[0] == 'S':
            current_line_data = process_line(characters)
            shape_number = current_line_data[0]
            v1 = current_line_data[1]
            v2 = current_line_data[2]
            v3 = current_line_data[3]
            v4 = current_line_data[4]
            shapes[int(shape_number)] = (int(v1), int(v2), int(v3), int(v4), int(v1))

    return vertices, shapes


def plot_vs3_shape(vertices, shapes):
    fig = plt.figure()

    # syntax for 3-D projection
    ax = plt.axes(projection='3d')

    for shape in shapes:
        x = []
        y = []
        z = []
        for i in range(len(shapes[shape])):
            v_num = shapes[shape][i]
            x_1, y_1, z_1 = vertices[v_num]
            x.append(x_1)
            y.append(y_1)
            z.append(z_1)
        ax.plot3D(x, y, z)
    for v in vertices:
        x,y,z=vertices[v]
        ax.text(x,y,z,str(v))

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    return fig, ax


def i_plot_vs3_shape(vertices, shapes):

    surfaces = []
    for shape in shapes:
        x = []
        y = []
        z = []
        names = []
        for i in range(len(shapes[shape])):
            v_num = shapes[shape][i]
            x_1, y_1, z_1 = vertices[v_num]
            x.append(x_1)
            y.append(y_1)
            z.append(z_1)
            names.append(str(v_num))
        surface = go.Scatter3d(x=x, y=y, z=z, mode='markers+lines+text',
                               text=names)
        surfaces.append(surface)
    fig = go.Figure(data=surfaces)

    return fig
