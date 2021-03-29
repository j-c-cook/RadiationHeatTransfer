# Jack C. Cook
# Sunday, March 28, 2021

import RadiationHeatTransfer as rht


def main():
    vs3file = 'cube.vs3'
    vertices, shapes = rht.viewgr.read_vs3_file(vs3file, delimiter=' ')  # this file is space delimited

    print(vertices)

    print(shapes)

    fig, ax = rht.viewgr.plot_vs3_shape(vertices, shapes)

    # save a static pdf
    fig.savefig('vs3_plot.pdf')

    # save an interactive html
    fig = rht.viewgr.i_plot_vs3_shape(vertices, shapes)
    fig.write_html('file.html')


if __name__ == '__main__':
    main()
