import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

from matplotlib.backends.backend_pdf import PdfPages

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
        
        
def add_arrow(line, direction='right', size=100, color=None):
    """
    add an arrow to the midpoint of each line object in a 2D plot.
    usuage:
    lines = plt.plot(x,y)[0]
    add_arrow(lines)

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    
    copied from: https://stackoverflow.com/a/34018322
    modified 2022/03/10 M. McCre for arrow to midpoint of all, and to all lines instead of first line.
    useful example gallery: https://matplotlib.org/3.1.0/gallery/text_labels_and_annotations/fancyarrow_demo.html
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()
    
    for i in range(len(xdata)-1):
      line.axes.annotate('',
          xytext=(xdata[i], ydata[i]),
          xy=((xdata[i]+xdata[i+1])/2, (ydata[i]+ydata[i+1])/2),
          arrowprops=dict(arrowstyle="-|>", color=color),
          size=size
    )

def add_poslabel(line, size=14, color=None):
    """
    add position label to start of each line
    usuage:
    lines = plt.plot(x,y)[0]
    add_arrow(line)

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    
    modified from: https://stackoverflow.com/a/34018322
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()
    
    for x,y in zip(xdata,ydata):
      line.axes.annotate('%.4f,%.4f'%(x,y), (x, y))

def SinglePagePdf(filename, figs=None, dpi=300):
    '''
    output a multipage pdf with the given filename, for all open plots if figs=none, or for a given list of figures.
    '''
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for ii, fig in figs:
        fig.savefig("%s-%i"%(filename,ii), format='pdf',dpi=dpi)

def MultiPagePdf(filename, figs=None, dpi=300):
    '''
    output a multipage pdf with the given filename, for all open plots if figs=none, or for a given list of figures.
    '''
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf',dpi=dpi)
    pp.close()