
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def plotter2d( x1,y1,x2,y2,x3,y3,param1,param2,param3,size,dpip):
    """ A helper function to make a graph

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    data1 : array
       The x data

    data2 : array
       The y data

    parameters set the following :
    pyplot.figure(num=None, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True)

    Returns
    -------
    out : list
        list of artists added
    """
#    plt.rcParams['font.family'] = 'serif'
#    plt.rcParams['font.serif'] = 'Bitstream Vera Serif'
#    plt.rcParams['font.monospace'] = 'DejaVu Sans Mono'
    plt.rcParams['font.size'] = 16
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['legend.fontsize'] = 15
    plt.rcParams['figure.titlesize'] = 18
    plt.rcParams['mathtext.default'] = 'regular'

#    print('stylefile to use :%s'%stylefile )
#    plt.style.use(stylefile)
    if size==1:
        plt.figure(figsize=(13.8889,9.72222),dpi=dpip)
    elif size==2:
        plt.figure(figsize=(6,4),dpi=dpip)
    plt.ylabel('Normalized intensity ')

#    plt.ylabel(r'Relative intensity $\alpha_i > \beta_i$')
    plt.xlabel(r'Ramanshift / $cm^{-1}$')
    plt.grid(True)

    out = plt.plot(x1,y1,param1,x2,y2,param2,x3,y3,param3)
    plt.xticks(np.arange(-1200,1800,200))
    plt.yticks(np.arange(-0.1,1.1,.1))

    filename=time.strftime("image-%m%d%Y_%H%M.png")
    print(filename)
    plt.savefig(filename, dpi=dpip)

    plt.clf()    # clear the figure

    # Bar graph , where bar represents the band intensity
    # common parameters:
    plt.figure(figsize=(15.8889,7.72222),dpi=dpip)
    barWidth=6
    plt.grid(True)

    total=len(x1) #+len(x2)+len(x3)
    colorList = ["0"]*total
    for i in range(len(x1)):
        colorList[i]="red"
    #print(len(colorList),total,colorList)
    plt.bar(x1,y1, width=barWidth ,color=colorList,edgecolor=colorList)

    total=len(x2) #+len(x2)+len(x3)
    colorList = ["0"]*total
    for i in range(len(x2)):
        colorList[i]="blue"
    #print(len(colorList),total,colorList)
    plt.bar(x2,y2, width=barWidth ,color=colorList,edgecolor=colorList )

    total=len(x3) #+len(x2)+len(x3)
    colorList = ["0"]*total
    for i in range(len(x3)):
        colorList[i]="green"
    #print(len(colorList),total,colorList)
    plt.bar(x3,y3, width=barWidth ,color=colorList,edgecolor=colorList )

    plt.xticks(np.arange(-1200, 1800,200))
    plt.xticks(rotation=90)

    # Annotation
    #for i in range(len(x1)):
    #    posn=x1[i]
    #    coord=((round(x1[i],2),round(y1[i],2))
        #plt.annotate('%s' %posn, xy=coord, textcoords='data')
        #plt.txt('%s' %posn, xy=coord,rotation=vertical)

    plt.ylabel('Normalized intensity ')
    plt.xlabel(r'Ramanshift / $cm^{-1}$')
    filename2=time.strftime("image-%m%d%Y_%H%M-bar.png")
    plt.savefig(filename2, dpi=dpip)
    return out

#********************************************************************