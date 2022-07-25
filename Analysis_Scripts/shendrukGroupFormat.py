# Make our own custom colours
from matplotlib.colors import LinearSegmentedColormap as lsc
import matplotlib.pyplot as plt

# To find where to put the shendrukGroupStyle.mplstyle file run the following two lines in python
# import matplotlib as mpl
# print mpl.get_configdir()

# Edinburgh colours
# Official Edinburgh colours (cyclic)
saphire = [0, 50./255., 95./255.]
crimson = [193./255., 0, 67./255.]
capri = [0, 196./255., 223./255.]
amber = [244./255., 170./255., 0]
plum = [129./255., 2./255., 98./255.]
cerulean = [0, 145./255., 181./255.]
ruby = [212./255., 0, 114./255.]
# Official Edinburgh colours (non-cyclic)
cardinal = [172./255., 0, 64./255.]
cinnamon = [205./255., 90./255., 19./255.]
limegreen = [41./255., 188./255., 41./255.]
gold = [141./255., 116./255., 74./255.]
taupe = [110./255., 80./255., 72./255.]
teal = [69./255., 126./255., 129./255.]
forestgreen = [0, 70./255., 49./255.]
mahogany = [106./255., 51./255., 40./255.]
silver = [194./255., 211./255., 223./255.]
oldrose = [184./255., 133./255., 141./255.]
curry = [156./255., 154./255., 0]
cobalt = [0, 80./255., 114./255.]
# Modified Edinburgh colours for colormaps
rubydarker = [197./255., 0, 99./255.]
purple = [56./255., 6./255., 92./255.]
cardinaldarker = [97./255., 0, 36./255.]
ceruleandarker = [0, 113./255., 140./255.]
amberlighter = [240./255., 191./255., 79./255.]
amberbrighter = [245./255., 242./255., 88./255.]
white = [1, 1, 1]
onyx = [15./255., 15./255., 15./255.]
bggrey = [0.95, 0.95, 0.95]

# Font sizes
fontsize=22
fontsize_title=32
fontsize_axes=28
fontsize_tick=fontsize
fontsizeSmall=18
fontsizeTiny=14

# Lists of colours
def merge_two_dicts(x, y):
  z = x.copy()   # start with x's keys and values
  z.update(y)    # modifies z with y's keys and values & returns None
  return z
class clist:
    def __init__(self):
        self.cyclic = { "saphire":saphire, "crimson":crimson, "capri":capri, "amber":amber, "plum":plum, "cerulean":cerulean, "ruby":ruby }
        self.noncyclic = { "cardinal":cardinal, "cinnamon":cinnamon, "limegreen":limegreen, "gold":gold, "taupe":taupe, "teal":teal, "forestgreen":forestgreen, "mahogany":mahogany, "silver":silver, "oldrose":oldrose, "curry":curry, "cobalt":cobalt }
        self.modified = { "rubydarker":rubydarker, "purple":purple, "cardinaldarker":cardinaldarker, "ceruleandarker":ceruleandarker, "amberlighter":amberlighter, "amberbrighter":amberbrighter, "white":white, "onyx":onyx, "bggrey":bggrey }
        self.official = merge_two_dicts(self.cyclic, self.noncyclic)
        self.all = merge_two_dicts(self.official, self.modified)

# GOOD custom colour maps
plasma = lsc.from_list("", [cobalt,plum,rubydarker,cinnamon,amber])
viridis = lsc.from_list("", [purple,ceruleandarker,limegreen,amberbrighter])
inferno = lsc.from_list("", [onyx,purple,plum,crimson,cinnamon,amber,white])
rain = lsc.from_list("", [onyx,saphire,cobalt,cerulean,capri,silver,white])
seismic = lsc.from_list("", [saphire,white,cardinal])
xmas = lsc.from_list("", [cardinaldarker,white,forestgreen])
cool = lsc.from_list("", [capri,ruby])
deepsea = lsc.from_list("", [purple,ceruleandarker,limegreen])

# Worse custom colour maps
bombpops = lsc.from_list("", [cerulean,white,ruby])
hot = lsc.from_list("", [onyx,cardinaldarker,crimson,cinnamon,amber,white])
cherry = lsc.from_list("", [onyx,cardinaldarker,crimson,white])
spring = lsc.from_list("", [crimson,ruby,amber,amberbrighter,white])
autumn = lsc.from_list("", [crimson,cinnamon,amberlighter])

def reverse_colourmap(cmap, name = 'my_cmap_r'):
    # return lsc(name, plt.cm.revcmap(cmap._segmentdata))
    return cmap.reversed()

plasma_r = reverse_colourmap(plasma)
inferno_r = reverse_colourmap(inferno)
seismic_r = reverse_colourmap(seismic)
bombpops_r = reverse_colourmap(bombpops)
hot_r = reverse_colourmap(hot)
deepsea_r = reverse_colourmap(deepsea)
viridis_r = reverse_colourmap(viridis)
rain_r = reverse_colourmap(rain)
cherry_r = reverse_colourmap(cherry)
spring_r = reverse_colourmap(spring)
cool_r = reverse_colourmap(cool)
xmas_r = reverse_colourmap(xmas)
autumn_r = reverse_colourmap(autumn)

def plot_colortable(colorsDict, title="" ):
    cell_width = 212
    cell_height = 52
    swatch_width = 48
    margin = 12
    topmargin = 40

    names = list(colorsDict)
    n = len(names)
    ncols = 3
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * ncols + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    # fig.subplots_adjust(margin/width, margin/height,(width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * ncols)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    # ax.set_title(title, fontsize=24, loc="left", pad=10)
    ax.set_title(title, fontsize=24, loc="left")

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=colorsDict[name], linewidth=18)
    #
    return fig

def errorbar_fill(x, y, yerr=0, alpha=0.25, **kwargs):
    ax = kwargs.pop('ax', plt.gca())
    base_line, = ax.plot(x, y, **kwargs)
    ax.fill_between(x, y+yerr, y-yerr, facecolor=base_line.get_color(), alpha=alpha)
