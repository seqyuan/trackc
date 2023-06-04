import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def hex2rgb(value):
    # convert hex to rgb
    value = value.lstrip('#')
    lv = len(value)
    rgb = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
    return np.array([rgb[0]/255, rgb[1]/255, rgb[2]/255, 1])


def colorC(cname="RdBu_r", 
           bottom_color="#ffffff", 
           bad_color="white", 
           over_color="white", 
           under_color="white", 
           alpha=0):
    #cname = ["twilight_shifted", "jet", "RdBu_r", "RdGy_r", "BrBG_r", "hot_r", "Spectral_r"]
    #cname = ["terrain_r", "ocean", "gist_earth", "gist_stern_r", "tab20b", "twilight"]
    from matplotlib.colors import ListedColormap
    if isinstance(cname, str):
        #cmap=plt.get_cmap(cname)
        cmap = get_cmap(cname)
    else:
        cmap=cname

    cmap.set_bad(color=bad_color, alpha=alpha)
    cmap.set_over(color=over_color, alpha=alpha)
    cmap.set_under(color=under_color, alpha=alpha)

    if bottom_color==None:
        return cmap

    bottom_color = hex2rgb(bottom_color)
    newcolors = cmap(np.linspace(0, 1, 256))
    #white = np.array([1, 1, 1, 1])
    newcolors[0, :] = bottom_color
    newcmap = ListedColormap(newcolors)

    return newcmap

fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)
#fruitpunch2 = sns.blend_palette(['white', 'blue'], as_cmap=True)
fruitpunch2 = sns.blend_palette(['white', 'purple'], as_cmap=True)
fruitpunch3 = LinearSegmentedColormap.from_list('fruitpunch3', 
                                             [(0, 'white'),
                                              (0.08, 'w'),
                                              (0.4, 'r'),
                                              (1, '#CF3F35')], N=100)




#-----------------------------------------------------------------

