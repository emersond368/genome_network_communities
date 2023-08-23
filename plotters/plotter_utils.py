
from typing import List

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from Bin import Bin
from algorithms.Community import Community


def plot_heatmap(array, outfile, bin_map: List[Bin], communities=None):

    array_copy = np.copy(array)

    # set the colorscale and color map
    pct = np.percentile(array, 98)
    colorscale = (0.0, pct)
    # colors in color map crediged to Thomas Gilgenast
    cmap = colors.LinearSegmentedColormap.from_list(
        'pvalue_obs', colors=['#666666', '#797979', '#858585', '#a6a6a6',
                              '#d9d9d9', 'white', '#ffb833', 'red', 'black'])

    array_copy[array_copy < colorscale[0]] = colorscale[0]
    array_copy[~np.isfinite(array_copy)] = colorscale[0] - 1
    heatmap = plt.subplot(1, 1, 1)
    heatmap.set_xlim([0, len(array_copy)])
    heatmap.set_ylim([len(array_copy), 0])

    heatmap.imshow(array_copy, interpolation='none',
                   vmin=colorscale[0],
                   vmax=colorscale[1], cmap=cmap,
                   extent=[0, len(array_copy), len(array_copy), 0])

    if communities:
        boundary = Community(
            bin_map[0].chrom, bin_map[0].start,  bin_map[-1].end)

        for i in range(len(communities)):

            community = communities[i]

            start = float(community.start - boundary.start) * \
                len(array_copy) / (boundary.end - boundary.start)
            end = float(community.end - boundary.start) * \
                len(array_copy) / (boundary.end - boundary.start)

            color = '#fc03f8'

            heatmap.plot([end, end], [start, end], color=color, linestyle='-',
                         linewidth=2)
            heatmap.plot([start, end], [start, start], color=color, linestyle='-',
                         linewidth=2)

    # clean ticks
    heatmap.get_xaxis().set_ticks([])
    heatmap.get_yaxis().set_ticks([])

    plt.savefig(outfile)
    plt.clf()
    plt.close()


def make_gif_across_community_resolutions(gammas, region, directory):

    frames = []
    for gamma in gammas:
        heatmap_fname = f"{directory}{region}_{gamma}.png"
        new_frame = Image.open(heatmap_fname)
        frames.append(new_frame)

    gif_fname = f"{directory}{region}.gif"
    frames[0].save(gif_fname, format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=300, loop=0)
    print(f"GIF of results generated {gif_fname}")
