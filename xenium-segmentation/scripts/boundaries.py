#!/usr/bin/env python
# coding: utf-8

# Cell boundaries

import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
from matplotlib.patches import Rectangle
import numpy as np
import warnings
import geopandas as gpd
import tifffile
import os
import json
from fast_alphashape import alphashape
from shapely.ops import transform
warnings.filterwarnings('ignore')


# Import file

def import_files() -> pd.DataFrame:
    df_nuclear = pd.read_csv("data/transcripts_cellpose.csv")
    df_baysor = pd.read_csv("data/transcripts.csv")
    df_nuclear = df_nuclear[['transcript_id', 'cell_id']]
    df = pd.merge(left=df_baysor, right=df_nuclear, how="left", left_on="transcript_id", right_on='transcript_id')
    return df

df = import_files()

def import_image(path: str):
    file = os.path.join(path, "morphology_mip.ome.tif")
    img = tifffile.imread(file)
    return img

img = import_image("data/xenium")

def get_pixel_size(path: str) -> float:
    file = open(os.path.join(path, "experiment.xenium"))
    experiment = json.load(file)
    pixel_size = experiment['pixel_size']
    return pixel_size

pixel_size = get_pixel_size("data/xenium")


# Subset to a smaller FOV

# Larger slides might take a long time to plot and the image would be too crowded to actually see the boundaries. Hence, we subset it to a smaller field of view (FOV)

max_width = int(os.getenv("WIDTH"))
max_height = int(os.getenv("HEIGHT"))

x_offset = int(os.getenv("X_OFFSET"))
y_offset = int(os.getenv("Y_OFFSET"))

# check boundaries

if max_width > img.shape[1]:
    max_width = img.shape[1]
    x_offset = 0
if max_height > img.shape[0]:
    max_height = img.shape[0]
    y_offset = 0

if (x_offset < 0) and (img.shape[1] > max_width):
    x_offset = round(img.shape[1] /2 - max_width /2)
if (max_width + x_offset) > img.shape[1]:
    x_offset = img.shape[1] - max_width

if (y_offset < 0) and (img.shape[0] > max_height):
    y_offset = round(img.shape[0] /2 - max_height /2)
if (max_height + y_offset) > img.shape[0]:
    y_offset = img.shape[0] - max_height

fig = plt.figure(figsize = (8,6))
ax = fig.add_subplot(111)
ax.imshow(
    img,
    vmin=np.percentile(img, 99)*0.1,
    vmax=np.percentile(img, 99)*1.1,
    cmap=sns.dark_palette("#bfcee3", reverse=False, as_cmap=True)
)
img_size = Rectangle((0,0),img.shape[1], img.shape[0], edgecolor='b', facecolor='none')
fov = Rectangle((x_offset, y_offset), max_width, max_height, edgecolor='r', facecolor='none')
ax.add_patch(img_size)
ax.add_patch(fov)
ax.set_xlim((0, img.shape[1]))
ax.set_ylim((img.shape[0], 0))
plt.gca().set_aspect('equal')

def subset_fov(img, df, max_width, max_height):
    img = img[y_offset:(y_offset+max_height), x_offset:(x_offset + max_width)]

    df = df[
        ((df.x / pixel_size) >= (x_offset - 20)) &
        ((df.x / pixel_size) <= (x_offset + max_width + 20)) &
        ((df.y / pixel_size) >= (y_offset - 20 )) &
        ((df.y / pixel_size) <= (y_offset + max_height + 20))
    ]

    df.x = df.x - (x_offset * pixel_size)
    df.y = df.y - (y_offset * pixel_size)
    return (img, df)

(img, df) = subset_fov(img, df, max_width, max_height)


# Create cell boundaries

# Using alphashapes

def make_alphashape(points: pd.DataFrame, alpha: float):
    points = np.array(points)
    shape = alphashape(points, alpha=alpha)
    return shape

shapes = df[~pd.isnull(df.cell)].groupby("cell")[['x', 'y']].apply(make_alphashape, alpha=0.05)
shapes = gpd.GeoSeries(shapes)

# Plot image

def scale_to_image(x, y):
    return(x/pixel_size, y/pixel_size)

fig = plt.figure(figsize = (15,15))
ax = fig.add_subplot(111)
ax.imshow(
    img,
    vmin=np.percentile(img, 99)*0.1,
    vmax=np.percentile(img, 99)*1.1,
    cmap=sns.dark_palette("#bfcee3", reverse=False, as_cmap=True)
)

ax.set_xlim((0, img.shape[1]))
ax.set_ylim((img.shape[0], 0))

colors = sns.color_palette()[3]
shapes.apply(lambda x: transform(scale_to_image, x)).plot(facecolor=colors, edgecolor='none', alpha=0.2, ax=ax)
shapes.apply(lambda x: transform(scale_to_image, x)).plot(facecolor="none", edgecolor=colors, alpha=0.7,  ax=ax)
