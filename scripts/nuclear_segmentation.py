#!/usr/bin/env python
# coding: utf-8

# Nuclear segmentation

# Using the DAPI staining and cellpose

import imageio as io
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from tqdm.notebook import tqdm
import pathlib
from cellpose import models, core
import json

xenium_path = 'data'


# Read in Xenium DAPI
# In this part we import the DAPI OME TIFF, create a max projection of the different layers.


def read_dapi_image(path: str, downscale_factor: int = 2) -> np.ndarray:
    img_fpath = pathlib.Path(os.path.join(path, 'morphology_mip.ome.tif'))
    tif = tifffile.TiffFile(img_fpath)
    img = tif.asarray()
    return downscale_image(img, downscale_factor=downscale_factor)

def downscale_image(img: np.ndarray, downscale_factor: int = 2) -> np.ndarray:
    # Calculate the amount of padding needed for each axis
    pad_height = (downscale_factor - img.shape[0] % downscale_factor) % downscale_factor
    pad_width = (downscale_factor - img.shape[1] % downscale_factor) % downscale_factor

    # Pad the array with zeros
    img = np.pad(img, ((0, pad_height), (0, pad_width)), mode='constant')
    return img


maxed_xenium = read_dapi_image(xenium_path, downscale_factor=1)


# Run cellpose 
# Here, we use the pretrained model to perform a nuclear segmentation with cellpose.

def run_cellpose(img: np.ndarray, model_path: str) -> (np.ndarray, np.ndarray, np.ndarray):
    use_GPU = core.use_gpu()
    model = models.CellposeModel(gpu=use_GPU, pretrained_model= model_path  )
    channels = [0,0]
    masks, flows, styles = model.eval([img], channels=channels, diameter=model.diam_labels,flow_threshold=0, cellprob_threshold=0)
    return (masks, flows, styles)

masks, flows, styles = run_cellpose(
    maxed_xenium,
    model_path = r'models/DAPI'
)


# Plot and save segmentation

plt.imshow(masks[0])


# Add the new segmentation to the transcripts.csv

detected_transcripts = pd.read_csv(os.path.join(xenium_path, 'transcripts.csv.gz'))

# Get the pixel to um conversion

def get_pixel_size(path: str) -> float:
    file = open(os.path.join(path, "experiment.xenium"))
    experiment = json.load(file)
    pixel_size = experiment['pixel_size']
    return pixel_size

pixel_size = get_pixel_size(xenium_path)

detected_transcripts['x_location_pixels'] = detected_transcripts.x_location.values*(1/pixel_size)
detected_transcripts['y_location_pixels'] = detected_transcripts.y_location.values*(1/pixel_size)

detected_cells = masks[0][detected_transcripts.y_location_pixels.values.astype(int), detected_transcripts.x_location_pixels.values.astype(int)]
detected_transcripts['cell_id'] = detected_cells
detected_transcripts['overlaps_nucleus'] = (detected_cells > 0).astype(int)
detected_transcripts.to_csv("transcripts_cellpose.csv")

