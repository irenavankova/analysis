#!/usr/bin/env python3

import os
import glob
import re
import numpy as np
from PIL import Image
import imageio.v3 as iio

# =========================================================================
# 1. CONFIGURATION
# =========================================================================
FILE_PREFIX = "ColSpeed"  # Matches your PLOT_VARIABLE choice ('SIConc' or 'SIVol')
FRIS_LOC = '/Users/ivankova/Desktop/Fris_hr/Fris_plots_prl'
OUTPUT_GIF_PATH = f'{FRIS_LOC}/{FILE_PREFIX}_mosaic_timeline.gif'
DURATION_MS = 250  # Display duration for each frame in milliseconds (250ms = 4 frames/sec)

simulations = {
    '8': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '4': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '2': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2')],
    '1': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')]
}


def extract_date(filename):
    match = re.search(r'_(\d{4}-\d{2}-\d{2})\.png$', filename)
    return match.group(1) if match else None


# =========================================================================
# 2. GATHER FILES BY DATE
# =========================================================================
data = {f'F{k}': {'Spin1': {}, 'Spin6': {}} for k in simulations.keys()}

for Fnum, cases in simulations.items():
    dx = f'F{Fnum}'
    for sec, subsec in cases:
        subsec_str = subsec if sec == 'Spin1' else ''
        snapshot_dir = f'{FRIS_LOC}/snapshots_{FILE_PREFIX}/{dx}_{sec}{subsec_str}'
        pattern = f'{snapshot_dir}/{FILE_PREFIX}_{dx}_{sec}{subsec_str}_*.png'

        for fpath in sorted(glob.glob(pattern)):
            date_str = extract_date(fpath)
            if date_str:
                data[dx][sec][date_str] = fpath

# Create timelines based on data availability
all_spin1_dates = sorted(list(set(d for dx in data for d in data[dx]['Spin1'])))
all_spin6_dates = sorted(list(set(d for dx in data for d in data[dx]['Spin6'])))

# --- Spin6 Only Fallback Condition ---
if not all_spin1_dates:
    print("⚠️ No Spin1 data found. Falling back to Spin6 only (First 5 years / 60 frames).")
    target_timeline = all_spin6_dates[:60]
    spin6_only_mode = True
else:
    # Standard 10-year master timeline (60 months Spin1 + 60 months Spin6)
    target_timeline = all_spin1_dates[:60] + all_spin6_dates[:60]
    spin6_only_mode = False

print(f"Timeline configured: {len(target_timeline)} frames total.")

# =========================================================================
# 3. BUILD MOSAIC PANELS AND SAVE TO GIF
# =========================================================================
gif_frames = []
blank_canvas = None

for idx, date_str in enumerate(target_timeline):
    # Determine phase based on mode and index
    if spin6_only_mode:
        phase = 'Spin6'
    else:
        phase = 'Spin1' if idx < len(all_spin1_dates[:60]) else 'Spin6'

    print(f"Stitching Frame {idx + 1}/{len(target_timeline)} ({phase} - {date_str})")


    # Safe image loader with fallback gray block
    def safe_load(res):
        global blank_canvas
        fpath = data[res][phase].get(date_str)
        if fpath and os.path.exists(fpath):
            return Image.open(fpath)

        # Build placeholder if missing
        if blank_canvas is None:
            for r in ['F1', 'F2', 'F4', 'F8']:
                valid_f = next(iter(data[r][phase].values()), None)
                if valid_f:
                    w, h = Image.open(valid_f).size
                    blank_canvas = Image.new('RGB', (w, h), color=(220, 220, 220))
                    break
        return blank_canvas.copy() if blank_canvas else Image.new('RGB', (800, 600), color=(220, 220, 220))


    img_f1 = safe_load('F1')
    img_f2 = safe_load('F2')
    img_f4 = safe_load('F4')
    img_f8 = safe_load('F8')

    # Normalize dimensions to F1 layout baseline
    w, h = img_f1.size
    img_f2 = img_f2.resize((w, h))
    img_f4 = img_f4.resize((w, h))
    img_f8 = img_f8.resize((w, h))

    # Paste onto a 2x2 Canvas
    grid_img = Image.new('RGB', (2 * w, 2 * h))
    grid_img.paste(img_f1, (0, 0))
    grid_img.paste(img_f2, (w, 0))
    grid_img.paste(img_f4, (0, h))
    grid_img.paste(img_f8, (w, h))

    # Keep a list of arrays in memory
    gif_frames.append(np.array(grid_img))

print("\nWriting out the GIF file... (this may take a minute depending on image sizes)")

# Save the list of images as an animated GIF
iio.imwrite(OUTPUT_GIF_PATH, gif_frames, extension=".gif", duration=DURATION_MS, loop=0)

print(f"\nSuccess! GIF animation successfully saved to:\n--> {OUTPUT_GIF_PATH}")