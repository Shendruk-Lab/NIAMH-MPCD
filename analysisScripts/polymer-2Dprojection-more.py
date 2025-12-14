"""
    NIAMH-MPCD
    2D Projection Rendering Script for Nematic MPCD Simulations

    This script generates high-quality 2D visualisations of either the velocity field
    or the nematic director field from 3D simulation data, projected along a chosen axis.
    It overlays scalar order or velocity magnitude maps, director/flow arrows, polymer chains,
    and optional topological defects.

    Projections are automatically handled based on user input (x, y, or z), with arrows
    rendered at physically accurate lengths (no quiver rescaling).

    - Uses `shendrukGroupStyle` for formatting (from https://github.com/Shendruk-Lab/).
      You must install it or remove calls to `shendrukGroupFormat`.

    Originally developed by Tyler N. Shendruk  
    Adapted by Zahra Valei 
"""

import os
import sys
import ast
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from subprocess import call
from pathlib import Path
import json
import shendrukGroupFormat as ed

plt.style.use('shendrukGroupStyle')

# === Constants ===
C = 0.5
KEEP_FRAMES = True
START_FRAME = 0
FINISH_FRAME = 999999999
BITRATE = 5000
FRAMERATE = 6
CODEC = "libx264"
MYLW = 1.0
FS = 25

# === Argument Parsing ===
if len(sys.argv) < 5:
    raise ValueError("Usage: script.py <data_path> <defect:0|1> <field:flowfield|directorfield> <projection:x|y|z>")

data_path = Path(sys.argv[1])
show_defects = int(sys.argv[2])
field_type = sys.argv[3]
projection = sys.argv[4].lower()

if projection not in ['x', 'y', 'z']:
    raise ValueError("Projection must be one of: x, y, z")

proj_dim = {'x': 0, 'y': 1, 'z': 2}[projection]
plot_dims = [0, 1, 2]
plot_dims.remove(proj_dim)
dimX, dimY = plot_dims
labX, labY = [ 'x', 'y', 'z' ][dimX], [ 'x', 'y', 'z' ][dimY]

# === Field-specific parameters ===
if field_type == "flowfield":
    cmap = ed.deepsea
    label = r"Velocity, $\left|\vec{v}\right|$"
    arrow_args = dict(headlength=5, headwidth=3, headaxislength=4.5)
else:
    cmap = ed.plasma
    label = r"Scalar order parameter, $S_c$"
    arrow_args = dict(headlength=0, headwidth=0, headaxislength=0)

# === Load input.json ===
with open(data_path / "input.json", 'r') as f:
    sim_config = json.load(f)

domain = sim_config["domain"]
xyz_size = domain + [1] * (3 - len(domain))

md_mpcd = sim_config.get("stepsMD", 50)
field_out = sim_config["flowOut"] if field_type == "flowfield" else sim_config["dirSOut"]

# === Load md.inp ===
with open(data_path / "md.inp", 'r') as file:
    for line in file:
        if "polyN" in line:
            mono_n = int(line.split()[-1].strip("()"))
        elif "polyM" in line:
            poly_m = int(line.split()[-1].strip("()"))
        elif "stepAvg" in line:
            out_md = int(ast.literal_eval(line.split('=')[1].strip())[-1])

md_per_mpcd = field_out / (out_md / md_mpcd)

atoms_per_frame = mono_n * poly_m

# === Load polymer file ===
print('\tSearching for VMD polymer file ...')

poly_name = next(
    (os.path.join(root, f) for root, _, files in os.walk(data_path)
     for f in files if '-vmd.vtf' in f),
    None
)

if not poly_name:
    raise FileNotFoundError("VMD polymer file not found.")


# === Prepare arrays ===
XYZ = np.zeros((3, *xyz_size))
DIR = np.zeros((3, *xyz_size))
S_full = np.zeros(xyz_size)

cm = []

# === Projection Logging ===
print(f"\tProjection along '{projection}' (axis {proj_dim})")
print(f"\tPlotting: {labX} vs {labY}")
print(f"\tPlot size: {xyz_size[dimX]} * {xyz_size[dimY]}")

# === Nematic Defect Drawing ===
def drawNematicDefect(position, charge, orientation, multiplier):
    """
    Draws a nematic defect with a given position, charge, and orientation.
    Positive charge (+1/2) is represented with one line and a circle,
    Negative charge (-1/2) with three lines and a circle.
    """
    def draw_line(coord, theta, dtheta, start_offset, length, colour, lw):
        x0 = coord[0] + start_offset * np.cos(theta + dtheta)
        y0 = coord[1] + start_offset * np.sin(theta + dtheta)
        x1 = coord[0] + (start_offset + length) * np.cos(theta + dtheta)
        y1 = coord[1] + (start_offset + length) * np.sin(theta + dtheta)
        return plt.plot([x0, x1], [y0, y1], color=colour, lw=lw, zorder=3)

    def draw_circle(coord, radius, colour, lw):
        angles = np.linspace(0, 2 * np.pi, 150)
        x = coord[0] + radius * np.cos(angles)
        y = coord[1] + radius * np.sin(angles)
        plt.fill(x, y, 'w', zorder=3)
        return plt.plot(x, y, color=colour, lw=lw, zorder=3)

    # Set default visual parameters
    length = 1.5 * multiplier
    radius = 0.5 * multiplier
    linewidth = 2 * multiplier

    if multiplier == 2:
        length = 1.2 * multiplier
        radius = 1.0 * multiplier
        linewidth = 1.5 * multiplier

    if np.sign(charge) > 0:
        # +1/2 defect
        colour = ed.limegreen
        draw_line(position, orientation, 0.0, radius, length, colour, linewidth)
        draw_circle(position, radius, colour, linewidth)

    elif np.sign(charge) < 0:
        # –1/2 defect
        colour = ed.saphire
        angles = [0.0, 2 * np.pi / 3, 4 * np.pi / 3]
        for angle in angles:
            draw_line(position, orientation, angle, radius, length, colour, linewidth)
        draw_circle(position, radius, colour, linewidth)

# === Open Data Files ===
field_file = open(data_path / f"{field_type}.dat", "r")
for _ in range(13): field_file.readline()
poly_file = open(poly_name, "r")
while poly_file.readline()[0] == "#": pass
poly_file.readline(); poly_file.readline()
POLY = np.zeros((3, mono_n))

if show_defects:
    defect_file = open(data_path / "defects.dat", "r")
    for _ in range(14): defect_file.readline()

i, j, frame_id = 0, 0, -1
while True:
    i += 1
    line = field_file.readline()
    if not line: break

    tokens = line.strip().split()
    if len(tokens) == 8:
        _, qx, qy, qz, vx, vy, vz, sval = tokens
    else:
        _, qx, qy, qz, vx, vy, vz = tokens
        sval = np.sqrt(float(vx)**2 + float(vy)**2 + float(vz)**2)

    x, y, z = map(int, (qx, qy, qz))
    XYZ[:, x, y, z] = [float(qx)+0.5, float(qy)+0.5, float(qz)+0.5]
    DIR[:, x, y, z] = list(map(float, (vx, vy, vz)))
    S_full[x, y, z] = float(sval)

    index = (x, y, z)
    plot_index = (index[dimX], index[dimY])


    if i == np.prod(xyz_size):
        # Read polymer positions
        poly_file.readline(); poly_file.readline()
        for k in range(mono_n):
            _, Qx, Qy, Qz = poly_file.readline().split()
            POLY[:, k] = float(Qx), float(Qy), float(Qz)

        for _ in range(int(md_per_mpcd - 1)):
            poly_file.readline(); poly_file.readline()
            for _ in range(mono_n): poly_file.readline()

        # Unwrap and compute CM
        pos = np.zeros((mono_n, 3))
        pos[0] = POLY[:, 0]
        for m in range(1, mono_n):
            for d in range(3):
                dx = POLY[d, m] - POLY[d, m - 1]
                if dx > 0.5 * xyz_size[d]: dx -= xyz_size[d]
                elif dx < -0.5 * xyz_size[d]: dx += xyz_size[d]
                pos[m, d] = pos[m - 1, d] + dx

        cm_val = np.mean(pos, axis=0)
        cm.append(cm_val)

        if j >= START_FRAME:
            fig, ax = plt.subplots()
            ax.set_xlim(0, xyz_size[dimX])
            ax.set_ylim(0, xyz_size[dimY])
            ax.set_aspect('equal')

            slice_idx = int(cm_val[proj_dim]) % xyz_size[proj_dim]

            qX = np.take(XYZ[dimX], slice_idx, axis=proj_dim)
            qY = np.take(XYZ[dimY], slice_idx, axis=proj_dim)
            dX = np.take(DIR[dimX], slice_idx, axis=proj_dim)
            dY = np.take(DIR[dimY], slice_idx, axis=proj_dim)
            S  = np.take(S_full,   slice_idx, axis=proj_dim)

            ax.quiver(qX, qY, dX, dY, color="black", zorder=2, scale=1, scale_units='xy', **arrow_args)
            sc = ax.imshow(S.T, cmap=cmap, origin='lower',
                           vmin=0, vmax=S.max(),
                           zorder=1, extent=[0, xyz_size[dimX], 0, xyz_size[dimY]])

            for k in range(mono_n):
                circle = Circle((POLY[dimX, k], POLY[dimY, k]), radius=0.5,
                                color="white", alpha=0.5, zorder=2)
                ax.add_patch(circle)

            if show_defects:
                buff = defect_file.readline().split()
                if not buff: break
                number = int(buff[1])
                def_pos = np.zeros((number, 2))
                def_cha = np.zeros(number)
                def_ori = np.zeros(number)
                for d in range(number):
                    buff = defect_file.readline().split()
                    def_pos[d] = float(buff[dimX]), float(buff[dimY])
                    def_cha[d] = float(buff[2])
                    def_ori[d] = float(buff[3])
                defect_file.readline()

                for p, c, o in zip(def_pos, def_cha, def_ori):
                    drawNematicDefect(p, c, o, multiplier=1)

            cb = fig.colorbar(sc, ax=ax, fraction=0.08, pad=0.04, shrink=0.8)
            cb.ax.set_position([0.8, 0.1, 0.05, 0.8])
            cb.ax.set_ylabel(label, fontsize=FS)
            ax.axis('off')
            plt.savefig(f"{data_path}/frame{j:04d}.png", format='png', dpi=600)
            plt.close()

        frame_id += 1
        j += 1
        i = 0
        if j > FINISH_FRAME: break

field_file.close()
poly_file.close()
if show_defects:
    defect_file.close()

# === Animate ===
output_name = data_path / (
    "flow.mp4" if field_type == "flowfield" else f"director_def{show_defects}.mp4"
)
call(f"rm -f '{output_name}'", shell=True)
call(
    f"ffmpeg -f image2 -r {FRAMERATE} -i {data_path}/frame%04d.png "
    f"-vcodec {CODEC} -b {BITRATE}k -r {FRAMERATE} '{output_name}'",
    shell=True
)

