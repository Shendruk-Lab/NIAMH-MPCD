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

    Developed by Zahra Valei 
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
import argparse

###########################################################
### Set up argsparse
###########################################################
parser = argparse.ArgumentParser(description='Flow field rendering script.')
parser.add_argument("datapath", type=str, help="Path to the data (should be directory with input.json)")
parser.add_argument('fieldType', type=str, help="Field type: 'vel', or 'nem' for velocity or director field, respectively. For a blank canvas, select 'none'.")
parser.add_argument("-a", "--avdim", type=str, help="Dimension to project along: 'x', 'y', or 'z'", default='z')
parser.add_argument("-s","--start", type=int, help="Starting timestep for averaging", default=1)
parser.add_argument("-f","--finish", type=int, help="Finishing timestep for averaging", default=9999999)
parser.add_argument("-l", "--length", type=float, help="Length of director lines", default=0.5)
parser.add_argument("-A", "--myAspect", type=str, help="'auto' or 'equal'", default='equal')
parser.add_argument("-k", "--keepFrames", type=int, help="0=don't keep (delete) frames; 1=keep frames", default=0)
parser.add_argument("-p", "--savePDF", type=int, help="1 saves transparent pdfs for papers, 0 for none", default=0)
parser.add_argument("-d", "--defectData", type=int, help="Show defects (0=False; 1=True)", default="0")
parser.add_argument("-c", "--colourbar", type=int, help="1 adds a colourbar to the plot, 0 for none", default=1)
parser.add_argument("-g", "--groupStyle", type=int,
                    help="0=don't use group style; 1=use group style",
                    default=1)
parser.add_argument("--qx", type=int, help="Only show every qx arrow in x",default=1)
parser.add_argument("--qy", type=int, help="Only show every qy arrow in y",default=1)

args = parser.parse_args()

###########################################################
### Read arguments
###########################################################
print("Arguments:")
for arg, value in vars(args).items():
	print(f"\t{arg}: {value}")
data_path = Path(args.datapath)
show_defects = args.defectData
fieldType = args.fieldType
projection = args.avdim
START_FRAME = args.start
FINISH_FRAME = args.finish
keepFrames = args.keepFrames
savePDF = args.savePDF
defectData = args.defectData
C = args.length
cbFlag = args.colourbar
myAspect = args.myAspect
groupStyle = args.groupStyle

makeTransparent = False # Transparent backgrounds make crappy videos, but look good on webpages

fieldType=fieldType.lower()
if(fieldType=="v" or fieldType=="vel" or fieldType=="velocity"):
  fieldType="vel"
elif(fieldType=="n" or fieldType=="nem" or fieldType=="nematic" or fieldType=="dir" or fieldType=="director"):
  fieldType="nem"
elif(fieldType=="none" or fieldType=="0"):
  fieldType="none"

###########################################################
### Style/formating stuff
###########################################################
if groupStyle:
    plt.style.use('shendrukGroupStyle')
    import shendrukGroupFormat as ed
BITRATE = 5000
FRAMERATE = 6
CODEC = "libx264"
MYLW = 1.0
FS = 25
myAlpha = 1.0

if groupStyle:
    lineColour = ed.silver
    monomerColour = ed.onyx
else:
    lineColour = 'white'
    monomerColour = 'black'

if projection not in ['x', 'y', 'z']:
    raise ValueError("Projection must be one of: x, y, z")

proj_dim = {'x': 0, 'y': 1, 'z': 2}[projection]
plot_dims = [0, 1, 2]
plot_dims.remove(proj_dim)
dimX, dimY = plot_dims
labX, labY = [ 'x', 'y', 'z' ][dimX], [ 'x', 'y', 'z' ][dimY]

# === Field-specific parameters ===
if fieldType == "vel":
    if groupStyle:
        cmap = ed.deepsea
    else:
        cmap = plt.cm.winter
    label = r"Velocity, $\left|\vec{v}\right|$"
    arrow_args = dict(headlength=5, headwidth=3, headaxislength=4.5)
elif fieldType == "nem":
    if groupStyle:
        cmap = ed.plasma
    else:
        cmap = plt.cm.plasma
    label = r"Scalar order, $S$"
    arrow_args = dict(headlength=0, headwidth=0, headaxislength=0)
elif fieldType == "none":
    pass
else:
    print( "Field type not recognized.")
    exit()


# === Load input.json ===
with open(data_path / "input.json", 'r') as f:
    sim_config = json.load(f)

domain = sim_config["domain"]
xyz_size = domain + [1] * (3 - len(domain))

# Figure
width=8
height=6
if not cbFlag:
	width=height
fig1, ax = plt.subplots(figsize=(width, height))
if myAspect == 'auto':
    shrink_factor = 1.0
elif float(xyz_size[dimX])<float(xyz_size[dimY]):
    shrink_factor = float(xyz_size[dimX])/float(xyz_size[dimY])
else:
    shrink_factor = 1.0

md_mpcd = sim_config.get("stepsMD", 50)
if fieldType == "vel":
    field_out = sim_config["flowOut"]
elif fieldType == "nem":
    field_out = sim_config["dirSOut"]
elif fieldType == "none":
    field_out = 1
else:
    print("Field type not recognized.")
    exit()

# === Load md.inp ===
with open(data_path / "md.inp", 'r') as file:
    for line in file:
        if "polyN" in line:
            mono_n = int(line.split()[-1].strip("()"))
        elif "polyM" in line:
            poly_m = int(line.split()[-1].strip("()"))
        elif "stepAvg" in line:
            out_md = int(ast.literal_eval(line.split('=')[1].strip())[-1])

atoms_per_frame = mono_n * poly_m

md_per_mpcd = field_out / (out_md / md_mpcd)

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

POLY = np.zeros((poly_m, mono_n, 3))
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
        if groupStyle:
            colour = ed.limegreen
        else:
            colour = 'g'
        draw_line(position, orientation, 0.0, radius, length, colour, linewidth)
        draw_circle(position, radius, colour, linewidth)

    elif np.sign(charge) < 0:
        # –1/2 defect
        if groupStyle:
            colour = ed.saphire
        else:
            colour = 'b'
        angles = [0.0, 2 * np.pi / 3, 4 * np.pi / 3]
        for angle in angles:
            draw_line(position, orientation, angle, radius, length, colour, linewidth)
        draw_circle(position, radius, colour, linewidth)

# === Open Data Files ===
if fieldType == "vel": field_file = open(data_path / f"flowfield.dat", "r")
elif fieldType == "nem": field_file = open(data_path / f"directorfield.dat", "r")
elif fieldType == "none": field_file = None
else: raise ValueError("Field type must be 'vel' or 'nem'.")

if fieldType != "none":
    for _ in range(13): field_file.readline()
poly_file = open(poly_name, "r")
# Skip header lines
while poly_file.readline()[0] == "#": pass
for _ in range( poly_m + 1 ): poly_file.readline()

if show_defects:
    defect_file = open(data_path / "defects.dat", "r")
    for _ in range(14): defect_file.readline()

i, j, frame_id = 0, 0, -1
while True:
    if fieldType != "none":
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
    else:
        i = np.prod(xyz_size)

    if i == np.prod(xyz_size):
        # Read polymer positions
        poly_file.readline(); poly_file.readline()
        for p in range(poly_m):
            for m in range(mono_n):
                line = poly_file.readline()
                if not line: break
                _, Qx, Qy, Qz = line.split()
                POLY[p, m] = float(Qx), float(Qy), float(Qz)

        for _ in range(int(md_per_mpcd - 1)):
            poly_file.readline(); poly_file.readline()
            for _ in range(atoms_per_frame):
                poly_file.readline()
        
        if not line: break

        # Unwrap and compute CM
        cm_val = np.zeros(3)
        for p in range(poly_m):
            pos = np.zeros((mono_n, 3))
            pos[0] = POLY[p, 0]
            for m in range(1, mono_n):
                for d in range(3):
                    dx = POLY[p, m, d] - POLY[p, m - 1, d]
                    if dx > 0.5 * xyz_size[d]: dx -= xyz_size[d]
                    elif dx < -0.5 * xyz_size[d]: dx += xyz_size[d]
                    pos[m, d] = pos[m - 1, d] + dx
            cm_val += np.mean(pos, axis=0)

        cm_val /= poly_m
        cm.append(cm_val)

        if j >= START_FRAME:
            fig,ax = plt.subplots(figsize=(width, height))
            ax.set_xlim(0, xyz_size[dimX])
            ax.set_ylim(0, xyz_size[dimY])
            ax.set_aspect(myAspect)

            slice_idx = int(cm_val[proj_dim]) % xyz_size[proj_dim]

            if fieldType != "none":
                qX = np.take(XYZ[dimX], slice_idx, axis=proj_dim)
                qY = np.take(XYZ[dimY], slice_idx, axis=proj_dim)
                dX = np.take(DIR[dimX], slice_idx, axis=proj_dim)
                dY = np.take(DIR[dimY], slice_idx, axis=proj_dim)
                S  = np.take(S_full,   slice_idx, axis=proj_dim)

                for x in range(xyz_size[dimX]):
                    for y in range(xyz_size[dimY]):
                        if S[x, y] > 0.05:
                            x0 = qX[x, y] - C * dX[x, y]
                            x1 = qX[x, y] + C * dX[x, y]
                            y0 = qY[x, y] - C * dY[x, y]
                            y1 = qY[x, y] + C * dY[x, y]
                            ax.plot([x0, x1], [y0, y1], color=lineColour, lw=MYLW, zorder=2)

                sc = ax.imshow(S.T, cmap=cmap, origin='lower',
                            vmin=0, vmax=S.max(),
                            zorder=1, extent=[0, xyz_size[dimX], 0, xyz_size[dimY]])
                # fig.canvas.draw()

            for p in range(poly_m):
                for k in range(mono_n):
                    circle = Circle((POLY[p, k, dimX], POLY[p, k, dimY]), radius=0.5,
                                    color=monomerColour, alpha=myAlpha, zorder=4)
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

            if cbFlag and fieldType != "none":
                cb=fig.colorbar(sc,shrink=shrink_factor,aspect=20*shrink_factor, pad=0.04)
                cb.ax.set_ylabel(label, fontsize=FS)
            ax.axis('off')
            plt.savefig(f"{data_path}/frame{j:06d}.png", format='png', bbox_inches='tight', transparent=makeTransparent)
            if savePDF:
                plt.savefig(f"{data_path}/frame{j:06d}.pdf", format='pdf', bbox_inches='tight', transparent=makeTransparent)
            plt.close()

        frame_id += 1
        j += 1
        i = 0
        if j > FINISH_FRAME: break

if fieldType != "none":
    field_file.close()
poly_file.close()
if show_defects:
    defect_file.close()

# === Animate ===
output_name = f"{data_path}/polymer"
if fieldType == "vel":
    output_name+="_flow"
elif fieldType == "nem":
    output_name+="_director"
elif fieldType == "none":
    pass
else:
    print("Field type not recognized.")
    error()
if defectData:
    output_name += "_defects"
output_name += f"_{projection}.mp4"
call(f"rm -f '{output_name}'", shell=True)
call(
    f"ffmpeg -f image2 -r {FRAMERATE} -i {data_path}/frame%06d.png "
    f"-vcodec {CODEC} -b {BITRATE}k -r {FRAMERATE} '{output_name}'",
    shell=True
)

if not keepFrames:
	call("rm frame*.png",shell=True)