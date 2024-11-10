'''Configs of all plots'''
### Name
EXPERIMENT_NAME = "lid-driven"

#### locations 
ROOT_LOCATION = "../../energy_results/"
MEAN_LOCATION = "/kinetic_energy/"
DET_LOCATION = "/deterministic/kinetic_energy/"
IND_LOCATION = "/individual/ind_kinetic_energy"
DATA_SOURCE = "refinement_9.csv"

#### output 
OUTPUT_LOCATION = f"output/{EXPERIMENT_NAME}/kinetic_energy/"

#### stochatic
NUMBER_SAMPLES = 100
NOISE_TYPES = ["lid-driven_exp1", "lid-driven_exp2", "lid-driven_exp3"]

#### matching experiment to its p-value
P_VALUE = {"lid-driven_exp1": 3/2.0, "lid-driven_exp2": 2, "lid-driven_exp3": 3}

#### stationary
STATIONARY_TIME = {"lid-driven_exp1": 0.4, "lid-driven_exp2": 0.2, "lid-driven_exp3": 0.1}
STATIONARY_ENERGY = {"lid-driven_exp1": 0.00864871866995075, "lid-driven_exp2": 0.013356505566619, "lid-driven_exp3": 0.0181430728419994}

#### plotting configs
# colours
COLOURS_MEAN = {"lid-driven_exp1": "#1b9e77", "lid-driven_exp2": "#d95f02", "lid-driven_exp3": "#7570b3"}
COLOURS_INDIVIDUAL  = {"lid-driven_exp1": "#66c2a5", "lid-driven_exp2": "#fc8d62", "lid-driven_exp3": "#8da0cb"}
BLACK = "#000000"

# histogram plot
HIST_DPI = 300
HIST_FILEFORMAT = "pdf"

YMAX =  {"lid-driven_exp1": 0.06, "lid-driven_exp2": 0.016, "lid-driven_exp3": 0.014}
LINEAR_PLOT = True
LOG_PLOT = False
HIST_XAXIS_LOG: bool = False

# trajectory plot
TRAJ_DPI = 300
TRAJ_FILEFORMAT = "png"

LINEWIDTH_MEAN = 2
LINEOPACITY_MEAN = 1

LINEWIDTH_SD = 1.5
LINEOPACITY_SD = 1
LINESTYLE_SD = (0, (1, 1))

LINEWIDTH_INDIVIDUAL = 0.1
LINEOPACITY_INDIVIDUAL = 1

LABEL_FONTSIZE = 20
TICK_FONTSIZE = 20
TRAJ_YAXIS_SCALE = "linear"

