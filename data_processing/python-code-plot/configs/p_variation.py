'''Configs of all plots'''
### Name
EXPERIMENT_NAME = "p-variation"
EXPERIMENTS = {1: f"{EXPERIMENT_NAME}_exp1", 2: f"{EXPERIMENT_NAME}_exp2", 3: f"{EXPERIMENT_NAME}_exp3" }

#### matching experiment to its p-value
P_VALUE = {1: 3/2.0, 2: 2, 3: 3}

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

#### stationary
STATIONARY_TIME = {1: 0.5, 2: 0.5, 3: 0.5}

#### plotting configs
# colours
COLOURS_MEAN = {1: "#1b9e77", 2: "#d95f02", 3: "#7570b3"}
COLOURS_INDIVIDUAL  = {1: "#66c2a5", 2: "#fc8d62", 3: "#8da0cb"}
BLACK = "#000000"

# histogram plot
HIST_DPI = 300
HIST_FILEFORMAT = "pdf"

#ymax should be dynamic
YMAX =  {1: 0.06, 2: 0.016, 3: 0.014}
LINEAR_PLOT = True
LOG_PLOT = False
HIST_XAXIS_LOG: bool = True

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

LABEL_FONTSIZE = 10
TICK_FONTSIZE = 10
TRAJ_YAXIS_SCALE = "log"

LINESTYLES_DET = {1: "dotted", 2: "dashed", 3: "dashdot"}

