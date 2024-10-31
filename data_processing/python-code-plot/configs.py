'''Configs of all plots'''

#### locations 
ROOT_LOCATION = "../../energy_results/"
MEAN_LOCATION = "/kinetic_energy/"
IND_LOCATION = "/individual/ind_kinetic_energy"
DATA_SOURCE = "refinement_9.csv"

#### stochatic
NUMBER_SAMPLES = 10
NOISE_TYPES = ["p-variation_exp1", "p-variation_exp2", "p-variation_exp3"]

#### stationary
STATIONARY_TIME = {"p-variation_exp1": 0.13, "p-variation_exp2": 0.32, "p-variation_exp3": 0.13}
STATIONARY_ENERGY = 0.04206492724892853

#### plotting configs
# colours
COLOURS_MEAN = {"p-variation_exp1": "#1b9e77", "p-variation_exp2": "#d95f02", "p-variation_exp3": "#7570b3"}
COLOURS_INDIVIDUAL  = {"p-variation_exp1": "#66c2a5", "p-variation_exp2": "#fc8d62", "p-variation_exp3": "#8da0cb"}
BLACK = "#000000"

# histogram plot
HIST_DPI = 300
HIST_FILEFORMAT = "pdf"

YMAX =  {"p-variation_exp1": 0.06, "p-variation_exp2": 0.016, "p-variation_exp3": 0.014}
LINEAR_PLOT = True
LOG_PLOT = False

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

