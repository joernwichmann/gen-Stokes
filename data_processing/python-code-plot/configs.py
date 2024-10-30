'''Configs of all plots'''

#### locations 
ROOT_LOCATION = "../data/energy_results/"
MEAN_LOCATION = "/kinetic_energy/"
IND_LOCATION = "/individual/ind_kinetic_energy"
DATA_SOURCE = "refinement_9.csv"

#### stochatic
NUMBER_SAMPLES = 1000
NOISE_TYPES = ["intro-additive", "intro-multiplicative", "intro-transport"]

#### stationary
STATIONARY_TIME = {"intro-additive": 0.13, "intro-multiplicative": 0.32, "intro-transport": 0.13}
STATIONARY_ENERGY = 0.04206492724892853

#### plotting configs
# colours
COLOURS_MEAN = {"intro-additive": "#1b9e77", "intro-multiplicative": "#d95f02", "intro-transport": "#7570b3"}
COLOURS_INDIVIDUAL  = {"intro-additive": "#66c2a5", "intro-multiplicative": "#fc8d62", "intro-transport": "#8da0cb"}
BLACK = "#000000"

# histogram plot
HIST_DPI = 300
HIST_FILEFORMAT = "pdf"

YMAX =  {"intro-additive": 0.06, "intro-multiplicative": 0.016, "intro-transport": 0.014}
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

