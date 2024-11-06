'''Configs of all plots'''

EXPERIMENT_NAME = "p-variation"
EXPERIMENTS = {1: f"{EXPERIMENT_NAME}_exp1", 2: f"{EXPERIMENT_NAME}_exp2", 3: f"{EXPERIMENT_NAME}_exp3" }
#### matching experiment to its p-value
P_VALUE = {1: 3/2.0, 2: 2, 3: 3}

#### locations 
ROOT_LOCATION = "../../point_statistic_results/"
DET_LOCATION = "/deterministic/p1/mean"
MEAN_LOCATION = "/p1/mean"
IND_LOCATION = "/p1/seed"
DATA_SOURCE = "/refinement_10.csv"

#### stochatic
NUMBER_SAMPLES = 100

#### stationary
STATIONARY_TIME = {1: 0.4, 2: 0.4, 3: 0.4}
STATIONARY_VAL_X = {1:0.044043617290565246,
                    2:-0.014916417092931724,
                    3:-0.010145255394873131}
STATIONARY_VAL_Y = {1:3.7317082618283726,
                    2:0.5703828845767361,
                    3:0.23141628595966934}

#### plotting configs
# colours
COLOURS_MEAN = {1: "#1b9e77", 2: "#d95f02", 3: "#7570b3"}
COLOURS_INDIVIDUAL  = {1: "#66c2a5", 2: "#fc8d62", 3: "#8da0cb"}
BLACK = "#000000"

# histogram plot
HIST_DPI = 300
HIST_FILEFORMAT = "pdf"

YMAX =  {1: 0.06, 2: 0.016, 3: 0.014}
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

