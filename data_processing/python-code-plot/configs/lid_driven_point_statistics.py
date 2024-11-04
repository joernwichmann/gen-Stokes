'''Configs of all plots'''

EXPERIMENT_NAME = "lid-driven"
EXPERIMENTS = {1: f"{EXPERIMENT_NAME}_exp1", 2: f"{EXPERIMENT_NAME}_exp2", 3: f"{EXPERIMENT_NAME}_exp3" }
#### matching experiment to its p-value
P_VALUE = {1: 3/2.0, 2: 2, 3: 3}

#### locations 
ROOT_LOCATION = "../../point_statistic_results/"
DET_LOCATION = "/deterministic/p1/mean"
MEAN_LOCATION = "/p1/mean"
IND_LOCATION = "/p1/seed"
DATA_SOURCE = "/refinement_9.csv"

#### stochatic
NUMBER_SAMPLES = 10

#### stationary
#STATIONARY_TIME = {1: 0.13, 2: 0.32, 3: 0.13}
STATIONARY_TIME = {1: 0.4, 2: 0.4, 3: 0.4}
STATIONARY_VAL_X = {1:0.027870811123365345,
                    2:0.0018390038770400882,
                    3:0.05266201306509131}
STATIONARY_VAL_Y = {1:0.0018067095598706987,
                    2:0.001796688611121015,
                    3:0.0027004089468416533}

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

