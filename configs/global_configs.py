"""Contains global parameter configuration."""
################               GLOBAL configs               ############################
LOG_LEVEL: str = "info"  #supported levels: debug, info, warning, error, critical 

### Dump-location
DUMP_LOCATION: str = "sample_dump"

################               GENERATE configs               ############################
### Model
MODEL_NAME: str = "Navier--Stokes" #see src.algorithms.select.py for available choices
P_VALUE: float = 2

### Discretisation
# Time
INITIAL_TIME: float = 0
END_TIME: float = 1
REFINEMENT_LEVELS: list[int] = list(range(2,10))
INITIAL_INTENSITY: float = 1

# Space
NUMBER_SPACE_POINTS: int = 12
MESH_NAME: str = "unit square"  #see 'src.discretisation.mesh' for available choices
NAME_BOUNDARY_CONDITION: str = "zero"  #see 'src.discretisation.mesh' for available choices

# Stochastic
MC_SAMPLES: int = 1000
NOISE_INCREMENTS: str = "classical" # see 'src.noise' for available choices
NOISE_INTENSITY: float = 1000

################               ANALYSE configs               ############################
#Convergence
TIME_CONVERGENCE: bool = True
TIME_COMPARISON_TYPE: str = "absolute"       ## "absolute" and "relative" are supported

#Stability
STABILITY_CHECK: bool = True

#Energy
ENERGY_CHECK: bool = True

#Statistics
STATISTICS_CHECK: bool = True
