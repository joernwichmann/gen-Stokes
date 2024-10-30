"""Contains global parameter configuration."""
################               GLOBAL configs               ############################
LOG_LEVEL: str = "info"  #supported levels: debug, info, warning, error, critical 

### Dump-location
DUMP_LOCATION: str = "sample_dump"

################               GENERATE configs               ############################
### Model
MODEL_NAME: str = "p-Stokes" #see src.algorithms.select.py for available choices
KAPPA_VALUE: float = 0.1

# Deterministic forcing
FORCING: str = "trigonometric"   #see 'src.predefined_data' for available choices
FORCING_FREQUENZY_X: int = 2
FORCING_FREQUENZY_Y: int = 4
FORCING_INTENSITY: float = 1

# Time
INITIAL_TIME: float = 0
END_TIME: float = 1
REFINEMENT_LEVELS: list[int] = list(range(9,10))

# Initial data
INITIAL_CONDITION_NAME: str = "polynomial - HL projected with BC"    #see 'src.predefined_data' for available choices
INITIAL_FREQUENZY_X: int = 2
INITIAL_FREQUENZY_Y: int = 4
INITIAL_INTENSITY: float = 1

# Elements
VELOCITY_ELEMENT: str = "CG"    #see firedrake doc for available spaces
VELOCITY_DEGREE: int = 2       

PRESSURE_ELEMENT: str = "CG"    #see firedrake doc for available spaces
PRESSURE_DEGREE: int = 1

# Mesh
NUMBER_SPACE_POINTS: int = 12
MESH_NAME: str = "unit square"  #see 'src.discretisation.mesh' for available choices
NAME_BOUNDARY_CONDITION: str = "zero"  #see 'src.discretisation.mesh' for available choices

# Monte Carlo
MC_SAMPLES: int = 10
NOISE_INCREMENTS: str = "classical" # see 'src.noise' for available choices

# Noise coefficient
NOISE_INTENSITY: float = 1
NOISE_COEFFICIENT_NAME: str = "polynomial" #see 'src.predefined_data' for available choices
NOISE_FREQUENZY_X: int = 2
NOISE_FREQUENZY_Y: int = 4

################               ANALYSE configs               ############################
#Convergence
TIME_CONVERGENCE: bool = False
TIME_COMPARISON_TYPE: str = "absolute"       ## "absolute" and "relative" are supported

#Stability
STABILITY_CHECK: bool = False

#Mean energy
ENERGY_CHECK: bool = True

#Individual energy
IND_ENERGY_CHECK: bool = True
IND_ENERGY_NUMBER: int = 1000

#Statistics
STATISTICS_CHECK: bool = False
