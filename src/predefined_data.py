from firedrake import *
from typing import TypeAlias, Callable, Any, Tuple

from src.discretisation.space import SpaceDiscretisation
from src.discretisation.projections import Stokes_projection, HL_projection, HL_projection_withBC


### Converter that maps string representation of functions to its implementation
def get_function(name_requested_function: str, space_disc: SpaceDiscretisation, 
                 index_x: int = 1, index_y: int = 1) -> Function | Tuple[Function,Function]:
    """
    Return a discrete velocity field. 

    Available fields: 'zero', 'x: hill, y: wave', 'non-solenoidal', 'solenoidal', "x: hill, y: wave - projected", "non-solenoidal - projected", "solenoidal - projected"

    Explanations with respect to unit square mesh:
        'x: hill, y: wave' -- first component is a hill, second component is a 2x2 checker board; the diagonal and off-diagonal contain hills and valleys, respectively
        'non-solenoidal' -- product of sin waves with frequency specified by index_x, index_y
        'solenoidal' -- product of sin and cos waves with frequency specified by index_x, index_y
    """
    match name_requested_function:
        ### unprojected functions
        ## returns velocity
        case "zero":
            return Function(space_disc.velocity_space)
        case "x: hill, y: wave":
            return _hill_wave(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "non-solenoidal":
            return _non_solenoidal(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "solenoidal":
            return _solenoidal(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "polynomial":
            return _polynomial(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "polynomial - no BC":
            return _polynomial_non_bc(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "polynomial - no div":
            return _polynomial_non_div(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "trigonometric":
            return _trig(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "trigonometric- no BC":
            return _trig_non_bc(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "trigonometric - no div":
            return _trig_non_div(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "gravity":
            return _gravity(velocity_space=space_disc.velocity_space)
        case "lid-driven-weak":
            return _lid_driven_weak(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        case "lid-driven-strong":
            return _lid_driven_strong(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space)
        
        ### Stokes projected functions
        ## returns velocity and pressure
        case "x: hill, y: wave - Stokes projected":
            return Stokes_projection(_hill_wave(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "non-solenoidal - Stokes projected":
            return Stokes_projection(_non_solenoidal(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc) 
        case "solenoidal - Stokes projected":
            return Stokes_projection(_solenoidal(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - Stokes projected":
            return Stokes_projection(_polynomial(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - no BC - Stokes projected":
            return Stokes_projection(_polynomial_non_bc(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - no div - Stokes projected":
            return Stokes_projection(_polynomial_non_div(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "zero - Stokes projected":
            return Stokes_projection(Function(space_disc.velocity_space),space_disc)
        case "trigonometric - Stokes projected":
            return Stokes_projection(_trig(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric- no BC - Stokes projected":
            return Stokes_projection(_trig_non_bc(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric - no div - Stokes projected":
            return Stokes_projection(_trig_non_div(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        
        ### HL projected functions
        ## returns velocity and pressure
        case "x: hill, y: wave - HL projected":
            return HL_projection(_hill_wave(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "non-solenoidal - HL projected":
            return HL_projection(_non_solenoidal(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc) 
        case "solenoidal - HL projected":
            return HL_projection(_solenoidal(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - HL projected":
            return HL_projection(_polynomial(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - no BC - HL projected":
            return HL_projection(_polynomial_non_bc(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - no div - HL projected":
            return HL_projection(_polynomial_non_div(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric - HL projected":
            return HL_projection(_trig(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric- no BC - HL projected":
            return HL_projection(_trig_non_bc(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric - no div - HL projected":
            return HL_projection(_trig_non_div(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        
        ### HL projected functions with BC
        ## returns velocity and pressure
        case "zero - HL projected with BC":
            return HL_projection_withBC(Function(space_disc.velocity_space),space_disc)
        case "polynomial - HL projected with BC":
            return HL_projection_withBC(_polynomial(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - no BC - HL projected with BC":
            return HL_projection_withBC(_polynomial_non_bc(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "polynomial - no div - HL projected with BC":
            return HL_projection_withBC(_polynomial_non_div(mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric - HL projected with BC":
            return HL_projection_withBC(_trig(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric- no BC - HL projected with BC":
            return HL_projection_withBC(_trig_non_bc(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        case "trigonometric - no div - HL projected with BC":
            return HL_projection_withBC(_trig_non_div(j=index_x,k=index_y,mesh=space_disc.mesh,velocity_space=space_disc.velocity_space),space_disc)
        
        ### others
        case other:
            raise NotImplementedError

#############################        Function generator   
### abstract concept
FunctionGenerator: TypeAlias = Callable[[Any],Function]

### implementation

#see Section 6 in 'Time-splitting Methods to solve the stochastic incrompressible Stokes equation'
def _non_solenoidal(j: int, k: int, mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        sin(j*pi*x)*sin(k*pi*y),
        -1*sin(j*pi*x)*sin(k*pi*y)
        ])
    return project(expr, velocity_space)

#see Section 6 in 'Time-splitting Methods to solve the stochastic incrompressible Stokes equation'
def _solenoidal(j: int, k: int, mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    #only solenoidal for j=k
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        -1.0*cos(j*pi*x - pi/2.0)*sin(k*pi*y - pi/2.0),
        sin(j*pi*x - pi/2.0)*cos(k*pi*y - pi/2.0)
        ])
    return project(expr, velocity_space)

def _trig(j: int, k: int, mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        sin(j*pi*x)*sin(k*pi*y),
        -1*sin(k*pi*x)*sin(j*pi*y)
        ])
    return project(expr, velocity_space)

def _trig_non_div(j: int, k: int, mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        sin(j*pi*x)*sin(k*pi*y),
        -1*sin(2*k*pi*x)*sin(2*j*pi*y)
        ])
    return project(expr, velocity_space)

def _trig_non_bc(j: int, k: int, mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        sin(j*pi*x)*cos(k*pi*y),
        -1*cos(k*pi*x)*sin(j*pi*y)
        ])
    return project(expr, velocity_space)

def _hill_wave(mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        sin(pi*x)*sin(pi*y),
        sin(2*pi*x)*sin(2*pi*y)
        ])
    return project(expr, velocity_space)

def _gravity(velocity_space: FunctionSpace) -> Function:
    expr = as_vector([0,-10])
    return project(expr, velocity_space)


#exactly divergence-free and vanishes on the boundary of the unit square
def _polynomial(mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    scaling: float = 1.0
    expr = as_vector([
        scaling*(x*x*(1-x)*(1-x)*(2-6*y+4*y*y)*y),
        scaling*(-y*y*(1-y)*(1-y)*(2-6*x+4*x*x)*x)
        ])
    return project(expr, velocity_space)

def _polynomial_non_bc(mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    scaling: float = 1.0
    expr = as_vector([
        scaling*(x*x*(1-x)*(1-x)*(2-6*y+4*y*y)*y + 1),
        scaling*(-y*y*(1-y)*(1-y)*(2-6*x+4*x*x)*x + 1)
        ])
    return project(expr, velocity_space)

def _polynomial_non_div(mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    scaling: float = 1.0
    expr = as_vector([
        scaling*(x*x*(1-x)*(1-x)*(2-6*y+4*y*y)*y + x*(1-x)*y*(1-y)),
        scaling*(-y*y*(1-y)*(1-y)*(2-6*x+4*x*x)*x + x*(1-x)*y*(1-y))
        ])
    return project(expr, velocity_space)

def _lid_driven_weak(mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        (4*y*x*(1-x)),
        (0)
        ])
    return project(expr, velocity_space)

def _lid_driven_strong(mesh: MeshGeometry, velocity_space: FunctionSpace) -> Function:
    x, y = SpatialCoordinate(mesh)
    expr = as_vector([
        (1*y**100),
        (0)
        ])
    return project(expr, velocity_space)
        
