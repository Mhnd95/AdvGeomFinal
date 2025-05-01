#boundary Import necessary libraries
using Gridap                # Main Gridap package for finite element analysis
using Gridap.Geometry       # For mesh and geometry handling
using Gridap.FESpaces       # For finite element spaces
using Gridap.MultiField     # For coupled multi-physics problems
using Gridap.Io             # For input/output operations
using Gridap.Fields         # For field operations
using Gridap.TensorValues   # For tensor operations
using Gridap.ODEs           # For time-dependent problems
using Gridap.CellData       # For cell data operations and projection
using WriteVTK              # For VTK file output (visualization)
using GridapGmsh            # For Gmsh mesh integration

# ============================================================================
# PROBLEM DESCRIPTION
# ============================================================================
# This is a plane strain poroelasticity formulation
# In plane strain, we assume εzz = 0 (no strain in z-direction)
# but σzz ≠ 0 (stress in z-direction can exist)
# This is appropriate for modeling soil/rock layers where the z-dimension is 
# constrained but stresses can develop in that direction

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

# Material properties
E = 20.0e6         # Young's modulus (Pa)
nu = 0.2          # Poisson's ratio
B = 0.8           # Biot coefficient (coupling between fluid pressure and solid stress)
M = 1.0e9         # Biot modulus (Pa) - related to fluid and solid compressibility
k = 1.0e-3        # Permeability (m^2) - how easily fluid flows through the medium
mu = 1.0e-3       # Fluid viscosity (Pa·s)

# Loading conditions
Pb = 31.5e6
p0 = 20.0e6

# Time stepping parameters
T = 0.0005          # Final time (s)
num_steps = 3   # Number of time steps
dt = T / num_steps # Time step size (s)

# ============================================================================
# DERIVED MATERIAL PROPERTIES
# ============================================================================
# Calculate Lamé parameters for plane strain formulation
# For plane strain, we use the same Lamé parameters as in 3D
lambda = E * nu / ((1 + nu) * (1 - 2 * nu))  # First Lamé parameter (Pa)
mu = E / (2 * (1 + nu))                      # Second Lamé parameter (shear modulus) (Pa)
k_mu = k / mu                                # Hydraulic conductivity (permeability/viscosity)

dirichlet_tags = ["top_bottom", "wellbore"]

# ============================================================================
# SETUP OUTPUT AND MESH
# ============================================================================
# Create output directory if it doesn't exist
output_dir = "results"
if !isdir(output_dir)
    mkdir(output_dir)
end

# Load the Gmsh mesh from file
# The mesh should be a square domain with properly tagged boundaries
model = GmshDiscreteModel("wellbore.msh")

# Define boundary tags for applying boundary conditions
# These tags should match the physical groups defined in the Gmsh file
dirichlet_tags = ["top_bottom", "wellbore"]

# Export the mesh for visualization
writevtk(model, "model")  # Save model for visualization in ParaView

# Print information about boundary entities for debugging
# This helps verify that boundary conditions will be applied to the correct entities
labels = get_face_labeling(model)
for tag in dirichlet_tags
    println("Entities tagged as $tag: ", findall(labels.tag_to_name .== tag))
end

# ============================================================================
# DOMAIN AND INTEGRATION SETUP
# ============================================================================
# Set up integration degree for numerical quadrature
degree = 2  # Quadrature order

# Create triangulation and integration measures
Ω = Triangulation(model)           # Domain triangulation
dΩ = Measure(Ω, degree)            # Volume integration measure
Γ = BoundaryTriangulation(model)   # Boundary triangulation
dΓ = Measure(Γ, degree)            # Boundary integration measure

# ============================================================================
# FINITE ELEMENT SPACES
# ============================================================================
# Define polynomial orders for the mixed formulation
order_u = 2  # P2 (quadratic) elements for displacement
order_p = 1  # P1 (linear) elements for pressure - satisfies LBB condition

# Create reference finite elements
reffe_u = ReferenceFE(lagrangian, VectorValue{2,Float64}, order_u)  # Vector-valued for displacement
reffe_p = ReferenceFE(lagrangian, Float64, order_p)                 # Scalar-valued for pressure

# ============================================================================
# BOUNDARY CONDITIONS
# ============================================================================
# Displacement space with Dirichlet BC on top and bottom (fixed in y-direction)
δu = TestFESpace(model, reffe_u, conformity=:H1, 
                 dirichlet_tags=["top_bottom"])
u = TrialFESpace(δu, x -> VectorValue(0.0, 0.0))  # Zero displacement at bottom boundary

# Pressure space with Dirichlet BC on left and right sides (drained boundaries)
δp = TestFESpace(model, reffe_p, conformity=:H1, dirichlet_tags=["wellbore"])
p = TrialFESpace(δp, Pb)  # Boundary pressure at wellbore

# Create multi-field space for the coupled problem
Y = MultiFieldFESpace([δu, δp])  # Combined test space for displacement and pressure