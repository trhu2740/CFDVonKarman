from mpi4py import MPI
from petsc4py import PETSc
import dolfinx as dfx
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.io import gmshio
import ufl as ufl
import numpy as np
import time
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------------
# update and implement the integrateFluidStress function
# -------------------------------------------------------------------------------
def integrateFuidStress(a_U, a_P, a_Mu, a_N, a_Mesh, a_GammaP):
    eps = 0.5 * (ufl.grad(a_U) + ufl.grad(a_U).T)
    sig = -a_P * ufl.Identity(2) + 2.0 * a_Mu * eps

    traction = ufl.dot(sig, a_N)

    forceX = traction[0] * a_GammaP
    forceY = traction[1] * a_GammaP

    fXVal = dfx.fem.assemble_scalar(dfx.fem.form(forceX))
    fYVal = dfx.fem.assemble_scalar(dfx.fem.form(forceY))

    return [fXVal, fYVal]


# -------------------------------------------------------------------------------
# Meshing: make sure to change the file names below
# -------------------------------------------------------------------------------
meshFile = "TwoWayVonKarman.msh"
outFileV = "NSE-V.xdmf"
outFileP = "NSE-P.xdmf"
forceFile = "forces.dat"

viscosity = 0.002
density = 1.0
U0 = 2.5
diam = 0.2
boxH = 4.1 * diam

dt = 0.005
t_start = 0.0
t_end = 7.0
t_theta = 0.5

Reynolds = density * (2.0 / 3.0) * U0 * diam / viscosity
print("The Problem Reynolds Number is:", Reynolds)


def noSlipBC(x):
    return np.stack((np.zeros(x.shape[1]), np.zeros(x.shape[1])))


def pressureBC(x):
    return np.zeros(x.shape[1])


def inletBC_LEFT(x):
    vals = np.zeros((mesh.geometry.dim, x.shape[1]))
    vals[0] = 6.0 * U0 * x[1] * (boxH - x[1]) / (boxH * boxH)
    vals[1] = 0.0
    return vals


def inletBC_RIGHT(x):
    vals = np.zeros((mesh.geometry.dim, x.shape[1]))
    vals[0] = -6.0 * U0 * x[1] * (boxH - x[1]) / (boxH * boxH)
    vals[1] = 0.0
    return vals


# -------------------------------------------------------------------------------
# Meshing: make sure to update/modify the ids as below for your mesh/geometry
# -------------------------------------------------------------------------------
ID_INLET_LEFT = 30
ID_TOP_LEFT = 31
ID_FUNNEL_TOP_LEFT = 32
ID_OUTLET_TOP = 33
ID_FUNNEL_TOP_RIGHT = 34
ID_TOP_RIGHT = 35
ID_INLET_RIGHT = 36
ID_BOTTOM_RIGHT = 37
ID_FUNNEL_BOTTOM_RIGHT = 38
ID_OUTLET_BOTTOM = 39
ID_FUNNEL_BOTTOM_LEFT = 40
ID_BOTTOM_LEFT = 41
ID_CYL_LEFT = 42
ID_FIN_LEFT = 43
ID_CYL_RIGHT = 46
ID_FIN_RIGHT = 44


startTime = time.time()

mesh, cell_markers, facet_markers = gmshio.read_from_msh(
    meshFile, MPI.COMM_WORLD, gdim=2
)

nVec = ufl.FacetNormal(mesh)

tdim = mesh.topology.dim
fdim = tdim - 1

V = ufl.VectorElement("Lagrange", mesh.ufl_cell(), 1)
P = ufl.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
M = ufl.MixedElement([V, P])
W = dfx.fem.FunctionSpace(mesh, M)

V_sub, V_submap = W.sub(0).collapse()
P_sub, P_submap = W.sub(1).collapse()

print("here")

# -------------------------------------------------------------------------------
# Meshing: make sure to update/modify the boundary dofs as below, for your
# mesh/geometry in case if needed
# -------------------------------------------------------------------------------
b_dofs_INLET_LEFT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_INLET_LEFT)
)
b_dofs_INLET_RIGHT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_INLET_RIGHT)
)
b_dofs_TOP_LEFT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_TOP_LEFT)
)
b_dofs_TOP_RIGHT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_TOP_RIGHT)
)
b_dofs_FUNNEL_TOP_LEFT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_FUNNEL_TOP_LEFT)
)
b_dofs_FUNNEL_TOP_RIGHT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_FUNNEL_TOP_RIGHT)
)
b_dofs_FUNNEL_BOTTOM_RIGHT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_FUNNEL_BOTTOM_RIGHT)
)
b_dofs_FUNNEL_BOTTOM_LEFT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_FUNNEL_BOTTOM_LEFT)
)
b_dofs_OUTLET_TOP = dfx.fem.locate_dofs_topological(
    (W.sub(1), P_sub), fdim, facet_markers.find(ID_OUTLET_TOP)
)
b_dofs_OUTLET_BOTTOM = dfx.fem.locate_dofs_topological(
    (W.sub(1), P_sub), fdim, facet_markers.find(ID_OUTLET_BOTTOM)
)
b_dofs_BOTTOM_RIGHT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_BOTTOM_RIGHT)
)
b_dofs_BOTTOM_LEFT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_BOTTOM_LEFT)
)
b_dofs_CYL_LEFT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_CYL_LEFT)
)
b_dofs_CYL_RIGHT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_CYL_RIGHT)
)
b_dofs_FIN_LEFT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_FIN_LEFT)
)
b_dofs_FIN_RIGHT = dfx.fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_markers.find(ID_FIN_RIGHT)
)

uD_Wall = dfx.fem.Function(V_sub)
uD_Wall.interpolate(noSlipBC)

uD_Inlet_LEFT = dfx.fem.Function(V_sub)
uD_Inlet_LEFT.interpolate(inletBC_LEFT)

uD_Inlet_RIGHT = dfx.fem.Function(V_sub)
uD_Inlet_RIGHT.interpolate(inletBC_RIGHT)

uD_Outlet = dfx.fem.Function(P_sub)
uD_Outlet.interpolate(pressureBC)

bc_INLET_LEFT = dfx.fem.dirichletbc(uD_Inlet_LEFT, b_dofs_INLET_LEFT, W.sub(0))
bc_INLET_RIGHT = dfx.fem.dirichletbc(uD_Inlet_RIGHT, b_dofs_INLET_RIGHT, W.sub(0))
bc_TOP_LEFT = dfx.fem.dirichletbc(uD_Wall, b_dofs_TOP_LEFT, W.sub(0))
bc_TOP_RIGHT = dfx.fem.dirichletbc(uD_Wall, b_dofs_TOP_RIGHT, W.sub(0))
bc_FUNNEL_TOP_RIGHT = dfx.fem.dirichletbc(uD_Wall, b_dofs_FUNNEL_TOP_RIGHT, W.sub(0))
bc_FUNNEL_TOP_LEFT = dfx.fem.dirichletbc(uD_Wall, b_dofs_FUNNEL_TOP_LEFT, W.sub(0))

bc_FUNNEL_BOTTOM_RIGHT = dfx.fem.dirichletbc(
    uD_Wall, b_dofs_FUNNEL_BOTTOM_RIGHT, W.sub(0)
)
bc_FUNNEL_BOTTOM_LEFT = dfx.fem.dirichletbc(
    uD_Wall, b_dofs_FUNNEL_BOTTOM_LEFT, W.sub(0)
)
bc_OUTLET_TOP = dfx.fem.dirichletbc(uD_Outlet, b_dofs_OUTLET_TOP, W.sub(1))
bc_OUTLET_BOTTOM = dfx.fem.dirichletbc(uD_Outlet, b_dofs_OUTLET_BOTTOM, W.sub(1))
bc_BOTTOM_RIGHT = dfx.fem.dirichletbc(uD_Wall, b_dofs_BOTTOM_RIGHT, W.sub(0))
bc_BOTTOM_LEFT = dfx.fem.dirichletbc(uD_Wall, b_dofs_BOTTOM_LEFT, W.sub(0))
bc_CYL_LEFT = dfx.fem.dirichletbc(uD_Wall, b_dofs_CYL_LEFT, W.sub(0))
bc_CYL_RIGHT = dfx.fem.dirichletbc(uD_Wall, b_dofs_CYL_RIGHT, W.sub(0))
bc_FIN_LEFT = dfx.fem.dirichletbc(uD_Wall, b_dofs_FIN_LEFT, W.sub(0))
bc_FIN_RIGHT = dfx.fem.dirichletbc(uD_Wall, b_dofs_FIN_RIGHT, W.sub(0))

bc = [
    bc_INLET_LEFT,
    bc_INLET_RIGHT,
    bc_TOP_LEFT,
    bc_TOP_RIGHT,
    bc_FUNNEL_TOP_RIGHT,
    bc_FUNNEL_TOP_LEFT,
    bc_FUNNEL_BOTTOM_RIGHT,
    bc_FUNNEL_BOTTOM_LEFT,
    bc_OUTLET_TOP,
    bc_OUTLET_BOTTOM,
    bc_BOTTOM_RIGHT,
    bc_BOTTOM_LEFT,
    bc_CYL_LEFT,
    bc_CYL_RIGHT,
    bc_FIN_LEFT,
    bc_FIN_RIGHT,
]

ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_markers)

# TODO: Get integrateFuidStress for both left and right cylinder#
Gamma_CYL = ds(ID_CYL_LEFT)

(v, q) = ufl.TestFunctions(W)

mu = dfx.fem.Constant(mesh, dfx.default_scalar_type(viscosity))
rho = dfx.fem.Constant(mesh, dfx.default_scalar_type(density))
idt = dfx.fem.Constant(mesh, dfx.default_scalar_type(1.0 / dt))
theta = dfx.fem.Constant(mesh, dfx.default_scalar_type(t_theta))
b = dfx.fem.Constant(mesh, PETSc.ScalarType((0.0, 0.0)))
I = ufl.Identity(2)

W1 = dfx.fem.Function(W)
(u, p) = ufl.split(W1)

T1_1 = rho * ufl.inner(v, ufl.grad(u) * u) * ufl.dx
T2_1 = mu * ufl.inner(ufl.grad(v), ufl.grad(u)) * ufl.dx
T3_1 = p * ufl.div(v) * ufl.dx
T4_1 = q * ufl.div(u) * ufl.dx
T5_1 = rho * ufl.dot(v, b) * ufl.dx
L_1 = T1_1 + T2_1 - T3_1 - T4_1 - T5_1

W0 = dfx.fem.Function(W)
(u0, p0) = ufl.split(W0)

T1_0 = rho * ufl.inner(v, ufl.grad(u0) * u0) * ufl.dx
T2_0 = mu * ufl.inner(ufl.grad(v), ufl.grad(u0)) * ufl.dx
T3_0 = p * ufl.div(v) * ufl.dx
T4_0 = q * ufl.div(u0) * ufl.dx

T5_0 = rho * ufl.dot(v, b) * ufl.dx
L_0 = T1_0 + T2_0 - T3_0 - T4_0 - T5_0

F = idt * ufl.inner((u - u0), v) * ufl.dx + (1.0 - theta) * L_0 + theta * L_1

uNorm = ufl.sqrt(ufl.inner(u0, u0))
h = ufl.CellDiameter(mesh)
tau = (
    (2.0 * theta * idt) ** 2 + (2.0 * uNorm / h) ** 2 + (4.0 * mu / h**2) ** 2
) ** (-0.5)

residual = (
    idt * rho * (u - u0)
    + theta
    * (rho * ufl.grad(u) * u - mu * ufl.div(ufl.grad(u)) + ufl.grad(p) - rho * b)
    + (1.0 - theta)
    * (rho * ufl.grad(u0) * u0 - mu * ufl.div(ufl.grad(u0)) + ufl.grad(p) - rho * b)
)

F_SUPG = tau * ufl.inner(ufl.grad(v) * u, residual) * ufl.dx
F_PSPG = -tau * ufl.inner(ufl.grad(q), residual) * ufl.dx

F = F + F_SUPG + F_PSPG

problem = NonlinearProblem(F, W1, bcs=bc)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-7
solver.report = True

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "gmres"
opts[f"{option_prefix}pc_type"] = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "umfpack"
ksp.setFromOptions()

t = t_start
tn = 0

if outFileV.endswith("pvd"):
    vFile = dfx.io.VTKFile(MPI.COMM_WORLD, outFileV, "w")
elif outFileV.endswith("xdmf"):
    vFile = dfx.io.XDMFFile(
        mesh.comm, outFileV, "w", encoding=dfx.io.XDMFFile.Encoding.ASCII
    )
    vFile.write_mesh(mesh)

if outFileP.endswith("pvd"):
    pFile = dfx.io.VTKFile(MPI.COMM_WORLD, outFileP, "w")
elif outFileP.endswith("xdmf"):
    pFile = dfx.io.XDMFFile(
        mesh.comm, outFileP, "w", encoding=dfx.io.XDMFFile.Encoding.ASCII
    )
    pFile.write_mesh(mesh)

time_Arr = []
fX_Arr = []
fY_Arr = []
cD_Arr = []
cL_Arr = []

while t < t_end:

    t_in = time.time()
    n, converged = solver.solve(W1)
    assert converged
    t_out = time.time()

    print(f"t = {t:.6f}; Number of iterations: {n:d}; compute time: {t_out-t_in:f}")

    uf = W1.split()[0].collapse()
    pf = W1.split()[1].collapse()

    uf.name = "vel"
    pf.name = "pres"

    if outFileV.endswith("pvd"):
        vFile.write_function(uf, tn)
    elif outFileV.endswith("xdmf"):
        vFile.write_function(uf, t)

    if outFileP.endswith("pvd"):
        pFile.write_function(pf, tn)
    elif outFileP.endswith("xdmf"):
        pFile.write_function(pf, t)

    # ---------------------------------------------------------------------------
    # FOR ASSIGNMENT:
    # implement the call to the stress integration function, and the
    # calculation of drag and lift coefficients
    # ---------------------------------------------------------------------------
    [fX, fY] = integrateFuidStress(uf, pf, mu, nVec, mesh, Gamma_CYL)
    cD = fX * 2.0 / ((2.0 * U0 / 3.0) ** 2 * diam)
    cL = fY * 2.0 / ((2.0 * U0 / 3.0) ** 2 * diam)
    time_Arr.append(t)
    fX_Arr.append(fX)
    fY_Arr.append(fY)
    cD_Arr.append(cD)
    cL_Arr.append(cL)

    W0.x.array[:] = W1.x.array

    t += dt
    tn += 1

vFile.close()
pFile.close()

forceData = np.column_stack([time_Arr, fX_Arr, fY_Arr, cD_Arr, cL_Arr])
np.savetxt(forceFile, forceData)

endTime = time.time()

print("Total simulation time:", endTime - startTime)

fig, ax = plt.subplots(1, 2)
ax[0].plot(time_Arr[-200:], cD_Arr[-200:], "r")
ax[1].plot(time_Arr[-200:], cL_Arr[-200:], "m")
ax[0].set_xlabel("time", fontweight="bold")
ax[1].set_xlabel("time", fontweight="bold")
ax[0].set_ylabel("streamwise force", fontweight="bold")
ax[1].set_ylabel("cross-stream force", fontweight="bold")
plt.show()
