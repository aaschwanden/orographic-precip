"""
Nonlinear 1-D transport equation with Dirichlet conditions
at the surface and an inhomogeneous Neumann (flux) condition
at the base. The domain is the interval from a to b.

u = u_0 at x=a, u = u_1 at x=b.
dudt-div(q(u)*nabla_grad(u) + velocity*nabla_grad(u)) = f,

Solution method: automatic, i.e., by a NonlinearVariationalProblem/Solver
(Newton method).
"""

from dolfin import *
import sys
import numpy as np
import pylab as plt
from argparse import ArgumentParser

set_log_level(ERROR)
tol = 1E-14

class DirichletBCTransientNonlinearSolver(object):

    def __init__(self, kappa, velocity, f, g, bcs, *args, **kwargs):
        super(DirichletBCTransientNonlinearSolver, self).__init__()

        # get time control parameters
        time_control = kwargs['time_control']

        t_a = time_control['t_a']
        t_e = time_control['t_e']
        dt = time_control['dt']
        theta = time_control['theta']

        print('--------------------------------------------------------')
        print('Running Transient Nonlinear Solver with Dirichlet BCs')
        print('--------------------------------------------------------\n')

        t = t_a
        print('time: {} (start)'.format(t))

        u_sol = []

        # get initial condition
        u_init = kwargs['u_init']
        u_init.t = t

        u_surf, u_base = bcs

        u_ = Function(V)  # most recently computed solution

        
        # We use the exact solution of the linear diffusion problem
        # as an initial condition at t=t_a
        u0 = interpolate(u_init, V)

        u_prev = u0
        u_sol.append(u0.vector().array())
        t += dt

        while t <= t_e:
            print('time: {}'.format(t))
            

            # u_(n+theta)
            u_mid = (1.0-theta)*u_prev + theta*E

            # necessary quantities for streamline upwinding :
            h      = 2 * CellSize(mesh)

            # skewed test function :
            psihat = psi + h/2 * sign(velocity) * psi.dx(0)

            F = ((u-u_prev)*psihat*dx + dt*(inner(kappa*nabla_grad(u_mid), nabla_grad(psi))*dx 
                + velocity*E.dx(0)*psihat*dx 
                + f*psihat*dx))

            F  = action(F, u_)
            J = derivative(F, u_, u)

            # Compute solution
            problem = NonlinearVariationalProblem(F, u_, bcs, J)
            solver  = NonlinearVariationalSolver(problem)
            prm = solver.parameters
            prm['newton_solver']['absolute_tolerance'] = 1E-8
            prm['newton_solver']['relative_tolerance'] = 1E-7
            prm['newton_solver']['maximum_iterations'] = 25
            prm['newton_solver']['relaxation_parameter'] = 1.0
            solver.solve()

            t += dt

            u_prev.assign(u_)
            u_sol.append(u_.vector().array())
        print('Time stepping done')

        # Fixme: What should we return??
        self.u_ = u_
        self.u_sol = u_sol


# Set up geometry and mesh
a = 30000
nx = 10000
mesh = IntervalMesh(nx, -a, a)
ele_order = 1
# Define function space
V = FunctionSpace(mesh, 'Lagrange', ele_order)


# Define boundary conditions

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0]+a) < tol

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0]-a) < tol

boundary_parts = FacetFunction("size_t", mesh, 1)
boundary_parts.set_all(0)

Gamma_d = SurfaceBoundary()
Gamma_d.mark(boundary_parts, 1)

Gamma_g = LowerBoundary()
Gamma_g.mark(boundary_parts, 2)

# Define variational problem
psi  = TestFunction(V)
u  = TrialFunction(V)
u_mid = E
f = Constant(0.)
ds = ds[boundary_parts]

# Update surface boundary condition
bcs = [DirichletBC(V, Constant(u_0), boundary_parts, 1), DirichletBC(V, Constant(u_1), boundary_parts, 2)]
