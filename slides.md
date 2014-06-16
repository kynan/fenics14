class: center, middle, inverse

# Firedrake: Re-imagining FEniCS by Composing Domain-specific Abstractions

## **Florian Rathgeber**<sup>1</sup>, Lawrence Mitchell<sup>1</sup>, David Ham<sup>1,2</sup>, Michael Lange<sup>3</sup>, Andrew McRae<sup>2</sup>, Fabio Luporini<sup>1</sup>, Gheorghe-teodor Bercea<sup>1</sup>, Paul Kelly<sup>1</sup>

.footnote[<sup>1</sup> Department of Computing, Imperial College London
<sup>2</sup> Department of Mathematics, Imperial College London  
<sup>3</sup> Department of Earth Science & Engineering, Imperial College London]

???

Good morning! My name is Florian and it is my pleasure to introduce Firedrake, a
re-imagination of the FEniCS concept by composing domain-specific abstraction.
Firedrake is developed by a group at Imperial, some which are in the room, so
I'll be presenting work of a number of people.

---

.scale[![FEniCS](http://fenicsproject.org/_static/fenics_banner.png)]

> The FEniCS Project is a collection of free software for automated, efficient
> solution of differential equations.
>
> .source[&mdash; fenicsproject.org]

???

When I say "re-imagining FEniCS", let's recap how the FEniCS project defines
itself on the website: ...

---

.scale[![Firedrake](http://firedrakeproject.org/_static/banner.png)]

> Firedrake is an automated system for the portable solution of partial
> differential equations using the finite element method (FEM).
>
> .source[&mdash; firedrakeproject.org]

???

Firedrake's mission statement is very similar: ...

--

Two-layer abstraction for FEM computation from high-level descriptions:
* Firedrake: a portable finite-element computation framework  
  *Drive FE computations from a high-level problem specification*
* PyOP2: a high-level interface to unstructured mesh based methods  
  *Efficiently execute kernels over an unstructured grid in parallel*

???

The key concept is a layering of abstractions, with Firedrake as a portable
finite-element computation framework on top, allowing to drive computations from
a very high-level problem specification. The lower layer is PyOP2, a framework
for unstructured mesh computations, whose role is to efficiently execute kernels
over an unstructured mesh in parallel.

---

## The Firedrake/PyOP2 tool chain

![Firedrake](images/firedrake_toolchain.svg)

???

This diagram shows an overview of the Firedrake / PyOP2 toolchain:

* Decoupling of Firedrake (FEM) and PyOP2 (parallelisation) layers
* (mostly) DOLFIN compatible API
* UFL to describe finite-element discretisation
* PETSc for linear/nonlinear solvers and mesh building / distribution
* Platform-specific runtime code generation and JIT compilation
* Portability for unstructured mesh applications: FEM, non-FEM or combinations
* Extensible framework beyond FEM computations (e.g. image processing)
* Revisit this diagram later

---

class: center, middle
# Parallel computations on unstructured meshes with PyOP2

???

Let's start at the bottom layer and look at parallel computations on
unstructured meshes with PyOP2.

---

.pull-left[
## Scientific computations on unstructured meshes

* Independent *local operations* for each element of the mesh described by a *kernel*.
* *Reductions* aggregate contributions from local operations to produce the final result.

### PyOP2

A domain-specific language embedded in Python for parallel computations on unstructured meshes or graphs.

### Unstructured mesh
.scale[![PyOP2 mesh](images/op2_mesh.svg)]
]

???

PyOP2 was born from the realisation that scientific computations on unstructured
meshes often share the same structure: there is an independent local operation
that needs to be performed for every entity of the mesh and can be described by
a computational kernel. The final result is obtained from a reduction which
aggregates contributions from these local operations.

PyOP2 is a domain-specific language embedded in Python which implements this
principle for parallel computations on unstructured meshes or graphs.

--

.pull-right[
## PyOP2 Data Model

### Mesh topology
* ``Sets`` – cells, vertices, etc
* ``Maps`` – connectivity between entities in different sets

### Data
* ``Dats`` – Defined on sets (hold pressure, temperature, etc)

### Kernels / parallel loops
* Executed in parallel on a set through a parallel loop
* Read / write / increment data accessed via maps

### Linear algebra
* Sparsities defined by mappings
* Matrix data on sparsities
* Kernels compute a local matrix – PyOP2 handles global assembly
]

???

PyOP2 uses a number of simple primitives to describe unstructured meshes and
data defined on them:
* Set: abstractly defines class of entities, only know how big it is
* Map: defines connectivity between elements of sets, lookup table (which
    entities of the target Set associated with an entitity of the source Set)
* Dat: abstracted array, defined on a Set, contains actual values, which can
  live in  CPU or GPU memory
* Kernels / parallel loops:
  * define operations/computations to be performed independently for
    every entity of a given Set.
  * executed over the entire Set (or subset) in parallel (parallel loop)
  * can access data associated via Maps with one level of indirection
  * have certain access modes for data they access
* Linear algebra:
  * Matrix defined on sparsity, which are defined by pairs of Maps
  * Parallel loops can assemble matrices with local assembly operation defined
    by kernel
  * PyOP2 takes care of global assembly
* take home: PyOP2 objects are bare/simple objects, but give powerful tools to
  express FE objects/constructs

---

## PyOP2 Architecture

.scale[![PyOP2 implementation](images/pyop2_architecture.svg)]

???

PyOP2 architecture shown in this diagram:
* provides an API to the user (which may be another program, e.g. Firedrake)
* declare data types, exectute parallel loops (dynamic construct)
* PyOP2 runtime schedules computation efficiently (colouring avoids data races)
* code generated and JIT compiled at runtime for different backends
* CPU JIT: shell out to compiler to compile generated kernel + marshalling code
* use ctypes to load the compiled shared object
* No OpenCL and CUDA results this time (focussing on FEM features)

---

class: center, middle

# Finite-element computations with Firedrake

???

How can Firedrake make use of that for FE computations?

---

## Firedrake vs. DOLFIN/FEniCS tool chains

.scale[![Firedrake architecture](images/firedrake_toolchain_dolfin.svg)]

???

Version of the diagram that compares Firedrake and DOLFIN/FEniCS tool chains:

* Key design decision:
  * Python as main language (c.f. language bar in the middle)
  * Cython to speed up performance-critical library code (mesh, sparsity
      building) and interface to 3rd party libraries (PETSc via petsc4py)
  * only lower to C for kernel execution
  * we control C kernels completely, can use ctypes for lightweight interfacing
  * DOLFIN: C++ library exposing a Python interface via SWIG (other way round)
  * decompose and selectively lower high-level Firedrake constructs instead of
    exposing functionality of C++ API to Python via SWIG
* FFC not responsible for optimisation of code (role is only to produce an
  abstract kernel loop nest suitable for optimisation by COFFEE)
  * DOLFIN: FFC generates C++ string conforming to UFC interface
  * Firedrake: FFC generates a kernel suitable for execution by PyOP2
* UFC vs. parallel loops
  * in UFC every kernel type has a special interface
  * have to extend UFC if you want to do anything that is not yet specified
  * parallel loop interface is completely flexible
* escape hatch for things not expressible in UFL: e.g. compute maximum of CG
  and DG functions for every CG DOF (slope limiters)
* PyOP2 as parallel execution layer for assembly kernels: responsible for
  storage, transfer and communication of data
* PETSc used for meshes (DMPlex), nonlinear solves (SNES), linear solves (KSP, PC)
* *No parallel code*: parallelism handled by PyOP2 + PETSc

---

.left30[
### Function
Field defined on a set of degrees of freedom (DoFs), data stored as PyOP2 `Dat`
### FunctionSpace
Characterized by a family and degree of FE basis functions, defines DOFs for function and relationship to mesh entities
### Mesh
Defines abstract topology by sets of entities and maps between them (PyOP2 data structures)
]

.right70[
## Firedrake concepts
![Firedrake types](images/firedrake_types.svg)
]

???

Firedrake uses same concepts of Functions defined on FunctionSpaces defined on a
Mesh as DOLFIN, however all of these are implemented in terms of the PyOP2 data
structures I showed earlier:
* Function's data stored as PyOP2 Dat (Firedrake does not directly touch field
    data!)
* Dat defined on the DOFs DataSet of the FunctionSpace (defined on the nodes Set)
* nodes Set related to Sets of cells, exterior/interior facets of the Mesh via Maps

---

## Driving Finite-element Computations in Firedrake

Solving the Helmholtz equation in Python using Firedrake:
`$$\int_\Omega \nabla v \cdot \nabla u - \lambda v u ~dV = \int_\Omega v f ~dV$$`

```python
from firedrake import *

# Read a mesh and define a function space
mesh = Mesh('filename')
V = FunctionSpace(mesh, "Lagrange", 1)

# Define forcing function for right-hand side
f = Expression("- (lmbda + 2*(n**2)*pi**2) * sin(X[0]*pi*n) * sin(X[1]*pi*n)",
               lmbda=1, n=8)

# Set up the Finite-element weak forms
u = TrialFunction(V)
v = TestFunction(V)

lmbda = 1
a = (dot(grad(v), grad(u)) - lmbda * v * u) * dx
L = v * f * dx

# Solve the resulting finite-element equation
p = Function(V)
*solve(a == L, p)
```

???

Now that we've discussed the concepts, let's look at how to drive FE
computations in Firedrake

This slide shows the code for solving the Helmholtz equation and should be
familiar to any DOLFIN user modulo the import from Firedrake.

---

## Behind the scenes of the solve call

* Firedrake always solves nonlinear problems in resdiual form `F(u;v) = 0`
* Transform linear problem into residual form:
  ```python
  J = a
  F = ufl.action(J, u) - L
  ```
  * Jacobian known to be `a`
  * **Always** solved in a single Newton (nonlinear) iteration
* Use Newton-like methods from PETSc SNES
* PETSc SNES requires two callbacks to evaluate residual and Jacobian:
  * evaluate residual by assembling residual form
    ```python
    assemble(F, tensor=F_tensor)
    ```
  * evaluate Jacobian by assembling Jacobian form
    ```python
    assemble(J, tensor=J_tensor, bcs=bcs)
    ```

???

* What happens when the final solve is called on the previous slide?
* Unified solving interface even behind the scenes

* If Jacobian not provided by the user, Firedrake uses automatic differentiation:
  ```python
  J = ufl.derivative(F, u)
  ```
* Kernels generated by FFC and executed as PyOP2 parallel loops
  * Firedrake builds PyOP2 parallel loop call, using FFC-generated kernel
  * iterate over cells (for cell integrals) or facets (interior/exterior facet integrals)
  * output tensor depends on rank of the form (PyOP2 `Mat`, `Dat` or `Global`)
  * input arguments: coordinate field and any coefficients present in the form

---
name: bcs

## Applying boundary conditions

* Always preserve symmetry of the operator
* Avoid costly search of CSR structure to zero rows/columns
* Zeroing during assembly, but requires boundary DOFs:
  * negative row/column indices for boundary DOFs during addto
  * instructs PETSc to drop entry, leaving 0 in assembled matrix

???

* How can we call assembly before knowing final BCs?
* BCs may change between point of assembly and solve
* assembly returns unassembled matrix with assembly "thunk" (recipe), called with BCs when solving
* Assembly is cached
  * pre-assembly not required in most circumstances
  * Matrices record BCs they have been assembled with, no need for reassembly
  * assembly cache has FIFO eviction strategy

---
template: bcs

## Preassembly

```python
A = assemble(a)
b = assemble(L)
solve(A, p, b, bcs=bcs)
```

---
template: bcs

## Preassembly

```python
*A = assemble(a)  # A unassembled, A.thunk(bcs) not yet called
b = assemble(L)
solve(A, p, b, bcs=bcs)
```

---
template: bcs

## Preassembly

```python
A = assemble(a)  # A unassembled, A.thunk(bcs) not yet called
b = assemble(L)
*solve(A, p, b, bcs=bcs)  # A.thunk(bcs) called, A assembled
```

---
template: bcs

## Preassembly

```python
A = assemble(a)  # A unassembled, A.thunk(bcs) not yet called
b = assemble(L)
solve(A, p, b, bcs=bcs)  # A.thunk(bcs) called, A assembled
# ...
*solve(A, p, b, bcs=bcs)  # bcs consistent, no need to reassemble
```

---
template: bcs

## Preassembly

```python
A = assemble(a)  # A unassembled, A.thunk(bcs) not yet called
b = assemble(L)
solve(A, p, b, bcs=bcs)  # A.thunk(bcs) called, A assembled
# ...
solve(A, p, b, bcs=bcs)  # bcs consistent, no need to reassemble
# ...
*solve(A, p, b, bcs=bcs2)  # bcs differ, reassemble, call A.thunk(bcs2)
```

---

## Distributed Parallel Computations with MPI

.scale[![Decomposed mesh](images/pyop2_mpi_mesh.svg)]

???

* Enforce constraint on local mesh numbering for efficient comp-comm overlap
* Local mesh entities partioned into four consecutive sections
  * **Core:** Entities owned by this processor which can be processed without
    accessing halo data.
  * **Owned:** Entities owned by this processor which access halo data when
    processed.
  * **Exec halo:** Off-processor entities redundantly executed over
    because they touch owned entities.
  * **Non-exec halo:** Off-processor entities which are not processed, but
    read when computing the exec halo.
* Computations on boundaries require up-to-date *halo* data
* Entities that do not touch the boundary (core entities, by construction) can
  be computed while halo data exchange is in flight
* Halo exchange is automatic and happens only if halo is "dirty"

---

## Benchmarks

### Hardware
* Intel Xeon E5-2620 @ 2.00GHz (Sandy Bridge)
* 16GB RAM

### Compilers
* Intel Compilers 14.0.1
* Intel MPI 3.1.038
* Compiler flags: -O3 -xAVX

### Software
* DOLFIN 389e0269 (April 4 2014)
* Firedrake 570d999 (May 13 2014)
* PyOP2 e775c5e (May 9 2014)

### Problem setup
* DOLFIN + Firedrake: RCM mesh reordering enabled
* DOLFIN: quadrature with optimisations enabled
* Firedrake: quadrature with COFFEE loop-invariant code motion enabled

???

* only quadrature supported in Firedrake

---

## Poisson benchmark

.left-column[
preassembled system

### Solver
CG
### Preconditioner
Hypre Boomeramg
]

.right-column[
```python
V = FunctionSpace(mesh, "Lagrange", degree)

# Dirichlet BC for x = 0 and x = 1
bc = DirichletBC(V, 0.0, [3, 4])

# Test, trial and coefficient functions
u = TrialFunction(V)
v = TestFunction(V)
f = Function(V).interpolate(Expression(
        "10*exp(-(pow(x[0] - 0.5, 2) + \
        pow(x[1] - 0.5, 2)) / 0.02)"))
g = Function(V).interpolate(Expression("sin(5*x[0])"))

# Bilinear and linear forms
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Pre-assemble and solve
u = Function(V)
A = assemble(a, bcs=bc)
b = assemble(L)
bc.apply(b)

solve(A, u, b, solver_parameters=params)
```
]

---

![Poisson single core](plots/Poisson_loglog_dim3_degree3_np1.svg)

### solid: Firedrake, dashed: DOLFIN

---

![Poisson scaling](plots/PoissonParallel_loglog_dim3_degree3_size35.svg)

### solid: Firedrake, dashed: DOLFIN

---

![Poisson speedup](plots/PoissonParallel_plot_dim3_degree3_size35_speedup.svg)

### solid: Firedrake, dashed: DOLFIN

---

### Incompressible Navier-Stokes benchmark (Chorin's method)

.left-column[
preassembled system

### Solver
* GMRES for tentative velocity + velocity correction
* CG for pressure correction

### Preconditioner
* block-Jacobi
* ILU block preconditioner
]

.right-column[
```python
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
u, p = TrialFunction(V), TrialFunction(Q)
v, q = TestFunction(V), TestFunction(Q)

dt = 0.01
nu = 0.01
p_in = Constant(0.0)

noslip = DirichletBC(V, Constant((0.0, 0.0)), (1, 3, 4, 6))
inflow = DirichletBC(Q, p_in, 5)
outflow = DirichletBC(Q, 0, 2)
bcu = [noslip]
bcp = [inflow, outflow]

u0, u1, p1 = Function(V), Function(V), Function(Q)
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
    nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1, L1 = lhs(F1), rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx
```
]

---

![Navier-Stokes single core RHS](plots/NavierStokes_RHS_loglog_np1.svg)

### solid: Firedrake, dashed: DOLFIN

---

![Navier-Stokes single core solve](plots/NavierStokes_solve_loglog_np1.svg)

### solid: Firedrake, dashed: DOLFIN

---

![Navier-Stokes scaling RHS](plots/NavierStokesParallel_RHS_loglog_scale0.2.svg)

### solid: Firedrake, dashed: DOLFIN

---

![Navier-Stokes scaling solve](plots/NavierStokesParallel_solve_loglog_scale0.2.svg)

### solid: Firedrake, dashed: DOLFIN

---

![Navier-Stokes speedup RHS](plots/NavierStokesParallel_RHS_plot_scale0.2_speedup.svg)

### solid: Firedrake, dashed: DOLFIN

---

![Navier-Stokes speedup solve](plots/NavierStokesParallel_solve_plot_scale0.2_speedup.svg)

### solid: Firedrake, dashed: DOLFIN

---

## Summary and additional features

### Summary
* Two-layer abstraction for FEM computation from high-level descriptions
  * Firedrake: a performance-portable finite-element computation framework  
    *Drive FE computations from a high-level problem specification*
  * PyOP2: a high-level interface to unstructured mesh based methods  
    *Efficiently execute kernels over an unstructured grid in parallel*
* Decoupling of Firedrake (FEM) and PyOP2 (parallelisation) layers
* Firedrake concepts implemented with PyOP2/PETSc constructs
* Portability for unstructured mesh applications: FEM, non-FEM or combinations
* Extensible framework beyond FEM computations (e.g. image processing)

--

### Preview: Firedrake features not covered
* Automatic optimization of generated assembly kernels with COFFEE (Fabio's
    talk)
* Solving PDEs on extruded (semi-structured) meshes (Doru + Andrew's talk)
* Building meshes using PETSc DMPlex
* Using fieldsplit preconditioners for mixed problems
* Solving PDEs on immersed manifolds
* ...

---

## Thank you!

Contact: Florian Rathgeber, [@frathgeber](https://twitter.com/frathgeber), <f.rathgeber@imperial.ac.uk>

### Resources

  * **PyOP2** https://github.com/OP2/PyOP2
    * *[PyOP2: A High-Level Framework for Performance-Portable Simulations on Unstructured Meshes](http://dx.doi.org/10.1109/SC.Companion.2012.134)*
      Florian Rathgeber, Graham R. Markall, Lawrence Mitchell, Nicholas Loriant, David A. Ham, Carlo Bertolli, Paul H.J. Kelly,
      WOLFHPC 2012
    * *[Performance-Portable Finite Element Assembly Using PyOP2 and FEniCS](http://link.springer.com/chapter/10.1007/978-3-642-38750-0_21)*
       Graham R. Markall, Florian Rathgeber, Lawrence Mitchell, Nicolas Loriant, Carlo Bertolli, David A. Ham, Paul H. J. Kelly ,
       ISC 2013
  * **Firedrake** https://github.com/firedrakeproject/firedrake
    * *COFFEE: an Optimizing Compiler for Finite Element Local Assembly*
      Fabio Luporini, Ana Lucia Varbanescu, Florian Rathgeber, Gheorghe-Teodor Bercea, J. Ramanujam, David A. Ham, Paul H. J. Kelly,
      submitted
  * **UFL** https://bitbucket.org/mapdes/ufl
  * **FFC** https://bitbucket.org/mapdes/ffc

**This talk** is available at http://kynan.github.io/fenics14 ([source](https://github.com/kynan/fenics14))

Slides created with [remark](http://remarkjs.com)
