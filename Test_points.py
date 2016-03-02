from dolfin import *

mesh = UnitSquareMesh(5,5)
Q = FunctionSpace(mesh, "CG",1)
xx = interpolate(Expression("x[0]"), Q)
yy = interpolate(Expression("x[1]"), Q)
print xx.vector().array()
print yy.vector().array()
plot(mesh,interactive=True)
