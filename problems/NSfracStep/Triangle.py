from dolfin import *
from ..NSfracStep import *

#mesh = UnitSquareMesh(20,20)
mesh = Mesh("Triangle_corner_nice.xml")


NS_parameters.update(
	nu = 0.01,
	T = 1.0,
	dt = 0.0001,
	use_krylov_solver = False,
	save_step = 10)

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 0.5)

class Nos(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and not near(x[1],0.5)

top = Top()
nos = Nos()
bound = FacetFunction("size_t", mesh)
bound.set_all(0)
top.mark(bound, 1)
nos.mark(bound,2)
plot(bound); interactive()


def create_bcs(V,** NS_namespace):
	bc00 = DirichletBC(V,30,bound,1)
	bc01 = DirichletBC(V,0,bound,2)
	bc02 = DirichletBC(V,0,bound,1)


	#bc01 = DirichletBC(V,0,top)
	#bc02 = DirichletBC(V,0,right_side)
	#bc03 = DirichletBC(V,0,left_side)
	#bc02 = DirichletBC(V,0,boundary)
	return dict(u0 =[bc01,bc00],u1 = [bc01,bc02],p = [])

