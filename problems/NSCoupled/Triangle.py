from dolfin import *
from ..NSCoupled import *
mesh = Mesh("Triangle_corner_nice.xml")


#plot(mesh) ;interactive()

NS_parameters.update(
    omega = 0.1,
	nu = 1.0/1000.0,
	#use_krylov_solver = False,
    max_error = 1e-13,
    max_iter = 400,
    plot_interval = 10,
    output_timeseries_as_vector = True)

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 1.0)

class Nos(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and not near(x[1], 1.0)

top = Top()
nos = Nos()
bound = FacetFunction("size_t", mesh)
bound.set_all(0)
top.mark(bound, 1)
nos.mark(bound, 2)
#plot(bound); interactive()


def create_bcs(VQ, **NS_namespace):        
    bc0 = DirichletBC(VQ.sub(0), (0, 0), nos)
    bc1 = DirichletBC(VQ.sub(0), (1.0, 0), top)
    return dict(up = [bc0, bc1])

"""
def theend_hook(u_, p_, mesh, **NS_namespace):
    file = File("Stream.xdmf");
    file << psi;
    plot(u_, title='Velocity')
    interactive()
    plot(p_, title='Pressure')
    interactive()
    try:
        from fenicstools import StreamFunction
        psi = StreamFunction(u_, [], mesh, use_strong_bc=True)
        plot(psi, title='Streamfunction', interactive=True)
        interactive()
    except:
        pass


"""

def theend_hook(u_, p_,V, VQ, Q, p, q, **kw):
    bc0 = DirichletBC(V, (0, 0), nos)
    bc1 = DirichletBC(V, (1.0, 0), top)
    uu = project(u_, V, bcs=[bc0,bc1])
    #plot(uu)
    interactive()
    psi = Function(Q)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    solve(inner(grad(p), grad(q))*dx == inner(curl(u_), q)*dx, psi, bcs=[DirichletBC(Q, 0, "on_boundary")])
    pa = psi.vector().array().argmin()
    sort = psi.vector().array().argsort()
    print pa,"--",sort[0], sort[1]
    xx = interpolate(Expression("x[0]"), Q)
    yy = interpolate(Expression("x[1]"), Q)
    xm_1 = xx.vector()[sort[0]]
    ym_1 = yy.vector()[sort[0]]
    xm_2 = xx.vector()[sort[-1]]
    ym_2 = yy.vector()[sort[-1]]
    print "Center main-eddy: x: %.4f, %.4f " %(xm_1, ym_1)
    print "Stream function value at main-eddy: %.4f " %(psi(xm_1,ym_1))
    mycurl = project(curl(u_), Q, bcs=[DirichletBC(Q, 0, "on_boundary")])
    v_value = mycurl(xm_1,ym_1)
    print "Vorticity value at main-eddy: %.4f "%(v_value)
    
    print "-----------Here comes second eddy:------------"

    print "Center second-eddy: x: %.4f, %.4f " %(xm_2, ym_2)
    print "Stream function value at second-eddy: %.4f " %(psi(xm_2,ym_2))
    v_value = mycurl(xm_2,ym_2)
    print "Vorticity value at second-eddy: %.4f "%(v_value)
    




    plot(psi, title='Streamfunction', interactive=True)
    #plot(mycurl, title='Vorticity', interactive=True)
    file = File("Stream.xdmf");
    file << psi;

"""            
def initialize(x_1, x_2, bcs, **NS_namespace):
    for ui in x_1:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
    for ui in x_2:    
        [bc.apply(x_2[ui]) for bc in bcs[ui]]

def pre_solve_hook(mesh, velocity_degree, **NS_namespace):
    Vv = VectorFunctionSpace(mesh, 'CG', velocity_degree)
    return dict(uv=Function(Vv))

def temporal_hook(q_, tstep, u_, uv, p_, plot_interval, **NS_namespace):
    if tstep % plot_interval == 0:
        assign(uv.sub(0), u_[0])
        assign(uv.sub(1), u_[1])
        plot(uv, title='Velocity')
        plot(p_, title='Pressure')
      
def theend_hook(p, q, u_, p_, uv, mesh, **NS_namespace):
    assign(uv.sub(0), u_[0])
    assign(uv.sub(1), u_[1])
    plot(uv, title='Velocity')
    plot(p_, title='Pressure')

    try:
        #from fenicstools import StreamFunction
        #psi = StreamFunction(uv, [], mesh, use_strong_bc=True)
        psi = Function(Q)
        solve(inner(grad(p), grad(q))*dx == inner(curl(u_), q)*dx, psi, bcs=[DirichletBC(Q, 0, "on_boundary")])
        plot(psi, title='Streamfunction', interactive=True)

        f = File("psi.xdmf")
        f << psi
    except:
        pass
"""
