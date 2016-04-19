from dolfin import *
from ..NSCoupled import *
import numpy as np
import matplotlib.pyplot as plt
mesh = Mesh("Triangle_corner_nice.xml")


#plot(mesh) ;interactive()

NS_parameters.update(
    omega = 1.0,
	nu = 1.0/1.0,
    max_error = 1e-13,
    max_iter = 1000,
    plot_interval = 10,
    use_krylov_solvers = True)
    #velocity_degree = -100)# leke med denne)
    #output_timeseries_as_vector = True)#,

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

#dofs mesh
#dmp = DofMapPlotter(V)
#dmp.plot()
#dmp.show()



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
    """
    x_1 = np.linspace(0,-4,1000)
    x_2 = np.linspace(-3.7, -4,1000)
    u_values = []
    val_1 = zeros(len(x_1))
    val_2 = zeros(len(x_2))
    for i in range(len(x_1)):
        #u_values.append(u_[0](array([0.0,x_1[i]])))
        val_1[i] = abs(u_[0](array([0.0,x_1[i]])))
        val_2[i] = abs(u_[0](array([0.0,x_2[i]])))
    #print val

    plt.plot(x_1,val_1); plt.yscale('log');plt.xlim((-4.1,0));plt.ylim((1e-28,1.0));\
    plt.gca().invert_xaxis();plt.grid();plt.xlabel("Distance down centre-line ");plt.ylabel("Absolute transverse velocity") ;plt.show()
    plt.plot(x_2,val_2); plt.yscale('log');plt.xlim((-4.01,-3.7));plt.ylim((1e-28,1e-8));\
    plt.gca().invert_xaxis();plt.grid();plt.xlabel("Distance down centre-line ");plt.ylabel("Absolute transverse velocity"); plt.show()


    """


    #bc0 = DirichletBC(V, (0, 0), nos)
    #bc1 = DirichletBC(V, (1.0, 0), top)
    #uu = project(u_, V)#, bcs=[bc0,bc1])
    #plot(uu)
    interactive()
    psi = Function(Q)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    solve(inner(grad(p), grad(q))*dx == inner(curl(u_), q)*dx, psi, bcs=[DirichletBC(Q, 0, "on_boundary")])
    pa = psi.vector().array().argmin()
    sort = psi.vector().array().argsort()
    #print pa,"--",sort[0], sort[1]
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
    #print "delta x = ", (0.3315-abs(xm_1))
    #print "delta y = ", (0.3555 - abs(ym_1))
    """
    print "-----------Here comes second eddy:------------"

    print "Center second-eddy: x: %.4f, %.4f " %(xm_2, ym_2)
    print "Stream function value at second-eddy: %.4f " %(psi(xm_2,ym_2))
    v_value = mycurl(xm_2,ym_2)
    print "Vorticity value at second-eddy: %.4f "%(v_value)
    """

    plot(psi, title='Streamfunction', interactive=True)
    #plot(mycurl, title='Vorticity', interactive=True)
    file = File("Psi_100_degree_1_isoceles.xdmf");
    file << psi;
