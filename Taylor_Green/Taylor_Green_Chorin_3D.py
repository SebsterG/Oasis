from dolfin import *
import matplotlib.pyplot as plt
import time
import numpy as np
set_log_active(False)
start_time = time.time()

N = 10
mesh = BoxMesh(Point(-pi, -pi, -pi), Point(pi, pi, pi), N, N, N)
#plot(mesh,interactive=True)

class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool((near(x[0], -pi) or near(x[1], -pi) or near(x[2], -pi)) and
                        (not (near(x[0], pi) or near(x[1], pi) or near(x[2], pi))) and on_boundary)

    def map(self, x, y):
        if near(x[0], pi) and near(x[1], pi) and near(x[2],pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1] - 2.0*pi
            y[2] = x[2] - 2.0*pi
        elif near(x[0], pi) and near(x[1], pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1] - 2.0*pi
            y[2] = x[2]
        elif near(x[1], pi) and near(x[2], pi):
            y[0] = x[0]
            y[1] = x[1] - 2.0*pi
            y[2] = x[2] - 2.0*pi
        elif near(x[1], pi):
            y[0] = x[0]
            y[1] = x[1] - 2.0*pi
            y[2] = x[2]
        elif near(x[0], pi) and near(x[2], pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1]
            y[2] = x[2] - 2.0*pi
        elif near(x[0], pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1]
            y[2] = x[2]
        else:
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - 2.0*pi
V = VectorFunctionSpace(mesh, "CG", 2, constrained_domain=PeriodicBoundary())
Q = FunctionSpace(mesh,"CG", 1,constrained_domain=PeriodicBoundary())
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

PB = PeriodicBoundary()
bound = FacetFunction("size_t", mesh)
bound.set_all(0)
PB.mark(bound,1)
#plot(bound,interactive=True)


nu = 20.0*pi/1600 # Re = 6000
p_0=Expression('1./16.*(cos(2*x[0])+cos(2*x[1]))*(cos(2*x[2])+2)')
u0 = interpolate(Expression(('sin(x[0])*cos(x[1])*cos(x[2])','-cos(x[0])*sin(x[1])*cos(x[2])',"0")),V)

#plot(u0)#,interactive=True)
u1 = Function(V)
p1 = Function(Q)

bcs=[]
bcp=[]

dt = 0.001

k = Constant(dt)
f = Constant((0.0, 0.0,0.0))
nu = Constant(nu)
# first without Pressure
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx +  nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# correction with Pressure
a2 = -k*inner(grad(p),grad(q))*dx
L2 = div(u1)* q *dx

# last step

a3 = inner(u,v)*dx
L3 = inner(u1,v)*dx - k*inner(grad(p1),v)*dx


#ufile = File("results/velocity.pvd")
#pfile = File("results/pressure.pvd")
#curlfile = File("results/curl.pvd")

T = 20.0
t = dt
counter = 0
dKdt = []
while t < T + DOLFIN_EPS:
    # Update pressure boundary condition
    solve(a1==L1,u1,bcs)

    #pressure correction
    solve(a2==L2,p1,bcp)
    #print norm(p1)

    #last step
    solve(a3==L3,u1,bcs)

    u0.assign(u1)
    print "Timestep: ", t
    if (counter%100==0 or counter%100 == 1):
        kinetic_e = assemble(0.5*dot(u1,u1)*dx)/(2*pi)**3
        if (counter%100)==0:
            kinetic_hold = kinetic_e
        if (counter%100)==1:
            dKdt.append((kinetic_e - kinetic_hold)/dt)
            print "kinetic energy: ",kinetic_e
            dissipation_e = assemble(nu*inner(grad(u1), grad(u1))*dx) / (2*pi)**3
            print "dissipation: ", dissipation_e
            #plot(u1,rescale=False)
    	    #ufile << u1
    	    #pfile << p1
            #curl_ = curl(u1)
            #curlfile << project(curl_[2],Q) #project(curl(u1),V)

    #plot(p1,rescale=True)
    counter+=1
    t += dt

print("--- %s seconds ---" % (time.time() - start_time))
#plt.plot(dKdt)
#plt.show()
np.savetxt('dKdt.txt', dKdt, delimiter=',')
