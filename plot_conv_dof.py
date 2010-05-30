# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.yscale("log")
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
axis('equal')
data = numpy.loadtxt("conv_dof.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="H1 error")
legend()


# finalize
show()
