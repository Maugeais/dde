import ctypes
import numpy as np

import matplotlib.pyplot as plt
import subprocess
 

class ReVal(ctypes.Structure):
    _fields_ = [("t", ctypes.POINTER(ctypes.c_double)),("x", ctypes.POINTER(ctypes.c_double))]
    

dde_c = ctypes.CDLL("./dde.so")
dde_c.rk4.restype = ReVal

dde_c.preDerPol(3, ctypes.c_float(1))

#subprocess.check_output(['ls', '-l'])  # All that is technically needed...

def compile(fname) :
    
    print('compilation de ', fname)
    
    root = fname.split('.')[0]
    subprocess.check_output(['gcc', '-shared',  '-Wl,-soname,'+root, '-o', root+'.so', '-fPIC', fname])
    
    fun_c = ctypes.CDLL("./"+root+".so").fun
    tau_c = ctypes.CDLL("./"+root+".so").tau

    return(fun_c, tau_c)

    

def rk4Delay(t0, X0, T, fname) :
    
    fun_c, tau_c = compile(fname);
    
    h = t0[1]-t0[0]
    N = int((T-t0[-1])/h)
    
    N0, dim = X0.shape
    
    t0 = np.array(t0).astype(np.double)
    X0 = np.array(X0).astype(np.double)


    r = dde_c.rk4(dim, N0, t0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), X0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.c_float(T), fun_c, tau_c, N)



    t = np.ctypeslib.as_array((ctypes.c_double * (N+N0)).from_address(ctypes.addressof(r.t.contents)))
    x = np.ctypeslib.as_array((ctypes.c_double * (N+N0)*dim).from_address(ctypes.addressof(r.x.contents)))

    x = x.reshape((N+N0, dim))
    
    print(t.shape, x.shape)
    
    return(t, x)




#t, x = rk4Delay(0, 100, np.array([1, 1]), 'test1.c', 100)

#plt.plot(t, x[:, 0])
#plt.plot(t, x[:, 1])

#plt.figure()
#plt.plot(x[:, 0], x[:, 1])

#plt.show()


def exampledde1() :    
    
    tau0 = 3
    t0 = np.linspace(-10, 0, 10000)
    X0 = np.zeros((len(t0), 2))    
    X0[:, 0] = np.cos(np.pi/2*t0/tau0)
    X0[:, 1] = 0*np.cos(3*np.pi/2*t0/tau0)
    
    #plt.plot(t0, X0[0, :])
    
    t, X = rk4Delay(t0, X0, 600, 'test1.c') #solve(fundde, taudde, 0., 30., h, X0)
    
    plt.plot(t, X[:, 0])
    plt.plot(t, X[:, 1])
    plt.show()
    
    
    
exampledde1()    
