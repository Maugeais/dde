import ctypes
import numpy as np

import matplotlib.pyplot as plt
import subprocess 
import os

class ReVal(ctypes.Structure):
    _fields_ = [
        ("t", ctypes.POINTER(ctypes.c_longdouble)),
        ("x", ctypes.POINTER(ctypes.c_longdouble))
    ]


dde_c = ctypes.CDLL(os.path.split(os.path.realpath(__file__))[0]+"/ddec.so")

#dde_c.rk4Neutral.restype = ReVal
dde_c.eulerNeutral.restype = ReVal
dde_c.eulerImpNeutral.restype = ReVal
dde_c.impTrNeutral.restype = ReVal
dde_c.rk4Neutral.restype = ReVal
dde_c.rk4.restype = ReVal

#dde_c.freeme.argtypes = c_void_p,
dde_c.freeme.restype = None

dll = None

dlclose_func = ctypes.CDLL(None).dlclose
dlclose_func.argtypes = [ctypes.c_void_p]
dlclose_func.restype = ctypes.c_int

def compile(fname) :
    
    global dll
    
    # Il faut effacer la dll si elle existe déjà sinon on garde l'ancienne
    if dll != None :
        dlclose_func(dll._handle)

        
    root = '.'.join(fname.split('.')[:-1])
    
    # Test si le fichier .so n'exite pas, ou bien que la date du .c est > à la date du .so
    
    if  not os.path.isfile("./"+root+".so") or (os.path.getmtime("./"+root+".so") < os.path.getmtime("./"+root+".c") ) :
                        
        print('Compilation de', fname+'...')
    
        subprocess.check_output(['gcc', '-O3', '-shared',  '-Wl,-soname,'+root, '-o', root+'.so', '-fPIC', fname])
        
        
    dll = ctypes.CDLL("./"+root+".so")
    fun_c = dll.fun
    tau_c = dll.tau
    
    return(fun_c, tau_c)

    

def rk4Delay(t0, X0, T, fname, params, alg = 'rk4Neutral', interpOrder = 2, free=False) :
    
    dde_c.preDerPol(interpOrder)

    
    fun_c, tau_c = compile(fname);
        
        
    print('Computation with '+alg+'...')

    h = t0[1]-t0[0]
    N = int((T-t0[-1])/h)
    
    N0, dim = X0.shape
    
    t0 = np.array(t0).astype(np.double)
    X0 = np.array(X0).astype(np.double)
    params = np.array(params).astype(np.double)
    
 
    
    if alg == 'rk4Neutral' :
        
        r = dde_c.rk4Neutral(dim, N0, t0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), 
                  X0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), ctypes.c_float(T), 
                  fun_c, tau_c, N, params.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)))
        
    elif alg == 'eulerNeutral' :
        
        
        r = dde_c.eulerNeutral(dim, N0, t0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), 
                  X0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), ctypes.c_float(T), 
                  fun_c, tau_c, N, params.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)))          
        
    elif alg == 'eulerImpNeutral' :

        r = dde_c.eulerImpNeutral(dim, N0, t0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), 
                  X0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), ctypes.c_float(T), 
                  fun_c, tau_c, N, params.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)))        
        
    elif alg == 'impTrNeutral' :
        r = dde_c.impTrNeutral(dim, N0, t0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), 
                  X0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), ctypes.c_float(T), 
                  fun_c, tau_c, N, params.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)))
        
    else :
        
        print("Unknown function")
        return([0], [0])
    
        # r = dde_c.rk4(dim, N0, t0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), 
        #             X0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), ctypes.c_float(T), 
        #             fun_c, tau_c, N, params.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)))  
    
    t = np.ctypeslib.as_array((ctypes.c_double * (N+N0)).from_address(ctypes.addressof(r.t.contents)))
    
    
    x = np.ctypeslib.as_array((ctypes.c_double * ((N+N0)*dim)).from_address(ctypes.addressof(r.x.contents)))
    
    x = x.reshape((N+N0, dim))

    return(t, x)

def rk4(t0, X0, h, T, fname, params) :
    
    fun_c, tau_c = compile(fname);
        
    N = int((T-t0)/h)
    X0 = np.array(X0).astype(np.double)
    
    dim = X0.shape[0]
            
    params = np.array(params).astype(np.double)
        
    r = dde_c.rk4(dim, N, ctypes.c_float(t0), X0.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)), ctypes.c_float(T), 
                    fun_c, params.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)))  

    try :
        t = np.ctypeslib.as_array((ctypes.c_double * (N)).from_address(ctypes.addressof(r.t.contents)))
        
        x = np.ctypeslib.as_array((ctypes.c_double * ((N)*dim)).from_address(ctypes.addressof(r.x.contents)))
        
        x = x.reshape((N, dim))
                
        return(t, x)
    except :
        print('erreur de convergence !!!')

    
def freeme(a) :
    dde_c.freeme(a.ctypes.data_as(ctypes.POINTER(ctypes.c_longdouble)))

def exampledde1() :    
    
    tau0 = 3
    t0 = np.linspace(-10, 0, 10000)
    X0 = np.zeros((len(t0), 2))    
    X0[:, 0] = np.cos(np.pi/2*t0/tau0)
    X0[:, 1] = 0*np.cos(3*np.pi/2*t0/tau0)
    
    #plt.plot(t0, X0[0, :])
    
    t, X = rk4Delay(t0, X0, 600, 'test1.c', [0.7]) 
    
    plt.plot(t, X[:, 0])
    plt.plot(t, X[:, 1])
    plt.show()
    
if __name__ == '__main__' :
    exampledde1()
    
