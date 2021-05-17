# Benchmark

## bellen.c

The function f is polynomial, but the derivative of solution is not continuous. As a result, all the proposed methods converge at the same speed


![error](images/error_bellen.svg)

Therefore the computation time is in favour of the most simple one: explicit Euler method with interpolation order 1.

![computation_time](images/computation_time_bellen.svg)

## stateArt.c

The function f is polynomial only C1, but the solutions are regular. The higher order methods are to be favored. The interpolation order (indicated by solid, dashed or dotted curves) is to be taken equal to 2  (2 and 3 give the same results !), instead of order 1 which gives poorer results.

The comparison is made with rk4 method, interpolation order = 2, and step size twice as small as the smallest.

![error](images/error_stateArt.svg)

Computation time is in favor of explicit Euler method, but owing to its slow convergence, rk4 is a better choice.

![computation_time](images/computation_time_stateArt.svg)

Remarks: 
* eulerNeutral has a tendency to add energy to the system, it can therfore diverge for small values of h
* eulerImpNeutral is very stable, and has a tendency to stay close to the zero solution for small values of h
* impTrNeutral looks bad only because of a delay in the destabilisation process, hence a shift in the phase of the solution.
* both implicit solutions are very stable, and need quite a big intinial condition to be destabilised and go to the non constant solution (1e-6 instead of 1e-10) 
