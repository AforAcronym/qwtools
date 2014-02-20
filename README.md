Quantum Wells Tools
=======

Some simple quantum mechanical calculations used for semiconductor heterostructures simulations. The problem is intended to be one-dimensional now and is likely to remain such for ever.

Here I have transfer matrix method (TMM, in tmm.jl) implementation that correctly calculates quantum-confined levels for a single rectangular quantum well (QW). It is not suited for multiple QWs and doesn't calculate wavefunctions (to be done in near future). Function called test there shows how to use the tool.

qwsymm.lj stands for a single rectangular QW computation taken from the direct Shrodinger equation solution and is used to check the initial TMM implementation.

Goals of utils.lj and constants.jl are obvious. One can see that calculations are made in CGS units system.
