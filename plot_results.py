#!/usr/bin/env python

from pylab import *
import sys

X = loadtxt(sys.argv[1])

epsvals = X[:,0]
Ces = X[:,1]
Des = X[:,2]
Ies = X[:,3]

grid(True)
plot(epsvals,Ces, 'bo-', ); 
xscale('log'); yscale('log'); 
ylim(0.8,1.1e4)
xlim(0.000009,1.1)
xlabel(r"$log(\epsilon)$") ; ylabel(r"$log(C(\epsilon))$")
savefig("ces.pdf")

clf()
grid(True)
plot(epsvals,Ies, 'bo-', ); xscale('log'); yscale('log')
xlim(0.5e-5,1)
ylim(1,1.2e4)
xlabel(r"$log(\epsilon)$") ; ylabel(r"$log(I(\epsilon))$")
savefig("ies.pdf")

clf()
grid(True)
plot(epsvals,Des, 'bo-', ); xscale('log'); yscale('log')
xlim(0.000005,1.2)
ylim(0.000005,2)
xlabel(r"$log(\epsilon)$") ; ylabel(r"$log(D(\epsilon))$")
savefig("des.pdf")
