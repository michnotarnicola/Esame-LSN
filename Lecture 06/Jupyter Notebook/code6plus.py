#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:39:57 2020

@author: acer
"""
from numpy import *
from math import *
from random import *
import pylab as plt
# #EQUILIBRAZIONE
# #-------------------------------------------------#
# graph1 = genfromtxt('E_inst.dat')
# graph2 = genfromtxt('M_inst.dat')
# x=arange(0,200001,1)
# plt.figure(figsize=(16.,6.))
# #plt.clf()
# plt.xlim(0000,20000)
# plt.plot(x,graph1[:], label='E')
# plt.plot(x,graph2[:], label='M')
# plt.xlabel('time (adimensional units)')
# plt.ylabel('(adimensional units)')
# plt.legend()
# plt.grid(True)
# plt.show()

##--------------grafici finali--------#
# ENERGIA INTERNA#

points=100
T = linspace(0.3,3.0,num=points)
beta = 1/T
J = 1.0
Ns = 50
th=zeros(points)
for i in range(points) :
    th[i] = tanh(J/T[i])
thN= th**Ns
ch = 1/th
e = -J*( th + ch*thN )/( 1 + thN )

graph = genfromtxt('TvsE.metro')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2], color='b', label='Metropolis', linestyle='None', marker='o', markersize=2)
graph = genfromtxt('TvsE.gibbs')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2], color='r', label='Gibbs', linestyle='None', marker='o', markersize=2)
plt.plot(T, e, color='C1')
plt.title('Ising 1D, internal energy')
plt.xlabel('T')
plt.ylabel('U/N')
plt.legend()
plt.grid(True)
plt.show()


#CAPACITÀ TERMICA#

points=100
T = linspace(0.3,3.0,num=points)
beta = 1/T
J = 1.0
Ns = 50
th=zeros(points)
for i in range(points) :
    th[i] = tanh(J/T[i])
thN= th**Ns
ch = 1/th
heat=((beta*J)**2)*(((1+thN+(Ns-1)*(th**2)+(Ns-1)*(ch**2)*thN)/(1+thN))-Ns*((th+ch*thN)/(1+thN))**2)

#plt.xlim(0.5,1.12)
#plt.ylim(0.41,0.45)
graph = genfromtxt('TvsC.metro')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2], color='b',label='Metropolis', linestyle='None', marker='o', markersize=2)
graph = genfromtxt('TvsC.gibbs')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2], color='r',label='Gibbs', linestyle='None', marker='o', markersize=2)
plt.plot(T, heat,color='C1')
plt.title('Ising 1D, Heat capacity')
plt.xlabel('T')
plt.ylabel('C_v')
#plt.xlim(0.49,0.51)
#plt.ylim(0.32,0.34)
plt.legend()
plt.grid(True)
#plt.savefig('Cv.png')
plt.show()


# #MAGNETIZZAZIONE#

points=100
T = linspace(0.3,3.0,num=points)
beta = 1/T
J = 1.0
Ns = 50
th=zeros(points)
for i in range(points) :
    th[i] = tanh(J/T[i])
thN= th**Ns
ch = 1/th
h=0.2 #external field
M=zeros(points)
for p in range(points):
    b = 1/T[p]
    l1 = exp(b*J)*cosh(b*h)+sqrt(exp(2*b*J)*cosh(b*h)*cosh(b*h)-2*sinh(2*b*J))
    l2 = exp(b*J)*cosh(b*h)-sqrt(exp(2*b*J)*cosh(b*h)*cosh(b*h)-2*sinh(2*b*J))
    Z = l1**Ns + l2**Ns
    M[p] = (exp(b*J)*sinh(b*h)*((l1**(Ns-1))*(1+exp(b*J)*cosh(b*h)/sqrt(exp(2*b*J)*cosh(b*h)*cosh(b*h)-2*sinh(2*b*J))) 
        + (l2**(Ns-1))*(1-exp(b*J)*cosh(b*h)/sqrt(exp(2*b*J)*cosh(b*h)*cosh(b*h)-2*sinh(2*b*J)))))/(Z)

#plt.xlim(0.5,0.66)
#plt.ylim(0.989,0.9895)
    
graph = genfromtxt('TvsM.metro')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2], color='b',label='Metropolis', linestyle='None', marker='o', markersize=2)
graph = genfromtxt('TvsM.gibbs')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2], color='r',label='Gibbs', linestyle='None', marker='o', markersize=2)
plt.plot(T, M,color='C1')
plt.title('Ising 1D, magnetization M with h = 0.02')
plt.xlabel('T')
plt.ylabel('M')
plt.grid(True)
plt.legend()
plt.show()

#SUSCETTIVITÀ#

points=100
T = linspace(0.3,3.0,num=points)
J = 1.0
Ns = 50
th=zeros(points)
for i in range(points) :
    th[i] = tanh(J/T[i])
thN= th**Ns
ch = 1/th
X=zeros(points)
for p in range(points):
    beta=1/T[p]
    X[p] = beta*exp(2*beta*J)*(1-thN[p])/(1+thN[p])
#plt.xlim(0.49,0.51)
#plt.ylim(68,76)
graph = genfromtxt('TvsX.metro')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2], color='b',label='Metropolis', linestyle='None', marker='o', markersize=2)
graph = genfromtxt('TvsX.gibbs')
plt.errorbar(graph[:,0],graph[:,1], yerr=graph[:,2],color='r', label='Gibbs', linestyle='None', marker='o', markersize=2)
plt.plot(T, X, color='C1')
plt.title('Ising 1D, Susceptibility')
plt.xlabel('T')
plt.ylabel('$\chi$')
plt.legend()
plt.grid(True)
plt.show()
