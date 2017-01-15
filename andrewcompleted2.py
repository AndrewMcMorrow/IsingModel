# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 20:27:44 2017

@author: fmcmo
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 20:10:15 2017

@author: fmcmo
"""

import numpy as np
import random
import matplotlib.pyplot as plt

n = 16

def coli(n):     #defining colinear matrix
    c = [[0 for x in range(n)] for y in range(n)] 
    j=0
    for i in range(n):    #x and y coordinates 10*10
      for j in range(n):   
        if j%2==1:    #setting it as colinear
          c[i][j] = 1
        else:
          c[i][j] = -1       #making every second column -1 and other column +1

    return np.array(c)

def randomlattice(n):     #defining random lattice
    randl=[]
    pos=1 #positive spins
    neg=-1 #negative
    for i in range(n):      #for loop to create random spin arrays
        row=[]
        for j in range(n):
            if random.random()<0.5: #50-50 chance of each lattice being + or -
                row.append(neg)
            else:
                row.append(pos)
        randl.append(row)     #appending column to each row
    return np.array(randl) 

def metro(lattice,temp):    #this decides whether a spin is flipped or not
    for i in range(200000):  #iterates 400000 times
        p = random.randint(0, n-1) #random row co ordinate      
        q = random.randint(0, n-1) #random column co ordinate
        s = lattice[p,q]    #checks whether co ordinate p,q is + or -
        bc=lattice[(p+1)%n,q]+lattice[p,(q+1)%n]+lattice[(p-1)%n,q] + lattice[p, (q-1)%n]  #adds sign of neighbours, periodic conditions included with %n
        hamo=-2*s*bc #change in hamiltonian
        if hamo >= 0: 
            s *=-1 #lattice flips if hamiltonian >= 0
        elif random.random() < np.exp(hamo/float(temp)):    
            s *=-1 #if hamo < ), check exp(hamo/temp), if a random number between 0 and 1 is less than this, then flip spin
        lattice[p,q]=s   
    return lattice   #shows new spin

#print coli(n)   #prints original lattice

temp=0.2 #initalizing temp  #if using plt.imshow, change temp here to change temp in contourf plot
#print metro(coli(n), temp)           #uncomment this if using plt.imshow
#plt.imshow(metro(coli(n),temp))   #uncomment this if using plt.imshow()
while (temp <5.0):   #runs loop from temp 0.1 to 5.0
    sumoflattices1=np.array(metro(coli(n),temp))  #puts updated 2d matrix into 1d array
    magnetization2=(np.sum(sumoflattices1)/float(n**2))   #calculates magnetization per site
    temp =temp +0.1      #increases temperature in steps of 0.1
    result3=abs(magnetization2)              #obtain absolute value of magnetization
    print "temp = %.1f Magnetization = %.10f" % (temp,result3)
    if (temp >= 0.1 and temp <= 2.0 and result3 >= 0.9 ) or (temp > 2.0) :  #removes extreme outliers at low temperatures 
        plt.plot(temp, result3, 'o')   #plots individual points on graph
    plt.xlabel("Temperate J/Kb")
    plt.ylabel("Magnetization per site")
    plt.xlim(0.0, 6.5)   #limit of x axis
    plt.xticks([0, 1, 2, 3, 4, 5, 6])  
    plt.ylim(0, 1.2)     #limit of y axis
    plt.title("Phase Transition of 2D Ising Model")
    plt.savefig("andrewplot.png", bbox_inches='tight')   #saves figure
plt.axvline(2.25, linestyle='--')  #plots vertical line at critical temperature
