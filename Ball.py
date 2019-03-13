# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 19:34:00 2017

@author: clara
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import linalg as LA
import math

 
class Ball:
    
    """
    This class creates a ball, with a specific position and velocity. 
    The mass, radius and color can also be decided.
    """
    K_b=1.38 #defines the Boltzmann constant
    
    def __init__(self,r=[0.0,0.0], v=[0.0,0.0],mass=16, rad=0.5,clr='r' ):
        self.__r=np.array(r, dtype='float') 
        self.v=np.array(v, dtype='float')
        self.radius=rad
        self.mass=mass
        self.__patch = plt.Circle(self.__r, self.radius, fc=clr)
        
        
    def pos(self):
        """Returns the position of the ball"""
        return self.__r
    
    def vel(self):
        """Returns the velocity of the ball"""
        return self.v
    
    def move(self,dt):
        """Moves the ball by a time step dt"""
        self.__r= self.__r+(self.v *dt)
        self.__patch.center = self.__r
        return self.__r
    
    def time_to_collision(self, other):
        """Solves the quadratic equation for the time to the next collision"""
        R= (self.__r-other.__r)
        V= (self.v-other.v)
        RA= self.radius+other.radius #This is for collision with another ball
        a_term= 2*np.dot(R,V)
        b_term= (((a_term)**2)-4*(LA.norm(R)**2-RA**2)*(LA.norm(V)**2)) 
        dt1=(-a_term+np.sqrt(b_term))/(2*(LA.norm(V)**2))
        dt2=(-a_term-np.sqrt(b_term))/(2*(LA.norm(V)**2))
        if dt1<=5e-10 or dt1!=dt1: 
            newdt1=np.inf 
        else:
            newdt1=dt1
        if dt2<=5e-10 or dt2!=dt2:
            newdt2=np.inf
        else:
            newdt2=dt2  #this removes unwanted solutions, and the nan case.
        return np.amin(np.array([newdt1,newdt2]))

    def collide(self,other):#takes the other ball as argument
        """By converting to the center of mass frame, changes the speed of the balls"""
        reduced_mass= (self.mass - other.mass)/(self.mass +other.mass)
        bigger_mass=2/(self.mass+other.mass)
        
        R=self.__r-other.__r
        R_mag= (np.dot(R,R))**0.5 #magnitude of R
        R_hat= R/R_mag
        
        parallel_v1= np.dot(self.v,R_hat)*R_hat
        perpendicular_v1=self.v-parallel_v1
        
        parallel_v2= np.dot(other.v,R_hat)*R_hat
        perpendicular_v2=other.v-parallel_v2
        
        new_parallel_v1= reduced_mass*parallel_v1+(other.mass *bigger_mass*parallel_v2)      
        new_parallel_v2=(self.mass*bigger_mass*parallel_v1)- (reduced_mass*parallel_v2)
        
        
        if other.mass==np.inf: #special case of the container
            new_parallel_v1=-parallel_v1
            new_parallel_v2=0
            mom_2_initial=0.0
            mom_2_final=0.0
            KE_2_initial=0
        
        
        
        #update speeds
        new_v1= new_parallel_v1+perpendicular_v1
        new_v2= new_parallel_v2+perpendicular_v2
        
        #momentum conservation check
        mom_1_initial=self.mass*self.v
        mom_2_initial=other.mass*other.v
        mom_initial= LA.norm(mom_1_initial+mom_2_initial)
        
        mom_1_final=self.mass*new_v1
        mom_2_final=other.mass*new_v2
        mom_final=LA.norm(mom_1_final+mom_2_final)
        
        if other.mass==np.inf:
            mom_initial=LA.norm(mom_1_initial)
            mom_final= LA.norm(mom_1_final)
            
        
        
        if np.abs(1-((mom_final)/(LA.norm(mom_initial))))>0.01: #conservation laws
            print ('No violation of conservation of momentum hypothesis rejected to 1%')
            
        #KE conservation check
        KE_1_initial=self.mass*0.5*(LA.norm(self.v)**2)
        KE_2_initial=other.mass*0.5*(LA.norm(other.v)**2)   
        
        KE_1_final=self.mass*0.5*(LA.norm(new_v1)**2)
        KE_2_final=other.mass*0.5*(LA.norm(new_v2)**2)
        
        if other.mass==np.inf:
            KE_2_initial=0
            KE_2_final=0
        
        if np.abs(1-((KE_1_initial+KE_2_initial)/(KE_1_final+KE_2_final)))>0.01:       
            print ('No violation of conservation of energy hypothesis rejected at 1%')
            
        self.v=new_v1
        other.v=new_v2
        
        return self.v, other.v
    
    
    def time_to_container(self,container_radius): 
        """Solves the quadratic equation for the case of the container"""
        
        R= (self.__r)
        V= (self.v)
        RA= container_radius-self.radius 
        a_term= 2*np.dot(R,V)
        b_term= (((a_term)**2)-4*(LA.norm(R)**2-RA**2)*(LA.norm(V)**2))
        dt1=(-a_term+np.sqrt(b_term))/(2*(LA.norm(V)**2))
        dt2=(-a_term-np.sqrt(b_term))/(2*(LA.norm(V)**2))
        if dt1<=5e-10 or dt1!=dt1:
            newdt1=np.inf
        else:
            newdt1=dt1
        if dt2<=5e-10 or dt2 !=dt2: #removes nan
            newdt2=np.inf
        else:
            newdt2=dt2 #ignores the negative solutions to allow us to find smallest dt
        
        return np.amin(np.array([newdt1,newdt2])) 
        
    def get_patch(self):
        """Used to access the hidden attribute patch"""
        return self.__patch
    
        
    def __repr__(self):
        """Ball representation"""
        return "A ball with (r=array[%g,%g], v=array[%g,%g], mass = %g, radius= %g)"%(self.pos()[0], self.pos()[1],self.v[0],self.v[1], self.mass,self.radius)    

    def __str__(self):
        return "%g,%g,%g,%g,%g,%g" %(self.pos()[0], self.pos()[1],self.v[0],self.v[1], self.mass,self.radius)
     