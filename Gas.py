# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 19:34:01 2017

@author: clara
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import linalg as LA
import math
import Ball as B


class Gas():     
    """
    This creates a gas, with N particles of a certain mass and radius. The container radius
    and the scaling temperature are also parameters. Use My_Gas.animation() to get the 
    animation. Use My_Gas.information() to just get the data.
    """
    K_b=1.38 #Boltzmann constant
    def __init__(self,N=20,T=293,radius_N=0.5, container_radius=42,mass_N=16):
        """Instantiate the class with the parameters wanted"""
       
        
        self.container_radius=container_radius
        self.container=B.Ball([0,0],[0,0],np.inf, self.container_radius)
        self.ball_list=[]
        self.radius=radius_N
        self.ideal_temperature=T #scaling factor for the temperature
        self.temperature=0
        self.pressure=0
        self.total_KE=0
        self.N=N
        self.mass=mass_N
        self.dt=0.003 #This timestep is set up for optimal animation speed, it can be changed using the method
        
        
        position_list=np.zeros(self.N)
        #random generation
        #This starts the particles off in a square
        for i in range(0,self.N): #Sets up the length of the square container
            number=math.ceil(np.sqrt(self.N))
            length= np.sqrt(2)*self.container_radius
            random_positions= np.linspace(-(length/2-self.radius),(length/2-self.radius),number)
            
            r=np.zeros(2) #position vector
            v=np.zeros(2) #speed vector
    
            
            for j in range (0,number):
                if i in range(j*number,(j+1)*number):
                    k=i-(j*number)
                    r[1]=random_positions[j]
                    r[0]=random_positions[k]
            

            #scales the particle speed according to the average speed
            # x velocity and y velocity have an average of 0
            velo=np.sqrt((3*self.K_b*T)/self.mass)
            number2=np.linspace(-velo,velo,self.N)
            number3=(np.linspace(-velo,velo,self.N))
            np.random.shuffle(number2)
            np.random.shuffle(number3)
        
            v[0]=number2[i]
            v[1]=number3[i]
       
        
            self.ball_list.append(B.Ball(r,v,self.mass,self.radius))
            np.append(position_list,r)
            
            
        #ensures the ball number/radius isn't too big for the container    
        if number*2*self.radius>(length-2*self.radius):  
            raise Exception ("Ball radius is too big for the container size!")
        else:
            pass
        
        
        
    def modify_timestep(self,new_timestep):
        """Used to access and change the timestep"""
        self.timestep=new_timestep
        return self.timestep
        
    def init_figure(self):
        """
        Initialise the container diagram and add it to the plot.
        This method is called once by the animation method.
        Returns a list or tuple of the 'patches' to be animated. 
        """ 
        # add the big circle to represent the container
        #QUOTE LAB SCRIPT
        ax = plt.axes(xlim=(-(self.container_radius+5), self.container_radius+5), ylim=(-(self.container_radius+5),self.container_radius+5))
        ax.axes.set_aspect('equal')
        BigCirc = plt.Circle([0.0,0.0], self.container.radius, ec = 'b', fill = False, ls = 'solid') #change syntax
        ax.add_artist(BigCirc)
        # initialise the text to be animated and add it to the plot
        self.__text0 = ax.text(-(self.container_radius-0.1),(self.container_radius-1),"f={:4d}".format(0,fontsize=12))
        patches = [self.__text0]
        # add the patches for the balls to the plot
        for b in self.ball_list:
            pch = b.get_patch()
            ax.add_patch(pch)
            patches.append(pch)
        return patches    
    
    
    
    
    def next_frame(self, framenumber):
        """
        Does the next frame of the animation.
        This method is called by the animation method with a single argument 
        representing the frame number.
        Returns a list of the 'patches' being animated.
        """
        
        #empty lists to initiate the gas.
        ball_times=[]
        container_times=[]
        
        
        
        self.__text0.set_text("f={:4d}".format(framenumber))
        patches = [self.__text0]
        
        
        for ball in self.ball_list: #keep all balls times to container
             a=ball.time_to_container(self.container_radius)
             container_times.append(a)
        min_container_time=np.amin(np.array(container_times))
        ball_number=np.argmin(np.array(container_times))
          
            
            
        ball_times=np.ones((len(self.ball_list),len(self.ball_list))) #initiate the matrix of collision times
        for i in range(0,len(self.ball_list)):
            for j in range(0,i+1): 
                ball_times[(i,j)]=np.inf
                
        for i in range (0,len(self.ball_list)): #adds all the times to collision to the array
            for j in range (i+1,len(self.ball_list)):           
                ball_times[(i,j)]=self.ball_list[i].time_to_collision(self.ball_list[j])
                
        closest_ball=np.argmin(ball_times,axis=1) #gives the closest ball to the ball collided
        min_collision_times_per_ball=np.amin(ball_times,axis=1) #min time to collision for each ball. This is a flat array
        min_collision_time=np.amin(min_collision_times_per_ball) #gets the minimum of the flat array
        ball_collided=np.argmin(min_collision_times_per_ball) #operates on flat array, returns ball number that collided.
        
        
        if min_container_time < min_collision_time:
            #the ball collides with container
            collision_time = min_container_time
            new_dt=self.dt
            
            if collision_time<new_dt:
                ball1=self.ball_list[ball_number] #this is the ball that collides with the container
                ball1.collide(self.container)
                ball1.move(new_dt)
                
                #momentum change
                r=ball1.pos() #accesses the position
                theta= np.arctan(r[1]/r[0]) #works out the angle (polar coordinates)
                R_hat= np.array([np.cos(theta),np.sin(theta)]) #calculates radial unit vector
                v_radial= np.abs(np.dot(ball1.v,R_hat)) #finds radial speed
                d_momentum=2*v_radial*ball1.mass #computes small change in momentum
                d_pressure=(d_momentum/new_dt)/(2*np.pi*self.container_radius) #Pressure=Force/Area
                self.pressure += d_pressure
                
                
                for i in self.ball_list: 
                    if i==self.ball_list[ball_number]: 
                        patches.append(i.get_patch()) #the ball that was just collided doesnt move
                    else :
                        if i.time_to_container(self.container_radius) < new_dt:
                            i.collide(self.container)  #prevents balls escaping the container
                            patches.append(i.get_patch())
                            
                        #can do the same with sticking balls but need a think
                        else:
                            i.move(new_dt)
                            patches.append(i.get_patch())
                        #collect dp
                
                    
            else:
                for i in self.ball_list:
                    i.move(new_dt)
                    patches.append(i.get_patch())
            
                
            return patches
            
        else:
            #the two balls collide
            collision_time=min_collision_time
            new_dt=self.dt
            
            if collision_time<new_dt:
                self.ball_list[ball_collided].collide(self.ball_list[closest_ball[ball_collided]]) #extract the ball number from previously
                self.ball_list[ball_collided].move(new_dt)
                self.ball_list[closest_ball[ball_collided]].move(new_dt) 
                
                for i in self.ball_list: #don't move the balls that have already been moved
                    if i==self.ball_list[closest_ball[ball_collided]] or i==self.ball_list[ball_collided]:
                        patches.append(i.get_patch())
                    else :
                        if i.time_to_container(self.container_radius) < new_dt: #same as before
                            i.collide(self.container)
                            patches.append(i.get_patch())
                        else:
                            i.move(new_dt)
                            patches.append(i.get_patch())
                    
            else :
                for i in self.ball_list:
                    i.move(new_dt)
                    patches.append(i.get_patch())
       
                
            return patches
        
            
       
                
    def animation(self,No_frames):
        """
        This method is used to start the animation, it asks the use to input the 
        number of frames for which the animation should run.
        """
        
        fig = plt.figure()
        ax = plt.axes(xlim=(-(self.container_radius+5), self.container_radius+5), ylim=(-(self.container_radius+5),self.container_radius+5))
        ax.axes.set_aspect('equal')  
        
        self.anim = animation.FuncAnimation( fig, 
                                    self.next_frame, 
                                    init_func = self.init_figure, 
                                    frames=No_frames,
                                    interval = 1, 
                                    repeat=False,
                                    blit = True)
        
        plt.show()
        
    def data_output(self, time): # essentially the same as next_frame, but only outputs data
        """
        This method uses the same loop structure as next_frame method, is used to extract
        data only. The user choses a time (in seconds) for how long to run the 
        animation for. A counter is set off and when it is reached, the 
        relevant data is returned in a list.
        """
        
        t=0
        while t<time:
            
            ball_times=[]
            container_times=[]
            KE_ball=[]
            velocity_squared=[]
            velocity=[]
           
            
            for ball in self.ball_list: #keep all balls time to container
                 a=ball.time_to_container(self.container_radius)
                 container_times.append(a)
            min_container_time=np.amin(np.array(container_times))
            ball_number=np.argmin(np.array(container_times))
         
            
            
            ball_times=np.ones((len(self.ball_list),len(self.ball_list))) #initiate the matrix
            for i in range(0,len(self.ball_list)):
                for j in range(0,i+1): 
                    ball_times[(i,j)]=np.inf
                
            for i in range (0,len(self.ball_list)): #adds all the times to collision to the array
                for j in range (i+1,len(self.ball_list)):           
                    ball_times[(i,j)]=self.ball_list[i].time_to_collision(self.ball_list[j])
                
            closest_ball=np.argmin(ball_times,axis=1) #gives the closest ball number 
            min_collision_times_per_ball=np.amin(ball_times,axis=1) #min time to collision for each ball. Flat array
            min_collision_time=np.amin(min_collision_times_per_ball) #compare to other minimum
            ball_collided=np.argmin(min_collision_times_per_ball) #operates on flat array, returns ball number(add 1!)
        
        
            if min_container_time < min_collision_time:
            #the ball collides with container
                collision_time=min_container_time
                new_dt=self.dt
               
            
                if collision_time<new_dt:
                    ball1=self.ball_list[ball_number] 
                    ball1.collide(self.container)
                
                #momentum change
                    r=ball1.pos() #position
                    theta= np.arctan(r[1]/r[0])
                    R_hat= np.array([np.cos(theta),np.sin(theta)])
                    v_radial= np.abs(np.dot(ball1.v,R_hat))
                    d_momentum=2*v_radial*ball1.mass
                    d_pressure=(d_momentum/time)/(2*np.pi*self.container_radius)
                    self.pressure += d_pressure #adds the small pressure change to the variable 
                    
                    
                    for i in self.ball_list: #only moves the balls that haven't collided
                        if i!= self.ball_list[ball_number]:
                            if i.time_to_container(self.container_radius) < new_dt:
                                i.collide(self.container)
                            
                            else:
                                i.move(new_dt)
                                
                                       
                
                    
                else:
                    for i in self.ball_list:
                        i.move(new_dt)
                        
                            
                for b in self.ball_list:
                    a=LA.norm(b.v)**2
                    
                    KE_ball.append(a*b.mass*0.5) #adds the kinetic energy of the ball to the list
                    self.total_KE =+ a*b.mass*0.5 #adds up to the total KE of the system
                    
                    velocity_squared.append(a)#adds the velocity squared to the list
                    velocity.append(a**0.5) #adds the velocity to a list of velocities
                    
                    average_v= sum(velocity_squared)/(len(velocity_squared)) #find average velocity squared
                    self.temperature=(average_v*self.mass)/(3*self.K_b) #use the thermodynamic relation
                    self.velocity=velocity
                
                t+= new_dt #a step has been taken
                
            
            else:
            #the two balls collide
                collision_time=min_collision_time
                new_dt=self.dt
            
                if collision_time<new_dt:
                    self.ball_list[ball_collided].collide(self.ball_list[closest_ball[ball_collided]])
                    self.ball_list[ball_collided].move(new_dt)
                    self.ball_list[closest_ball[ball_collided]].move(new_dt)
                
                    for i in self.ball_list: #don't move the balls that have already been moved
                        if i==self.ball_list[closest_ball[ball_collided]] or i==self.ball_list[ball_collided]:
                            pass    
                        else :
                            if i.time_to_container(self.container_radius) < new_dt:
                                i.collide(self.container)
                                
                            else:
                                i.move(new_dt)
                                
                    
                else :
                    for i in self.ball_list:
                        i.move(new_dt)
                        
        
        
                for b in self.ball_list:
                    a=LA.norm(b.v)**2
                    
                    KE_ball.append(a*b.mass*0.5) #adds the kinetic energy of the ball to the list
                    self.total_KE =+ a*b.mass*0.5 #adds up to the total KE of the system
                    
                    velocity_squared.append(a)#adds the velocity squared to the list
                    velocity.append(a**0.5) #adds the velocity to a list of velocities
                    
                    average_v= sum(velocity_squared)/(len(velocity_squared)) #find average velocity squared
                    self.temperature=(average_v*self.mass)/(3*self.K_b) #use the thermodynamic relation
                    self.velocity=velocity
                
                
                t += new_dt #a step has been taken
                
                
        result=[self.temperature,self.pressure,self.N] #output list                       
        return result
        
    
    def area(self):
        """Returns the area of the container"""
        return ((self.container_radius)**2*(np.pi))

    
    def ideal(self): 
        """Uses the ideal gas law to give the pressure expected"""
        pressure=(self.N*self.K_b*self.ideal_temperature)/(((self.container_radius)**2)*np.pi)
        return pressure

    def information(self,run_time):
        """
        Uses the data_output method to print the statements. Takes for argument
        how long the simulation should run, in seconds.
        """
        result=self.data_output(run_time)
        print('Pressure is', result[1], 'Nm^-1')
        print ('Ideal Pressure is', self.ideal(), 'Nm^-1')
        print('Temperature is', result[0], 'K')
        print('Number of particles is', result[2])
        print('Area of the container is', self.area(), 'm^2')
        
    def __repr__(self):
        """Gas representation"""
        return "This gas has particles N=%g with radius %g and mass %g. Model temperature is %g K and container radius %g"%(self.N,self.radius,self.mass,self.ideal_temperature,self.container_radius)
    
    def __str__(self):
        return "%g,%g,%g,%g,%g"%(self.N,self.radius,self.mass,self.ideal_temperature,self.container_radius)
        
 
