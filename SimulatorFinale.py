#from vpython import canvas,sphere,vector,rate,vec,pi,cos
from vpython import *
def Simulate(a1,a2,e,i,T):
    
    print('The Simulation depicts the Motion of the Binary according to the data inputted.')
    print('The inclination ',i,' is added to the fame in order to display the motion as percieved from Earth. ')
    
    a1 *= 3.086e+16
    a2 *= 3.086e+16
    a1,a2 = sorted((a1,a2))
    print(a1,a2)
    w = 2*pi/T


      
    scene = canvas()




        
        
        



    
    r_1 = a1*(1 - e**2)
    r_2 = a2*(1 - e**2)
    Star_1,Star_2 = sphere(radius = (a2/a1)*5e10 ,make_trail = True,emissive = True,color=color.orange),sphere(radius = 5e10 , make_trail = True, emissive = True)
    light_1 , light_2 = local_light(pos = vec(0,0,0)),local_light(pos = vec(0,0,0),color=color.orange)
    
    r1 = lambda a,t:a/(1-e*cos(t))
    r2 = lambda a,t:a/(1-e*cos(t))
    th2 = lambda t:pi-w*t
    th1 = lambda t:pi + w*t
    h = 1/(100)
    t = 0 
    scene.camera.rotate(angle = i,axis = vec(1,0,0),origin = vec(0,0,0))
    scene.center = vec(0,0,0)

    while True:
        rate(2000)
       
        
        
        
        
        Star_1.pos,Star_2.pos = vec(-2*a1*e + r1(r_1,th1(t))*cos(th1(t)),r1(r_1,th1(t))*sin(th1(t)),0),vec(2*a2*e-r2(r_2,th2(t))*cos(th2(t)),r2(r_2,th2(t))*sin(th2(t)),0)
        
        light_1.pos,light_2.pos = Star_1.pos,Star_2.pos

        t+=h
    

    
    
    
    

    
