import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation # animaciones
import matplotlib.animation

def to_grid(func, N):
    x_values = np.linspace(-1, 1, N)
    y_values = np.linspace(-1, 1, N)

    # Create a grid of coordinates
    X, Y = np.meshgrid(x_values, y_values)

    # Evaluate the function at each point in the grid
    Z = func(X, Y)

    return Z, x_values[1]-x_values[0]

def anima_array(hist,title,interva=10,s_min = -1, s_max = 1):
    fig, ax = plt.subplots()
    n = hist[0].shape[0]
    # Funcion para animar
    def animate(i):
        #Imitamos el codigo del lab 0
        ax.clear() # limpiar ejes
        ax.contour(hist[i], levels=[0], colors='black')
        ax.contourf(hist[i])
        #ax.imshow(awa2[i], cmap='hot')
        #ax.colorbar()
        ticks = np.linspace(s_min,s_max,5)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        plt.yticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        ax.set_title(title)
        #ax.text(0.85, 0.9, f"t = "+str(i*dt), horizontalalignment ='left', verticalalignment ='center', transform = ax.transAxes)
        #ax.legend(loc = "upper left")

    animation = matplotlib.animation.FuncAnimation(fig, animate, frames = len(hist) , interval = interva , repeat = True)

    return animation

def anima_array_2(hist,title,interva=10,s_min = -1, s_max = 1):
    fig, ax = plt.subplots()
    n = hist[0].shape[0]
    # Funcion para animar
    def animate(i):
        #Imitamos el codigo del lab 0
        ax.clear() # limpiar ejes
        ax.contour(hist[i], levels=[0], colors='black')
        #ax.contourf(hist[i])
        #ax.imshow(awa2[i], cmap='hot')
        #ax.colorbar()
        ticks = np.linspace(s_min,s_max,5)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        plt.yticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        ax.set_title(title)
        #ax.text(0.85, 0.9, f"t = "+str(i*dt), horizontalalignment ='left', verticalalignment ='center', transform = ax.transAxes)
        #ax.legend(loc = "upper left")

    animation = matplotlib.animation.FuncAnimation(fig, animate, frames = len(hist) , interval = interva , repeat = True)

    return animation


#Anima con imagen de fondo
def anima_array_imagen(hist,title,img,interva=10,s_min = -1, s_max = 1):
    fig, ax = plt.subplots()
    n = hist[0].shape[0]
    # Funcion para animar
    def animate(i):
        #Imitamos el codigo del lab 0
        ax.clear() # limpiar ejes
        ax.imshow(img,cmap='gray')
        ax.contour(hist[i], levels=[0], colors='red')


        #ax.contourf(hist[i])
        #ax.imshow(awa2[i], cmap='hot')
        #ax.colorbar()
        ticks = np.linspace(s_min,s_max,5)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        plt.yticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        ax.set_title(title)
        #ax.text(0.85, 0.9, f"t = "+str(i*dt), horizontalalignment ='left', verticalalignment ='center', transform = ax.transAxes)
        #ax.legend(loc = "upper left")

    animation = matplotlib.animation.FuncAnimation(fig, animate, frames = len(hist) , interval = interva , repeat = True)

    return animation
