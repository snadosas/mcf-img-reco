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

#Anima campo vectorial
def anima_vector_field(hist_x,hist_y,title,interva=10,s_min = -1, s_max = 1,stride=2):
    fig, ax = plt.subplots()
    n = hist_x[0].shape[0]

    # Aided by ChatGPT
    # Define the grid
    #x = np.linspace(-1, 1, n)  # Adjust the number of points to control density
    #y = np.linspace(-1, 1, n)
    x = np.arange(0, n, 1)
    y = np.arange(0, n, 1)
    X, Y = np.meshgrid(x, y,indexing='xy')

    # Funcion para animar
    def animate(i):
        #Imitamos el codigo del lab 0
        ax.clear() # limpiar ejes
        #ax.contour(hist[i], levels=[0], colors='red')
        # Define the vector field components
        f_X = hist_x[i]  # Your x-component array
        f_Y = hist_y[i]  # Your y-component array

        # Subsample the grid for visibility
        X_sub = X[::stride, ::stride]
        Y_sub = Y[::stride, ::stride]
        f_X_sub = f_X[::stride, ::stride]
        f_Y_sub = f_Y[::stride, ::stride]

        # Plot the quiver plot
        plt.quiver(X_sub, Y_sub, f_X_sub, f_Y_sub)
        # Add a reference arrow of length 1
        plt.arrow(n//2, n//2, 0, 1, head_width=0.05, head_length=0.1, color = 'r')  # Upward arrow

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

    animation = matplotlib.animation.FuncAnimation(fig, animate, frames = len(hist_x) , interval = interva , repeat = True)

    return animation

def anima_VF_arr(hist,hist_x,hist_y,title,title2,interva=10,s_min = -1, s_max = 1,stride=2):
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16, 8))
    n = hist_x[0].shape[0]

    # Aided by ChatGPT
    # Define the grid
    #x = np.linspace(-1, 1, n)  # Adjust the number of points to control density
    #y = np.linspace(-1, 1, n)
    x = np.arange(0, n, 1)
    y = np.arange(0, n, 1)
    X, Y = np.meshgrid(x, y,indexing='xy')

    # Funcion para animar
    def animate(i):
        #Imitamos el codigo del lab 0
        ax1.clear() # limpiar ejes
        ax2.clear()
        ###ax1
        #Imitamos el codigo del lab 0
        ax1.clear() # limpiar ejes
        ax1.contour(hist[i], levels=[0], colors='black')
        ax1.contourf(hist[i])
        #ax.imshow(awa2[i], cmap='hot')
        #ax.colorbar()
        ticks = np.linspace(s_min,s_max,5)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_xticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        ax1.set_yticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        ax1.set_title(title)
        #ax.contour(hist[i], levels=[0], colors='red')
        # Define the vector field components
        f_X = hist_x[i]  # Your x-component array
        f_Y = hist_y[i]  # Your y-component array

        # Subsample the grid for visibility
        X_sub = X[::stride, ::stride]
        Y_sub = Y[::stride, ::stride]
        f_X_sub = f_X[::stride, ::stride]
        f_Y_sub = f_Y[::stride, ::stride]

        # Plot the quiver plot
        ax2.quiver(X_sub, Y_sub, f_X_sub, f_Y_sub)
        # Add a reference arrow of length 1
        ax2.arrow(n//2, n//2, 0, 1, head_width=0.05, head_length=0.1, color = 'r')  # Upward arrow

        #ax.contourf(hist[i])
        #ax.imshow(awa2[i], cmap='hot')
        #ax.colorbar()
        ticks = np.linspace(s_min,s_max,5)
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')


        ax2.set_xticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        ax1.set_yticks(labels=ticks,ticks=np.linspace(0,n-1,5))
        ax2.set_title(title2)
        #ax.text(0.85, 0.9, f"t = "+str(i*dt), horizontalalignment ='left', verticalalignment ='center', transform = ax.transAxes)
        #ax.legend(loc = "upper left")

    animation = matplotlib.animation.FuncAnimation(fig, animate, frames = len(hist_x) , interval = interva , repeat = True)

    return animation