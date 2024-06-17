
import sys
import os

from  aux_functions import *
from finite_difference import *
from reinitialization import *
import skfmm
import time


def MCF_1(phi_0,n_iter,dx,dt,b):

    phi = phi_0.copy()
    phi_array = [phi]

    for i in range(n_iter):
        grad_x,grad_y = fd_gradient(phi,dx)
        hess_x,hess_y = fd_hessian(phi,dx)
        _, hess_xy = fd_gradient(grad_x,dx)

        norm_grad = np.sqrt( grad_x * grad_x + grad_y * grad_y )
        curv = grad_x * grad_x * hess_y - 2 * grad_x * grad_y * hess_xy  + grad_y * grad_y * hess_x


        phi = phi + (dt * curv * b) / (norm_grad * norm_grad)

        phi_array.append(phi.copy())

    return phi, phi_array

def MCF_2(phi_0, max_iter , dx , dt, b ): #Usa nuestra signed distance function

    #phi = phi_0.copy()
    phi_hist = []
    phi, _, _, _ = reinit_fd(phi_0,max_iter*10,dx,dt,b, smoothed_sgn_3, thresold = 0.001)
    for _ in range(max_iter):
        phi = np.pad(phi, 1, mode='maximum')

        # Compute the Laplacian using central difference scheme
        laplacian = ((np.roll(phi, 1, axis=0) + np.roll(phi, -1, axis=0) + \
                    np.roll(phi, 1, axis=1) + np.roll(phi, -1, axis=1) - \
                    4 * phi) / (dx * dx))

        phi = phi + dt*b*laplacian

        phi, _, _, _ = reinit_fd(phi[1:-1, 1:-1],max_iter*10,dx,dt,b, smoothed_sgn_3, thresold = 0.001)

        phi_hist.append(phi.copy())

    return phi, phi_hist

def MCF_3(phi_0, max_iter , dx , dt, b ): #Codigo modificado de chatGPT y snapshot anterior

    #phi = phi_0.copy()
    phi_hist = []
    phi = skfmm.distance(phi_0, dx=dx, order=1)
    for _ in range(max_iter):
        phi = np.pad(phi, 1, mode='maximum')

        # Compute the Laplacian using central difference scheme
        laplacian = ((np.roll(phi, 1, axis=0) + np.roll(phi, -1, axis=0) + \
                    np.roll(phi, 1, axis=1) + np.roll(phi, -1, axis=1) - \
                    4 * phi) / (dx * dx))

        phi = phi + dt*b*laplacian

        phi = skfmm.distance(phi[1:-1, 1:-1], dx=dx, order=1)

        phi_hist.append(phi.copy())

    return phi, phi_hist

if __name__ == '__main__':
    method = str(sys.argv[1])
    shape = sys.argv[2]
    N = int(sys.argv[3])
    max_iter = int(sys.argv[4])
    dt = float(sys.argv[5])
    b = float(sys.argv[6])
    if len(sys.argv) == 8:
        gif_title = sys.argv[7]
    else:
        gif_title = 'mcf.gif'

    if shape == 'circle':
        shape_f = lambda x,y : x**2 + y**2 - 0.5**2
        title = "Mean Curvature Flow on the circle"
    elif shape == 'heart':
        shape_f = lambda x,y : np.abs( y - 0.75*np.cbrt(x**2) ) - 0.5 * np.sqrt(1-x**2)
        title = "Mean Curvature Flow on the heart"
    else:
        def shape_f(z,w):
            x = z
            y = w
            return eval(shape)
        #sys.exit(0)
        title = 'Mean Curvature Flow'

    phi_0, dx = to_grid(shape_f, N)

    start_time = time.time()

    if method == '1':
        tmp_phi, tmp_phi_array = MCF_1(phi_0,max_iter,dx,dt,b)
    elif method == '2':
        tmp_phi, tmp_phi_array = MCF_2(phi_0, max_iter, dx, dt, b)
    elif method == '3':
        tmp_phi, tmp_phi_array = MCF_3(phi_0,max_iter,float(dx) , dt, b)
    else:
        sys.exit(0)

    print("El proceso se ha demorado --- %s segundos ---" % (time.time() - start_time))

    anima = anima_array_2(tmp_phi_array, title)

    if not os.path.exists('anims'):
        os.makedirs('anims')

    while os.path.exists('anims/' + gif_title):
        from datetime import datetime
        gif_title = unique_filename(gif_title)


    anima.save('anims/' + gif_title, writer='pillow')

