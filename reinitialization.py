import numpy as np

from finite_difference import *

#Calcula la funcion de signo discontinua
def discont_sgn(phi,dx):
    return np.sign(phi)

#Calcula la función (en grilla) de la funcion signo sugerida arriba
def smoothed_sgn(phi_0,dx):
    phi = phi_0.copy()

    grad_x,grad_y = fd_gradient(phi,dx)
    norm_grad_squared = grad_x * grad_x + grad_y * grad_y

    return (phi) / (phi * phi + dx * dx * norm_grad_squared)

def smoothed_sgn_2(phi_0,dx):
    phi = phi_0.copy()


    return (phi) / (phi * phi + dx * dx )

# Funcion con signo con argumento identico a smth_sgn_2, pero con un thresold que controle que no explota
def smoothed_sgn_3(phi_0,dx):
    phi = phi_0.copy()

    phi = (phi) / (phi * phi + dx * dx )

    n = phi.shape[0]
    #Control de que no explote
    for i in range(n):
        for j in range(n):
            if phi[i,j] > 1:
                phi[i,j] = 1
            elif phi[i,j] < -1:
                phi[i,j] = -1

    from scipy.ndimage import gaussian_filter

    return gaussian_filter(phi, sigma=0.1)

#Evoluciona la funcion a traves de la EDP de reinicializacion
def reinit_fd(phi_0,n_iter,dx,dt,b, sgn_fn, thresold = 0.001):

    phi = phi_0.copy()
    phi_array = [phi]
    phi_x, phi_y = fd_gradient(phi,dx)
    phi_x_array = [phi_x]
    phi_y_array = [phi_y]
    smth_sgn = sgn_fn(phi_0,dx)

    for i in range(n_iter):
        #print("A.MAX: ", np.max(phi))
        x_f,x_b,y_f,y_b = fd_gradient_nonsym(phi,dx)
        #x_chos, y_chos = gudonov_choice(smth_sgn,x_f,x_b,y_f,y_b)
        x_chos, y_chos = RouyTourin_Formula(smth_sgn,x_f,x_b,y_f,y_b)
        norm_grad = np.sqrt( x_chos + y_chos )
        #norm_grad = np.sqrt(x_chos * x_chos + y_chos * y_chos )

        ### Se aplica Runge Kutta
        phi_1 = phi - dt * b * smth_sgn * ( norm_grad - 1)

        x_f,x_b,y_f,y_b = fd_gradient_nonsym(phi_1,dx)
        #x_chos, y_chos = gudonov_choice(smth_sgn,x_f,x_b,y_f,y_b)
        x_chos, y_chos = RouyTourin_Formula(smth_sgn,x_f,x_b,y_f,y_b)
        norm_grad = np.sqrt( x_chos + y_chos )

        phi_2 = phi_1 - dt * b* smth_sgn * ( norm_grad - 1)

        phi_1_2 = 0.75*phi + 0.25*phi_2

        x_f,x_b,y_f,y_b = fd_gradient_nonsym(phi_1_2,dx)
        #x_chos, y_chos = gudonov_choice(smth_sgn,x_f,x_b,y_f,y_b)
        x_chos, y_chos = RouyTourin_Formula(smth_sgn,x_f,x_b,y_f,y_b)
        norm_grad = np.sqrt( x_chos + y_chos )

        phi_3_2 = phi_1_2 - dt* b * smth_sgn * ( norm_grad - 1)

        phi = (1/3)*phi+(2/3)*phi_3_2
        phi_x, phi_y = fd_gradient(phi, dx)


        phi_array.append(phi.copy())
        phi_x_array.append(phi_x.copy())
        phi_y_array.append(phi_y.copy())


        g_x, g_y = fd_gradient(phi,dx)
        slope = dt* b * smth_sgn * ( np.sqrt(g_x*g_x+g_y*g_y) - 1)

        print(np.sum( np.abs(slope) < thresold), phi_0.shape[0]*phi_0.shape[0]*0.8)
        print('m' , np.mean(np.abs(slope)))
        if np.sum(np.abs(slope) < thresold ) > phi_0.shape[0]*phi_0.shape[0]*0.8 :
            print("Breaking on :", i)
            break
        #print('a',np.mean(n[1:-1, 1:-1]), '\t',np.max(n[1:-1, 1:-1]), '\t',np.min(n[1:-1, 1:-1]))
        #print('s',np.mean(s[1:-1, 1:-1]), '\t',np.max(s[1:-1, 1:-1]), '\t',np.min(s[1:-1, 1:-1]))

    if i == n_iter-1:
        print("Max Iteration achieved")
        sys.exit(0)

    return phi,phi_array, phi_x_array, phi_y_array

if __name__ == '__main__':
    import sys
    import os

    from aux_functions import *

    shape = sys.argv[1]
    sgn_f = sys.argv[2]
    N = int(sys.argv[3])
    max_iter = int(sys.argv[4])
    dt = float(sys.argv[5])
    b = float(sys.argv[6])
    if len(sys.argv) == 8:
        gif_title = sys.argv[7]
    else:
        gif_title = 'reinit1.gif'

    if shape == 'circle':
        shape_f = lambda x,y : x**2 + y**2 - 0.5**2
        title = "Reinitialization of the circle"
    elif shape == 'heart':
        shape_f = lambda x,y : np.abs( y - 0.75*np.cbrt(x**2) ) - 0.5 * np.sqrt(1-x**2)
        title = "Reinitialization of the heart"
    else:
        def shape_f(z,w):
            x = z
            y = w
            return eval(shape)
        #sys.exit(0)
        title = 'Reinitialization of the 0 level set'

    phi_0, dx = to_grid(shape_f, N)

    if sgn_f == 'discontinuous':
        tmp_phi, tmp_phi_array, tmp_phi_x_array, tmp_phi_y_array = reinit_fd(phi_0, max_iter, dx, dt, b, discont_sgn)
    elif sgn_f == 'continuous':
        tmp_phi, tmp_phi_array, tmp_phi_x_array, tmp_phi_y_array = reinit_fd(phi_0, max_iter, dx, dt, b, smoothed_sgn_3)
    else:
        sys.exit(0)


    #anima = anima_array(tmp_phi_array, title)
    #anima = anima_vector_field(tmp_phi_x_array,tmp_phi_y_array,title)
    anima = anima_VF_arr(tmp_phi_array,tmp_phi_x_array,tmp_phi_y_array,'Reinitialized Function','Gradient of reinitialized function')

    if not os.path.exists('anims'):
        os.makedirs('anims')

    while os.path.exists('anims/' + gif_title):
        from datetime import datetime
        gif_title = gif_title + '-' + datetime.today().strftime('%Y-%m-%d')

    print("Writing...")
    anima.save('anims/' + gif_title, writer='pillow')

    #import skfmm

    #print(np.max(np.abs(tmp_phi_array[-1] - skfmm.distance(tmp_phi_array[-1], dx=float(dx), order=1))))