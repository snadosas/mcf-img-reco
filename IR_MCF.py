
import sys
import os

from PIL import Image
from  aux_functions import *
from finite_difference import *
from reinitialization import *
import skfmm
import time

# Square
def square(x, y):
    return np.maximum(np.abs(x), np.abs(y)) - 0.9

#def square(x,y):
#    return (x + 0.5)**2 + (y - 0)**2 - 0.1**2

#La funcion g del enunciado
def g(x):
    return 1/(1+x*x)

# Flujo segmentador #Utiliza signed distance function
def SF_1(Img_0,max_iter,dt,b):
    #Discretización de la curva inicial (cuadrado)
    N = Img_0.shape[0]
    phi, dx = to_grid(square,N)
    dx = float(dx)
    phi = skfmm.distance(phi, dx=dx, order=1)
    phi_hist = [phi]

    #Se calcula el gradiente de g compuesto con la norma del gradiente de Img_0
    grad_I_x,grad_I_y = fd_gradient(Img_0,dx)
    norm_grad_I  = np.sqrt( grad_I_x*grad_I_x + grad_I_y*grad_I_y)
    g_I= g(norm_grad_I)
    grad_gI_x,grad_gI_y = fd_gradient(g_I,dx)

    for _ in range(max_iter):
        # Compute the gradient of phi
        grad_x,grad_y = fd_gradient(phi,dx)

        #phi = np.pad(phi, 1, mode='maximum')
        #phi = np.pad(phi, 1, mode='constant', constant_values=(0))
        #phi = np.pad(phi, 1, mode='edge')
        phi = np.pad(phi, 1, mode='constant', constant_values=(1))

        # Compute the Laplacian using central difference scheme
        laplacian = ((np.roll(phi, 1, axis=0) + np.roll(phi, -1, axis=0) + \
                    np.roll(phi, 1, axis=1) + np.roll(phi, -1, axis=1) - \
                    4 * phi) / (dx * dx))

        #euler forward (central differencing
        slope1 = dt*b*laplacian[1:-1, 1:-1]*g_I

        #La idea es calcular slope2 con upwind differencing
        up_grad_x, up_grad_y = upwind_differencing(grad_gI_x, grad_gI_y, phi[1:-1, 1:-1], dx)
        slope2 = dt*(grad_gI_x*up_grad_x + grad_gI_y*up_grad_y )
        #slope2 =  dt * b * (grad_gI_x * grad_x + grad_gI_y * grad_y)
        #Caso localmente color constante
        for i in range(slope2.shape[0]):
            for j in range(slope2.shape[1]):
                if grad_I_x[i,j] == 0 and grad_I_y[i,j] == 0:
                    slope2[i,j] = 0
                    if laplacian[i+1,j+1] < 0.01: #Caso Curvatura Plana en color plano
                        slope1[i,j] = 0.001


        phi = phi[1:-1, 1:-1] + slope1 + slope2

        #phi = phi[1:-1, 1:-1] + dt*b*laplacian[1:-1, 1:-1]*g_I + dt*b*(grad_gI_x*grad_x + grad_gI_y*grad_y)


        phi = skfmm.distance(phi, dx=dx, order=1)

        phi_hist.append(phi.copy())


    return phi, phi_hist

# Flujo segmentador
def SF_2(Img_0,max_iter,dt,b, sug_timestep = False):
    N = Img_0.shape[0]
    phi, dx = to_grid(square,N)
    dx = float(dx)
    #phi = skfmm.distance(phi, dx=dx, order=1)
    phi_hist = [phi]


    grad_I_x,grad_I_y = fd_gradient(Img_0,dx, 'edge')
    norm_grad_I  = np.sqrt( grad_I_x*grad_I_x + grad_I_y*grad_I_y)
    g_I= g(norm_grad_I)
    grad_gI_x,grad_gI_y = fd_gradient(g_I,dx)

    if sug_timestep:
        suggestion = 1 / ( np.max(np.abs(grad_gI_x))/dx +   np.max(np.abs(grad_gI_y))/dx + (4*b)/(dx*dx) )
        print("dt suggestion: ", suggestion)
        dt = 0.9 * suggestion

    for _ in range(max_iter):

        grad_x,grad_y = fd_gradient(phi,dx,'edge')
        hess_x,hess_y = fd_hessian(phi,dx,'edge')
        _, hess_xy = fd_gradient(grad_x,dx,'edge')

        #print('gx: ', grad_x)

        norm_grad =  grad_x * grad_x + grad_y * grad_y
        curv = grad_x * grad_x * hess_y - 2 * grad_x * grad_y * hess_xy  + grad_y * grad_y * hess_x

        #print('c: ', curv)

        deno = norm_grad + 0.01
        #print(np.mean(np.abs(deno)))

        #print('deno: ',deno)
        #print("------------")
        slope1 = dt*b*curv*g_I/deno
        # La idea es calcular slope2 con upwind differencing
        up_grad_x, up_grad_y = upwind_differencing(grad_gI_x, grad_gI_y, phi, dx)
        slope2 = dt * (grad_gI_x * up_grad_x + grad_gI_y * up_grad_y)

        for i in range(phi.shape[0]):
            for j in range(phi.shape[1]):
                if grad_I_x[i,j] < 0.001 and grad_I_y[i,j] < 0.001:
                    slope2[i,j] = 0

        phi = phi + slope1 + slope2


        #phi = skfmm.distance(phi, dx=dx, order=1)

        # normalization routine
        '''for i in range(phi.shape[0]):
            for j in range(phi.shape[0]):
                if phi[i,j] > 1:
                    phi[i,j] = 1
                elif phi[i,j] < -1 :
                    phi[i, j] = -1

        from scipy.ndimage import gaussian_filter

        phi =  gaussian_filter(phi, sigma=0.1)
        #print('max: ', np.max(np.abs(phi)))'''

        # second normalization routine
        #print('max:', np.max(np.abs(phi)))
        #phi = phi / np.max(np.abs(phi))

        phi_hist.append(phi.copy())

    return phi, phi_hist

# Runge Kutta con  scikit
def SF_3(Img_0,max_iter,dt,b):
    #Discretización de la curva inicial (cuadrado)
    N = Img_0.shape[0]
    phi, dx = to_grid(square,N)
    dx = float(dx)
    phi = skfmm.distance(phi, dx=dx, order=1)
    phi_hist = [phi]

    #Se calcula el gradiente de g compuesto con la norma del gradiente de Img_0
    grad_I_x,grad_I_y = fd_gradient(Img_0,dx)
    norm_grad_I  = np.sqrt( grad_I_x*grad_I_x + grad_I_y*grad_I_y)
    g_I= g(norm_grad_I)
    grad_gI_x,grad_gI_y = fd_gradient(g_I,dx)

    for _ in range(max_iter):
        # Compute the gradient of phi
        grad_x,grad_y = fd_gradient(phi,dx)

        #phi = np.pad(phi, 1, mode='maximum')
        #phi = np.pad(phi, 1, mode='constant', constant_values=(0))
        #phi = np.pad(phi, 1, mode='edge')
        phi = np.pad(phi, 1, mode='constant', constant_values=(1))

        # Compute the Laplacian using central difference scheme
        laplacian = ((np.roll(phi, 1, axis=0) + np.roll(phi, -1, axis=0) + \
                    np.roll(phi, 1, axis=1) + np.roll(phi, -1, axis=1) - \
                    4 * phi) / (dx * dx))

        #RK
        phi_1 = phi[1:-1, 1:-1] + dt*b*laplacian[1:-1, 1:-1]*g_I + dt*b*(grad_gI_x*grad_x + grad_gI_y*grad_y)
        phi_1 = skfmm.distance(phi_1, dx=dx, order=1)
        grad_x, grad_y = fd_gradient(phi_1, dx)
        phi_1 = np.pad(phi_1, 1, mode='constant', constant_values=(1))

        # Compute the Laplacian using central difference scheme
        laplacian = ((np.roll(phi_1, 1, axis=0) + np.roll(phi_1, -1, axis=0) + \
                      np.roll(phi_1, 1, axis=1) + np.roll(phi_1, -1, axis=1) - \
                      4 * phi_1) / (dx * dx))

        # RK


        phi_2 = phi_1[1:-1, 1:-1] + dt * b * laplacian[1:-1, 1:-1] * g_I + dt * b * (
                    grad_gI_x * grad_x + grad_gI_y * grad_y)

        phi_2= skfmm.distance(phi_2, dx=dx, order=1)

        phi_1_2 = 0.75 * phi[1:-1, 1:-1] + 0.25 * phi_2
        phi_1_2 = skfmm.distance(phi_1_2, dx=dx, order=1)
        grad_x, grad_y = fd_gradient(phi_1_2, dx)

        phi_1_2 = np.pad(phi_1_2, 1, mode='constant', constant_values=(1))

        laplacian = ((np.roll(phi_1_2, 1, axis=0) + np.roll(phi_1_2, -1, axis=0) + \
                      np.roll(phi_1_2, 1, axis=1) + np.roll(phi_1_2, -1, axis=1) - \
                      4 * phi_1_2) / (dx * dx))



        phi_3_2 = phi_1_2[1:-1, 1:-1] + dt * b * laplacian[1:-1, 1:-1] * g_I + dt * b * (
                grad_gI_x * grad_x + grad_gI_y * grad_y)

        phi = (1 / 3) * phi[1:-1, 1:-1] + (2 / 3) * phi_3_2
        phi = skfmm.distance(phi, dx=dx, order=1)

        phi_hist.append(phi.copy())

    return phi, phi_hist

# reinitizaled con la reinitializaicon creada por mi
def SF_4(Img_0,max_iter,dt,b):
    #Discretización de la curva inicial (cuadrado)
    N = Img_0.shape[0]
    phi, dx = to_grid(square,N)
    dx = float(dx)
    #phi = skfmm.distance(phi, dx=dx, order=1)
    phi, _, _, _ = reinit_fd(phi, max_iter * 10, dx, dt, b, smoothed_sgn_3, thresold=0.001)
    phi_hist = [phi]

    #Se calcula el gradiente de g compuesto con la norma del gradiente de Img_0
    grad_I_x,grad_I_y = fd_gradient(Img_0,dx)
    norm_grad_I  = np.sqrt( grad_I_x*grad_I_x + grad_I_y*grad_I_y)
    g_I= g(norm_grad_I)
    grad_gI_x,grad_gI_y = fd_gradient(g_I,dx)

    for _ in range(max_iter):
        # Compute the gradient of phi
        grad_x,grad_y = fd_gradient(phi,dx)

        #phi = np.pad(phi, 1, mode='maximum')
        #phi = np.pad(phi, 1, mode='constant', constant_values=(0))
        #phi = np.pad(phi, 1, mode='edge')
        phi = np.pad(phi, 1, mode='constant', constant_values=(1))

        # Compute the Laplacian using central difference scheme
        laplacian = ((np.roll(phi, 1, axis=0) + np.roll(phi, -1, axis=0) + \
                    np.roll(phi, 1, axis=1) + np.roll(phi, -1, axis=1) - \
                    4 * phi) / (dx * dx))

        #euler forward
        slope1 = dt*b*laplacian[1:-1, 1:-1]*g_I
        #La idea es calcular slope2 con upwind differencing
        up_grad_x, up_grad_y = upwind_differencing(grad_gI_x, grad_gI_y, phi[1:-1, 1:-1], dx)
        slope2 = dt*b*(grad_gI_x*up_grad_x + grad_gI_y*up_grad_y )

        #Caso localmente color constante
        for i in range(slope2.shape[0]):
            for j in range(slope2.shape[1]):
                if grad_I_x[i,j] == 0 and grad_I_y[i,j] == 0:
                    slope2[i,j] = 0
                    if laplacian[i+1,j+1] < 0.01: #Caso Curvatura Plana en color plano
                        slope1[i,j] = 0.001


        phi = phi[1:-1, 1:-1] + slope1 + slope2

        #phi = phi[1:-1, 1:-1] + dt*b*laplacian[1:-1, 1:-1]*g_I + dt*b*(grad_gI_x*grad_x + grad_gI_y*grad_y)


        #phi = skfmm.distance(phi, dx=dx, order=1)

        phi, _, _, _ = reinit_fd(phi, max_iter, dx, dt, b, smoothed_sgn_3, thresold=0.001)

        phi_hist.append(phi.copy())


    return phi, phi_hist

#Primitivo con RK
def SF_5(Img_0,max_iter,dt,b, sug_timestep = False):
    N = Img_0.shape[0]
    phi, dx = to_grid(square,N)
    dx = float(dx)
    #phi = skfmm.distance(phi, dx=dx, order=1)
    phi_hist = [phi]


    grad_I_x,grad_I_y = fd_gradient(Img_0,dx, 'edge')
    norm_grad_I  = np.sqrt( grad_I_x*grad_I_x + grad_I_y*grad_I_y)
    g_I= g(norm_grad_I)
    grad_gI_x,grad_gI_y = fd_gradient(g_I,dx)

    if sug_timestep:
        suggestion = 1 / ( np.max(np.abs(grad_gI_x))/dx +   np.max(np.abs(grad_gI_y))/dx + (4*b)/(dx*dx) )
        print("dt suggestion: ", suggestion)
        dt = 0.9 * suggestion

    for _ in range(max_iter):

        grad_x,grad_y = fd_gradient(phi,dx,'edge')
        hess_x,hess_y = fd_hessian(phi,dx,'edge')
        _, hess_xy = fd_gradient(grad_x,dx,'edge')

        #print('gx: ', grad_x)

        norm_grad =  grad_x * grad_x + grad_y * grad_y
        curv = grad_x * grad_x * hess_y - 2 * grad_x * grad_y * hess_xy  + grad_y * grad_y * hess_x

        #print('c: ', curv)

        deno = norm_grad + 0.01
        #print(np.mean(np.abs(deno)))

        #print('deno: ',deno)
        #print("------------")
        slope1 = dt*b*curv*g_I/deno
        # La idea es calcular slope2 con upwind differencing
        up_grad_x, up_grad_y = upwind_differencing(grad_gI_x, grad_gI_y, phi, dx)
        slope2 = dt * (grad_gI_x * up_grad_x + grad_gI_y * up_grad_y)

        for i in range(phi.shape[0]):
            for j in range(phi.shape[1]):
                if grad_I_x[i,j] < 0.001 and grad_I_y[i,j] < 0.001:
                    slope2[i,j] = 0

        phi1 = phi + slope1 + slope2

        #phi2
        grad_x,grad_y = fd_gradient(phi1,dx,'edge')
        hess_x,hess_y = fd_hessian(phi1,dx,'edge')
        _, hess_xy = fd_gradient(grad_x,dx,'edge')

        #print('gx: ', grad_x)

        norm_grad =  grad_x * grad_x + grad_y * grad_y
        curv = grad_x * grad_x * hess_y - 2 * grad_x * grad_y * hess_xy  + grad_y * grad_y * hess_x

        #print('c: ', curv)

        deno = norm_grad + 0.01
        #print(np.mean(np.abs(deno)))

        #print('deno: ',deno)
        #print("------------")
        slope1 = dt*b*curv*g_I/deno
        # La idea es calcular slope2 con upwind differencing
        up_grad_x, up_grad_y = upwind_differencing(grad_gI_x, grad_gI_y, phi1, dx)
        slope2 = dt * (grad_gI_x * up_grad_x + grad_gI_y * up_grad_y)

        for i in range(phi.shape[0]):
            for j in range(phi.shape[1]):
                if grad_I_x[i,j] < 0.001 and grad_I_y[i,j] < 0.001:
                    slope2[i,j] = 0

        phi2 = phi1 + slope1 + slope2


        phi12 = 0.75*phi + 0.25*phi2

        #phi32
        grad_x, grad_y = fd_gradient(phi12, dx, 'edge')
        hess_x, hess_y = fd_hessian(phi12, dx, 'edge')
        _, hess_xy = fd_gradient(grad_x, dx, 'edge')

        # print('gx: ', grad_x)

        norm_grad = grad_x * grad_x + grad_y * grad_y
        curv = grad_x * grad_x * hess_y - 2 * grad_x * grad_y * hess_xy + grad_y * grad_y * hess_x

        # print('c: ', curv)

        deno = norm_grad + 0.01
        # print(np.mean(np.abs(deno)))

        # print('deno: ',deno)
        # print("------------")
        slope1 = dt * b * curv * g_I / deno
        # La idea es calcular slope2 con upwind differencing
        up_grad_x, up_grad_y = upwind_differencing(grad_gI_x, grad_gI_y, phi12, dx)
        slope2 = dt * (grad_gI_x * up_grad_x + grad_gI_y * up_grad_y)

        for i in range(phi.shape[0]):
            for j in range(phi.shape[1]):
                if grad_I_x[i, j] < 0.001 and grad_I_y[i, j] < 0.001:
                    slope2[i, j] = 0

        phi32 = phi12 + slope1 + slope2
        #phi = skfmm.distance(phi, dx=dx, order=1)

        phi = phi/3 + 2*phi32/3
        # normalization routine
        '''for i in range(phi.shape[0]):
            for j in range(phi.shape[0]):
                if phi[i,j] > 1:
                    phi[i,j] = 1
                elif phi[i,j] < -1 :
                    phi[i, j] = -1

        from scipy.ndimage import gaussian_filter

        phi =  gaussian_filter(phi, sigma=0.1)
        #print('max: ', np.max(np.abs(phi)))'''

        # second normalization routine
        #print('max:', np.max(np.abs(phi)))
        #phi = phi / np.max(np.abs(phi))

        phi_hist.append(phi.copy())

    return phi, phi_hist

if __name__ == '__main__':
    method = str(sys.argv[1])
    image_dir = str(sys.argv[2])
    desired_res = int(sys.argv[3])
    max_iter = int(sys.argv[4])
    dt = float(sys.argv[5])
    b = float(sys.argv[6])
    if len(sys.argv) == 8:
        gif_title = sys.argv[7]
    elif len(sys.argv) == 9:
        gif_title = sys.argv[7]
        anim_flag = int(sys.argv[8])
    else:
        gif_title = 'SegmentFlow'
        anim_flag = True

    title = 'Segmenting Flow'

    #Rutina para importar imagen usando pillow sugerido por CHATGPT
    image = Image.open(image_dir)

    # Convert the image to grayscale
    gray_image = image.convert('L')

    # Convert the image to a NumPy array
    image_array = np.array(gray_image)

    #Se hace cuadrada la imagen y se reduce de acuerdo a N_factor
    #max_asp = max(image_array.shape) // N_factor
    max_asp = desired_res
    print('Resolution has been reduced to: ', max_asp)

    image_array = np.array(gray_image.resize((max_asp, max_asp))) / 255

    start_time = time.time()

    if method == '1':
        test_rec, test_rec_hist = SF_1(image_array, max_iter, dt, b)
    elif method == '2':
        test_rec, test_rec_hist = SF_5(image_array, max_iter, dt, b, True)
    elif method == '3':
        test_rec, test_rec_hist = SF_3(image_array, max_iter, dt, b)
    elif method == '4':
        test_rec, test_rec_hist = SF_4(image_array, max_iter, dt, b)
    else:
        sys.exit(0)

    print("El proceso se ha demorado --- %s segundos ---" % (time.time() - start_time))


    if not os.path.exists('anims'):
        os.makedirs('anims')

    if not os.path.exists('img_out'):
        os.makedirs('img_out')

    #title detail appending
    gif_title = gif_title + "_met" + str(method) + "_res" + str(desired_res) + "_iter" + str(max_iter) + "_dt" + str(dt) + "_b" + str(b)

    while os.path.exists('anims/' + gif_title + '.gif') or os.path.exists('img_out/' + gif_title + '.png'):
        from datetime import datetime
        gif_title = unique_filename(gif_title + '.gif')[:-4]

    if anim_flag:
        start_time = time.time()
        rec_anim = anima_array_imagen(test_rec_hist, title, image_array)
        rec_anim.save('anims/' + gif_title, writer='pillow')
        print("La creacion de la animacion se ha demorado --- %s segundos ---" % (time.time() - start_time))
    # Image saving

    fig = compound_image(image,test_rec_hist[-1])
    gif_title = gif_title[:-4] + '.png'
    fig.savefig('img_out/' + gif_title, bbox_inches='tight', pad_inches=0)

