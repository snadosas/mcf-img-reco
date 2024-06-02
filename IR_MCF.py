
import sys
import os

from PIL import Image
from  aux_functions import *
from finite_difference import *
import skfmm

# Square
def square(x, y):
    return np.maximum(np.abs(x), np.abs(y)) - 0.9
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

        phi = np.pad(phi, 1, mode='maximum')

        # Compute the Laplacian using central difference scheme
        laplacian = ((np.roll(phi, 1, axis=0) + np.roll(phi, -1, axis=0) + \
                    np.roll(phi, 1, axis=1) + np.roll(phi, -1, axis=1) - \
                    4 * phi) / (dx * dx))

        #euler forward
        phi = phi[1:-1, 1:-1] + dt*b*laplacian[1:-1, 1:-1]*g_I + dt*b*(grad_gI_x*grad_x + grad_gI_y*grad_y)


        phi = skfmm.distance(phi, dx=dx, order=1)

        phi_hist.append(phi.copy())

    return phi, phi_hist

# Flujo segmentador
def SF_2(Img_0,max_iter,dt,b):
    N = Img_0.shape[0]
    phi, dx = to_grid(square,N)
    dx = float(dx)
    #phi = skfmm.distance(phi, dx=dx, order=1)
    phi_hist = [phi]


    grad_I_x,grad_I_y = fd_gradient(Img_0,dx)
    norm_grad_I  = np.sqrt( grad_I_x*grad_I_x + grad_I_y*grad_I_y)
    g_I= g(norm_grad_I)
    grad_gI_x,grad_gI_y = fd_gradient(g_I,dx)

    for _ in range(max_iter):

        grad_x,grad_y = fd_gradient(phi,dx)
        hess_x,hess_y = fd_hessian(phi,dx)
        _, hess_xy = fd_gradient(grad_x,dx)

        norm_grad = np.sqrt( grad_x * grad_x + grad_y * grad_y )
        curv = grad_x * grad_x * hess_y - 2 * grad_x * grad_y * hess_xy  + grad_y * grad_y * hess_x

        phi = np.pad(phi, 1, mode='maximum')

        phi = phi[1:-1, 1:-1] + dt*b*curv*norm_grad*g_I + dt*b*(grad_gI_x*grad_x + grad_gI_y*grad_y)


        #phi = skfmm.distance(phi, dx=dx, order=1)

        phi_hist.append(phi.copy())

    return phi, phi_hist


if __name__ == '__main__':
    method = str(sys.argv[1])
    image_dir = str(sys.argv[2])
    N_factor = int(sys.argv[3])
    max_iter = int(sys.argv[4])
    dt = float(sys.argv[5])
    b = float(sys.argv[6])
    if len(sys.argv) == 8:
        gif_title = sys.argv[7]
    else:
        gif_title = 'SegmentFlow.gif'

    title = 'Segmenting Flow'

    #Rutina para importar imagen usando pillow sugerido por CHATGPT
    image = Image.open(image_dir)

    # Convert the image to grayscale
    gray_image = image.convert('L')

    # Convert the image to a NumPy array
    image_array = np.array(gray_image)

    #Se hace cuadrada la imagen y se reduce de acuerdo a N_factor
    max_asp = max(image_array.shape) // N_factor
    print('Resolution has been reduced to: ', max_asp)

    image_array = np.array(gray_image.resize((max_asp, max_asp)))

    if method == '1':
        test_rec, test_rec_hist = SF_1(image_array, max_iter, dt, b)
    elif method == '2':
        test_rec, test_rec_hist = SF_2(image_array, max_iter, dt, b)
    else:
        sys.exit(0)

    rec_anim = anima_array_imagen(test_rec_hist, title, image_array)

    if not os.path.exists('anims'):
        os.makedirs('anims')


    while os.path.exists('anims/' + gif_title):
        from datetime import datetime
        gif_title = gif_title + '-' + datetime.today().strftime('%Y-%m-%d')

    rec_anim.save('anims/' + gif_title, writer='pillow')
