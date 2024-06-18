import numpy as np


#Diferencias finitas para primeras derivadas
# (second order central difference)
def fd_gradient(phi_0,h, pad_val=1):
    if pad_val == 'edge':
        phi_pad = np.pad(phi_0, 1, mode='edge',)
    else:
        phi_pad = np.pad(phi_0, 1, mode='constant', constant_values=(1))

    dif_x = (np.roll(phi_pad, -1, axis=1) - np.roll(phi_pad, 1, axis=1))/(2*h)
    dif_y = (np.roll(phi_pad, -1, axis=0) - np.roll(phi_pad, 1, axis=0))/(2*h)

    return dif_x[1:-1, 1:-1], dif_y[1:-1, 1:-1]

#Diferencias finitas para segundas derivadas
def fd_hessian(phi_0,h, pad_val):
    if pad_val == 'edge':
        phi_pad = np.pad(phi_0, 1, mode='edge',)
    else:
        phi_pad = np.pad(phi_0, 1, mode='constant', constant_values=(1))

    dif_x = (np.roll(phi_pad, -1, axis=1) + np.roll(phi_pad, 1, axis=1) - 2*phi_pad)/(h*h)
    dif_y = (np.roll(phi_pad, -1, axis=0) + np.roll(phi_pad, 1, axis=0) - 2*phi_pad)/(h*h)

    return dif_x[1:-1, 1:-1], dif_y[1:-1, 1:-1]

#Calcula el gradiente usando esquemas no simetricos (forward y backward)
def fd_gradient_nonsym(phi_0,h):
    phi_pad = np.pad(phi_0, 1, mode='constant', constant_values=(1))
    dif_x_forw = (np.roll(phi_pad, -1, axis=1) - phi_pad)/(h)
    dif_x_back = (phi_pad - np.roll(phi_pad, 1, axis=1))/(h)
    dif_y_forw = (np.roll(phi_pad, -1, axis=0) - phi_pad)/(h)
    dif_y_back = (phi_pad - np.roll(phi_pad, 1, axis=0))/(h)

    return dif_x_forw[1:-1, 1:-1],dif_x_back[1:-1, 1:-1], dif_y_forw[1:-1, 1:-1],dif_y_back[1:-1, 1:-1]

def upwind_differencing(a_x,a_y,phi_0,h):
    x_forw, x_back, y_forw, y_back = fd_gradient_nonsym(phi_0, h)

    x_grad = np.zeros_like(x_forw)
    y_grad = np.zeros_like(x_forw)

    for i in range(x_forw.shape[0]):
        for j in range(x_forw.shape[0]):
            if a_x[i,j] > 0:
                x_grad[i,j] = x_back[i,j]
            else:
                x_grad[i,j] = x_forw[i,j]

            if a_y[i,j] > 0:
                y_grad[i,j] = y_back[i,j]
            else:
                y_grad[i,j] = y_forw[i,j]

    return x_grad, y_grad
#Elige que aproximacion usar de acuerdo al esquema de Gudonov
def gudonov_choice(a,x_forw,x_back,y_forw,y_back):
    x_chosen = np.zeros(x_forw.shape)
    y_chosen = np.zeros(y_forw.shape)
    for i in range( x_forw.shape[0] ):
        for j in range( x_forw.shape[1] ):
            if (a[i,j]*x_forw[i,j] > 0) and (a[i,j]*x_back[i,j] > 0):
                x_chosen[i,j] = x_back[i,j]
            elif  (a[i,j]*x_forw[i,j] < 0) and (a[i,j]*x_back[i,j] < 0):
                x_chosen[i,j] = x_forw[i,j]
            elif (a[i,j]*x_forw[i,j] >= 0) and (a[i,j]*x_back[i,j] <= 0):
                x_chosen[i,j] = 0
            elif (a[i,j]*x_forw[i,j] <= 0) and (a[i,j]*x_back[i,j] >= 0):
                if np.abs(a[i,j]*x_forw[i,j]) < np.abs(a[i,j]*x_back[i,j]):
                    x_chosen[i,j] = x_back[i,j]
                else:
                    x_chosen[i,j] = x_forw[i,j]

        for i in range( y_forw.shape[0] ):
            for j in range( y_forw.shape[1] ):
                if (a[i,j]*y_forw[i,j] > 0) and (a[i,j]*y_back[i,j] > 0):
                    y_chosen[i,j] = y_back[i,j]
                elif  (a[i,j]*y_forw[i,j] < 0) and (a[i,j]*y_back[i,j] < 0):
                    y_chosen[i,j] = x_forw[i,j]
                elif (a[i,j]*y_forw[i,j] >= 0) and (a[i,j]*y_back[i,j] <= 0):
                    y_chosen[i,j] = 0
                elif (a[i,j]*x_forw[i,j] <= 0) and (a[i,j]*y_back[i,j] >= 0):
                    if np.abs(a[i,j]*y_forw[i,j]) < np.abs(a[i,j]*y_back[i,j]):
                        y_chosen[i,j] = y_back[i,j]
                    else:
                        y_chosen[i,j] = y_forw[i,j]
    #En realidad creo que puedo hacerlo too en un loop
    return x_chosen, y_chosen

#Formula de Rouy Tourin: Dan una aproximacion para los terminos phi_x ^2, phi_y^2
def RouyTourin_Formula(a,x_forw,x_back,y_forw,y_back):
    x_chosen = np.zeros(x_forw.shape)
    y_chosen = np.zeros(y_forw.shape)


    for i in range( x_forw.shape[0] ):
        for j in range( x_forw.shape[1] ):
            if a[i,j] > 0:
                a_x = x_back[i,j]; b_x = x_forw[i,j]
                a_y = y_back[i,j]; b_y = y_forw[i,j]
            elif a[i,j] < 0:
                a_x = x_forw[i,j]; b_x = x_back[i,j]
                a_y = y_forw[i,j]; b_y = y_back[i,j]
            else:
                a_x = 0;           b_x = 0
                a_y = 0;           b_y = 0

            x_chosen[i,j] = max(max(a_x,0)**2,min(b_x,0)**2)
            y_chosen[i,j] = max(max(a_y,0)**2,min(b_y,0)**2)

    return x_chosen, y_chosen