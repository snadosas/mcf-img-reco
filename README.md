# Mean Curvature Flow for Image Recognition

This repository contains the code of our project for the **Numerical Analysis of PDEs** course in 2024 fall semester at the FCFM, Universidad de Chile.

The outline of the project is to first simulate mean curvature flow using finite differences, and then modify the flow so that it can recognize images.

## State of the repository

At the time of this commit, we don't plan to keep updating the repository as the project has finished, the state of the methods that shall later be mentioned is:
- MCF.py:
  - All the three methods work fine, but the best is method 3 that uses scikitfmm, method 2 may have slow times under unlucky parameters and may need manual adjusting
- IR_MCF.py
  - Method 2 does not work in almost all cases, a programming error is suspected.
  - Method 4 work in some cases, it's heavily dependent on manual adjusting so that the reinitialization PDEs can be solved quickly.
  - Method 1 works well
  - Method 3 does not currently work as it does not have a heuristic, this one should be considered deprecated and not useable.
## Usage

Clone the repo and open a terminal in the folder of the repo and run 

`pip install requirements.txt`

 ``MCF.py`` contains the code related to the simulation of mean curvature flow using an implicit function.
It should be called like this inside the project folder 

` python MCF.py method shape N n_iter dt b gif_title='mcf.gif' `

Where:
- **method** : An integer that refers to the method used to calculate the flow
  - 1 : Uses finite differences the approximate the curvature and the norm of the gradient directly.
  - 2 : Uses signed distance function obtained from finite differences and the simplified flow equation.
  - 3 : Uses signed distance function from the scikit.fmm library and the simplified flow equation.
    (More on those methods can be found in the **context** section)
- **shape**: Refers to the function whose level set 0 is the closed curve where the flow will be applied, the only built-in functions are 'circle' and 'heart'. This input can also be given as the argument of a function where the first coordinate is x and the second is y, later the function is defined using python's `eval`.
- **N** : Refers to the resolution of the image
- **n_iter** : Refers to the number of steps the curve will be evolved due to the flow.
- **dt** : Size of the time step.
- **b** : Size of the velocity of the flow.
- **gif_title** : Title of the gif to be created.

Such line will produce a gif in `/anims/gif_title`.

``IR_MCF.py`` contains the code related to the simulation of a square on the corners on an image, such square will be under the motion of a ''mean curvature''-like flow that stops it's movement when there is a great change in the gradient of the image.

`` python IR_MCF.py method image_dir N_factor n_iter dt b gif_title anim_flag``

Where:
- **method** : An integer that refers to the method used to calculate the flow
  - 1 : Uses signed distance function from the scikit.fmm library and the simplified flow equation.
  - 2 : Uses finite differences the approximate the curvature and the norm of the gradient directly.
  - 3 : A modification of method 1 that uses Runge Kutta (deprecated)
  - 4: Uses signed distance function obtained from finite differences and the simplified flow equation.
    (More on those methods can be found in the **context** section)
- **image_dir**: Directory of the image.
- **desired_res** : Refers to the resolution of the image that will be used in the code. 
- **n_iter** : Refers to the number of steps the curve will be evolved due to the flow.
- **dt** : Size of the time step.
- **b** : Size of the velocity of the flow.
- **gif_title** : Title of the gif to be created.
- **anim_flag** : If set to 1, then an animation will be created.
## Examples

` python MCF.py 1 'circle'  25 200 0.1 0.01 'circle_example.gif'`

![circle example](/anims/circle_example.gif)

`python MCF.py 1 '(x**2+y**2)**3 - 4*(x**2)*(y**2)'  50 150 0.1 0.01 'flower_example.gif'`

![flower example](/anims/flower_example_met2.gif)

`python MCF.py 3 'heart' 30 150 0.1 0.01 'heart_example.gif'`

![heart example](/anims/heart_example.gif)


`python .\IR_MCF.py 1 'img/twoapples.jpg' 80 1000 0.05 0.01 'twoapples' 1`

![apple example](/anims/twoapples_met1_res80_iter1000_dt0.05_b0.01_eef.gif)

`python .\IR_MCF.py 4 'img/cat_bg.png' 50 600 0.001 0.01 'catto' 1`

![apple example](/anims/catto_met4_res50_iter600_dt0.001_b0.01.gif)


cat_bg_met4_res50_iter425_dt0.01_b0.01
## Context

### Mean Curvature Flow

Let $\Omega = [-1,1]^2$, $\varphi: \Omega \to \mathbb{R}$ be a smooth function, let $`\Gamma = \{ x \in \Omega, \varphi(x) = 0 \}`$ the level set $0$ of $\varphi$, suppose that $\Gamma$ is a closed curve that is a subset of $(0,1)^2$.

In terms of the function $\varphi$, the following geometrical objects can be calculated:
- The unit normal field $\vec{N} = \frac{\nabla \varphi}{| \nabla \varphi|}$
- The curvature $\kappa = \nabla \cdot \vec{N} =  (\varphi_x^2 \varphi_{yy} - 2\varphi_x\varphi_y \varphi_{xy}+\varphi_y^2\varphi_{xx} )/|\nabla \varphi|^3$

The effect of the flow of a mean curvature normal vector field $\vec{V} = -b\kappa \vec{N}$ on $\Gamma$ where $b > 0$ is some constant, can be written in as a PDE on $\varphi$

$$
  \varphi_t = b\kappa | \nabla \varphi | 
$$

In the case where $| \nabla \varphi | = 1 $, the PDE becomes

$$
  \varphi_t = b \Delta \varphi
$$

this is the case when $\varphi$ is a signed distance function, that is, $|\varphi(x)| = d(x, \Gamma)$ and, for points inside the region defined by $\Gamma$, $\varphi < 0$, and for points outside $\varphi > 0$.

### Finite Differences

Method 1 in MCF.py uses centralized differences to calculate the curvature, this seems to be enough due to the dissipative nature of the equation.

$$
D^{+}_{x} \varphi_{ij} = \frac{\varphi-\varphi}{\Delta x}
$$

$$
D^{-}_{x} \varphi_{i,j} = \frac{\varphi_{i,j}-\varphi_{i-1,j}}{\Delta x}
$$

$$
D^{0}_{x} \varphi_{i,j} = \frac{\varphi_{i+1,j}-\varphi_{i-1,j}}{2\Delta x}
$$

Where the second order derivatives are calculated as follows: $\varphi_{xy} = D^0_y D^0_x \varphi, \varphi_{xx} =D^-_x D^+_x \varphi$. 

Method 2 and 3 in MCF.py a five-point stencil finite-difference is used for the laplacian.
$$
\Delta \varphi_{i,j} = \frac{\varphi_{i,j-1}+\varphi_{i,j+1}+\varphi_{i-1,j}+\varphi_{i+1,j}-4\varphi_{i,j}}{(\Delta x)^2}
$$

Method 2 uses the Rouy-Touring Formula to calculate the norm of the gradient.

If $S(\varphi_0) > 0$ set $\varphi_x ^2 = \max(\max(\Delta^{-}_{x} \varphi, 0)^2,\min(\Delta^{+}_{x} \varphi,0)^2)$. 

If $S(\varphi_{0})<0$ set $\varphi_{x}^{2} = \max(\min (\Delta^{-}_{x} \varphi,\max(\Delta^{+}_{x} \varphi,0)^{2})$.

### Image Recognition

Consider a grayscale image $A$, such image can be discretized into a matrix $I$ where each entry represents the brightness of a pixel, one wishes for the flow to stop itself when there's a jump on $|\nabla I|$, that is, a big change of brightness, this is accomplished by solving the following PDE.

$$
\varphi_t = g(| \nabla I|)|\nabla \varphi|  \text{div} (\frac{\nabla \varphi}{|\nabla \varphi|}) + \nabla g(|\nabla I|) \cdot \nabla \varphi
$$

Where $g$ is a positive strictly decreasing function such that $g(0) = 1$ and $\lim_{s \to + \infty}g(s) = 0$.

One such $g$ is $g(x) = (1+x^2)^{-1}$.



## References

[1] Osher, S., & Fedkiw, R. (2003). Level set methods and dynamic implicit surfaces. In Applied mathematical sciences. 