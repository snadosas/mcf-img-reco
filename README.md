# Mean Curvature Flow for Image Recognition

This repository contains the code of our project for the **Numerical Analysis of PDEs** course in 2024 fall semester at the FCFM, Universidad de Chile.

The outline of the project is to first simulate mean curvature flow using finite differences, and then modify the flow so that it can recognize images.

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

`` python IR_MCF.py method image_dir N_factor n_iter dt b gif_title``

Where:
- **method** : An integer that refers to the method used to calculate the flow
  - 1 : Uses signed distance function from the scikit.fmm library and the simplified flow equation.
  - 2 : Uses finite differences the approximate the curvature and the norm of the gradient directly.
    (More on those methods can be found in the **context** section)
- **image_dir**: Directory of the image.
- **N_factor** : Refers to the factor the resolution of the image will be rescaled. 
- **n_iter** : Refers to the number of steps the curve will be evolved due to the flow.
- **dt** : Size of the time step.
- **b** : Size of the velocity of the flow.
- **gif_title** : Title of the gif to be created.

## Examples

` python MCF.py 1 'circle'  25 200 0.1 0.01 'circle_example.gif'`

![circle example](/anims/circle_example.gif)

`python MCF.py 3 'heart' 30 150 0.1 0.01 'heart_example.gif'`

![heart example](/anims/heart_example.gif)

## Context

### Mean Curvature Flow

Let $\Omega = [-1,1]^2$, $\varphi: \Omega \to \mathbb{R}$ be a smooth function, let $\Gamma = \{ x \in \Omega, \varphi(x) = 0 \} the level set $0$ of $\varphi$, suppose that $\Gamma$ is a closed curve that is a subset of $(0,1)^2$.

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

this is the case when $\varphi$ is a signed distance function, that is, $|\varphi(x)| = d(x, \Gamma)$ and, for points inside the region defined by $\Gamma$, $\varphi < 0$, and for points outside $\varphi > 0$

### Finite Differences

### Image Recognition

Consider a grayscale image $A$, such image can be discretized into a matrix $I$ where each entry represents the brightness of a pixel, one wishes for the flow to stop itself when there's a jump on $|\nabla I|$, that is, a big change of brightness, this is accomplished by solving the following PDE.

$$
\varphi_t = g(| \nabla I|)|\nabla \varphi|  \text{div} (\frac{\nabla \varphi}{|\nabla \varphi|}) + \nabla g(|\nabla I|) \cdot \nabla \varphi
$$

Where $g$ is a positive strictly decreasing function such that $g(0) = 1$ and $\lim_{s \to + \infty}g(s) = 0$.

One such $g$ is $g(x) = (1+x^2)^{-1}$.

## To Do List
- Fix the method 2 in `MCF.py`
- Fix `IR_MCF.py`

## References

[1] Osher, S., & Fedkiw, R. (2003). Level set methods and dynamic implicit surfaces. In Applied mathematical sciences. 