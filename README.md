# Mean Curvature Flow for Image Recognition

This repository contains the code of our proyect for the course **Numerical Analysis of PDEs** imparted at the FCFM, Universidad de Chile.

The outline of the project is to first simulate mean curvature flow using finite differences, and then modify the flow so that it can recognize images.


## Use

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

## To Do List
- Fix the method 2 in `MCF.py`
- Fix `IR_MCF.py`

## References

[1] Osher, S., & Fedkiw, R. (2003). Level set methods and dynamic implicit surfaces. In Applied mathematical sciences. 