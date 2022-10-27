# cell_competition_2D
The model is mimicking a Moran process [1] or SP model [2] in 2D

![dense19](https://user-images.githubusercontent.com/61786710/197974874-4587394c-be42-4c61-b6ee-6af143b13e1e.png)

https://user-images.githubusercontent.com/61786710/197976053-24e20499-0c2e-418a-a0f9-197d76221ab1.mp4


**How to Run the Model**
Opening the model brings you in the interface.

On the interface are:

1. A large black square, the canvas, this is the “world” on which patches (containing “cells”) will be displayed when running the model.

3. A number of widgets (interface icons).
    a. To the left of the canvas: to run the model (the grey buttons “setup” and “go” and “go-once”) and to adjust parameters (green sliders and on/off switches).
    b. below the canvas (green choosers and switches): to control the visualisation of the patches on the canvas.
    
3. At the right of the canvas: a set of plots and monitors and an output-window showing the performance of  
    the model.
    
    ![image](https://user-images.githubusercontent.com/61786710/198335981-ee5fec5f-ed68-47cf-a5d9-005015deeefc.png)



**Before running the model: Decide on the Setting**
1. Set the size of the “world”.
     The default settings are: minimal x, y coordinates -20, maximum x, y coordinates +20, i.e. a size of 41 x 41  (= 1681 patches). 
     To change the world-size: click on the “Settings” icon and change min-/max- px/ycor. Note that the Patch size and hence the size of the canvas will change in accordance. To re-adjust the size of the canvas, change Patch size or click with the left-hand mouse button anywhere on the canvas. A drop-down menu appears: click “select”. As a result, the canvas is surrounded by a grey border with black squares at the corners. Next, put the mouse on one of the black squares, click the right-hand mouse button while dragging till the desired canvas size has been obtained.
     
     ![image](https://user-images.githubusercontent.com/61786710/198336113-a022742b-cd90-43a6-954e-a64b69b3af51.png)
 



**Refrences**

1. Moran Random processes
2. Piedrafita et al. 2020
