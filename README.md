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
 
**2. Set the initial parameter values:**

     a. Choose the number of cells that are of type M1 (“mutant 1”) and of type M2 (“mutant 2”) by clicking with the left-hand mouse button on the sliders #M1 and #M2 and choose “edit” from the drop-down list that appears. Type the desired number of cells of the chosen cell type in the box “Value” and press OK.

![image](https://user-images.githubusercontent.com/61786710/198336897-60aef4bd-94bd-45da-a70d-d14e2c457063.png)

This can also be done by shifting the wedges along the sliders #M1 and #M2, but because the maximum is large ( 1000 by default, but you may change this) it is hard to get the desired number exactly.

      The number of Wild Type cells is the total number of patches - #M1 - #M2 by default.
      
     b. Similarly, select values (between 0 and 1) for the parameters CellType_mean_strength and CellType  
         _sd_strength (with CellType = WT, M1 or M2).  The default values are 0.5 and 0.1 and it is advised to keep the standard deviation (sd) at 0.1.
         

3. Set the **stopping condition** by clicking with the right-hand mouse button on the red selection arrow of 
    the chooser “Stopping-Condition” and select one of the options:
    ![image](https://user-images.githubusercontent.com/61786710/198337410-55f9c2d3-195a-43ac-a881-0111fc421591.png)
                                  
•	“No Wildtype” runs the model until all wild type cells have been displaced by one of the mutant cell types.

•	“Dominant Lineage” until one linage has overtaken the world.

•	“Dominant CellType” until one cell type has displaced the othe r cell types.

•	“STOP-time” until the number of ticks, set by the slider “STOP-timer”, has been reached.



**Refrences**

1. Moran Random processes
2. Piedrafita et al. 2020
