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

This can also be done by shifting the wedges along the sliders #M1 and #M2, but because the maximum is large ( 1000 by default, but you may change this) it is hard to get the desired number exactly. The number of Wild Type cells is the total number of patches - #M1 - #M2 by default.
      
b. Similarly, select values (between 0 and 1) for the parameters CellType_mean_strength and CellType  
         _sd_strength (with CellType = WT, M1 or M2).  The default values are 0.5 and 0.1 and it is advised to keep the standard deviation (sd) at 0.1.
         

3. Set the **stopping condition** by clicking with the right-hand mouse button on the red selection arrow of 
    the chooser “Stopping-Condition” and select one of the options:
    ![image](https://user-images.githubusercontent.com/61786710/198337410-55f9c2d3-195a-43ac-a881-0111fc421591.png)
                                  
•	“No Wildtype” runs the model until all wild type cells have been displaced by one of the mutant cell types.

•	“Dominant Lineage” until one linage has overtaken the world.

•	“Dominant CellType” until one cell type has displaced the othe r cell types.

•	“STOP-time” until the number of ticks, set by the slider “STOP-timer”, has been reached.

4. Choose how the patches containing the cells will be displayed on the canvas.
    
    This can be done by clicking with the right-hand mouse button on the red selection arrow of the choosers
    
    “PatchColouring” and “Labelling”.
    
•	“By Lineage” colours the cells by their lineage identification number. 

At the start each cell has a unique lineage-ID. When the model runs, at each tick each cell – if “strong” enough (= stronger than the weakest of  its eight neighbours) and not displaced by others – generates a clone which shares its lineage-ID. If a cell of a given lineage reproduces, the number of patches coloured in accordance to that lineage increase. 

   The lineage-ID of cells can be made visible by choosing the option “By Lineage” from the chooser “Labelling”.

   The type of cell is shown when the option “By Celltype” is selected from the chooser “Labelling”.

   If “high-light-mutated-cells?” is switched ON only “mutant” cell types are shown (wild type cells are coloured white).

•	“By Type and Lineage” colours wild type cells grey, cells of type M1 orange and cells of type M2 violet. The intensity of shading reflects the lineageID of the cell.

   The lineageID of a cell can be made visible by choosing the option “By Lineage” from the chooser “Labelling”. 

   If “high-light-mutated-cells?” is switched ON only “mutant” cell types are shown (wild type cells are coloured white).


•	“By Type and Strength” does the same, but now the intensity of shading is proportional to the strength of the cell (darker = stronger).

   The strength of a cell can be made visible by choosing the option “By Strength” from the chooser “Labelling”. 

   If “high-light-mutated-cells?” is switched ON only “mutant” cell types are shown (wild type cells are coloured white).

•	All “Labelling” options can be combined with each choice of “PatchColouring” and undone by the option “No Label”. 


The number of cells of a lineage (clone size) and mean clone size for the cell types show up in the output-window (at the far right of the interface) if  “list-output?” is switched ON.


The Competition Model “Charlie” has the additional possibility of setting the initial spatial arrangement of the “mutant” cell types:

•	Restricted_Area = “ON”

Mutants are going to be placed at random in two circular areas (area A for cell type M1, area B for cell type M2) of which the sizes (radius) can be set by respectively the sliders Ra and Rb. The sizes are normalised to unit length (i.e. they range from 0 to 1). The position of the circular areas is controlled by the sliders xcorAreaA, ycorAreaA and xcorAreaB, ycorAreaB which range from -1 to +1 ( set at 0,0 would put an area in the centre of the world).

The choices are activated by pressing the “setup” bottom and should be made before running the model (they cannot be changed while the model is running).

The areas can be made visual by “Labelling” =  “By Area” (to undo this, choose the option “No Label”).

The monitors below the sliders inform about the size of the areas (in terms of number of patches)

and the density of the cell types put in the patches. Note that the number of mutant cells to be implemented is limited by the size of the areas.

•	Restricted_Area = “OFF” 

IF random? = “ON”
              The default setting of cell types distributed at random over the whole world will be applied.
              
ELSE (random? = “OFF”)

Cells will be placed at regular distances from each other. For this only one mutant (#M1 > 0) should be implemented and the number of cells should be a of power two (a quadrate; if any other number is chosen it will be rounded to the nearest quadrate).

IF maxspacing =”ON”

 The cells of a mutant cell type will be spaced out as much as possible while putting   
 them at regular distances from each other.
 
ELSE

  The number of patches that separates the cells is zero. 
  
  ( I am working on a version where the between-cell patch distance can be chosen by     
     the user).

After the choice for the settings has been made and implemented, press the “setup” button at the top of the Interface.


**Running the model**

Press the “go-once” button for a single run or the “go” button if the model is to run continuously for a time period defined by the “Stopping-Condition”. The model can be run either “on ticks” or “continuous” and the speed can be regulated by  the at the top of the interface.

Besides the number of ticks (displayed below the speed slider ), the actual run time is presented (below Figure 5) as duration indicator. 

Both choosers and the switch can be changed while the model is running. 




**Refrences**

1. Moran Random processes
2. Piedrafita et al. 2020
