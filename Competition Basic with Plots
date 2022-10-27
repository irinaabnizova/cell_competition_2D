extensions [palette ]

globals
 [#patches
  strength-list   WTstrength-list    M1strength-list    M2strength-list
  lineage-list    WTlineage-list     M1lineage-list     M2lineage-list
  unique-lineages unique-WTlineages  unique-M1lineages  unique-M2lineages
  #lineages       #WTlineages        #M1lineages        #M2lineages
  clonesize-list  WTclonesize-list   M1clonesize-list   M2clonesize-list  CList
  meanclonesize   meanWTclonesize    meanM1clonesize    meanM2clonesize
  #Displaced
  runtime
 ]

patches-own
 [PID#
  cellType
  lineageID mlineageID M2lineageID
  strength
  dw ds
  state
 ]

to setup
 ca
 set #patches count patches
 set_lists

 ask n-of #M1 patches [set cellType 1]                                                   ;;In "Un-Restricted" mode, chosen numbers of cells with mutant M1 (1) or M2 (2) are distributed at random
 ask n-of #M2 patches with [cellType != 1] [set cellType 2]

 ;;TYPIFYING PATCHES BY PATCHNUMBER (PID#), LINEAGE AND "STRENGTH"
 ask patches
  [
   set PID# (pxcor + max-pxcor) + (max-pycor - pycor) * (2 * max-pycor + 1)                    ;;each patch is numbered from the first (PID# = 0) to the last, going from top left to bottom right
   set lineageID PID#                                                                          ;;each patch is allocated a lineageID, initially identical to its patchID                                                                                                ;;each patch is allocated an initial strength, based on its mutation status

   if cellType = 0 [set strength random-normal WT_mean_strength WT_sd_strength]
   if cellType = 1 [set strength random-normal M1_mean_strength M1_sd_strength]
   if cellType = 2 [set strength random-normal M2_mean_strength M2_sd_strength]
  ]

 ask patches
  [
   update_lists
   update_patch_colors
  ]

 reset-ticks
 reset-timer
end

to go
;;SETTINGS
  set #Displaced 0                                                                           ;;important for tracking stratification events
  set_lists

  ;;DYNAMICS
  ask patches                                                                                ;;the main core (dynamics) of the model
   [
    if CellType = 0 [divide WT_mean_strength  WT_sd_strength]
    if CellType = 1 [divide M1_mean_Strength  M1_sd_Strength]                                ;;Calls for "division" (= creating potential daughter cells), specified for genotype of cells
    if CellType = 2 [divide M2_mean_strength  M2_sd_strength]
   ]

  ask patches
    [
      if state = "DV" [grow]                                                                 ;;Call for reproduction and replacing of neighbours by daughters. Note: "grow" is in ask patches separated from "divide" to avoid re-division
    ]

;;UPDATING
 ask patches
  [
   update_lists
   update_patch_colors
  ]


;;COMPUTATION OF PARAMETERS AND -OUTPUT
 compute_clonesizes                                                                ;;the computation of the parameter values (mean, stdev) is outlined in the section COMPUTATION PARAMETER VALUES below
 produce_output                                                                              ;;output options are specified in the section OUTPUT below

;;STOPPING CONDITIONS
if Stopping-Condition = "Dominant Lineage"
  [
   if #lineages = 1 [set runtime timer stop]
  ]
 if Stopping-Condition = "No Wildtype"
  [
   if count patches with [cellType = 0] = 0 [set runtime timer stop]
  ]
 if Stopping-Condition = "Dominant CellType"
  [
   if (count patches with [cellType = 0 ]) = #patches
        or
      (count patches with [cellType = 1 ]) = #patches
        or
      (count patches with [cellType = 2 ]) = #patches
     [set runtime timer stop]
  ]
 if Stopping-Condition = "STOP-time"
  [
    if ticks = STOP-Time [set runtime timer stop ]
  ]
  tick
end

;;DYNAMICS: LOCATING, DIVIDING & GROWING
to divide [m sd]                                                                             ;;Dividing = Make potential daughter-cells but only after weakest neighbour has been identified
   ifelse strength >= [strength] of weakest-neighbor                                         ;;and if progenitor (= maternal cell) is not the weakest in the neighbourhood
    [
     set state "DV"                                                                          ;;DV = "Dividing"
     set ds 0
     repeat 2
      [
       let strengthd (random-normal m sd) * ( 1 + (strength - m))                            ;;Assessing strength of strongest (ds) and weakest (dw) daughter
       ;let strengthd (random-normal m sd) * ( 1 + (strength - mean [strength] of neighbors))
       if strengthd > 1 [set strengthd 1]
       ifelse strengthd > ds [set dw ds set ds strengthd][set dw strengthd]
      ]
    ]
    [set state "ND"]
end

to-report weakest-neighbor
  report first sort-on [ strength ] neighbors
end

to grow                                                                                      ;;Clone grows by reproducing and daughter cell taking place of nearest neighbour
  ifelse dw >= [strength] of weakest-neighbor
   [
    let StrengthDw dw
    let CellTypeDw cellType
    let LineageDw  lineageID
    set state "RP"                                                                           ;;RP = "Replacing"
    ask weakest-neighbor                                                                     ;;Weakest neighbour becomes dw, or: dw takes position of ("conquers") weakest neighbor
     [
      set strength  StrengthDw
      set lineageID LineageDw
      set cellType  CellTypeDw
      set state "ST"                                                                         ;;ST = "Stratified"
      set #Displaced #Displaced + 1                                                          ;;Taking position by dw implies weakest neighbour is "stratified", i.e. a displacement occurs
     ]
   ]
  [
   set state "NR"                                                                            ;;NR = "Not Replacing"
   set strength ds                                                                           ;;daughter becomes "mother" cell to reproduce daughters in the next GO (=next generation)
  ]
end

;;COMPLETING LISTS
to set_lists
  set WTstrength-list [] set M1strength-list []  set M2strength-list []
  set WTlineage-list  [] set M1lineage-list  []  set M2lineage-list  []
  set lineage-list    [] set  CList []
end

to update_lists
  set lineage-list lput lineageID lineage-list                                                ;;lineageID of each patch is put in a list, one after another
  if CellType = 1 [set M1lineage-list  lput lineageID M1lineage-list]
  if CellType = 2 [set M2lineage-list  lput lineageID M2lineage-list]
  if CellType = 0 [set WTlineage-list  lput lineageID WTlineage-list]

  set unique-lineages sort remove-duplicates lineage-list
  if CellType = 1 [set unique-M1lineages  remove-duplicates M1lineage-list]             ;;duplicates of lineageID are removed to obtain unique lineageIDs
  if CellType = 2 [set unique-M2lineages  remove-duplicates M2lineage-list]
  if CellType = 0 [set unique-WTlineages  remove-duplicates WTlineage-list]

  set #lineages length unique-lineages
  if CellType = 1 [set #M1lineages  length unique-M1lineages]
  if CellType = 2 [set #M2lineages  length unique-M2lineages]
  if CellType = 0 [set #WTlineages  length unique-WTlineages]

  if CellType = 1 [set M1strength-list  lput (precision strength 2) M1strength-list ]
  if CellType = 2 [set M2strength-list  lput (precision strength 2) M2strength-list ]
  if CellType = 0 [set WTstrength-list  lput (precision strength 2) WTstrength-list ]
end

;COIMPUTING CLONE SIZES
to compute_clonesizes
 compute-clonesize-lists #WTlineages unique-WTlineages
 set WTclonesize-list CList
 set meanWTclonesize mean CList
 set CList []

 if any? patches with [CellType = 1]
   [
    compute-clonesize-lists #M1lineages unique-M1lineages
    set M1clonesize-list CList
    set meanM1clonesize mean CList
    set CList []
   ]

 if any? patches with [CellType = 2]
    [
     compute-clonesize-lists #M2lineages unique-M2lineages
     set M2clonesize-list CList
     set meanM2clonesize mean CList
     set CList []
    ]
end

to compute-clonesize-lists [#L UL]
 let j 0
  repeat #L
   [
    let clonesize count patches with [lineageID = item j UL]
    set Clist lput clonesize Clist
    set j j + 1
   ]
end

;;OUTPUT
to produce_output
   if any? patches with [cellType = 1 or cellType = 2]
  [
    if list-output?
    [output-print (word precision (meanWTclonesize) 2 "    "
                        count patches with [cellType = 1]"   "
                        precision (meanM1clonesize) 2 "    "
                        count patches with [cellType = 2]"   "
                        precision meanM2clonesize 2)
    ]
  ]
 print""
end

to update_patch_colors
  if PatchColouring = "By Lineage"
   [
    ifelse high-light-mutated-cells
     [
      ifelse any? patches with [CellType !=  0]
       [
        ifelse CellType = 1 or CellType = 2
         [set pcolor palette:scale-gradient palette:scheme-colors "Divergent" "Spectral" 10 lineageID 0 #patches]
         [set pcolor white]
       ]
       [set pcolor palette:scale-gradient palette:scheme-colors "Divergent" "Spectral" 10 lineageID 0 #patches]
     ]
    [set pcolor palette:scale-gradient palette:scheme-colors "Divergent" "Spectral" 10 lineageID 0 #patches]
  ]

  if PatchColouring = "By Strength"
   [
    ifelse any? patches with [CellType = 1 or CellType = 2]
     [
      ifelse high-light-mutated-cells = TRUE
        [set pcolor white]
        [set pcolor scale-color 7 strength 1 0.45]
      if CellType = 1 [set pcolor scale-color orange strength 1 0.35]
      if CellType = 2 [set pcolor scale-color violet strength 1 0.35]
     ]
     [set pcolor scale-color 7 strength 1 0.5]
   ]
  label-patches
end

to label-patches
  if Labelling = "By Lineage"
    [set plabel lineageID]
  if Labelling = "By Strength"
    [ifelse strength > 0.85 [set plabel-color yellow][set plabel-color black]
     set plabel precision strength 2]
  if Labelling = "By Cell Type"
    [
     ifelse CellType = 0
      [set plabel ""]
      [set plabel celltype]
     ]
end
@#$#@#$#@
GRAPHICS-WINDOW
168
12
647
492
-1
-1
11.5
1
6
1
1
1
0
1
1
1
-20
20
-20
20
0
0
1
ticks
10.0

BUTTON
13
14
76
47
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
11
119
160
152
WT_mean_strength
WT_mean_strength
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
11
160
159
193
WT_sd_strength
WT_sd_strength
0
10
0.1
0.01
1
NIL
HORIZONTAL

BUTTON
81
15
144
48
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
81
48
144
81
go-once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
994
68
1165
248
Stratification All Cells
Time
% Stratification Events
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"all cells" 1.0 0 -16777216 true "" "plot 100 * (#Displaced / #patches); ((max-pxcor * 2) + 1) ^ 2)"

MONITOR
673
14
760
59
Number of Cells
count patches
17
1
11

MONITOR
994
14
1160
59
Percentage of stratifiel cells
100 * (#Displaced / #patches)
2
1
11

TEXTBOX
14
81
164
99
STRENGTH
14
0.0
1

PLOT
991
358
1174
506
Number of Clones (Lineages)
tick
#lineages
0.0
10.0
0.0
441.0
true
true
"" ""
PENS
"M1" 1.0 0 -955883 true "" "if ticks > 0 [plot (length (unique-M1lineages))]"
"M2" 1.0 0 -8630108 true "" "if ticks > 0 [plot (length (unique-M2lineages))]"
"WT" 1.0 0 -7500403 true "" "if ticks > 0 [plot (length (unique-WTlineages))]"

PLOT
661
357
987
507
CloneSizes for Each Lineage
LineageID
CloneSize
0.0
441.0
0.0
441.0
false
true
"" ""
PENS
"WT" 1.0 1 -11053225 true "set-histogram-num-bars 500\nset-plot-x-range 0 #patches\nset-plot-y-range 0 #patches" "if ticks > 0 [histogram wtlineage-list]"
"M1" 1.0 1 -955883 true "set-histogram-num-bars 500\nset-plot-x-range 0 #patches\nset-plot-y-range 0 #patches" "if ticks > 0 [histogram m1lineage-list]"
"M2" 1.0 1 -8630108 true "set-histogram-num-bars 500\nset-plot-x-range 0 #patches\nset-plot-y-range 0 #patches" "if ticks > 0 [histogram M2lineage-list]"

MONITOR
660
307
734
352
# Lineages
;length unique-lineages\n#lineages
17
1
11

PLOT
663
514
989
686
Time series Mean Clone Size
tick
mean clonesize
0.0
10.0
0.0
25.0
true
true
"" ""
PENS
"M1" 1.0 0 -955883 true "" "if ticks > 0 and any? patches with [cellType = 1]\n[plot mean m1clonesize-list]"
"M2" 1.0 0 -8630108 true "" "if ticks > 0 and any? patches with [cellType = 2]\n[plot mean M2clonesize-list]"
"WT" 1.0 0 -9276814 true "" "if ticks > 0 and any? patches with [cellType = 0]\n[plot mean WTclonesize-list]"

CHOOSER
172
511
313
556
PatchColouring
PatchColouring
"By Lineage" "By Strength"
1

CHOOSER
318
511
459
556
Labelling
Labelling
"By Lineage" "By Strength" "No Label" "By Cell Type" "By Area"
2

SLIDER
12
268
155
301
M1_Mean_Strength
M1_Mean_Strength
0
10
0.52
0.01
1
NIL
HORIZONTAL

SLIDER
14
308
155
341
M1_sd_strength
M1_sd_strength
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
13
229
156
262
#M1
#M1
0
1000
100.0
1
1
NIL
HORIZONTAL

MONITOR
831
13
891
58
#M1 cells
count patches with [cellType = 1]
17
1
11

SWITCH
463
517
646
550
high-light-mutated-cells
high-light-mutated-cells
1
1
-1000

PLOT
660
67
983
249
Strength Distributions
strength
count
0.0
1.5
0.0
10.0
true
true
"" ""
PENS
"M1" 1.0 1 -955883 true "set-histogram-num-bars 30" "histogram m1strength-list"
"WT" 1.0 1 -11053225 true "set-histogram-num-bars 30" "histogram wtstrength-list"
"M2" 1.0 1 -8630108 true "set-histogram-num-bars 30" "histogram M2strength-list"

MONITOR
818
307
898
352
# M1lineages
#M1lineages
17
1
11

OUTPUT
1194
58
1478
665
11

BUTTON
1488
97
1588
130
NIL
clear-output\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1194
18
1486
60
 MEAN MUTANT CLONESIZE\n  WT            N(m1)         M1         N(m2)           M2\n
11
0.0
1

SLIDER
11
407
153
440
M2_Mean_Strength
M2_Mean_Strength
0
2
0.51
0.01
1
NIL
HORIZONTAL

SLIDER
12
444
154
477
M2_sd_Strength
M2_sd_Strength
0
2
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
11
371
153
404
#M2
#M2
0
1000
100.0
1
1
NIL
HORIZONTAL

MONITOR
900
307
984
352
NIL
#M2lineages
17
1
11

TEXTBOX
13
104
163
122
WildType\n
11
0.0
1

TEXTBOX
13
209
163
227
Mutant 1\n
11
25.0
1

TEXTBOX
12
355
162
373
Mutant 2
11
115.0
1

MONITOR
735
307
825
352
#WTLineages
min list #WTLineages (count patches with [celltype = 0])
17
1
11

MONITOR
763
13
828
58
# WT cells
count patches with [cellType = 0]
17
1
11

SWITCH
1488
59
1588
92
list-output?
list-output?
0
1
-1000

CHOOSER
10
512
161
557
Stopping-Condition
Stopping-Condition
"Dominant Lineage" "Dominant CellType" "No Wildtype" "STOP-time"
1

SLIDER
9
559
151
592
STOP-Time
STOP-Time
0
100
100.0
1
1
NIL
HORIZONTAL

PLOT
993
515
1175
662
LogPlot #Clones
ticks
ln #clones
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"M1" 1.0 0 -955883 true "" "if ticks > 0 [plotxy  ticks (ln length (unique-M1lineages))]"
"M2" 1.0 0 -8630108 true "" "if ticks > 0 [plotxy ticks (ln length (unique-M2lineages))]"
"WT" 1.0 0 -7500403 true "" "if ticks > 0 [plotxy  ticks (ln length (unique-WTlineages))]"

MONITOR
896
14
961
59
#M2 cells
count patches with [celltype = 2]
17
1
11

MONITOR
665
690
722
735
NIL
runtime
17
1
11

MONITOR
661
255
763
300
mean strength WT
mean WTstrength-list
2
1
11

MONITOR
767
255
865
300
mean strength M1
mean M1strength-list
2
1
11

MONITOR
871
256
984
301
mean strength M2
mean M2strength-list
2
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experimentIrinaModel-ABM-Patches2StopConditionOFF" repetitions="500" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>CloneSize</metric>
    <metric>#DisplacedTracedCells</metric>
    <metric>#Displaced</metric>
  </experiment>
  <experiment name="AllStrengths&amp;Divisions10by10" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>pstrength-list</metric>
    <metric>pdivision-list</metric>
  </experiment>
  <experiment name="virulance" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <enumeratedValueSet variable="M1_sd_strength">
      <value value="0.11"/>
      <value value="0.12"/>
      <value value="0.13"/>
      <value value="0.14"/>
      <value value="0.15"/>
      <value value="0.16"/>
      <value value="0.17"/>
      <value value="0.18"/>
      <value value="0.19"/>
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="M1_Mean_Strength">
      <value value="0.41"/>
      <value value="0.42"/>
      <value value="0.43"/>
      <value value="0.44"/>
      <value value="0.45"/>
      <value value="0.46"/>
      <value value="0.47"/>
      <value value="0.48"/>
      <value value="0.49"/>
      <value value="0.5"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
