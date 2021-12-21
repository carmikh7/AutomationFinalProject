# AutomationFinalProject
 
 This code does double layer mask assignment for a user defined layout.

 The user may define the layout in the "layout.txt" file that is in this folder.
 Areas where there is a figure should be marked with X; empty areas should be marked with _ .

 When the program is run, the user is propted to pick a new mininum distance between features or to keep the default distance of 2.80, or two blocks.

 A conflict graph representation is displayed before the node splitting. If the node splitting is successful, a conflict graph will be displayed at the end of the process. A mask assignment layout will be outputted in the terminal similar to the text layout defined by the user, but the mask assignments will be displayed with 0 or 1.

 The program goes through the following functions:

 sectionShapes
 Identifies the different blocks in the design and stores their coordinates in vectors

 iterCol and iterRow
 Are used to iterate through columns recursively when identifying blocks in sectionShapes

 findEdges
 Measure the smallest distance between blocks and identifies which one have edges that are below the minimum distance

 findConflict
 Conflict cycle detector, which finds loops with an odd number of vertices

 findNodesToSplit
 Figures out which nodes should be split to resolve conflicts using the shoftest distance between blocks

 nodeSplitting
 splits nodes at location identified by findNodesToSplit

 plotConflicts
 plots the conflict graph using centers from findCenters

 findCenters
 finds the center points of each block to be used in the graph construction in plotConflicts

 assignMask
 if all conflicts are resolved, assigns each block a mask

 iterateLoop
 iterative loop used when splitting nodes and checking for new conflicts

