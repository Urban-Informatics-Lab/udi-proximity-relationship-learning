## udi-proximity-relationship-learning
Code for the paper: [Urban Data Integration Using Proximity Relationship Learning for Design, Management, and Operations of Sustainable Urban Systems.](https://doi.org/10.1061/(ASCE)CP.1943-5487.0000806)

## supporting_functions.py

This file contains the functions that are called by the relationship_learning.py. Most of the functions in this file are supporting the geometric and relationship learning algorithms written in the main code. This file shold be present in the same directory as 

## relationship_learning.py
This file contains code for preprosessing the data and learning relationships. First the data from the element files( .xlsx) is loaded into memory and then the relationship learning algorithms are executed in the following order:

-between polygon and linear elements <br>
-between polygon and polygon elements<br>
-between polygon and point elements<br>

-between point and linear elements
-between point and point elements
-between linear and linear elements

## min_bounding_rect.py 
This file contains the code to find the minimum-area bounding box of a set of 2D points. This is a third-party supporting file which was not developed by the authors. Please see the file for more information.

## qhull_2d.py 
This file contains the code to compute the convex hull of a set of 2D points. This is a third-party supporting file which was not developed by the authors. Please see the file for more information.

## Data Files

This folder contains the data files on which we tests our framework. The folder contains the following 4 sub-folders:

-Polygon --> contains all polygon element files
-Linear --> contains all linear element files
-Point --> contains all point element files
-Block --> contains one file each (containing segments of corresponding linear elements) for all linear element files
 
The 4 sub-folders inside this folder should be stored in the same directory as relationship_learning.py

The program can only read .xlsx files right now.
