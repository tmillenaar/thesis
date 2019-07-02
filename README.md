# Imperfect Sorting Model of Linear Diffusion #

## Introduction ##

This model simulates the transport and deposition of sediment in an alluvial fan. 
The sediment is divided into two grain sizes, each of which is independently transported through linear diffusion. 
The diffusion equation used is based on the <a href="https://en.wikipedia.org/wiki/FTCS_scheme">FTCS</a> finite difference method.
The two grain sizes are allowed to occupy the same location in the model. Thich is referred to as *imperfect sorting* as opposed to *perfect sorting*, which is a method where each grain size is assigned its own seperate regime.

### Boundary conditions ###
The location of sediment is based on a mesh. This mesh is build up of several columns, each of which contains vertically stacked nodes.
The boundary at the left end of the mesh is the most proximal side of the model and starts with column 0. 
This side is assigned a sediment input (q0). The right end of the mesh represents the distal end of the fan/basin. 
This end of the model is set so that it cannot contain any sediment, so its height is always 0. 
Any sediment that reaches the distal end leaves the model and is considered to be the sediment output of the fan.
The basin can be subject to subsidence, which is treated as a beam rotating about a pivot point, located at the distal end of the model.
<br>
The diffusivity constant (i.e. sediment transport rate) and sediment input can be set per grain size. 
They can be defined to vary through time, cyclic or otherwise, as can the subsidence rate. 
