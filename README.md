# A minimal example to get the chaste serialization/deserialization crash due to NaN values

The app `MinimalSample` in this project randomly generates a periodic mesh with 3 by 3 cells by means of a Voronoi Tessellation. 
It then proceeds to calculate the elongation shape factor and serializes the simulation. 
On deserialization, the app crashes for my version of chaste and on my machine
