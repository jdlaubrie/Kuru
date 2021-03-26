// Gmsh project created on Wed Feb 13 16:44:14 2019
// Points
p1=newp; Point(p1) = {10., 0, 0, 1.0};
p2=newp; Point(p2) = {11.41, 0, 0, 1.0};
// Lines
l1=newl; Line(l1) = {p1, p2};
// Surfaces
thick[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {Line{l1};Layers{15};Recombine;};
// Volumes
wall1[] = Extrude {0, 0, 45.} {Surface{thick[1]};Layers{30};Recombine;};
wall2[] = Extrude {0, 0, 45.} {Surface{wall1[0]};Layers{30};Recombine;};

Transfinite Line {l1} = 3;

// Symmetry X
Physical Surface(91) = {22,44};
// Symmetry Y
Physical Surface(92) = {14,36};
// Top Z
Physical Surface(93) = {49};
// Bottom Z
Physical Surface(94) = {5};
// Inner wall
Physical Surface(95) = {26,48};
// Outer wall
Physical Surface(96) = {18,40};
// Volume of Cylinder
Physical Volume(97) = {1,2};

// Generate 3D mesh
Mesh 3;
// Remove duplicate entities
Coherence Mesh;

// Save the mesh
Mesh.MshFileVersion = 2.2;
Save "Cylinder.msh";
