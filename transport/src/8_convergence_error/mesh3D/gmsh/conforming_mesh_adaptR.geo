////////////////////////////////////////////////////////////////////
//Geometrical parameters

// Define the coordinate range of the mesh


DefineConstant
[
x0 = {-0.25, Name "xmin tissue"},
x1 = {+0.25, Name "xmax tissue"},
y0 = {-1.0, Name "ymin tissue"},
y1 = {+1.0,Name  "ymax tissue"},
z0 = {-1.0, Name "zmin tissue"},
z1 = {+1.0,Name  "zmax tissue"},
yc= {0.00, Name "y - center of vessel"},
zc= {0.00, Name "z - center of vessel"},
R = {0.25, Name "Radius of vessel"},
Num_el= {16.0, Name "Number of elements for unit size"}
adaptation={0.2, Name "Ratio of size element of vessel wall/ size element of tissue"}
];

// mesh element size at each point
cellsize=1./Num_el;



////////////////////////////////////////////////////////////////
//Define the 3d tissue slab
p=newp;

Point(p) = {x0, y0, z0, cellsize};
Point(p+1) = {x0, y1, z0, cellsize};
Point(p+2) = {x0, y1, z1, cellsize};
Point(p+3) = {x0, y0, z1, cellsize};

l= newl;
Line(l+0) = {p+0,p+1};
Line(l+1) = {p+1,p+2};
Line(l+2) = {p+2,p+3};
Line(l+3) = {p+3,p+0};

ll1=newll;
Line Loop(ll1) = {l+0, l+1,l+2,l+3};

//////////////////////////////////////////////////////////////////
// Define the vessel

p=newp;
Point(p)  ={x0,yc,zc, cellsize*adaptation};
Point(p+1)={x0,yc,zc+R, cellsize*adaptation};
Point(p+2)={x0,yc+R,zc, cellsize*adaptation};
Point(p+3)={x0,yc,zc-R, cellsize*adaptation};
Point(p+4)={x0,yc-R,zc, cellsize*adaptation};

c=newl;
Circle(c) = {p+1,p,p+2};
Circle(c+1) = {p+2,p,p+3};
Circle(c+2) = {p+3,p,p+4};
Circle(c+3) = {p+4,p,p+1};

ll2 = newll;
Line Loop(ll2) = {c,c+1,c+2,c+3};

////////////////////////////////////////////////////////////////
//build volume

// Define the surfaces
s = news;
Plane Surface(s+0) = {ll1,ll2};
Plane Surface (s+1)= {ll2};


// Extrude
tmp1[]=Extrude {x1-x0, 0, 0} { Surface{s+0}; };
tmp2[]=Extrude {x1-x0, 0, 0} { Surface{s+1}; };

//Define physical volume
v=newv;
Physical Volume(v) = {1,2}; 


