%% Test 1: test symmetry of the vector 
% % Test that this vector is symmetric
tol = 1e-14;
Ne=4;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 1;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,Ne,order,parameters);

%global vector and flipped with values above
GV = GlobalSourceVector(Ne,msh,GQ,order);
flippedGV = flip(GV);

diff = GV - flippedGV;
diffnorm = sum(sum(diff.*diff)); 
assert(abs(diffnorm) <= tol) 

%% Test 2: test that the Global Vector is evaluated correctly.
tol = 1e-14;
Ne=4;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 5;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,Ne,order,parameters);

%global vector with values above
GV = GlobalSourceVector(Ne,msh,GQ,order);

%analytical global vector 
analyticalGV = [0.625; 1.25; 1.25; 1.25; 0.625];

diff = GV - analyticalGV;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol) 

%% Test 3: test Global Vector Size
% Test that the size is equal to that of the nodes in the mesh
Ne=4;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 5;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,Ne,order,parameters);

%global vector with values above and size
GV = GlobalSourceVector(Ne,msh,GQ,order);
GVsize = size(GV);

assert(GVsize(1)==(msh.ngn)); %columns are same as node number
assert(GVsize(2)==1);         %rows are 1