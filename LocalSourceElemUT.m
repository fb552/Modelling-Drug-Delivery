%% Test 1: test symmetry of the vector
% % Test that this vector is symmetric using manual integration
tol = 1e-14;
eN=1;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 1;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,10,order,parameters);

elemat = LocalSourceElem(eN,msh,GQ,order);

assert(abs(elemat(1) - elemat(2)) <= tol)

%% Test 2: test 2 different elements of the same size produce same matrix
% % Test that for two elements of an equispaced mesh, the element matrices 
% % calculated are the same using manual integration
tol = 1e-14;
eN=1;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 1;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,10,order,parameters);

elemat1 = LocalSourceElem(eN,msh,GQ,order);

eN=2;

elemat2 = LocalSourceElem(eN,msh,GQ,order);

diff = elemat1 - elemat2;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 3: test that one vector is evaluted correctly
% % Test that the element vector is evaluated correctly using manual integration
tol = 1e-14;
eN=1;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 1;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,4,order,parameters);

elemat1 = LocalSourceElem(eN,msh,GQ,order);

elemat2 = [ 0.1250; 0.1250 ];
diff = elemat1 - elemat2;  %calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); %calculates the total squared error between the matrices
assert(abs(diffnorm) <= tol)

%% Test 4: test symmetry of the vector
% % Test that this vector is symmetric using Gaussian Quadrature
tol = 1e-14;
eN=1;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 1;
GQ.switch = '1';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,10,order,parameters);

elemat = LocalSourceElem(eN,msh,GQ,order);

assert(abs(elemat(1) - elemat(2)) <= tol)

%% Test 5: test 2 different elements of the same size produce same matrix
% % Test that for two elements of an equispaced mesh, the element matrices 
% % calculated are the same using Gaussian Quadrature
tol = 1e-14;
eN=1;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 1;
GQ.switch = '1';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,10,order,parameters);

elemat1 = LocalSourceElem(eN,msh,GQ,order);

eN=2;

elemat2 = LocalSourceElem(eN,msh,GQ,order);

diff = elemat1 - elemat2;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 6: test that one vector is evaluted correctly
% % Test that the element vector is evaluated correctly using Gaussian Quadrature
tol = 1e-14;
eN=1;
parameters.selection = '1';
parameters.D = 0;
parameters.lambda = 0;
parameters.f = 1;
GQ.switch = '1';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,4,order,parameters);

elemat1 = LocalSourceElem(eN,msh,GQ,order);

elemat2 = [ 0.1250; 0.1250 ];
diff = elemat1 - elemat2;  %calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); %calculates the total squared error between the matrices
assert(abs(diffnorm) <= tol)
