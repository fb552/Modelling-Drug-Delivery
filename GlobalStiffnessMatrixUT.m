%% Test 1: test symmetry of the matrix
% % Test that this matrix is symmetric
tol = 1e-14;
Ne=4;
parameters.selection = '1';
parameters.D = 1;
parameters.lambda = 1;
parameters.f = 0;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,Ne,order,parameters);

%global matrix and transpose with values above
GM = GlobalStiffnessMatrix(Ne,msh,GQ,order);
transposedGM = transpose(GM);

diff = GM - transposedGM;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 2: test that the Global Matrix is evaluated correctly
tol = 1e-14;
Ne=4;
parameters.selection = '1';
parameters.D = 5;
parameters.lambda = 0;
parameters.f = 0;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,Ne,order,parameters);

%global matrix with values above
GM = GlobalStiffnessMatrix(Ne,msh,GQ,order);

%analytical global matrix 
analyticalGM = [20,-20,0,0,0;-20,40,-20,0,0;0,-20,40,-20,0;0,0,-20,40,-20;0,0,0,-20,20];

diff = GM - analyticalGM;
diffnorm = sum(sum(diff.*diff)); 
assert(abs(diffnorm) <= tol)

%% Test 3: test Global Matrix Size
% Test that the size is equal to that of the nodes in the mesh
Ne=4;
parameters.selection = '1';
parameters.D = 5;
parameters.lambda = 5;
parameters.f = 0;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,Ne,order,parameters);

%global matrix with values above and size
GM = GlobalStiffnessMatrix(Ne,msh,GQ,order);
sizeGM = size(GM);

assert(sizeGM(1)==(msh.ngn)); %columns are same as node number
assert(sizeGM(2)==(msh.ngn)); %rows are same as node number

%% Test 4: test that values are only in the diagonals 
% Tests that all values are only in the diagonals, and the rest are zeros.
tol = 1e-14;
Ne=4;
parameters.selection = '1';
parameters.D = 1;
parameters.lambda = 1;
parameters.f = 0;
GQ.switch = '0';
GQ.npts = 2;
order = 1;

msh = OneDimLinearMeshGen(0,1,Ne,order,parameters);

%global matrix with values above
GM = GlobalStiffnessMatrix(Ne,msh,GQ,order);
%zeros same size as GM
zeroMatrix = zeros(msh.ngn);

mainGMDiagonal = diag(diag(GM));            %main diagonal
upperGMDiagonal = diag(diag(GM, 1), 1);     %upper diagonal
lowerGMDiagonal = diag(diag(GM, -1), -1);   %lower diagonal

%remove above diagonals from GM
emptyGM = GM - mainGMDiagonal - lowerGMDiagonal - upperGMDiagonal;

diff = emptyGM - zeroMatrix;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)