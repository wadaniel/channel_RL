fin = './case_folder/channel8000_AMD_RL.init';
Nx = 130; Ny = 41; Nz = 66;


fid = fopen(fin,'r','b');

n   = fread(fid,1,'int32','ieee-le');
x   = fread(fid,n,'float64','ieee-le');
assert(n == Nx)

n   = fread(fid,1,'int32','ieee-le');
y   = fread(fid,n,'float64','ieee-le');
assert(n == Ny)

n   = fread(fid,1,'int32','ieee-le');
z   = fread(fid,n,'float64','ieee-le');
assert(n == Nz)

n   = fread(fid,1,'int32','ieee-le');
xm  = fread(fid,n,'float64','ieee-le');
assert(n == Nx-1)

n   = fread(fid,1,'int32','ieee-le');
ym  = fread(fid,n,'float64','ieee-le');
assert(n == Ny-1)

n   = fread(fid,1,'int32','ieee-le');
zm  = fread(fid,n,'float64','ieee-le');
assert(n == Nz-1)

n   = fread(fid,3,'int32','ieee-le');
U   = fread(fid,n(1)*n(2)*n(3),'float64','ieee-le');
U   = reshape(U,n(1),n(2),n(3));
assert((n(1) == Nx) & (n(2) == Ny+1) & (n(3) == Nz+1))

n   = fread(fid,3,'int32','ieee-le');
V   = fread(fid,n(1)*n(2)*n(3),'float64','ieee-le');
V   = reshape(V,n(1),n(2),n(3));
assert((n(1) == Nx+1) & (n(2) == Ny) & (n(3) == Nz+1))

n   = fread(fid,3,'int32','ieee-le');
W   = fread(fid,n(1)*n(2)*n(3),'float64','ieee-le');
W   = reshape(W,n(1),n(2),n(3));
assert((n(1) == Nx+1) & (n(2) == Ny+1) & (n(3) == Nz))

fclose(fid);


disp('everything works!')
