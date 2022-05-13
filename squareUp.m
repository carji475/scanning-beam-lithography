function mrSquare = squareUp(vector,indices,Nx,Ny)
    mrSquare          = zeros(Nx,Ny);
    mrSquare(indices) = vector;
end