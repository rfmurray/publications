function ind = mesh2ind(nRows,nCols)

% mesh2ind  Determine which vertices form faces in the mesh
% 
%     usage:  mesh2ind( nRows, nCols )
% 
%     input arguments
%         nRows -- number of rows in mesh
%         nCols -- number of columns in mesh
%
%     output variables
%         ind   -- a matrix that shows which vertices in the mesh combine
%                  to form faces

m = nRows;
n = nCols;

% Get the indices for making faces out of a vertex mesh.
ind1 = reshape( 1:(m*n), [ m n ] );
ind1 = ind1(1:end-1,1:end-1);
ind2 = ind1+1;
ind3 = ind1+m;
ind4 = ind1+m+1;
ind123 = [ ind1(:)' ; ind2(:)' ; ind3(:)' ];
ind234 = [ ind2(:)' ; ind4(:)' ; ind3(:)' ];
ind = [ ind123(:) ; ind234(:) ];

ind = uint32(ind-1);

end