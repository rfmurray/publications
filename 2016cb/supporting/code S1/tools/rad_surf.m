function rad_surf( name, x, y, z, reflectance )

% RAD_SURF  Create a RADIANCE surface from MATLAB coordinate matrices
%
%     usage:  rad_surf( name, x, y, z, reflectance )
% 
%     name        -- the RADIANCE name that will be assigned to the surface
%     x, y, z     -- coordinate matrices that define the surface
%     reflectance -- the reflectance of the surface

% open .rad file
fid = fopen([ name '.rad' ],'w');

% write dot texture
fprintf(fid, 'void brightfunc dots\n');
fprintf(fid, '2 dots dots.cal\n');
fprintf(fid, '0\n');
fprintf(fid, '3 0.003000 0.180000 0.300000\n\n');

% write plastic material with dot texture
fprintf(fid, 'dots plastic greytextured\n');
fprintf(fid, '0 0 5 %.6f %.6f %.6f 0 0\n\n',reflectance,reflectance,reflectance);

% write patches
[m,n] = size(x);
for i = m:-1:2
    for j = 1:(n-1)
        
        % make tag for this patch
        tag = sprintf('%d_%d',m-i+1,j);
        
        % write triangles
        tri1name = sprintf('%s_tri1_%s',name,tag);
        tri2name = sprintf('%s_tri2_%s',name,tag);
        rad_polygon(fid,tri1name,'greytextured',[ [ x(i,j)     y(i,j)     z(i,j) ]'     [ x(i,j+1) y(i,j+1) z(i,j+1) ]' [ x(i-1,j) y(i-1,j) z(i-1,j) ]' ]);
        rad_polygon(fid,tri2name,'greytextured',[ [ x(i-1,j+1) y(i-1,j+1) z(i-1,j+1) ]' [ x(i,j+1) y(i,j+1) z(i,j+1) ]' [ x(i-1,j) y(i-1,j) z(i-1,j) ]' ]);
        
    end
end

% close file
fclose(fid);

end


function rad_polygon( fid, name, mod, vlist )

% RAD_POLYGON  Create a RADIANCE polygon
%
% rad_polygon( fid, name, mod, vlist )

% write polygon
fprintf(fid,'%s polygon %s\n',mod,name);
fprintf(fid,'0 0 %3d %12.6f %12.6f %12.6f\n',numel(vlist),vlist(:,1));
if size(vlist,2)>1
    fprintf(fid,'        %12.6f %12.6f %12.6f\n',vlist(:,2:end));
end
fprintf(fid,'\n');

end
