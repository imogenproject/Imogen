function saveObjMesh(name,varargin)

%=== SAVEOBJMESH ===================================================================================
% This function saves a 3D mesh as a Wavefront/Alias Obj file
%===================================================================================================
%
%-CALLS---------------------------------------------------------------------------------------------
%   saveObjMesh(name,x,y,z,nx,ny,nz)
%   saveObjMesh(name,facesAndVerts)
%-INPUT---------------------------------------------------------------------------------------------
% name      name of the obj file to save                                        str     *
% x,y,z     equally sized matrices with coordinates                             double  [nx ny nz]
% nx,ny,nz  are normal directions (optional)                                    double  [nx ny nz]
% facesAndVerts     structure with faces, vertices fields                       struct  matlab
%---------------------------------------------------------------------------------------------------
    
    %--- Prepare file name ---%
    if isempty(name) 
        answer = inputdlg('Enter file name for saving','Save isosurface as...');
        name = answer{1};
    end
    if isempty(strfind(name, '.obj')), name = [name '.obj']; end
    
    dataType = 'struct';
    normals = false;
    if (length(varargin) > 2 && isa(varargin{3}, 'double'))
        dataType = 'arrays';
        x = varargin{1}; y = varargin{2}; z = varargin{3};
        if (length(varargin) == 6)
            normals = true;
            nx = varargin{4}; ny = varargin{5}; nz = varargin{6};
        end
    elseif (length(varargin) == 1 && isa(varargin{1}, 'struct'))
       facesAndVerts = varargin{1}; 
    elseif (length(varargin) == 2)
        facesAndVerts = isosurface(varargin{1},varargin{2}(1));
    else
        error('ImogenUtil:SaveObjMeshArgsError','Unable to parse arugments to saveObjMesh.');
    end

    fid=fopen(name,'w');
    
    %------------------------------------------------------------
    if strcmpi(dataType,'arrays') % ARRAY CASE
    %------------------------------------------------------------
        l = size(x,1); h = size(x,2);  
        if ((l==1) || (h==1)) 
            [x,y] = meshgrid(x,y); 
            l=size(x,1); h=size(x,2);  
        end

        n=zeros(l,h);
        nn=1;
        for i=1:l
            for j=1:h
              n(i,j)=nn; 
              fprintf(fid, 'v %f %f %f\n',x(i,j),y(i,j),z(i,j)); 
              fprintf(fid, 'vt %f %f\n',(i-1)/(l-1),(j-1)/(h-1)); 
              if (normals), fprintf(fid, 'vn %f %f %f\n', nx(i,j),ny(i,j),nz(i,j)); end
              nn=nn+1;
            end
        end
        fprintf(fid,'g mesh\n');

        for i=1:(l-1)
            for j=1:(h-1)
                  if (normals) 
                    fprintf(fid,'f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n',n(i,j),n(i,j),n(i,j),n(i+1,j),n(i+1,j),n(i+1,j),n(i+1,j+1),n(i+1,j+1),n(i+1,j+1),n(i,j+1),n(i,j+1),n(i,j+1));
                  else
                    fprintf(fid,'f %d/%d %d/%d %d/%d %d/%d\n',n(i,j),n(i,j),n(i+1,j),n(i+1,j),n(i+1,j+1),n(i+1,j+1),n(i,j+1),n(i,j+1));
                  end
            end
        end
        fprintf(fid,'g\n\n');
    %------------------------------------------------
    else % STRUCTURE CASE
    %------------------------------------------------
        for i=1:size(facesAndVerts.vertices,1)
            v = facesAndVerts.vertices(i,:);
            fprintf(fid, 'v %f %f %f\n', v(1), v(2), v(3));
        end
        fprintf(fid,['g ', strrep(strrep(name,'.obj','_geo'),' ','_'), '\n']);
        for i=1:size(facesAndVerts.faces,1)
            f = facesAndVerts.faces(i,:);
            fprintf(fid, 'f %d %d %d\n', f(1), f(2), f(3));
        end
        fprintf(fid,'g\n\n');
    end
   
    fclose(fid);
end
  