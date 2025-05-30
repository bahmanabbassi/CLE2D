function yashow_spheric(mat, varargin)

% yashow_spheric - Display a matrix onto a sphere
%
% Syntax:
%    yashow_spheric(mat, ['cmap', cmap] ...
%                   [,'faces', faces] ...
%                   [,'mollweide'] ...
%                   [,'surf'])
%
% Description:
%    This function displays a matrix mapped onto a sphere, with options for
%    various display modes including surface plotting and Mollweide projection.
%
% Inputs:
%    mat - MATRIX: the matrix to be displayed.
%    'cmap' - STRING (optional): colormap to use for display.
%    'faces' - INTEGER (optional): number of faces for the sphere.
%    'mollweide' - BOOLEAN (optional): display Mollweide projection.
%    'surf' - BOOLEAN (optional): display surface plot.
%
% Example:
%    % Example usage:
%    yashow_spheric(rand(10), 'cmap', 'jet', 'surf');
%
% See Also:
%    yashow
%
% This file is part of YAW Toolbox (Yet Another Wavelet Toolbox).

%% Two cases: 'surf' display or mapping
%% First, the 'surf' case:

fig = gcf;

%% Set the mode of mat
if isreal(mat)
  defmode = 'real';
else
  defmode = 'abs';
end
[mode,varargin] = getopts(varargin,'mode',defmode);

mode = lower(mode);
modelist = {'abs', 'angle', 'real', 'imag'};
if all(~strcmp(modelist, mode)) 				  
  error(['The mode ''' mode ''' is undefined in yashow']);
end

nmat = feval(mode, double(mat));

if getopts(varargin,'surf',[],1)
  
  if (any(nmat(:) < 0))
    error(['The ''surf'' display canot be applied on matrix with negative ' ...
           'values.']);
  end
  
  if (size(nmat,3) ~= 1)
    error('No ''surf'' display for volume of data');
  end
  
  %% This flipping is required by the mapping
  %nmat = fliplr(nmat);
  [nth, nph] = size(nmat);
  
  %% The spherical grid
  if (rem(nth,2) == 1)
    [phi, theta] = sphgrid(nth - 1, nph, 'withpoles');
  else
    [phi, theta] = sphgrid(nth, nph);
    
    phi = phi([1 1:end end], :);
    
    theta = [zeros(1,nph); ...
	     theta; ...
	     ones(1,nph)*pi];
    
    nmat = nmat([1 1:end end],:,:);
    
    
    nmat(1,:) = mean(nmat(1,:));
    nmat(end,:) = mean(nmat(end,:));

  end
  
  phi = phi(:, [end 1:end]);
  theta = theta(:, [end 1:end]);
  nmat = nmat(:, [end 1:end]);
  
  %% Determing the 3D data.
  Radius       = nmat;

  X            = Radius .* sin(theta) .* cos(phi); 
  Y            = Radius .* sin(theta) .* sin(phi);
  Z            = Radius .* cos(theta);
  
  %% Display the corresponding surface
  surf(X,Y,Z);
  
  %% Mapping data on this surface
  h = findobj(fig,'Type','surface');
  set(h,'CData',nmat,'FaceColor','texturemap');
 
  axis('equal')
  axis('tight');
  shading flat;
  
  %% Adding some special effects 
  %camlight
  %lighting phong
  view(140,30);
  xlabel('X');
  ylabel('Y');
  zlabel('Z');

%% Second, the mollweide projection
elseif getopts(varargin,'mollweide',[],1) 
  
  [nth,nph] = size(nmat);

  %% The spherical grid
  if (rem(nth,2) == 1)
    [phi, theta] = sphgrid(nth - 1, nph, 'withpoles');
  else
    [phi, theta] = sphgrid(nth, nph);
  end

  X = (2^1.5/pi) .* (phi - pi) .* cos(pi/2 - theta);
  Y =  2^.5 .* sin(pi/2 - theta);
  pcolor(X,Y,nmat(:,[(end-nph/2+1):end 1:(nph/2)]));
  shading flat;
  axis equal tight;

elseif getopts(varargin,'contour',[],1)   
  
  [X,Y,Z] = sphere(40);
  surf(X,Y,Z,Z*0);
  shading flat;
  colormap(0.9*[1 1 1]);
  axis equal;
  
  hold on;
  
  [nth,nph] = size(nmat);

  %% The spherical grid
  if (rem(nth,2) == 1)
    phi = vect(0, 2*pi, nph, 'open');
    theta = vect(0, pi, nth);
  else
    phi = vect(0, 2*pi, nph, 'open');
    theta = vect(0, pi, nth, 'rlopen');
  end
  
  nb_contour = getopts(varargin,'contour',5);
  
  if (~isinteger(nb_contour))
    nb_contour = 5;
  end
  
  c = contourc(phi, theta, nmat, nb_contour);
  
  pos = 1;
  nbpos = size(c,2);
  r = 1.01;

  while (pos <= nbpos)
    sub_lg = c(2,pos);
    sub_val = c(1,pos);
    
    sub_c = c(:,(pos+1):(pos+sub_lg));
    
    th = sub_c(2,:);
    ph = sub_c(1,:);
    
    xc = r*sin(th).*cos(ph);
    yc = r*sin(th).*sin(ph);
    zc = r*cos(th);
    
    plot3(xc, yc, zc);
    
    pos = pos + sub_lg + 1;
  end
  
  hold off;
  
  return;
else %% Third, The mapping case:
    
  [nth, nph, ncolor] = size(nmat);
  
  nmat = nmat(end:-1:1,:,:);
  nph_2 = ceil(nph/2);
  nmat = nmat(:, [(nph_2+1):end 1:nph_2], :);
  
  %% Creating the sphere 
  [N,varargin] = getopts(varargin,'faces',40);
  sphere(N);

  %% Catching the surface object
  h = findobj(fig,'Type','surface');

  %% Mapping the result
  if (ncolor > 1)
    set(h,'CData',mat(end:-1:1,fftshift(1:end),:),'FaceColor','texturemap');
  else
    set(h,'CData',nmat,'FaceColor','texturemap');
  end    

  %% Miscellaneous modifications
  axis('equal')
  axis('tight');
  shading flat;
  view(140,30);
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
end



%% ++++++++ Another useful options +++++++++

%% Determining the colormap 
[cmap,varargin] = getopts(varargin, 'cmap', []);

if isempty(cmap)
  if islogical(mat)
    cmap = 'winter';
  else
    cmap = 'jet';
  end
  colormap(cmap);
elseif ~isnumeric(cmap)
  pos_par = yastrfind(cmap,'(');
  if isempty(pos_par)
    cmap_base = cmap;
  else
    cmap_base = cmap(1:(pos_par-1));
  end
  if (exist(cmap_base) ~= 2)
    error(['The colormap ''' cmap ''' doesn''t exist']);
  end
  
  eval(['colormap(' cmap ')']);
else
  colormap(cmap);
end