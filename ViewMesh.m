function p = ViewMesh(varargin)
%%%%%    ViewMesh(pt,trg) or ViewMesh(fpath,fname)
if nargin<1
    disp('Parameter : p = ViewMesh(pt,trg,\f\pos)  or p = ViewMesh(fname,\f\pos)');
    return
end
door=0;
pos=[-330.1630   95.0787  -93.5028];
pos=[27.6311  -27.1292  353.9745];
% ----------------------------------------------
%cow camera settings
% ----------------------------------------------
pos = [   10.5327    0.9505   -0.9002];
target =[   -0.0606   -0.0796    0  ];
upvector = [ 0 1 0];
upvectormode = 'manual';
viewangle = 8.0905;
if ischar(varargin{1})
    %fpath=varargin{1};
    fname=varargin{1};
    [ptVec trgnormal trgVec]= ReadObjShape([fname '.obj']);
    if nargin==2
        if length(varargin{2})==3
            pos = varargin{2};
        else f=varargin{2}; door=1;
        end
    end
    if nargin>2
        f = varargin{2}; door=1;
        pos = varargin{3};
    end
else ptVec=varargin{1};
    trgVec=varargin{2};
    if nargin==3
        if length(varargin{3})==3
            pos = varargin{3};
        else f=varargin{3}; door=1;
        end
    end
    if nargin>3
        f = varargin{3}; door=1;
        pos = varargin{4};
    end
end
color = [.99 .6863 0];
%color = [18 219 200]./255;

%figure

% pt=ptVec;
% trg=trgVec;

% for i=1:size(trgVec,1)
%     if trgnormal(i)>0
%         trg(i,:)=trgVec(i,[1 3 2]);
%     end
% end




p = patch('Vertices',ptVec, 'Faces',trgVec);
set(p,'FaceColor',color);
set(p,'edgecolor','none');
daspect([1 1 1])
view(3)

axis tight
axis off
axis equal
set(gca,'CameraPosition', pos);

material dull
camlight
lighting phong


if door==1
    set(p,'FaceVertexCData',f);
    set(p,'facecolor','interp');
    set(gcf,'renderer','zbuffer');
end

% camlight
% lighting gouraud
% % set(p,'edgecolor','k');



