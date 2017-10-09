clear;
clc;
load inside_ILT_nodes.mat
inside_ILT_nodes = lumen_nodes;

load outside_ILT_nodes.mat
outside_ILT_nodes = wall_nodes;

%read output file frome gmsh
original_file = fopen('010_12_INP.inp','r');
wall_stl_filename  = 'wall_P.stl';
lumen_stl_filename = 'lumen_P.stl';
planes_name        = 'PLANES.stl';
abaqus_inp_name = '010_12.inp';

A = fscanf(original_file,'%c'); %scan in file as a long butt string
start_nodes = strfind(A, '*NODE');
end_nodes = strfind(A, '******* E');

start_surface_elements = strfind(A,'Surface1');  
start_volume_elements  = strfind(A, 'Volume1');
where_elements         = strfind(A, '*ELEMENT');

all_nodes = A(start_nodes+5:end_nodes -1);
all_tri_elements = A(start_surface_elements + 8: where_elements(2)-1);
all_tet_elements = A(start_volume_elements + 7:end);

%Create temporary files
fid3 = fopen('temp_nodes.dat','w+');
fwrite(fid3,all_nodes);

fid1 = fopen('temp_tri_elements.dat','w+');
fwrite(fid1, all_tri_elements);

fid2 = fopen('temp_tet_elements.dat', 'w+');
fwrite(fid2, all_tet_elements);

fclose('all');

%read temporary files
all_nodes = dlmread('temp_nodes.dat');
nodes = all_nodes(:,1);
ax    = all_nodes(:,2);
ay    = all_nodes(:,3);
az    = all_nodes(:,4);

all_tri_elements = dlmread('temp_tri_elements.dat');
all_tet_elements = dlmread('temp_tet_elements.dat');
ILT_E = 1:length(all_tet_elements);
ILT_E = ILT_E';

%Get nodes and surfaces from original geometry (.stl) use stlread for stlread.
[wall_tri, wall_nodes] = stlread(wall_stl_filename);
[~, lumen_nodes ] = stlread(lumen_stl_filename);

we = 1:length(wall_tri); we = we'; %wall element connectivity
wn = 1:length(wall_nodes); wn = wn';
wx = wall_nodes(:,1);
wy = wall_nodes(:,2);
wz = wall_nodes(:,3);

lx = lumen_nodes(:,1);
ly = lumen_nodes(:,2);
lz = lumen_nodes(:,3);

%%Get boundaries
[bottom_boundary_nodes, top_boundary_nodes] = get_boundaries('PLANES.stl',all_nodes);

%Finds outside of ILT
logis=ismembertol([ax,ay,az],[wx,wy,wz],1e-8);
[badrows,~]=find(logis==0);
outside_tet_nodes=nodes; outside_tet_nodes(badrows,:)=[];

%Finds inside of ILT
logis=ismembertol([ax,ay,az],[lx,ly,lz],1e-8);
[badrows,~]=find(logis==0);
inside_tet_nodes=nodes; inside_tet_nodes(badrows,:)=[];

%Redundant filtering
bads = ismember(inside_tet_nodes,outside_tet_nodes);
keepnodes = ~bads;
inside_tet_nodes = inside_tet_nodes(keepnodes);

%Create lumen surface
axyz = [ax ay az];
a_nodes = 1:length(ax);
a_nodes = a_nodes';

real_lumen_xyz = axyz(inside_ILT_nodes,:);
lumen_tri = MyCrustOpen(real_lumen_xyz);
new_lumen_tri = [inside_ILT_nodes(lumen_tri(:,1));inside_ILT_nodes(lumen_tri(:,2));...
inside_ILT_nodes(lumen_tri(:,3))]';
lumen_tri_elem_nums = 1:length(new_lumen_tri);

%Create wall surface
real_wall_xyz = axyz(outside_ILT_nodes,:);
wall_tri = MyCrustOpen(real_wall_xyz);
new_wall_tri = [outside_ILT_nodes(wall_tri(:,1)),outside_ILT_nodes(wall_tri(:,2)),...
    outside_ILT_nodes(wall_tri(:,3))];

%Fixes element connectivity numbers
wall_tri_elem_nums = 1:length(new_wall_tri);
lumen_tri_elem_nums = lumen_tri_elem_nums+length(wall_tri_elem_nums)+length(ILT_E);
ILT_E = ILT_E+length(wall_tri_elem_nums);

%Write abaqus %INP file
fid3=fopen(abaqus_inp_name,'w');

%File header & nodes
fprintf(fid3,'*HEADING\n*PREPRINT,  ECHO=NO,  MODEL=NO,  HISTORY=NO, CONTACT=NO\n');
fprintf(fid3,'*Part, name=AAA\n*Node\n');
fprintf(fid3,'%d, %f, %f, %f\n', [all_nodes(:,1) 10*all_nodes(:,2:end)]');

%Print wall shell elements
fprintf(fid3,'*ELEMENT, type=S3R, ELSET=WALL\n');
fprintf(fid3,'%d, %d, %d, %d\n', [wall_tri_elem_nums', new_wall_tri]');

%Print ILT solid elements
fprintf(fid3, '*ELEMENT, type = C3D4H, ELSET=ILT\n');
fprintf(fid3,'%d, %d, %d, %d, %d\n', [ILT_E all_tet_elements(:,2:5)]');

%Lumen elements
fprintf(fid3,'*ELEMENT, type=S3R, ELSET=LUMEN\n');
fprintf(fid3,'%d, %d, %d, %d\n', [lumen_tri_elem_nums', new_lumen_tri]');

%Section Assignments
fprintf(fid3,'%s\n','*Solid Section, elset=ILT, material=ILT');
fprintf(fid3,'%s\n','1.,');
fprintf(fid3,'%s\n','*Shell Section, elset=WALL, material=WALL');
fprintf(fid3,'1.3, 5\n');
fprintf(fid3,'%s\n','*Shell Section, elset=LUMEN, material=ILT');
fprintf(fid3,'1, 5\n');

%Sets up top node and bottom node nsets for printing
toprem=mod(length(top_boundary_nodes),10); botrem=mod(length(bottom_boundary_nodes),10);
topheight=(length(top_boundary_nodes)-toprem)/10; botheight=(length(bottom_boundary_nodes)-botrem)/10;
botnodes=reshape(bottom_boundary_nodes(1:end-botrem),[botheight, 10]);
topnodes=reshape(top_boundary_nodes(1:end-toprem),[topheight, 10]);
botend=bottom_boundary_nodes(end-botrem:end); topend=top_boundary_nodes(end-toprem:end);

%Bottom & top boundaries
fprintf(fid3,'*Nset, nset=botnodes\n');
fprintf(fid3,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',botnodes);
fprintf(fid3,'%d, ',botend);
fprintf(fid3,'\n*Nset, nset=topnodes\n');
fprintf(fid3,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',topnodes);
fprintf(fid3,'%d, ',topend);

%Prepares nsets for printing (must print in rows of <15)
outrem=mod(length(outside_tet_nodes),10); inrem=mod(length(inside_tet_nodes),10);
outheight=(length(outside_tet_nodes)-outrem)/10; inheight=(length(inside_tet_nodes)-inrem)/10;
outnodes=reshape(outside_tet_nodes(1:end-outrem),[outheight, 10]); innodes=reshape(inside_tet_nodes(1:end-inrem),[inheight, 10]);
outend=outside_tet_nodes(end-outrem:end); inend=inside_tet_nodes(end-inrem:end);

%Remove top and bottom boundary from outnodes and innodes
remove_common_boundary_outnodes = ismember(outnodes, [topnodes; botnodes])
outnodes(remove_common_boundary_outnodes) = [];

remove_common_boundary_innodes  = ismember(innodes, [topnodes; botnodes]);
innodes(remove_common_boundary_innodes) = [];

%SURFACES FOR TIE CONTACT & LOADS:

%Lumen
fprintf(fid3,'\n*Surface, type=ELEMENT, name=lumensurf\nLUMEN, SNEG\n');

%Wall
fprintf(fid3,'*Surface, type=ELEMENT, name=wallsurf\nWALL, SNEG\n');

%Outside of ilt
fprintf(fid3,'*Nset, nset=outsideilt\n');
fprintf(fid3,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',outnodes);
fprintf(fid3,'%d, ',outend);
fprintf(fid3,'\n%s\n%s\n','*Surface, type=NODE, name=outsideiltsurf',...
    'outsideilt');

%Inside of ILT
fprintf(fid3,'*Nset, nset=insideilt\n');
fprintf(fid3,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',innodes);
fprintf(fid3,'%d, ',inend);
fprintf(fid3,'\n%s\n%s\n','*Surface, type=NODE, name=insideiltsurf',...
    'insideilt');
fprintf(fid3,'*End Part\n');

%Assembly
fprintf(fid3,'%s\n%s\n%s\n','*Assembly, name=Assembly',...
    '*Instance, name=AAA, part=AAA','*End Instance');

%Tie Contacts
fprintf(fid3,'%s\n%s\n','*Tie, name=tiecontact, adjust=yes, position tolerance=0',...
    'AAA.wallsurf, AAA.outsideiltsurf');
fprintf(fid3,'%s\n%s\n','*Tie, name=tiecontact2, adjust=yes, position tolerance=0',...
    'AAA.lumensurf, AAA.insideiltsurf');
fprintf(fid3,'*End Assembly\n');

%Material definitions
fprintf(fid3,'*Material, name=WALL\n');
fprintf(fid3,'*Hyperelastic, n=2\n');
fprintf(fid3,'1.74,0, 1.881,0,0,0,0\n'); %Material parameters for Raghavan-Vorp model from 2000 paper
fprintf(fid3,'%s\n','*Material, name=ILT');
fprintf(fid3,'%s\n','*Hyperelastic, n=2');
fprintf(fid3,'%s\n','0., 0.02804,    0.,    0., 0.02858,    0.,    0.');

%Step definitions
fprintf(fid3,'%s\n','*Step, name=Step1, nlgeom=YES');
fprintf(fid3,'%s\n','*Static');
fprintf(fid3,'%s\n','0.01, 1., 1e-5, 0.4');
fprintf(fid3,'%s\n','***********************************************');

%BCs
fprintf(fid3,'\n*Boundary\n');
fprintf(fid3,'AAA.topnodes, 1, 3\n');
fprintf(fid3,'AAA.botnodes, 1, 3\n');

%Loads
fprintf(fid3,'%s\n%s\n','*Dsload','AAA.lumensurf, P, 0.0157');

%Output requests
fprintf(fid3,'*Restart, write, frequency=1\n');
fprintf(fid3,'*Output, field\n');
fprintf(fid3,'*Node Output\n');
fprintf(fid3,'COORD, U\n');
fprintf(fid3,'*Element Output\n');
fprintf(fid3,'3\n');
fprintf(fid3,'EE, S\n');
fprintf(fid3,'*Output, history, variable=PRESELECT\n');
fprintf(fid3,'*El Print, freq=999999\n');
fprintf(fid3,'*Node Print, freq=999999\n');
fprintf(fid3,'*CONTROLS, PARAMETERS=TIME INCREMENTATION\n');
fprintf(fid3,'7, 10, 9, 25, 10, 7, 12, 8, 6, 3\n');
fprintf(fid3,'0.10, 0.5, 0.75, 0.85, 0.25, 0.75, 1.75, 0.75\n');
fprintf(fid3,'0.8, 1.5, 1.25, 2, 0.95, 0.1, 1, 0.95\n');
fprintf(fid3,'*End Step');

%End
fclose all;
disp('DONE!');
