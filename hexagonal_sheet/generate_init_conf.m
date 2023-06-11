close all
clear all

sigma_eq = 2^(1/6);
b = 1; %edge length
q = 1; %unit charge
zh = b;

X0 = [0,0,0; -b*cos(pi/3),b*sin(pi/3),0; 0,2*b*sin(pi/3),0; b,2*b*sin(pi/3),0; b+b*cos(pi/3),b*sin(pi/3),0; b,0,0; ...
      0,0,b; -b*cos(pi/3),b*sin(pi/3),b; 0,2*b*sin(pi/3),b; b,2*b*sin(pi/3),b; b+b*cos(pi/3),b*sin(pi/3),b; b,0,b; ...
      0,0,2*b; -b*cos(pi/3),b*sin(pi/3),2*b; 0,2*b*sin(pi/3),2*b; b,2*b*sin(pi/3),2*b; b+b*cos(pi/3),b*sin(pi/3),2*b; b,0,2*b;...
      b/2*cos(pi/3),-b/2*sin(pi/3),0; -b*cos(pi/3)-b/2*sin(pi/6),b*sin(pi/3)-b/2*cos(pi/6),0; -b/2,2*b*sin(pi/3),0;  b-b/2*cos(pi/3),2*b*sin(pi/3)+b/2*sin(pi/3),0; b+b*cos(pi/3)+b/2*sin(pi/6),b*sin(pi/3)+b/2*cos(pi/6),0; b+b/2,0,0; ...
      b/2*cos(pi/3),-b/2*sin(pi/3),2*b; -b*cos(pi/3)-b/2*sin(pi/6),b*sin(pi/3)-b/2*cos(pi/6),2*b; -b/2,2*b*sin(pi/3),2*b;  b-b/2*cos(pi/3),2*b*sin(pi/3)+b/2*sin(pi/3),2*b; b+b*cos(pi/3)+b/2*sin(pi/6),b*sin(pi/3)+b/2*cos(pi/6),2*b; b+b/2,0,2*b; ...
      -b/2,0,0; -b*cos(pi/3)-b/2*sin(pi/6),b*sin(pi/3)+b/2*cos(pi/6),0; b/2*sin(pi/6),2*b*sin(pi/3)+b/2*cos(pi/6),0; b+b/2,2*b*sin(pi/3),0; b+b*cos(pi/3)+b/2*sin(pi/6),b*sin(pi/3)-b/2*cos(pi/6),0; b-b/2*sin(pi/6),-b/2*cos(pi/6),0; ...
      -b/2,0,2*b; -b*cos(pi/3)-b/2*sin(pi/6),b*sin(pi/3)+b/2*cos(pi/6),2*b; b/2*sin(pi/6),2*b*sin(pi/3)+b/2*cos(pi/6),2*b; b+b/2,2*b*sin(pi/3),2*b; b+b*cos(pi/3)+b/2*sin(pi/6),b*sin(pi/3)-b/2*cos(pi/6),2*b; b-b/2*sin(pi/6),-b/2*cos(pi/6),2*b; ...
    ];
X0 = [X0; ...
      (X0(19,1:2)+X0(36,1:2))/2,zh; (X0(20,1:2)+X0(31,1:2))/2,zh; (X0(21,1:2)+X0(32,1:2))/2,zh; (X0(22,1:2)+X0(33,1:2))/2,zh; (X0(23,1:2)+X0(34,1:2))/2,zh; (X0(24,1:2)+X0(35,1:2))/2,zh ];



% plot3(X0(1:18,1),X0(1:18,2), X0(1:18,3),'o'); hold all
% plot3(X0(19:30,1),X0(19:30,2), X0(19:30,3),'*')
% plot3(X0(31:42,1),X0(31:42,2), X0(31:42,3),'s')
% plot3(X0(43:end,1),X0(43:end,2), X0(43:end,3),'d')

X = [];
mol = [];
ptype = [];
charge = [];
id = [];
m = 1;
atom = 1;

for lz=-2:2
for ly=-2:2
for lx=-2:2
    mol = [mol; m*ones(48,1)];
    ptype = [ptype; ones(18,1); 2*ones(12,1); 3*ones(12,1); 4*ones(6,1)];
    charge = [charge; zeros(18,1); q*ones(12,1); -q*ones(12,1); zeros(6,1)]; 
    R = [20.0*lx 20*ly 20*(lz)];  %2.732050807568878   cos(pi/6)*2+1
    for i=1: length(X0)
        X=[X; R + X0(i,:)];
        %plot3(X(end,1),X(end,2),X(end,3),'o'); hold all
    end
plot3(X(end,1),X(end,2),X(end,3),'o'); hold all 
id = [id; [atom:atom+47]'];
m = m+1;
atom = atom + 48;
end
end
end

axis equal
grid on

box2=max(max(X))+10;
box1=-box2; %min(min(X))-10;
out=sprintf('LAMMPS data file via MATLAB script\n\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n0 impropers\n\n4 atom types\n0 bond types\n0 angle types\n0 dihedral types\n\n', length(id),0,0,0);
out=[out, sprintf('%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\nAtoms\n\n',box1,box2,box1,box2,box1,box2)];
d=[id mol ptype charge X];
for j=1:length(id)
    out=[out,sprintf('%d %d %d %d %20.10f %20.10f %20.10f\n', d(j,1),d(j,2),d(j,3),d(j,4),d(j,5),d(j,6),d(j,7) )];
end
dlmwrite('conf', out, '');



