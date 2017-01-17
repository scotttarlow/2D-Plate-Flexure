function twopf = twoDpf()
%Two Dimensional Plate flexure Model, Indept Study Fall 2013,SIUC
%Scott Tarlow

%Initial Parameters
nx = 50; %number of nodes in the x,y directions
ny = 50;
xDist = 1000; %x/y km 
yDist = 1000; 
dx = xDist/nx;
dy = yDist/ny;
g = 9.8;
pc = 2700; % density of crust
pm = 3300; % density of mantle
w0 = 80000; %initial flexure
v = .3; %poissons ratio

% Horizontal Forces in GPa * Km^2
Hx = 0;
Hy = 0;

% Vertical Forces/Load
hf = zeros(nx,ny);
q = 7;
hf(20:23,10-q:15-q)= 7000; % height of seamount in meters 
%hf(20:23,17-q:20-q) = 4000; 
V = pc*g*hf;

%Rigidity in GPa * km^2
D = ones(nx*ny,ny*nx);
D(1:nx*ny,1:nx*ny) = 500*D;
for i = 1:nx
    for j=1:ny
        if i == j
            D(i,j) = 5;
            %D(i+1,j+1) = 5;
            %D(i+2,j+2) = 5;
        end
    end
end




%Building A matrix // Setting Boundary conditions
b = reshape(V,nx*ny,1);



%Flexure Portion 
W1 = zeros(nx*ny,nx*ny); % Western Side w=w0
for i = 1:nx
    for j= 1:ny
        if i == j;
            W1(i,j) = 1;
        end
    end
end
b(1:ny,1) = w0;


W2 = zeros(nx*ny,nx*ny); % Eastern Side w = 0
for i = (nx*ny)-nx:(nx*ny)
    for j = (nx*ny)-nx:(nx*ny)
        if i == j;
            W2(i,j) = 1;
        end
    end
end

%First derivative 
Onex = zeros(nx*ny,nx*ny); 
for i = 2+nx:(nx*ny)-(nx+1)
    Onex(i,i-ny) = -1/dx;
    Onex(i,i) = 1/dx;
end


Oney = zeros(nx*ny,nx*ny);
for i = 1+nx:(nx*ny)-(nx+1)
    Oney(i-1,i)= -1/dy;
    Oney(i,i) = 1/dy;
end

%second derivative

Twox = zeros(nx*ny,nx*ny);
for i = 2+nx:(nx*ny)-(nx+1);
    Twox(i,i-(ny+1)) = 1/(dx^2);
    Twox(i,i) = -2/(dx^2);
    Twox(i,i-(ny-1)) = -2/(dy^2);
end

Twoy = zeros(nx*ny,nx*ny);
for i = 2+nx:(nx*ny)-(nx+1);
    Twoy(i-1,i) = 1/(dy^2);
    Twoy(i,i) = -2/(dy^2);
    Twoy(i+1,1)= 1/(dy^2);
end



A = W1 + W2 + (Twox*(D*(Twox + v*Twoy)) + Twoy*(D*(Twoy + v*Twox)) + ...
    2*(Onex*Oney)*(D*(1-v)*(Onex*Oney)) + Hx*Twox + Hy*Twoy);



%restoring force

for i = 2+nx:(nx*ny)-(nx+1);
    A(i,i) = A(i,i) + (((pm-pc)/(pm))*g);
end

w = A\b;
w = reshape(w,nx,ny);
figure()
contourf(1:dx:xDist,1:dy:yDist,w);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
clabel(w,h);
axis ij


%twopf = D;

















