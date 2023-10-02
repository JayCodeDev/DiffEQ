%DiffEQ library of funcitons

function [PHI,S,X,Y,soln_time,A_mem] = poisson_2d(src_fun_name,L,H,num_pts)
Nx = num_pts;Ny = num_pts;Ntot = Nx*Ny;
x = linspace(0,L,Nx); dx = x(2)-x(1);y = linspace(0,H,Ny); dy = y(2)-y(1);
[X,Y] = meshgrid(x,y);
S = feval(src_fun_name,X,Y,L,H);
A = spalloc(Ntot,Ntot,5*Ntot);b = zeros(Ntot,1);
A_whos = whos('A');A_mem = A_whos.bytes;
i = 1;
for j=1:Ny
    n = get_n(i,j,Nx);A(n,n) = 1;b(n) = feval(src_fun_name, 0, y(j), L, H);
end
i = Nx;
for j=1:Ny
    n = get_n(i,j,Nx);A(n,n) = 1;b(n) = feval(src_fun_name, L, y(j), L, H);
end
j = 1;
for i=1:Nx
    n = get_n(i,j,Nx);A(n,n) = 1;b(n) = feval(src_fun_name, x(i), 0, L, H);
end
j = Ny;
for i=1:Nx
    n = get_n(i,j,Nx);A(n,n) = 1;b(n) = feval(src_fun_name, x(i), H, L, H);
end
factor_x = 1/(dx^2); factor_y = 1/(dy^2);factor_center = 2*(factor_x + factor_y);
for i=2:(Nx-1)
    for j=2:(Ny-1)
        n = get_n(i,j,Nx);A(n,n) = factor_center;A(n,get_n(i-1,j,Nx)) = factor_x;
        A(n,get_n(i+1,j,Nx)) = factor_x;A(n,get_n(i,j-1,Nx)) = factor_y;
        A(n,get_n(i,j+1,Nx)) = factor_y;b(n) = S(i,j);
    end
end
tic
phi = A\b;
soln_time = toc;
PHI = zeros(size(X));
for i=1:Nx
    for j=1:Ny
        n = get_n(i,j,Nx);PHI(i,j) = phi(n);
    end
end
end
function n = get_n(i,j,Nx)
n = (j-1)*Nx + i;
end
function test_poisson_2d_uniform()
src_fun_name = 'sourceuniform';
L = 1;H = 1;num_pts = 8;
[PHI, S, X, Y, soln_time, A_mem] = poisson_2d(src_fun_name, L, H, num_pts);
figure(1);
subplot(2,5,1);contourf(X,Y,S,10);colorbar;
xlabel('x');ylabel('y');title('Uniform-Source - s(x,y)');
subplot(2,5,6);contourf(X,Y,PHI,10);colorbar;
xlabel('x');ylabel('y');title('Uniform-Source - Poisson EQN SOLN \phi(x,y)');
end
function test_poisson_2d_zero()
src_fun_name = 'sourcezero';
L = 1;H = 1;num_pts = 8;
[PHI, S, X, Y, soln_time, A_mem] = poisson_2d(src_fun_name, L, H, num_pts);
figure(1);
subplot(2,5,2);contourf(X,Y,S,10); colorbar;
xlabel('x');ylabel('y');title('Zero-Source - s(x,y)');
subplot(2,5,7);contourf(X,Y,PHI,10);colorbar;
xlabel('x');ylabel('y');title('Zero-Source - Poisson EQN SOLN \phi(x,y)');
end
function test_poisson_2d_linearX()
src_fun_name = 'sourcelinearX';
L = 1;H = 1;num_pts = 8;
[PHI, S, X, Y, soln_time, A_mem] = poisson_2d(src_fun_name, L, H, num_pts);
figure(1);
subplot(2,5,3);contourf(X,Y,S,10);colorbar;
xlabel('x');ylabel('y');title('LinearX - s(x,y)');
subplot(2,5,8);contourf(X,Y,PHI,10);colorbar;
xlabel('x');ylabel('y');title('LinearX - Poisson EQN SOLN \phi(x,y)');
end
function test_poisson_2d_linearXY()
src_fun_name = 'sourcelinearXY';
L = 1;H = 1;num_pts = 8;
[PHI, S, X, Y, soln_time, A_mem] = poisson_2d(src_fun_name, L, H, num_pts);
figure(1);
subplot(2,5,4);contourf(X,Y,S,10);colorbar;
xlabel('x');ylabel('y');title('LinearXY - s(x,y)');
subplot(2,5,9);contourf(X,Y,PHI,10);colorbar;
xlabel('x');ylabel('y');title('LinearXY - Poisson EQN SOLN \phi(x,y)');
end
function test_poisson_2d_nonlinearXY()
src_fun_name = 'sourcenonlinearXY';
L = 1;H = 1;num_pts = 8;
[PHI, S, X, Y, soln_time, A_mem] = poisson_2d(src_fun_name, L, H, num_pts);
figure(1);
subplot(2,5,5);contourf(X,Y,S,10);colorbar;
xlabel('x');ylabel('y');title('Non-LinearXY - s(x,y)');
subplot(2,5,10);contourf(X,Y,PHI,10);colorbar;
xlabel('x');ylabel('y');title('Non-LinearXY - Poisson EQN SOLN \phi(x,y)');
end
function S = sourceuniform(X,Y,L,H)
S = ones(size(X));
end
function S = sourcezero(X,Y,L,H)
so = 0;S = so*ones(size(X));
end
function S = sourcelinearX(X,Y,L,H)
S = X/H;
end
function S = sourcelinearXY(X,Y,L,H)
S = X/L +2*Y/H;
end
function S = sourcenonlinearXY(X,Y,L,H)
S = (X/L)^2 + sin(10*Y/H);
end
function Misc()

M1=(2-1)/(1-(1/1.5))
M2=(2.5-1)/(1-(1/3))

%Finite-Difference Method
%Problem 1)

% Problem 1a & 1b
L = 1; % Length of the polymer strip (cm)
W = 1; % Width of the polymer strip (cm)
H = 0.1; % Thickness of the polymer strip (cm)
sig = 1000; % Conductivity of the polymer (S/cm)
V0 = 0; % Voltage at one end of the polymer (grounded)
V1 = 1; % Voltage at the other end of the polymer (1 V)

fprintf("\na)\n\nRefer to Figure(1)\n");

% Problem 1.a
nlist = [floor(logspace(1, 3, 10)) 1500];

tlist = zeros(size(nlist));
mlist = zeros(size(nlist));

for i = 1:length(nlist)
    N = nlist(i);
    dx = L/(N-1); % Distance between points
    x = linspace(0, L, N); % Array of x values
    A = zeros(N, N);
    b = zeros(N, 1);
    for j = 2:N-1
        A(j, j-1) = sig/dx^2;
        A(j, j) = -2*sig/dx^2;
        A(j, j+1) = sig/dx^2;
    end
    A(1, 1) = 1;
    b(1) = V0;
    A(N, N) = 1;
    b(N) = V1/sig*dx;
    tic;
    V = A\b;
    tlist(i) = toc;
    I = sig*(V1 - V0)/dx;
    mlist(i) = whos('A', 'b', 'V', 'I').bytes;
end

figure(1)
plot(x, V);
xlabel('Pos (cm)');
ylabel('Volt (V)');
title('Voltage Distribution Along a Conductive Polymer Strip');
grid on

fprintf("\nb)\n\nRefer to Figure(2)\n");

% Problem 1.b
figure(2)
yyaxis left;
plot(nlist, tlist);
xlabel('# of grid points');
ylabel('Comp time (s)');
yyaxis right;
plot(nlist, mlist);
ylabel('Memory utilization (bytes)');
title('Computation time and memory utilization as a function of grid points');
legend('Time', 'Memory');
grid on

fprintf("\nc)\n\nRefer to Figure(3)\n");

% Problem 1.c
nlist = [floor(logspace(1, 3, 10)) 1500 5000 10000 50000 100000 500000 1000000];
tlist = zeros(size(nlist));
mlist = zeros(size(nlist));

for i = 1:length(nlist)
    N = nlist(i);
    dx = L/(N-1); % Distance between points
    x = linspace(0, L, N); % Array of x values
    A = spdiags([ones(N,1), -2*ones(N,1), ones(N,1)], [-1, 0, 1], N, N);
    A(1,1) = 1;
    A(N,N) = 1;
    A = (sig/dx^2)*A;
    b = zeros(N, 1);
    b(1) = V0;
    b(N) = V1/sig*dx;
    tic;
    V = A\b;
    tlist(i) = toc;
    I = sig*(V1 - V0)/dx;
    mlist(i) = whos('A', 'b', 'V', 'I').bytes;
end

% comptime
figure(3)
yyaxis left;
loglog(nlist, tlist);
xlabel('# of grid points');
ylabel('Comp time (s)');
yyaxis right;
loglog(nlist, mlist);
ylabel('Memory utilization (bytes)');
title('Computation time and memory utilization vs grid points');
legend('Time', 'Memory');
grid on

fprintf("\nd)\n\nConductivity is not part of the equation because\nwe are calculating the rate of change of voltage\nand that is not greatly affected by the conductivity.\n");

fprintf("\ne)\n\nConductivity would need to be included in our\nequation as it would change our voltage over the\nstrip.\n");

fprintf("\nf)\n\nWe would need to consider the width of the line\nand consider a lower step size to account for the\nadditonal complexity. The calculation would also take\nmuch longer to run.\n");

end
function [xout,yout,zout,youtend] = Shooter(d2y,x0,y0,x1,y1,del,dy0,steps);
syms x y z
% %%Given
% d2y(x,y,z)=2*y+8*x*(9-x);
% x0=0;
% y0=0;
% x1=9;
% y1=0;
% 
% %%Pick
% del=3;
% dy0=-24;
% steps="true";
% [xout,yout,zout,youtend] = Shooter(d2y,x0,y0,x1,y1,del,dy0,steps)

%function
dy(x,y,z)=z;
dz(x,y,z)=d2y(x,y,z);

f1(x,y,z)=dy(x,y,z);
f2(x,y,z)=dz(x,y,z);

%Euler Method

xi=x0;
yi=y0;
zi=dy0;
xout(1,1)=xi;
yout(1,1)=yi;
zout(1,1)=zi;

I=(x1-x0)/del;

if steps=="true"
    fprintf("I\t\t\txi\t\t\tyi\t\t\tzi\n");
end

i=0;
while i<=I
    xii=xi+del;
    yii=vpa(yi+del*f1(xi,yi,zi),4);
    zii=vpa(zi+del*f2(xi,yi,zi),4);
    
    xout(1,i+1)=xii;
    yout(1,i+1)=yii;
    zout(1,i+1)=zii;
    
    if steps=="true"
        fprintf("%f\t%f\t%f\t%f\n",i,xi,yi,zi);
    end
    
    yi=yii;
    xi=xii;
    zi=zii;
    
    i=i+1;
end
youtend=yout(1,i);
if steps=="true"
fprintf("\n");
end

end
function [yout,xout] = RK4(xi,del,yi,df,I,steps)
%Inputs
% xi=0;
% yi=1;
% f(x,y)=x*y;
% I=10;
% steps="true";   %"true" or "false"
% del=0.1;
% [yout,xout] = RK4(xi,del,yi,f,I,steps);

xout(1,1)=xi;
yout(1,1)=yi;

if steps=="true"
    fprintf("I\txi\t\t\tyi\t\t\tyii\n");
end

i=1;
while i<=I
    k1 = del*df(xi,yi);
    k2 = del*df(xi+del/2,yi+k1/2);
    k3 = del*df(xi+del/2,yi+k2/2);
    k4 = del*df(xi+del,yi+k3);

    yii = yi + (1/6)*(k1+2*k2+2*k3+k4);
    
    if steps=="true"
        fprintf("%i\t%f\t%f\t%f\n",i,xi,yi,yii);
    end
    
    xi=xi+del;
    yi=yii;
    
    xout(1,i+1)=xi;
    yout(1,i+1)=yii;
    
    i=i+1;
end
if steps=="true"
    fprintf("\n");
end
end
function [yout,xout] = EulerMeth(xi,del,yi,df,I,steps)
%Inputs
% xi=0;
% del=0.1;
% yi=1;
% df(x,y)=x*y;
% I=10;
% steps="true";   %"true" or "false"
% [yout,xout] = EulerMeth(xi,del,yi,df,I,steps);

xout(1,1)=xi;
yout(1,1)=yi;

%Euler Method

if steps=="true"
    fprintf("I\t\t\txi\t\t\txii\t\t\tyi\t\t\tyii\n");
end

i=1;
while i<=I
    xii=xi+del;
    yii=vpa(yi+(xii-xi)*df(xi,yi),4);
    xout(1,i+1)=xii;
    yout(1,i+1)=yii;
    if steps=="true"
        fprintf("%f\t%f\t%f\t%f\t%f\n",i,xi,xii,yi,yii);
    end
    
    yi=yii;
    xi=xii;
    
    i=i+1;
end

if steps=="true"
fprintf("\n");
end

end
function Newton()

eps = 0.001;
x1 = 1;
y1 = 1;
next = [x1;y1]
for i=1:5
    x = next(1)
    y = next(2)
    J_i = [sinh(x) -1;y x]
    f_i = [cosh(x)-y; x*y-1]
    del_x = inv(J_i)*(-f_i)
    next(1) = [x+del_x(1)]
    next(2) = [y+del_x(2)]
    f_next = [cosh(next(1))-next(2);next(1)*next(2)-1]
    magnitude = sqrt(f_next(1)^2+f_next(2)^2)
    if magnitude <= eps
        solution = next
        i
        break
    end
end

end
function Broyden()

eps = 0.001;
x1 = 1;
y1 = 1;
next = [x1;y1]
B_i = [sinh(x1) -1;y1 x1]
for i=1:30
    i
    x = next(1)
    y = next(2)
    f_i = [cosh(x)-y;x*y-1]
    del_x = inv(B_i)*(-f_i)
    next(1) = [x+del_x(1)];
    next(2) = [y+del_x(2)]
    f_next = [cosh(next(1))-next(2);next(1)*next(2)-1]
    B_i = B_i+(f_next*del_x')/(del_x(1)^2+del_x(2)^2)
    magnitude = sqrt(f_next(1)^2+f_next(2)^2)
    if magnitude <= eps
        solution = next
        i
    break
    end
end

end
function SecantMethod(f,x0,x1,I)
% syms x
% f(x)=sqrt((4*x^2)-0.04)+0.24*asin(0.1/x)-10.055752;
% f(x)=-1-0.5*cotd(x)+1/sind(x);
% x0=60;   %xn-1
% x1=10;   %xn
% I=12;
% SecantMethod(f,x0,x1,I);

x0=vpa(x0);
x0_rep=vpa(x0,4);
x1=vpa(x1);
x1_rep=vpa(x1,4);

count=1;
while count<=I
    
    X=x1-f(x1)*((x1-x0)/(f(x1)-f(x0)));
    X_rep=vpa(X,4);
    
    fprintf("Iteration %i)\tXn-1 = %f & Xn = %f\t>>\tX = %f\n",count,x0_rep,x1_rep,X_rep);
    
    x0=x1;
    x0=vpa(x0);
    x0_rep=vpa(x0,4);
    x1=X;
    x1=vpa(x1);
    x1_rep=vpa(x1,4);
    
    count=count+1;
end

end
function NewtRaph(x0,f,I,e,Upper,Lower)
%%Inputs
% f(x)=x^2-15;
% I=10;
% e=0.001;
% x0=1;
% Upper=0;
% Lower=-100;
% NewtRaph(x0,f,I,e,Upper,Lower)
syms x
%%functions
df(x)=diff(f,1,x);  %Take the derivative of f(x)
Icount=1;   %Start the Iteration counter 
tol=999;    %Set tol value so the while loop will begin without error
OoBH=0;     %Set OoBH flag to 0 so the while loop will begin without error
OoBL=0;     %Set OoBL flag to 0 so the while loop will begin without error
while tol>e && Icount<=I && OoBH==0 && OoBL==0
    Icount_Rep=Icount;  %Save Icount value for display
    if df(x0)==0    %Determine if we hit a peak or trough
        x0=x0+0.00001;  %increase the guess slightly to try and resolve the error
    end
    X=x0-f(x0)/df(x0);  %Calculate the next guess
    tol=vpa(abs(X-x0),4);   %Calculate the current tolerance
    fprintf("Iteration %i \t x0 = %f \t X = %f\n",Icount,x0,X); %Display the calculations
    x0=X;   %Store the next guess as the new current guess
    Icount=Icount+1;    %Increment the Iteration counter
    if X>Upper  %Check if the next guess is above the Upper limit
        OoBH=1; %If the next guess is above the Upper limit set the OoBH flag to 1
    elseif X<Lower  %Check if the next guess is below the Lower limit
        OoBL=1; %If the next guess is below the Lower limit set the OoBH flag to 1
    end
end
if OoBH==1  %If OoBH Flag is set to 1, then display the Upper limit fail message
    fprintf("\nX hit the Upper limit of %f\n\n",Upper);
elseif OoBL==1  %If OoBL Flag is set to 1, then display the Lower limit fail message
    fprintf("\nX hit the Lower limit of %f\n\n",Lower);
elseif Icount>I  %If the maximum number of iterations is achieved, then display the Max iterations achieved message
    fprintf("\nMaximum iterations achieve without solution\n\n");
elseif OoBH==0 && OoBL==0 && Icount<=I  %If the acceptance critera is achieved, then display the found value
    fprintf("\nX = %f after %d iterations with a tolerance of %f\n\n",vpa(X,4),Icount_Rep,tol);
else    %If an unknown error has occured, then display the ID-10-T User Error Message
    fprintf("\nERROR: How did you possibly mess this thing up???\n\n");
end
end
function DiffEQ_HW2_1_Succ_Subs()
syms x

disp("Problem 1)");
f=x==cos(x);
x_aprx=vpasolve(f);
x_aprx=round(x_aprx,3)
disp("My solution was x = 0.739");
disp(" ");
disp("Problem 2)");
f=x^3+5*x-3==0;
x_aprx=vpasolve(f);
x_aprx=round(x_aprx(1,1),3)
disp("My solution was x = 0.564");
% x_val=0.56:0.001:0.57;
% f=x_val.^3+5.*x_val-3;
% plot(x_val,f)


end
function [L U]=my_lu(G,steps)
% G=[1 3 1 5;2 1 0 3;4 2 2 1;-3 1 3 2];
% steps="false";  %"true" or "false"
% [L U]=my_lu(G,steps)

[M,N]=size(G);
L=eye(M);

%Upper Triangle
col=1;
while col<M+1
    row=1;
    while row<M+1
        if(row<col)
            
        elseif(row==col)
            rowref=row;
            val=G(row,col);
        else
            val2=G(row,col);
            if(steps=="true")
                G(row,:)=G(row,:)-(val2/val)*G(rowref,:)
                L(row,col)=val2/val
            else
                G(row,:)=G(row,:)-(val2/val)*G(rowref,:);
                L(row,col)=val2/val;
            end
        end
        row=row+1;
    end
    col=col+1;
end

U=G;

end
function [G]=Gauss_Elim(A,b,steps)
% A=[12 -2 0 9;2 -10 5 0; 0 5 -22 7;9 0 7 -24];
% b=[100;0;0;200];
% steps="false";  %"true" or "false"
% [G]=Gauss_Elim(A,b,steps)

G=[A b]
[M,N]=size(G);

%Upper Triangle
col=1;
while col<M+1
    row=1;
    while row<M+1
        if(row<col)
            
        elseif(row==col)
            rowref=row;
            val=G(row,col);
        else
            val2=G(row,col);
            if(steps=="true")
                G(row,:)=G(row,:)-(val2/val)*G(rowref,:)
            else
                G(row,:)=G(row,:)-(val2/val)*G(rowref,:);
            end
        end
        row=row+1;
    end
    col=col+1;
end

%ones
col=1;
while col<M+1
    row=1;
    while row<M+1
        if(row<col)
            
        elseif(row==col)
            val3=G(row,col);
            if(steps=="true")
                G(row,:)=(1/val3)*G(row,:)
            else
                G(row,:)=(1/val3)*G(row,:);
            end
        else
            
        end
        row=row+1;
    end
    col=col+1;
end

%Identity
col=M;
while col>0
    row=M;
    while row>0
        if(row<col)
            val5=G(row,col);
            if(steps=="true")
                G(row,:)=G(row,:)-(val5/val4)*G(rowref,:)
            else
                G(row,:)=G(row,:)-(val5/val4)*G(rowref,:);
            end
        elseif(row==col)
            rowref=row;
            val4=G(row,col);
        else
            
        end
        row=row-1;
    end
    col=col-1;
end

end
function DiffEQ_HW1_2()
disp("*******************************************************************");
disp("1) (2 points) Unit vector");
    disp(" ");
    u=[-2;5;2;4]
    u_mag=norm(u)
    u_unit=u/u_mag
    %Check
    u_unit_mag=norm(u_unit);
    disp(" ");
disp("*******************************************************************");
disp("2) (2 points) Mag of vector");
    disp(" ");
    v=[1;-2;6]
    v_mag=norm(v)
    disp(" ");
disp("*******************************************************************");
disp("3) (2 points) Dot Product of vectors");
    disp(" ");
    a=[1;0;4]
    b=[1;2;-1]
    ab_dot=dot(a,b)
    disp(" ");
disp("*******************************************************************");
disp("4) (2 points) Angle between vectors");
    disp(" ");
    a=[1;0;4]
    b=[1;2;-1]
    disp("The angle between a and b can be found by taking the inverse");
    disp("cosine of the dot product of a and b and dividing it by the");
    disp("magnitude of a and the magnitude of b.");
    disp("acos(dot(a,b)/(|a||b|))");
    ab_dot=dot(a,b)
    a_mag=norm(a)
    b_mag=norm(b)
    ab_angle=acosd(ab_dot/(a_mag*b_mag))
    disp(" ");
disp("*******************************************************************");
disp("5) (2 points) Geometric significance of dot product between vectors");
    disp(" ");
    x=[1;2;3;4]
    y=[-2;1;4;-3]
    xy_dot=dot(x,y)
    fprintf("The significance of a zero dot product is that the angle between\nthe two vectors is 90 degrees so they are orthogonal.");
    disp(" ");
disp("*******************************************************************");
disp("6) (2 points) Product of AB");
    disp(" ");
    A=[1 -1 2;0 -3 1]
    B=[2;1;0]
    AB=A*B
    disp(" ");
disp("*******************************************************************");
disp("7) (2 points) Product of AB");
    disp(" ");
    A=[0 4 -2;-4 -3 0]
    B=[0 1;1 -1;2 3]
    AB=A*B
    disp(" ");
disp("*******************************************************************");
disp("8) (2 points) Problem 3.11");
    disp(" ");
    A=[-2 8 3; 2 -4 -6;5 0 7]
    B=[-3 7 3;4 -6 8;9 2 -5]
    a=A+B
    b=A-B
    c=2*A-3*B
    d=A'-2*B'
    disp(" ");
disp("*******************************************************************");
disp("9) (2 points) Problem 3.13");
    disp(" ");
    A=[2 -3 4;3 -1 5]
    B=[4 0;-1 6;3 -2]
    a=A*B
    b=B*A
    disp(" ");
disp("*******************************************************************");
end
function DiffEQ_HW0_4()
a_bin="10110";
disp("a_bin is "+a_bin);
b_bin="110011";
disp("b_bin is "+b_bin);
disp(" ");
disp("1) (5 points) Problem 1.22");
    disp(" ");
    a_dec=bin2dec(a_bin);
    disp(">> a_dec is "+a_dec);
    b_dec=bin2dec(b_bin);
    disp(">> b_dec is "+b_dec);
    disp(" ");
disp("2) (5 points) Problem 1.22, but interpret the binary numbers ");
disp("as signed integers versus unsigned integers");
    disp(" ");
    a_dec_signed=int8(a_dec);
    disp(">> a_dec_signed is "+a_dec_signed);
    a_dec_unsigned=uint8(a_dec);
    disp(">> a_dec_unsigned is "+a_dec_unsigned);
    b_dec_signed=int8(b_dec);
    disp(">> b_dec_signed is "+b_dec_signed);
    b_dec_unsigned=uint8(b_dec);
    disp(">> b_dec_unsigned is "+b_dec_unsigned);
    disp(" ");
disp("3) (5 points) Write code to determine machine epsilon for ");
disp("double precision floating point numbers. Machine epsilon ");
disp("is the smallest numberfor which the computer cannot tell ");
disp("the difference between 1 and 1 +eps");
    disp(" ");
    Epsilon=eps('double');
    disp(">> Machine epsilon for a double is "+Epsilon);
    disp(" ");
disp("4) (5 points) Write code to determine machine epsilon for ");
disp("single precision floating point numbers. Machine epsilon is ");
disp("the smallest number for which the computer cannot tell the ");
disp("difference between 1 and1 +eps. Explain the difference ");
disp("between single versus double precision floating point numbers.");
    disp(" ");
    Epsilon=eps('single');
    disp(">> Machine epsilon for a single is "+Epsilon);
    disp(">> In single precision, 32 bits are used to represent ");
    disp("floating-point number. In double precision, 64 bits are ");
    disp("used to represent floating-point number.");
end
function my_sqrt(A,x,ConvCriteria)
% Num Meth & Partial Diff Eq
% A=1500;    %SquareRoot(A)
% x(1,1)=1;   %Guess
% ConvCriteria=0.01;
% my_sqrt(A,x,ConvCriteria);

tol=999;
n=1;
while tol>ConvCriteria && n<100
    y(1,n)=(x(1,n)^2)/A;
    x(1,n+1)=(x(1,n)/8)*(15-10*y(1,n)+3*y(1,n)^2);
    tol=abs(x(1,n+1)-x(1,n));
    n=n+1;
end
plot(x);
SqrA=x(1,n);

disp("Square root of "+A+" is about "+SqrA+" after "+n+" iterations");

end
function mywierdfcn(x)
% Num Meth & Partial Diff Eq
% x=[-20:0.1:20];
% mywierdfcn(x);

len=length(x);
y=zeros(1,len);
i=1;
counter=x(1,i);
stop=x(1,len);
while counter<stop
    if i==1
        if counter<=2
            y(1,1)=x(1,1)*2;
        elseif counter>2
            y(1,1)=x(1,1)^2;
        end
    end
    i=i+1;
    counter=x(1,i);
    if counter<=2
        y(1,i)=x(1,i)*2;
    elseif counter>2
        y(1,i)=x(1,i)^2;
    end
end
plot(x,y);
end
function fib_golden_ratio(Num1,Num2,tolerance)
% Num Meth & Partial Diff Eq
% Num1=1;
% Num2=2;
% tolerance=0.001;
% fib_golden_ratio(Num1,Num2,tolerance);

n=100;
string=Num1+", "+Num2+", ";
count=1;
tol=inf;
while tol > tolerance
    FibNum=Num1+Num2;
    tol=abs(FibNum/Num2-Num2/Num1);
    if tol < tolerance
        string=string+FibNum;
    else
        string=string+FibNum+", ";
    end
    Num1=Num2;
    Num2=FibNum;
    count=count+1;
end
disp("Fibonacci sequence: "+string);
disp("Tolerance: "+tol);
end
function fib_while(Num1,Num2,n)
% Num Meth & Partial Diff Eq
% Num1=1;
% Num2=2;
% n=10;
% fib_while(Num1,Num2,n);

string=Num1+", "+Num2+", ";
count=1;
while count < n-1
    FibNum=Num1+Num2;
    if count == n-2
        string=string+FibNum;
    else
        string=string+FibNum+", ";
    end
    Num1=Num2;
    Num2=FibNum;
    count=count+1;
end
disp("Fibonacci sequence: "+string);
end
function SumOfSeriesCalc()
% Num Meth & Partial Diff Eq
Num=1;
n=32;

string=Num+", ";
SumOfSeries=Num;
count=2;
for count = 2:n
    Num=Num*2;
    SeriesNum=(1/Num)*(-1)^(count-1);
    SumOfSeries=SumOfSeries+SeriesNum;
    if count == n
        string=string+SeriesNum;
    else
        string=string+SeriesNum+", ";
    end
    count=count+1;
end
disp("Sum of series: "+SumOfSeries+" = "+string);

end
function noisy_dataplottool(noisy_data,x,n)
% Num Meth & Partial Diff Eq
% n=100;
% x=[1:n];
% noisy_data=10+2.0*randn(n,1);
% t = linspace(0, 10, n);
% x=t;
% noise = 1 + 0.25* randn (1 ,n);
% noisy_data = sin(t) + noise;
% noisy_data = transpose(noisy_data);
% noisy_dataplottool(noisy_data,x,n);

ynoise(1,1)=noisy_data(1,1);
ysmooth(1,1)=noisy_data(1,1);
count=2;
for count=2:n
    ynoise(1,count)=noisy_data(count,1);
    ysmooth(1,count)=smooth(noisy_data(count,1),noisy_data(count-1,1));
    
    count=count+1;
end

plot(x,ynoise,'b',x,ysmooth,'r');
legend('noisy data','smooth data','Location','NorthEast')
end
function xavg=smooth(x1,x2)
% x1=0;
% x2=1;
% xavg=smooth_data(x1,x2)

xavg=(x1+x2)/2;
end
function fib(Num1,Num2,n)
% Num Meth & Partial Diff Eq
% Num1=1;
% Num2=2;
% n=10;
% fib(Num1,Num2,n);

string=Num1+", "+Num2+", ";
count=0;
for count = 0:n-1
    FibNum=Num1+Num2;
    if count == n-1
        string=string+FibNum;
    else
        string=string+FibNum+", ";
    end
    Num1=Num2;
    Num2=FibNum;
    count=count+1;
end
disp("Fibonacci sequence: "+string);
end
function PosVSTimePlot(A,Alpha,StartTime,EndTime,Points,f,phi)
% Num Meth & Partial Diff Eq
% A=10;    %Initial Amplitude
% Alpha=0.25;    %damping constant
% StartTime=0;
% EndTime=12;
% Points=200;
% f=0.25;  %frequency
% phi=pi/4;  %phase shift
% PosVSTimePlot(A,Alpha,StartTime,EndTime,Points,f,phi);

%Functions
W=2*pi*f;   %angular frequency
deltaT=(EndTime-StartTime)/Points;
t=[StartTime:deltaT:EndTime];  %time
T=2*pi/W  %Period
y=A.*exp(-Alpha.*t).*cos(W.*t+phi);

%Outputs
disp("The period of the wave is " + T + " seconds");
plot(t,y);
title('Position vs Time Graph');
xlabel('Time (Seconds)');
ylabel('Position (Meters)');

end