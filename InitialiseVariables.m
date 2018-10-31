function [A,b,T,T_inf_max, T_inf_min] = InitialiseVariables()
%GIVEN VARIABLES----------------------------------------------------
k = 3; %thermal conductivity
delta_x = 0.01; %change in x 
delta_y = 0.01; %change in y
T_inf = 20; %ambient temperature
h = 20; %heat transfer coefficient

%MATRIX A (COEFFICIENT MATRIX)---------------------------------------
n = 33;
A = zeros(n);
m = ones(n,1);

A = diag(4*m,0);
A(1,1)=2;A(1,2)=-1;A(1,7)=-1;                     % row 1
A(2,2)=3;A(2,1)=-1;A(2,3)=-1;A(2,8)=-1;           % row 2
A(3,3)=3;A(3,2)=-1;A(3,4)=-1;A(3,9)=-1;           % row 3
A(4,4)=3;A(4,3)=-1;A(4,5)=-1;A(4,10)=-1;          % row 4
A(5,5)=3;A(5,4)=-1;A(5,6)=-1;A(5,11)=-1;          % row 5
A(6,6)=3;A(6,5)=-1;A(6,12)=-1;                    % row 6
A(7,7)=3;A(7,8)=-1;A(7,1)=-1;A(7,13)=-1;          % row 7
A(8,2)=-1;A(8,7)=-1;A(8,9)=-1;A(8,14)=-1;         % row 8
A(9,3)=-1;A(9,8)=-1;A(9,10)=-1;A(9,15)=-1;        % row 9
A(10,4)=-1;A(10,9)=-1;A(10,11)=-1;A(10,16)=-1;    % row 10
A(11,5)=-1;A(11,10)=-1;A(11,12)=-1;A(11,17)=-1;   % row 11
A(12,6)=-1;A(12,11)=-1;A(12,18)=-1;               % row 12
A(13,13)=3;A(13,7)=-1;A(13,14)=-1;A(13,19)=-1;    % row 13
A(14,8)=-1;A(14,13)=-1;A(14,15)=-1;A(14,20)=-1;   % row 14
A(15,9)=-1;A(15,14)=-1;A(15,16)=-1;A(15,21)=-1;   % row 15
A(16,10)=-1;A(16,15)=-1;A(16,17)=-1;A(16,22)=-1;  % row 16
A(17,11)=-1;A(17,16)=-1;A(17,18)=-1;A(17,23)=-1;  % row 17
A(18,12)=-1;A(18,17)=-1;A(18,24)=-1;              % row 18
A(19,19)=3;A(19,13)=-1;A(19,20)=-1;A(19,25)=-1;   % row 19
A(20,14)=-1;A(20,19)=-1;A(20,21)=-1;A(20,26)=-1;  % row 20
A(21,15)=-1;A(21,20)=-1;A(21,22)=-1;A(21,27)=-1;  % row 21
A(22,16)=-1;A(22,21)=-1;A(22,23)=-1;A(22,28)=-1;  % row 22
A(23,17)=-1;A(23,22)=-1;A(23,24)=-1;A(23,29)=-1;  % row 23

A(24,24)=4-(1-delta_y*h/k); %%%%%ROW 24 BOUNDRY 1 
A(24,18)=-1;A(24,23)=-1;

A(25,25)=3;A(25,19)=-1;A(25,26)=-1;A(25,30)=-1;   % row 25
A(26,20)=-1;A(26,25)=-1;A(26,27)=-1;A(26,31)=-1;  % row 26
A(27,21)=-1;A(27,26)=-1;A(27,28)=-1;A(27,32)=-1;  % row 27
A(28,22)=-1;A(28,27)=-1;A(28,29)=-1;A(28,33)=-1;  % row 28

A(29,29)=(4-2*(1-((sqrt(2)*delta_y*h)/(2*k)))); %%%%%ROW 29 BOUNDRY 2
A(29,28)=-1;A(29,23)=-1;

A(30,30)=3;A(30,25)=-1;A(30,31)=-1;               % row 30
A(31,26)=-1;A(31,30)=-1;A(31,32)=-1;              % row 31
A(32,27)=-1;A(32,31)=-1;A(32,33)=-1;              % row 32

A(33,33)=4-(1-delta_x*h/k); %%%%%ROW 33 BOUNDRY 3
A(33,32)=-1;A(33,28)=-1;


%MATRIX b-------------------------------------------------------
b=zeros(n,1);
b(6)=40;
b(12)=40;
b(18)=40;
b(24)=40+(delta_y*h*T_inf/k); %ROW 24
b(29)=sqrt(2)*h*delta_y*T_inf/k; %ROW 29
b(30)=70;
b(31)=70;
b(32)=70;
b(33)=70+(delta_x*h*T_inf/k); %ROW 33

%Find the numerical solution to Ax=b
x=A\b;

%TEMPERATURE VALUE MATRIX--------------------------------------------
T = zeros(6);
%Input boundary conditions
T(4:7,7)= 40;
T(1,1:4)= 70;
T(1:2,5:7)=nan;
T(3,6:7)=nan;
%Input calculated temperatures
T(7,1:6)=x(1:6)';
T(6,1:6)=x(7:12)';
T(5,1:6)=x(13:18)';
T(4,1:6)=x(19:24)';
T(3,1:5)=x(25:29)';
T(2,1:4)=x(30:33)';

%TO FIND AMBIENT TEMPERATURE RANGE
T_inf_max = 20;
T_inf_min = 20;
%to find max temp
x_amb = A\b;
while(x_amb(22)< 55)
    T_inf_max = T_inf_max + 0.01;
    
    b(24)=40+(delta_y*h*T_inf_max/k); %ROW 24
    b(29)=sqrt(2)*h*delta_y*T_inf_max/k; %ROW 29
    b(33)=70+(delta_x*h*T_inf_max/k); %ROW 33  
    
    x_amb = A\b;
end
T_inf_max = T_inf_max - 0.01;

%to find min temp
while(x_amb(22)> 50)
    T_inf_min = T_inf_min - 0.01;
    
    b(24)=40+(delta_y*h*T_inf_min/k); %ROW 24
    b(29)=sqrt(2)*h*delta_y*T_inf_min/k; %ROW 29
    b(33)=70+(delta_x*h*T_inf_min/k); %ROW 33  
    
    x_amb = A\b;
end

T_inf_min = T_inf_min + 0.01;

