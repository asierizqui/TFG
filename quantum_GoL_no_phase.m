n_row= 10;      % Dimensions of the 2-D lattice                     
n_column= 12;
steps=30;       % "steps" = number of steps,

% In the semiquantum GoL cell states are not a number, but a two-level vector: |\psi>=a\ket{1}+b\ket{0}=a\ket{alive}+b\ket{dead}
% Then, in the {|1>,|0>} basis, |\psi>=[a;b] : |a|^2 + |b|^2 =1.

for i= 1:n_row        % Fill the grid with death cells.
  for j = 1:n_column
    grid(i,j,1)=0;    % grid(i,j,1) = a (|1> state coef.)
    grid(i,j,2)=1;    % grid(i,j,2) = a (|0> state coef.)
  endfor
endfor

% Initial condition Example 1: Glider. Useful to check PBCs.
grid(4,2,1)=grid(5,2,1)=grid(5,4,1)=grid(6,2,1)=grid(6,3,1)=1;
grid(4,2,2)=grid(5,2,2)=grid(5,4,2)=grid(6,2,2)=grid(6,3,2)=0;

% Initial condition Example 2: "Tetromino" that becomes a "beehive" (stable figure)
%grid(4,4,1)=grid(5,4,1)=grid(5,5,1)=grid(5,6,1)=1;
%grid(4,4,2)=grid(5,4,2)=grid(5,5,2)=grid(5,6,2)=0;

% Initial condition Example 3: "Tetromino" that after the tenth generation becomes gour isolated (...)
  % (...) periodical blinkers. Each of them constitutes the so-called "traffic lights" periodical shape.
%grid(4,5,1)=grid(5,4,1)=grid(5,5,1)=grid(5,6,1)=1;
%grid(4,5,1)=grid(5,4,1)=grid(5,5,1)=grid(5,6,1)=1;


disp(['Grid, in the following  ', num2str(steps),'  steps'])

% --------------------------  INITIAL CONDITION DISPLAY  ---------------------------------                            
  colormap gray;       % Image traversing black to white in shades of gray.
  imagesc(grid(:,:,1));   % Only diplay |1> coef.) 
  axis equal tight     % To keep the square shape
  drawnow              
  pause(0.1)           % Pauses the GIF for (argument) seconds

% ----------------------------------  MAIN PROGRAM  -------------------------------------  
for t=1:steps
% grid at time t + 1 = "grid_next"; grid at the previous generation or time t = "grid".
 grid_next=grid;   % No-cloning theorem forbids this. Hence, synchronous application of GoL rules is NOT realizable in a QC!

 for i= 1:size(grid,1)    % size(grid,1)=n_row
   for j= 1:size(grid,2)  % size(grid,2)=n_column
      % x = counter of Moore neighbors's "health" = \sum_{i=1}^8 a_i = a_M/e^{i\phi}, with PBCs
      x= grid([mod((i+1)-1,size(grid,1))+1],j,1) + grid([mod((i-1)-1,size(grid,1))+1],j,1) ...         
        +grid(i,[mod((j+1)-1,size(grid,2))+1],1) + grid(i,[mod((j-1)-1,size(grid,2))+1],1) ...
        +grid([mod((i+1)-1,size(grid,1))+1],[mod((j+1)-1,size(grid,2))+1],1) ...
        +grid([mod((i+1)-1,size(grid,1))+1],[mod((j-1)-1,size(grid,2))+1],1) ...
        +grid([mod((i-1)-1,size(grid,1))+1],[mod((j+1)-1,size(grid,2))+1],1) ...
        +grid([mod((i-1)-1,size(grid,1))+1],[mod((j-1)-1,size(grid,2))+1],1);
    
     B=[1,1; 0,0]/[grid(i,j,1)+grid(i,j,2)]; % Birth operator
     D=[0,0; 1,1]/[grid(i,j,1)+grid(i,j,2)]; % Death operator
     S=[1,0; 0,1];                           % Survival op.
    
     if (0<=x) && (x<=1)
      G=D;
     endif
     if (1<x) && (x<=2)
      G=(sqrt(2)+1)*(2-x)*D+(x-1)*S;
     endif
     if (2<x) && (x<=3)
      G=(sqrt(2)+1)*(3-x)*S+(x-2)*B;
     endif
     if (3<x) && (x<=4)
      G=(sqrt(2)+1)*(4-x)*B+(x-3)*D;
     endif
     if (x>4)
      G=D;
     endif
     
     cell_next=G*[grid(i,j,1); grid(i,j,2)]; % Apply the rules.
     % Resulting grid, where cell_next=[cell_next(1,1); cell_next(2,1)]
     grid_next(i,j,1)=cell_next(1,1); grid_next(i,j,2)=cell_next(2,1);
        
   endfor % end the loop of column index "j".                 
 endfor  % end the loop of row index "i".

% ----------------------------------  DISPLAY -------------------------------------------                 
  colormap gray;       
  imagesc(grid(:,:,1));
  axis equal tight     
  drawnow                    
  pause(0.1)  

 grid=grid_next;   % Resulting grid becomes the reference one in the next time step.
                    
endfor  % end the cycle relative to the number of time steps ("t").