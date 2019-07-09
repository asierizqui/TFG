% Definition of useful expressions:
ket_0=[1 0]; % = |0>,   in the {|0>,|1>} computational basis.
ket_1=[0 1]; % = |1>
rho_A= kron(ket_0,ket_0'); % = |0>\otimes<0|= |0><0|,  where kron(A,B)= A \otimes B, and  ket_0'= (|0>)^\dagger = <0|.
sigma_x=[0 1; 1 0];
sigma_z=[1 0; 0 -1];
% CNOT = |0><0|\otimes 1_2 + |1><1|\otimes sigma_x :
CNOT= kron( kron(ket_0,ket_0'), eye(2) ) +  kron( kron(ket_1,ket_1'), sigma_x );

% Variable parameters.
n_indiv=3;    % Number of individuals.
epsilon=0.3; % , for instance. <sigma_z>_p(t_d) = 1-epsilon = "death age",  where t_d = t_death.
harm=0.1; % Harm= %10, e.g., of <H_int>_genotypes*prey's life.


alive_labels=[];    
for i=1:n_indiv;
 alive_labels(i)=[i];
endfor;


%  -------------------Initial "age" of n individuals ---------------------------- 
% Each initial (t=0) individual's pheno. will be described by a density matrix of the form, for instance, of Eq. (4.7): 
for i=1:n_indiv   % per individual
  a=(-1:0.1:(2-epsilon)/2)([randi([1 size(-1:0.1:(2-epsilon)/2, 2)])]);   % Eq. (4.15b) = 2a-1 = 1-epsilon -> to obtain the maximum value of a.
  % where, (-1:0.1:1)=[-1 -0.9 ... 0.9 1];     randi([min max])= random integer in [min,max];     size(A,2) = number of columns of the matrix A.
  % That is, we randomly choose an element. For example, (-1:0.1:1)(2)=-0.9;  (-1:0.1:1)(3)=-0.8;   (-1:0.1:1)(4)=-0.7;   ...
  b=0.7;  % For example.  For the moment where are interested just in a, because lifetime depends on it. See Eq. (4.15b)
  c=0.5;  % For example.
  
  % 3D matrix = matrix(row, column, page), the i th page, for the i th individual. On each page, there is a conventional 2D matrix.
  rho_g0(:,:,i)=[a, b-j*c; b+j*c, 1-a];    % j = imaginary number i.      Eq. (4.7)
  rho_p0(:,:,i)= CNOT*kron(rho_g0(:,:,i),rho_A)*CNOT';   % Eq. (4.8). It would be rho_g1(:,:,i) as well, but we are not implementing self-replication, for now.
  expectation_sigma_z_p0(i)=trace(rho_p0(:,:,i)*kron(eye(2),sigma_z));  % Eq. (4.5(2)). Instead, for the rho_g1 subspace it would be kron(sigma_z,eye(2))).
  expectation_sigma_z_g0(i)=expectation_sigma_z_p0(i);  % At t=0, <sigma_z>_p0=<sigma_z>_g0.
endfor


% For the interaction Hamiltonian H_int:
for i=1:n_indiv     
 a_k(i)=(-1:0.01:1)([randi([1 size(-1:0.01:1, 2)])]);   % a_k coef. in Eq. (4.24).
 b_k(i)=(-1:0.01:1)([randi([1 size(-1:0.01:1, 2)])]);   % b_k  "    "     "      .
endfor



 % Dimensions of the 2-D lattice.       
n_row= 3;             
n_column=4;
% Fill the grid with 0s. 0 = no indiv.
for i=1:n_row;
  for j=1:n_column;
     grid(i,j)=0;
  endfor
endfor


% Individuals get born in a random location of the grid.
% x_pos= [x position of the 1 st indiv., ...., x position of the n th indiv.)]
for i=1:n_indiv
   y_pos(i)=randi(n_row);      % randi(max) = random integer from 1 to max.
   x_pos(i)=randi(n_column);  
  while grid(y_pos(i),x_pos(i))!=0   % 2 indiv. or more cannot be born in the same cell.
   y_pos(i)=randi(n_row); x_pos(i)=randi(n_column);    
  endwhile  
  grid(y_pos(i),x_pos(i))=i; % We represent the i th individual in the grid by the integer i.
endfor 

disp(['________________________________GRID at t=0______________________________________']) %disp() = display string.
disp(grid); disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); disp(' ');



for t=0:0.1:1   % per time step
  disp(['-------------------------------------At t=',num2str(t),'--------------------------------------'])
  disp(' '); disp(' ')
  
  
  % CHECK if there is any DEAD individual.
  dead_labels=[];  % Individuals that die at this t step. 
  for i=1:size(alive_labels,2)  % per alive individual.
   if (expectation_sigma_z_p0(alive_labels(i))>= 1-epsilon)
     disp(['Individual ', num2str(alive_labels(i)),' has DIED!    (Death age=',num2str(1-epsilon),')'])
     grid(y_pos(alive_labels(i)),x_pos(alive_labels(i)))=0; % The indiv. dissapears from the grid.
     disp(['Then, individual ', num2str(alive_labels(i)),' disappears from the grid:']); disp(grid); disp(' ')
     % To avoid problems ("index out of bounds") erasing labels from "alive_labels":
     dead_labels=cat(2,dead_labels,alive_labels(i)); % cat(dim,A,B) concatenates B to the end of A along dimension dim. dim=2=column. 
   endif
  endfor  
  % Remove common elements from arrays "alive_labels" and "dead_labels":
  alive_labels=setdiff(alive_labels,dead_labels);  % setdiff(A,B) returns the data in A that is not in B.
  
  if size(alive_labels,1)==0%=size([],1), where [] is an empty array.
    disp('ALL individuals are dead, program ENDS!')
    break % Program ends.
  endif
  
  disp(' ')
  disp(['____________________________________AGE___________________________________'])
  eta=1-exp(-t); % Time in \lambda*t units: \lambda*t=t
  M_0=[1 0; 0 sqrt(1-eta)]; M_1=[0 sqrt(eta); 0 0];  % Eq. (B.22)
  for i=1:size(alive_labels,2)  % per alive individual.
    % Careful,  alive_labels(i) is the label of the alive individual, not i!
    rho_p0(:,:,alive_labels(i))=kron(eye(2),M_0)*rho_p0(:,:,alive_labels(i))*kron(eye(2),M_0)'+ kron(eye(2),M_1)*rho_p0(:,:,alive_labels(i))*kron(eye(2),M_1)';  % 1_2 otimes M_1 to act on pheno. subspace.
    expectation_sigma_z_p0(alive_labels(i))=trace(rho_p0(:,:,alive_labels(i))*kron(eye(2),sigma_z));
    disp(['     ',num2str(alive_labels(i)),':  <sigma_z>_p0=',num2str(expectation_sigma_z_p0(alive_labels(i)))]) 
  endfor     
  disp(' '); disp(' ');

  
  disp(['____________________________Random MOVEMENT_______________________________'])
  for i=1:size(alive_labels,2)  % per alive individual.
    if t==0
      disp('     NO movement at t=0.')
      break %  At t=0 individuals interact, but do not move until the next t step.
    endif
    x=randi([-1 1]); % movement in x axis.
    y=randi([-1 1]); %    "     "  y  "  .
    disp(['If individual ', num2str(alive_labels(i)),' moves in direction (x,y)=(',num2str(x),',',num2str(y),'):'])
    if (x!=0) || (y!=0) % = if the individual moves:
      % There is a - in ...(y_pos(alive_labels(i))-y)..., because row index y_pos() increases towards the negative direction of y:
      if (grid(mod((y_pos(alive_labels(i))-y)-1, n_row)+1, mod((x_pos(i)+x)-1,n_column)+1)!=0) % Cannot move to a cell where there already is an individual.
        disp(['Individual ', num2str(alive_labels(i)),' canNOT move there, that cell is already occupied by individual ', num2str(grid(mod((y_pos(alive_labels(i))-y)-1, n_row)+1, mod((x_pos(alive_labels(i))+x)-1,n_column)+1)),'.'])
      else
        grid(mod((y_pos(alive_labels(i))-y)-1,n_row)+1,mod((x_pos(alive_labels(i))+x)-1,n_column)+1)=grid(y_pos(alive_labels(i)),x_pos(alive_labels(i)));
        grid(y_pos(alive_labels(i)),x_pos(alive_labels(i)))=0; % The indiv. dissapears from the previous position.
        y_pos(alive_labels(i))=mod((y_pos(alive_labels(i))-y)-1, n_row)+1; % new position's row index.
        x_pos(alive_labels(i))=mod((x_pos(alive_labels(i))+x)-1, n_column)+1; % new position's column index.
      endif
    endif
    disp(grid)
  endfor      
  disp(' '); disp(' ')
  
  disp(['____________________________ INTERACTION__________________________________'])
  disp('(<sigma_z>_p(t=0)=<sigma_z>_g=trophic role)'); disp(' ');
  pred_prey_total=[];  % = [predator_label, prey_label], to clean later the repeated interactions. We will add a row per interaction.
  for i=1:size(alive_labels,2)  % per alive individual.
    for x_interac=-1:1    % Interaction range (Moore) in the x axis.
      for y_interac=-1:1  %     "         "     "      "  "   y   " .  
        if (x_interac!=0) || (y_interac!=0)  % Not to include the current indiv., where (x_interac,y_interac)=(0,0).
         %disp(['The cell at the (x,y)=(',num2str(x_interac),',',num2str(y_interac),') part of the ', num2str(i),' th indiv is: ', num2str(grid(mod((y_pos(i)-y_interac)-1, n_row)+1, mod((x_pos(i)+x_interac)-1,n_column)+1))])
          % There is a - in ...(y_pos(alive_labels(i))-y_interac)..., because row index y_pos() increases towards the negative direction of y:
          if (grid(mod((y_pos(alive_labels(i))-y_interac)-1, n_row)+1, mod((x_pos(alive_labels(i))+x_interac)-1,n_column)+1)!=0)
            interac_indiv=grid(mod((y_pos(alive_labels(i))-y_interac)-1, n_row)+1, mod((x_pos(alive_labels(i))+x_interac)-1,n_column)+1); % label of the individual that interacts with the i th one.
            disp(['Individual ', num2str(alive_labels(i)),' interacts with individual ', num2str(interac_indiv),'.'])
            disp(['     ',num2str(alive_labels(i)),':  <sigma_z>_p0(t=0)=',num2str(expectation_sigma_z_g0(alive_labels(i)))]) 
            disp(['     ', num2str(interac_indiv),':  <sigma_z>_p0(t=0)=',num2str(expectation_sigma_z_g0(interac_indiv))])  
            if ( expectation_sigma_z_g0(alive_labels(i)) >  expectation_sigma_z_g0(interac_indiv)  )
              disp(['     Thus, here: ', num2str(alive_labels(i)),'=predator; ', num2str(interac_indiv),'=prey.']); disp(' ')
              pred_prey_total=[pred_prey_total; [alive_labels(i), interac_indiv]];   
            elseif  ( expectation_sigma_z_g0(alive_labels(i)) <  expectation_sigma_z_g0(interac_indiv)  )
              disp(['     Thus, here: ', num2str(interac_indiv),'=predator; ', num2str(alive_labels(i)),'=prey.']); disp(' ')
              pred_prey_total=[pred_prey_total; [interac_indiv, alive_labels(i)]];  
            else % same level in the trophic pyramid. 
              disp(['     Thus, here, both are predators or both preys, they do NOT interact!.']); disp(' ')
              % Thus, here, we do not append a row to the matrix pred_prey_total. 
            endif
          endif
        endif
      endfor
    endfor
  endfor     
  
  % Clean the repeated interactions.
  pred_prey_not_repeated=unique(pred_prey_total,'rows');    % Return the unique rows.
  disp(' '); disp(' '); disp(' '); disp('EFFECTS:'); disp(' ');
  for i=1:size(pred_prey_not_repeated,1)  % Each row is an (unrepeated per t step) interaction. 
    disp(['If  ',num2str(pred_prey_not_repeated(i,1)),'=predator; ',num2str(pred_prey_not_repeated(i,2)),'=prey:'])
     
    % Interaction Hamiltonian H_int
    H_pred=a_k(pred_prey_not_repeated(i,1))*kron(sigma_z,sigma_z)+b_k(pred_prey_not_repeated(i,1))*kron(sigma_z,eye(2));
    H_prey=a_k(pred_prey_not_repeated(i,2))*kron(sigma_z,sigma_z)+b_k(pred_prey_not_repeated(i,2))*kron(eye(2),sigma_z);
    H_int=H_pred+H_prey;
    % <H_int>_kron(rho_g_pred,rho_g_prey):
    expectation_H_int_genotypes=trace( kron(rho_g0(:,:,pred_prey_not_repeated(i,1)),rho_g0(:,:,pred_prey_not_repeated(i,2)) )  * H_int);
    disp(['     <H_int>_genotypes = ',num2str(expectation_H_int_genotypes)]) 
    if (expectation_H_int_genotypes>0) % The predator HUNTs the prey. Harm is proportional to <H_int>_genotypes and prey's life:
     % Predator's age decreases, beneficial for it.
     expectation_sigma_z_p0(pred_prey_not_repeated(i,1))= expectation_sigma_z_p0(pred_prey_not_repeated(i,1)) - harm*expectation_H_int_genotypes* abs(expectation_sigma_z_p0(pred_prey_not_repeated(i,2)));
     % Prey's age increases by the same amount (feeding) predator's age decreases. Harmful for the prey.
     expectation_sigma_z_p0(pred_prey_not_repeated(i,2))= expectation_sigma_z_p0(pred_prey_not_repeated(i,2)) + harm*expectation_H_int_genotypes* abs(expectation_sigma_z_p0(pred_prey_not_repeated(i,2)));
     disp(['     The predator (',num2str(pred_prey_not_repeated(i,1)),') HUNTS the prey (',num2str(pred_prey_not_repeated(i,2)),'):      <sigma_z>_p0^pred=',num2str(expectation_sigma_z_p0(pred_prey_not_repeated(i,1))),';   <sigma_z>_p0^prey=',num2str(expectation_sigma_z_p0(pred_prey_not_repeated(i,2))),'.'])
    else % prey escapes 
       disp(['     The prey (',num2str(pred_prey_not_repeated(i,2)),') ESCAPES from the predator (',num2str(pred_prey_not_repeated(i,1)),'), so it does NOT get any harm!'])
    endif
    disp(' ')
  endfor 
  disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); disp(' '); 
  
endfor   % End of time t iterations.