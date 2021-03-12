function [] = kuramoto_1d_podgen(p_method,L_range)
%---------------------------------------------------------------------
%  kuramoto_1d_podgen.m - computes POD basis using different 
%                         methodologies
%
%  Version: 1.0
%
%  Variables:
%    p_method = 
%               0 - implements the standard POD basis construction
%                   for one value of the parameter (Re_range = value)
%               1 - standard snapshot stacking approach
%
%               2 - optimization interpretation
%
%               3 - generate using PID ideas
%
%    L_range =
%               for p_method = 0,
%                 integer (must have the file snapshots_at_L####.mat)
%               for p_method > 0,
%                 an array of lengths, eg., 
%                 L_range = [10:1:50]
%---------------------------------------------------------------------

  n_modes = 20;
  
%---------------------------------------------------------------------
%  Read, or generate the mass matrix for the finite element data
%
%  The points and connectivity should match those used in the
%  simulations
%---------------------------------------------------------------------
  L = L_range(1);
  if L<10
    number = strcat('0',int2str(100*L))
  else
    number = int2str(100*L);
  end
  
  solution_file = strcat('snapshots_at_L',number);
  load( solution_file, 'w_save' )
  
  [n_nodes,k_steps] = size(w_save);

  n_elements = (n_nodes-1);

  %  Generate node coordinates
  dx = 1/n_elements;
  for i=1:n_nodes
    x(i,1) = dx*(i-1);
  end

  %  Generate element connectivity
  ie = 0;
  for i=1:n_elements
    ie = ie + 1;
    e_conn(ie,:) = [i, i+1];
  end

  [M] = build_mass_1d(x, e_conn);
  [R] = chol(M);

%---------------------------------------------------------------------
%  P_METHOD == 0  CASE   (pod for one Re)
%---------------------------------------------------------------------
if (p_method==0)
  % read in snapshots
  L = L_range(1);
  if L<10
    number = strcat('0',int2str(100*L))
  else
    number = int2str(100*L);
  end
  
  solution_file = strcat('snapshots_at_L',number);
  load( solution_file, 'w_save' )

  % assemble the weighted snapshot matrix
  SN  = R*w_save;
  
  [U,S,V] = svd(SN,0);
  
  POD = R\U(:,1:n_modes);
end

%---------------------------------------------------------------------
%  P_METHOD == 1  CASE    (concatenate all data)
%---------------------------------------------------------------------
if (p_method==1)
  % read in the first set of snapshots to get sizes
  L = L_range(1)
  
  if L<10
    number = strcat('0',int2str(100*L));
  else
    number = int2str(100*L);
  end
  
  solution_file = strcat('snapshots_at_L',number);
  load( solution_file, 'w_save' )
  
  [n_states,n_snapshots] = size(w_save);
  n_ensemble = length(L_range);
  
  SN = zeros(n_states,n_ensemble*n_snapshots);
  SN(:,1:n_snapshots) = R*w_save;
  
  % read in the rest of the snapshots
  for i=2:n_ensemble
    L = L_range(i);
    
    if L<10
      number = strcat('0',int2str(100*L));
    else
      number = int2str(100*L);
    end

    solution_file = strcat('snapshots_at_L',number);
    load( solution_file, 'w_save' )

    SN(:,(i-1)*n_snapshots+1:i*n_snapshots) = R*w_save;
  end
  
  [U,S,V] = svd(SN,0);
  
  POD = R\U(:,1:n_modes);
end

%---------------------------------------------------------------------
%  P_METHOD == 2  CASE    (optimization interpretation)
%---------------------------------------------------------------------
if (p_method==2)
  % read in the first set of snapshots to get sizes
  L = L_range(1);
  solution_file = strcat('snapshots_at_L',int2str(100*L));
  load( solution_file, 'w_save' )
  
  [n_states,n_snapshots] = size(w_save);
  n_ensemble = length(L_range);
  
  SN = R*w_save * (L_range(2)-L_range(1))/2;
  
  % read in the rest of the snapshots
  for i = 2:n_ensemble-1
    L = L_range(i);
    solution_file = strcat('snapshots_at_L',int2str(100*L));
    load( solution_file, 'w_save' )

    SN = SN + R*w_save * (L_range(i+1)-L_range(i-1))/2;
  end
  
  % read in the last snapshot set
  L = L_range(end);
  solution_file = strcat('snapshots_at_L',int2str(100*L));
  load( solution_file, 'w_save' )

  SN = SN + R*w_save * (L_range(end)-L_range(end-1))/2;
  
  [U,S,V] = svd(SN,0);
  
  POD = R\U(:,1:n_modes);
end


%---------------------------------------------------------------------
%  P_METHOD == 3  CASE    (PID in parameter space)
%---------------------------------------------------------------------
if (p_method==3)

end


%---------------------------------------------------------------------
%  Plot/Output the POD Basis
%---------------------------------------------------------------------
  save kuramoto_1d_podgen POD L_range p_method
  
  % plot the first six POD modes
  figure(5)
  for i=1:10
    subplot(5,2,i)
    plot(x,POD(:,i))
    axis([0 1 -1.7 1.7])
    title_string = strcat('POD Mode ',int2str(i));
    title(title_string)
  end

  figure(6)
  semilogy(diag(S),'*')
