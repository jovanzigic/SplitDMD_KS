  function [] = kuramoto_1db(epsilon)
%-------------------------------------------------------------------------------
%  kuramoto_1db.m - solves the 1D Kuramoto-Sivashinsky equation in a
%                   unit periodic domain.  A scaling parameter, epsilon,
%                   has been introduced that is related to a domain 
%                   length:   epsilon = 1/L^2
%
%  Copyright (c) 2007, Jeff Borggaard, Virginia Tech
%  Version: 1.0
%
%  Usage:    [] = kuramoto_1db(epsilon)
%
%  Variables:  epsilon  -  scaling parameter
%                          (epsilon = 0.00633257, bifurcation)
%                          (epsilon = 0.00603102, heteroclinic bifurc.)
%                          (epsilon = 0.00579148, hopf bifurcation)
%                          (epsilon = 0.00000618, "chaos")
%% -----------------------------------------------------------------------------
  addpath('/Volumes/borggaard/Software/FEM/fem_functions')
  
  %-----------------------------------------------------------------------------
  %  Define "Input" Parameters
  %-----------------------------------------------------------------------------
%   if (nargin==0)
%     epsilon = 0.00701;
%   end
  
  n_quadrature     = 5;        % number of quadrature points

  n_nodes          = 161;      % spatial discretization information

  t_initial        = 0;        % time step information for heat equation
  t_step           = 0.02;
  t_save           = 0.02; %change to 0.1 for 4k, 0.02 for 400;
  t_final          = 400 ; %120.0;    % 400

  resid_tol        = 1e-8;     % Newton solver parameters
  step_tol         = 1e-8;
  max_iterations   = 40;
  
  %-----------------------------------------------------------------------------
  %  Geometry Module
  %-----------------------------------------------------------------------------
  n_elements = (n_nodes-1);

  %  Generate node coordinates
  dx = 1/n_elements;
  x = zeros(n_nodes,1);
  for i=1:n_nodes
    x(i,1) = dx*(i-1);
  end

  %  Generate element connectivity
  e_conn = zeros(n_elements,2);
  ie = 0;
  for i=1:n_elements
    ie = ie + 1;
    e_conn(ie,:) = [i, i+1];
  end

  %-----------------------------------------------------------------------------
  % Determine equation numbers and set up boundary condition information
  %-----------------------------------------------------------------------------
  ide = zeros(n_nodes,2);
  for n_nd=1:n_nodes-1    % this numbering is explicitly used elsewhere!
    ide(n_nd,1) = 2*n_nd-1;
    ide(n_nd,2) = 2*n_nd;
  end
  
  ide(n_nodes,1) = 1;
  ide(n_nodes,2) = 2;
  
  n_equations = 2*(n_nodes-1);  

  %-----------------------------------------------------------------------------
  %  Solver Module
  %-----------------------------------------------------------------------------
 
  %  Set Initial Conditions   (should do Galerkin projection)
  %-----------------------------------------------------------------------------
  w_current  =      sin(4*pi*x)/sqrt(epsilon) + 0.05*rand(size(x)); 
  wp_current = 4*pi*cos(4*pi*x)/sqrt(epsilon); 
  
%   w_current  =    400* sin(4*pi*x) + 0.05*rand(size(x)); 
%   wp_current = 400*4*pi*cos(4*pi*x); 

% special initial conditions for the L=400 case
%   w_current  =      sin(6*pi*x) + 0.05*rand(size(x)); 
%   wp_current = 6*pi*cos(6*pi*x); 

%   plot_hermite(2,x,e_conn,w_current,wp_current)    
% 
%   pause(0.001)

  %  Create storage
  n_steps = floor( (t_final-t_initial)/t_save + 0.1 );
  w_save  = zeros(length(x),n_steps);  % changed Dzeros to zeros
  wp_save = w_save;
  
  n_store = 1;
  w_save (:,n_store) = w_current;
  wp_save(:,n_store) = wp_current;

  %  Time Integration
  %-----------------------------------------------------------------------------
  for time = t_initial+t_step:t_step:t_final
    
    %---------------------------------------------------------------------------
    %  Newton Iterations
    %---------------------------------------------------------------------------
    w  = w_current;  % set initial guess for the iteration
    wp = wp_current;
        
    iteration = 0;
    converged = 0;
    
    while (~converged)
      iteration = iteration + 1;
      
      %-------------------------------------------------------------------------
      %  Compute Newton Step
      %-------------------------------------------------------------------------
      [r,J] = weak_resjac(x,e_conn,w,wp,w_current,wp_current,t_step,...
                          epsilon,ide,n_equations,n_elements,n_quadrature);

%       r0 = r;
%       J0 = J;
      
%       wp(1) = wp(1) + 0.00001;
%       wp(end) = wp(end) + 0.00001;

%       wp(2) = wp(2) + 0.00001;
%       [r,J] = weak_resjac(x,e_conn,w,wp,w_current,wp_current,t_step,...
%                           epsilon,ide,n_equations,n_elements,n_quadrature);
% 
%       sparse((r-r0)/0.00001)
%       J(:,4)
%       stop
      newton_step =-J\r;
      
      if ( iteration<5 && time< t_initial+3*t_step)
        lambda = 0.5;  % 0.5
      else
        lambda = 1;
      end

      %  Implicit update
      w (1:end-1,1) = w (1:end-1,1) + lambda*newton_step(1:2:end-1);
      wp(1:end-1,1) = wp(1:end-1,1) + lambda*newton_step(2:2:end  );
      
      %  Periodic boundary conditions
      w (n_nodes) = w (1);
      wp(n_nodes) = wp(1);
      
      norm_step   = norm(newton_step,2);
      norm_w      = norm(w,2);     % no wp for now!
      norm_resid  = norm(r);

%       fprintf('   norm_w = %g,  norm_r = %g \n',norm_w,norm_resid)
      
      if ( iteration==1 )
        norm_resid0 = max(norm_resid,1e-20);
      end

      converged = ( (norm_resid/norm_resid0 < resid_tol | ...
                     norm_resid < resid_tol           ) & ...
                     norm_step /norm_w      < step_tol );

      converged = converged | (iteration>max_iterations);
      if ( iteration>max_iterations )
        save kuramotos_1db_snap w_save wp_save epsilon x t_step
        error('Did not converge at time %g\n',time)
      end

    end
    
    % At this point, w_current is w(t-dt) and w is w(t)
                     
    w_current  = w;
    wp_current = wp;
                                              
    % Store solution in snapshot matrix
    if ( abs(mod(time-t_initial,t_save))<1e-12 )
      n_store = n_store+1;  % increment this for the next pass
      w_save (:,n_store) = w;
      wp_save(:,n_store) = wp;
      
      fprintf('\n At time=%g\n',time);
      
    end

  end % time iteration
    
  time = t_initial:t_save:t_final;
  save kuramoto_1db_snap w_save wp_save epsilon time x e_conn t_step
end
  
%% ---------------------------------------------------------------------                      
function [r,J] = weak_resjac(x,e_conn,w,wp,w_current,wp_current,t_step,...
                             epsilon,ide,n_equations,n_elements,n_quadrature)

%% ---------------------------------------------------------------------                      
  delta_w = 1e-7;         % stepsize for finite difference Jacobian
  
  theta   = 0.5;
  omt     = 1-theta;
  
  [rr,wt] = oned_quadrature(n_quadrature);
  one     = ones(size(rr));
  
  II      = zeros(10*10*n_elements,1);
  JJ      = zeros(10*10*n_elements,1);
  XX      = zeros(10*10*n_elements,1);  % changed Dzeros to zeros
  ntriplets = 0;

  % Dzeros changes to zeros
  r = zeros (n_equations,          1);
  
  for n_el=1:n_elements
    % compute value of each test function and their spatial derivaties
    % at the integration points  (1) linear, 3 point rule
    nodes_local           = e_conn(n_el,:);
    x_local               = x(nodes_local,:);
    [x_g,wt_g,phi0,phi1,p0_x,p1_x,p0_xx,p1_xx] = ...
                                           oned_shapeherm(x_local,rr,wt);
      
    % compute the value of each variable at the quadrature points
    w_local   = w(nodes_local,:);
    wp_local  = wp(nodes_local,:);

    wc_local  = w_current(nodes_local,:);    
    wpc_local = wp_current(nodes_local,:);
    
    w_g       = phi0 *w_local  + phi1 *wp_local;
    wx_g      = p0_x *w_local  + p1_x *wp_local;
    wxx_g     = p0_xx*w_local  + p1_xx*wp_local;
    
    wc_g      = phi0 *wc_local + phi1 *wpc_local;
    wcx_g     = p0_x *wc_local + p1_x *wpc_local;
    wcxx_g    = p0_xx*wc_local + p1_xx*wpc_local;
  
    %----------------------------------------------------------------------
    %  Integrate the weak form of the equations
    %----------------------------------------------------------------------
    r0_loc =          oned_f_int( (w_g-wc_g)/t_step       , phi0 , wt_g) ...
          - epsilon  *oned_f_int( theta*wx_g + omt*wcx_g  , p0_x , wt_g) ...
          + epsilon^2*oned_f_int( theta*wxx_g + omt*wcxx_g, p0_xx, wt_g) ...
          + epsilon*2*oned_f_int( theta*w_g.*wx_g                      ...
                                  + omt*wc_g.*wcx_g       , phi0 , wt_g);

    r1_loc =          oned_f_int( (w_g-wc_g)/t_step       , phi1 , wt_g) ...
          - epsilon  *oned_f_int( theta*wx_g + omt*wcx_g  , p1_x , wt_g) ...
          + epsilon^2*oned_f_int( theta*wxx_g + omt*wcxx_g, p1_xx, wt_g) ...
          + epsilon*2*oned_f_int( theta*w_g.*wx_g                      ...
                                  + omt*wc_g.*wcx_g       , phi1 , wt_g);

    %-------------------------------------------------------------------
    %  Assemble the global residual
    %-------------------------------------------------------------------
    unk_0    = [ide(n_el,1); ide(n_el+1,1)];
    r(unk_0) = r(unk_0) + r0_loc;

    unk_1    = [ide(n_el,2); ide(n_el+1,2)];
    r(unk_1) = r(unk_1) + r1_loc;
      
    % calculate the columns of J corresponding to the local unknowns
    % increment phi0 for the first node
    wplus_local   = w_local + [delta_w; 0];
        
    wplus_g       = phi0 *wplus_local  + phi1 *wp_local;
    wplusx_g      = p0_x *wplus_local  + p1_x *wp_local;
    wplusxx_g     = p0_xx*wplus_local  + p1_xx*wp_local;
 
    r0plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p0_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi0 , wt_g);

    r1plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p1_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi1 , wt_g);
      
  %  J(unk_0, unk_0(1)) = J(unk_0, unk_0(1)) + (r0plus_loc - r0_loc)/delta_w;
    ntriplets = ntriplets + 1;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_0(1);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;
    
    ntriplets = ntriplets + 2;
  %  J(unk_1, unk_0(1)) = J(unk_1, unk_0(1)) + (r1plus_loc - r1_loc)/delta_w;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_0(1);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;
  
    % increment phi0 for the second node
    wplus_local   = w_local + [0; delta_w];
        
    wplus_g       = phi0 *wplus_local  + phi1 *wp_local;
    wplusx_g      = p0_x *wplus_local  + p1_x *wp_local;
    wplusxx_g     = p0_xx*wplus_local  + p1_xx*wp_local;
 
    r0plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p0_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi0 , wt_g);

    r1plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p1_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi1 , wt_g);

    ntriplets = ntriplets + 2;
%    J(unk_0, unk_0(2)) = J(unk_0, unk_0(2)) + (r0plus_loc - r0_loc)/delta_w;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_0(2);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;
    
    ntriplets = ntriplets + 2;    
%    J(unk_1, unk_0(2)) = J(unk_1, unk_0(2)) + (r1plus_loc - r1_loc)/delta_w;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_0(2);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;

    % increment phi1 for the first node
    wpplus_local   = wp_local + [delta_w; 0];
        
    wplus_g       = phi0 *w_local  + phi1 *wpplus_local;
    wplusx_g      = p0_x *w_local  + p1_x *wpplus_local;
    wplusxx_g     = p0_xx*w_local  + p1_xx*wpplus_local;
 
    r0plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p0_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi0 , wt_g);

    r1plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p1_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi1 , wt_g);
    
    ntriplets = ntriplets + 2;
%    J(unk_0, unk_1(1)) = J(unk_0, unk_1(1)) + (r0plus_loc - r0_loc)/delta_w;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_1(1);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;

    ntriplets = ntriplets + 2;
%    J(unk_1, unk_1(1)) = J(unk_1, unk_1(1)) + (r1plus_loc - r1_loc)/delta_w;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_1(1);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;

    % increment phi1 for the second node
    wpplus_local   = wp_local + [0; delta_w];
        
    wplus_g       = phi0 *w_local  + phi1 *wpplus_local;
    wplusx_g      = p0_x *w_local  + p1_x *wpplus_local;
    wplusxx_g     = p0_xx*w_local  + p1_xx*wpplus_local;
 
    r0plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p0_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi0 , wt_g);

    r1plus_loc =     oned_f_int( (wplus_g-wc_g)/t_step   , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*wplusx_g + omt*wcx_g  , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*wplusxx_g + omt*wcxx_g, p1_xx, wt_g) ...
     + epsilon*2*oned_f_int( theta*wplus_g.*wplusx_g                    ...
                                  + omt*wc_g.*wcx_g      , phi1 , wt_g);
     
    ntriplets = ntriplets + 2;
%    J(unk_0, unk_1(2)) = J(unk_0, unk_1(2)) + (r0plus_loc - r0_loc)/delta_w;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_1(2);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;
    
    ntriplets = ntriplets + 2;
%    J(unk_1, unk_1(2)) = J(unk_1, unk_1(2)) + (r1plus_loc - r1_loc)/delta_w;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_1(2);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;
    ntriplets = ntriplets + 1;

  end  %  n_el=1:n_elements

  J = sparse( II(1:ntriplets), JJ(1:ntriplets), XX(1:ntriplets), ...
              n_equations, n_equations);

end % function weak_resjac

%-----------------------------------------------------------------------                      
function [r,J] = weak_sens(x,e_conn,s,sp,s_current,sp_current,...
                                    w,wp,w_current,wp_current,...
                           t_step,epsilon,ide,n_equations,n_elements,n_quadrature)

  delta_w = 1e-7;         % stepsize for finite difference Jacobian
  
  theta   = 0.5;
  omt     = 1-theta;
  
  [rr,wt] = oned_quadrature(n_quadrature);
  one     = ones(size(rr));
  
  II      = zeros(10*10*n_elements,1);
  JJ      = zeros(10*10*n_elements,1);
  XX      = zeros(10*10*n_elements,1); % changed Dzeros to zeros
  ntriplets = 0;

  r = zeros (n_equations,          1);  % changed Dzeros to zeros
  
  for n_el=1:n_elements
    % compute value of each test function and their spatial derivaties
    % at the integration points  (1) linear, 3 point rule
    nodes_local           = e_conn(n_el,:);
    x_local               = x(nodes_local,:);
    [x_g,wt_g,phi0,phi1,p0_x,p1_x,p0_xx,p1_xx] = ...
                                           oned_shapeherm(x_local,rr,wt);
      
    % compute the value of each variable at the quadrature points
    w_local   = w (nodes_local,:);
    wp_local  = wp(nodes_local,:);

    wc_local  = w_current (nodes_local,:);    
    wpc_local = wp_current(nodes_local,:);

    s_local   = s (nodes_local,:);
    sp_local  = sp(nodes_local,:);

    sc_local  = s_current (nodes_local,:);    
    spc_local = sp_current(nodes_local,:);
    
    w_g       = phi0 *w_local  + phi1 *wp_local;
    wx_g      = p0_x *w_local  + p1_x *wp_local;
    wxx_g     = p0_xx*w_local  + p1_xx*wp_local;
    
    wc_g      = phi0 *wc_local + phi1 *wpc_local;
    wcx_g     = p0_x *wc_local + p1_x *wpc_local;
    wcxx_g    = p0_xx*wc_local + p1_xx*wpc_local;
    
    s_g       = phi0 *s_local  + phi1 *sp_local;
    sx_g      = p0_x *s_local  + p1_x *sp_local;
    sxx_g     = p0_xx*s_local  + p1_xx*sp_local;
    
    sc_g      = phi0 *sc_local + phi1 *spc_local;
    scx_g     = p0_x *sc_local + p1_x *spc_local;
    scxx_g    = p0_xx*sc_local + p1_xx*spc_local;
  
    %----------------------------------------------------------------------
    %  Integrate the weak form of the equations
    %----------------------------------------------------------------------
    r0_loc =          oned_f_int( (s_g-sc_g)/t_step       , phi0 , wt_g) ...
          - epsilon  *oned_f_int( theta*sx_g + omt*scx_g  , p0_x , wt_g) ...
          -           oned_f_int( theta*wx_g + omt*wcx_g  , p0_x , wt_g) ...
          + epsilon^2*oned_f_int( theta*sxx_g + omt*scxx_g, p0_xx, wt_g) ...
          + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g, p0_xx, wt_g) ...
          + 2*epsilon*oned_f_int( theta*(s_g.*wx_g + w_g.*sx_g)          ...
                        + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi0 , wt_g) ...
          +         2*oned_f_int( theta*w_g.*wx_g                        ...
                        + omt*wc_g.*wcx_g                 , phi0 , wt_g);

    r1_loc =          oned_f_int( (s_g-sc_g)/t_step       , phi1 , wt_g) ...
          - epsilon  *oned_f_int( theta*sx_g + omt*scx_g  , p1_x , wt_g) ...
          -           oned_f_int( theta*wx_g + omt*wcx_g  , p1_x , wt_g) ...
          + epsilon^2*oned_f_int( theta*sxx_g + omt*scxx_g, p1_xx, wt_g) ...
          + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g, p1_xx, wt_g) ...
          + 2*epsilon*oned_f_int( theta*(s_g.*wx_g + w_g.*sx_g)          ...
                        + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi1 , wt_g) ...
          +         2*oned_f_int( theta*w_g.*wx_g                        ...
                        + omt*wc_g.*wcx_g                 , phi1 , wt_g);


    %-------------------------------------------------------------------
    %  Assemble the global residual
    %-------------------------------------------------------------------
    unk_0    = [ide(n_el,1); ide(n_el+1,1)];
    r(unk_0) = r(unk_0) + r0_loc;

    unk_1    = [ide(n_el,2); ide(n_el+1,2)];
    r(unk_1) = r(unk_1) + r1_loc;
      
    % calculate the columns of J corresponding to the local unknowns
    % increment phi0 for the first node
    splus_local   = s_local + [delta_w; 0];
        
    splus_g       = phi0 *splus_local  + phi1 *sp_local;
    splusx_g      = p0_x *splus_local  + p1_x *sp_local;
    splusxx_g     = p0_xx*splus_local  + p1_xx*sp_local;

    r0plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p0_x , wt_g) ...
     -           oned_f_int( theta*wx_g     + omt*wcx_g  , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ...
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi0 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                 , phi0 , wt_g);

    r1plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p1_x , wt_g) ...
     -           oned_f_int( theta*wx_g + omt*wcx_g      , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ... 
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi1 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                 , phi1 , wt_g);
      
    ntriplets = ntriplets + 1;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_0(1);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;
 %   J(unk_0, unk_0(1)) = J(unk_0, unk_0(1)) + (r0plus_loc - r0_loc)/delta_w;
 
 
    ntriplets = ntriplets + 2;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_0(1);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;
%    J(unk_1, unk_0(1)) = J(unk_1, unk_0(1)) + (r1plus_loc - r1_loc)/delta_w;
    
 
    % increment phi0 for the second node
    splus_local   = s_local + [0; delta_w];
        
    splus_g       = phi0 *splus_local  + phi1 *sp_local;
    splusx_g      = p0_x *splus_local  + p1_x *sp_local;
    splusxx_g     = p0_xx*splus_local  + p1_xx*sp_local;
 
    r0plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p0_x , wt_g) ...
     -           oned_f_int( theta*wx_g     + omt*wcx_g  , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ...
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi0 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                 , phi0 , wt_g);

    r1plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p1_x , wt_g) ...
     -           oned_f_int( theta*wx_g + omt*wcx_g      , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ...   
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi1 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                 , phi1 , wt_g);

    ntriplets = ntriplets + 2;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_0(2);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;
%    J(unk_0, unk_0(2)) = J(unk_0, unk_0(2)) + (r0plus_loc - r0_loc)/delta_w;

    ntriplets = ntriplets + 2;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_0(2);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;
%    J(unk_1, unk_0(2)) = J(unk_1, unk_0(2)) + (r1plus_loc - r1_loc)/delta_w;

    % increment phi1 for the first node
    spplus_local   = sp_local + [delta_w; 0];
        
    splus_g       = phi0 *s_local  + phi1 *spplus_local;
    splusx_g      = p0_x *s_local  + p1_x *spplus_local;
    splusxx_g     = p0_xx*s_local  + p1_xx*spplus_local;
 
    r0plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p0_x , wt_g) ...
     -           oned_f_int( theta*wx_g     + omt*wcx_g  , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ...
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi0 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                  , phi0 , wt_g);

    r1plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p1_x , wt_g) ...
     -           oned_f_int( theta*wx_g + omt*wcx_g      , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ...   
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi1 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                  , phi1 , wt_g);
    
    ntriplets = ntriplets + 2;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_1(1);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;
%    J(unk_0, unk_1(1)) = J(unk_0, unk_1(1)) + (r0plus_loc - r0_loc)/delta_w;

    ntriplets = ntriplets + 2;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_1(1);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;
%    J(unk_1, unk_1(1)) = J(unk_1, unk_1(1)) + (r1plus_loc - r1_loc)/delta_w;

    % increment phi1 for the second node
    spplus_local   = sp_local + [0; delta_w];
        
    splus_g       = phi0 *s_local  + phi1 *spplus_local;
    splusx_g      = p0_x *s_local  + p1_x *spplus_local;
    splusxx_g     = p0_xx*s_local  + p1_xx*spplus_local;
 
    r0plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi0 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p0_x , wt_g) ...
     -           oned_f_int( theta*wx_g + omt*wcx_g      , p0_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p0_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ...
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi0 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                 , phi0 , wt_g);

    r1plus_loc =      oned_f_int( (splus_g-sc_g)/t_step  , phi1 , wt_g) ...
     - epsilon  *oned_f_int( theta*splusx_g + omt*scx_g  , p1_x , wt_g) ...
     -           oned_f_int( theta*wx_g + omt*wcx_g      , p1_x , wt_g) ...
     + epsilon^2*oned_f_int( theta*splusxx_g + omt*scxx_g, p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*wxx_g + omt*wcxx_g    , p1_xx, wt_g) ...
     + 2*epsilon*oned_f_int( theta*(splus_g.*wx_g + w_g.*splusx_g)      ...  
                       + omt*(sc_g.*wcx_g + wc_g.*scx_g) , phi1 , wt_g) ...
     +         2*oned_f_int( theta*w_g.*wx_g                            ...
                       + omt*wc_g.*wcx_g                 , phi1 , wt_g);
     
    ntriplets = ntriplets + 2;
    II( ntriplets:ntriplets+1 ) = unk_0;
    JJ( ntriplets:ntriplets+1 ) = unk_1(2);
    XX( ntriplets:ntriplets+1 ) = ( r0plus_loc - r0_loc )/delta_w;
%    J(unk_0, unk_1(2)) = J(unk_0, unk_1(2)) + (r0plus_loc - r0_loc)/delta_w;

    ntriplets = ntriplets + 2;
    II( ntriplets:ntriplets+1 ) = unk_1;
    JJ( ntriplets:ntriplets+1 ) = unk_1(2);
    XX( ntriplets:ntriplets+1 ) = ( r1plus_loc - r1_loc )/delta_w;
%    J(unk_1, unk_1(2)) = J(unk_1, unk_1(2)) + (r1plus_loc - r1_loc)/delta_w;
    ntriplets = ntriplets + 1;

  end  %  n_el=1:n_elements

  J = sparse( II(1:ntriplets), JJ(1:ntriplets), XX(1:ntriplets), ...
              n_equations, n_equations);

end % function weak_resjac


%-----------------------------------------------------------------------                      
%-----------------------------------------------------------------------
function [] = plot_hermite(fig_num,x,e_conn,w,wp)
  figure(fig_num)
  hold on
  [n_elements, tmp] = size(e_conn);
  
  r = linspace(-1,1,11);  wt = zeros(size(r));  % changed Dzeros to zeros
  for n_el=1:n_elements
    nodes_local            = e_conn(n_el,:);
    x_local                = x(nodes_local,:);
    [x_g,w_g,phi0,phi1,p0_x,p1_x,p0_xx,p1_xx] = ...
                                           oned_shapeherm(x_local,r,wt);
    
    w_plot = phi0*w(nodes_local) + phi1*wp(nodes_local);
    plot(x_g,w_plot)
  end
  hold off
  %         up_plot = p0_x*u(nodes_local)' + p1_x*up(nodes_local)';
%         plot(x_g,up_plot,'r')
%     
%         upp_plot = p0_xx*u(nodes_local)' + p1_xx*up(nodes_local)';
%         plot(x_g,upp_plot,'g')

end % function plot_hermite
