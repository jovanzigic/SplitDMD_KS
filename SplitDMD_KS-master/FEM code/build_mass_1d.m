  function [M] = build_mass_1d(x, e_conn)
%-----------------------------------------------------------------------
%  build_mass_1d.m - builds the finite element "Mass" matrix
%                    assuming H^1(\Omega)  (no Dirichlet)
%
%  Copyright (c) 2007, Jeff Borggaard, Virginia Tech
%  Version: 1.0
%
%  Usage:    [M] = build_mass_1d
%
%  Variables:  x
%                        Nodal coordinates
%              e_conn
%                        1D element connectivity
%
%              M
%                        M_{ij} = \int h_i h_j dx
%-----------------------------------------------------------------------

  n_gauss          = 3;  % number of points used in integration
                         % 3 should be sufficient for low-order elems.

  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);
    
  %---------------------------------------------------------------------
  %  Build the finite element matrices
  %---------------------------------------------------------------------
  M = sparse(n_nodes,n_nodes);
  b = zeros (n_nodes,1      );

  [r,w] = oned_gauss(n_gauss);
  
  for n_el=1:n_elements
    % compute value of each test function and their spatial derivaties
    % at the integration points
    nodes_local           = e_conn(n_el,:);
    x_local               = x(nodes_local,:);
    [x_g,w_g,phi,p_x] = oned_shape(x_local,r,w);
      
    % compute the value of weighting function at the Gauss points
    p_g = p_function(x_g);
    
    %-------------------------------------------------------------------
    %  Integrate the element contributions to M
    %-------------------------------------------------------------------
    M_loc = oned_bilinear(p_g, phi, phi, w_g);

    %-----------------------------------------------------------------
    %  Assemble contributions into the global system matrix
    %-----------------------------------------------------------------
    M(nodes_local,nodes_local) = M(nodes_local,nodes_local) + M_loc;
    
  end
  
%-----------------------------------------------------------------------
%  Supporting Functions
%-----------------------------------------------------------------------

% twod_backm2d
% twod_gauss
% twod_shape
% twod_bilinear
% twod_f_int
% twod_plotc
% twod_plotm1
% twod_plotm2
% twod_projectd

function p = p_function(x)
  p = 1;

