clear, close all
for L = [402.35:0.01:402.35]
  kuramoto_1db(1/L^2)
  solution_file = strcat('snapshotsx_at_L',int2str(100*L),'_400_sens1');
  load kuramoto_1db_snap 
  save( solution_file, 'w_save', 'wp_save', 'epsilon' )
end

% for Re = [15:10:495]
%   burgers_2d(Re)
%   solution_file = strcat('snapshots_at_Re',int2str(Re));
%   load burgers_2d_snap w_save
%   save( solution_file, 'w_save', 'Re' )
% end

% i=0
% load burgers_2d_steady
% for Re = [10:10:230]
%   solution_file = strcat('steady_sol_at',int2str(Re));
%   load( solution_file )
%   twod_plotc(20, x, e_conn, w )
%   i = i + 1;
%   M(i) = getframe;
% end
% 
% movie(M)