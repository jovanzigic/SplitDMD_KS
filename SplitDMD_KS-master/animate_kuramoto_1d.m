function [M] = animate_kuramoto_1d(L)


 if (L<10)
   number = strcat('0',int2str(100*L));
 else
   number = int2str(100*L);
 end
 snapshot_file = strcat('snapshots_at_L',number);

 load(snapshot_file)
 [n_nodes,kmax] = size(w_save)
 w_min = min(min(w_save));
 w_max = max(max(w_save));

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

 for k=1:kmax
   oned_plot_hermite(20,x,e_conn,w_save(:,k),wp_save(:,k))
   caption = strcat(strcat(strcat('L = ',num2str(L)),', time = '),num2str((k-1)*0.05,'%6.2f'));
   title(caption)
   axis([0 1 w_min w_max])
  
   M(k) = getframe;
   close(20)
 end

 movie(M)

 movie_file = strcat('movie_L',int2str(100*L))
 save(movie_file,'M')
