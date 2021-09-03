% plotallsimulations.m
%
% Matlab script to plot all simulation results.

load grid_data

iz = find(inpolygon(X,Y,xb,yb) == 0);

for(i=1:1)
%for(i=1:9)
   
   figure(i);
   
   load(strcat('Simulation',num2str(i)));
   
   total = zeros(nx,ny);
   
   for(j=1:ns)
      
      total = total + datamat(:,:,j)./ns;
      
   end
   
   total1 = total';
   
   total1(iz) = 0;
   
   xa = xmin:dx:xmax;
   ya = ymin:dy:ymax;

   imagesc(xa,ya,log10(total1),[0 4.2]);
   
   a = gca;set(a,'ydir','normal');axis equal;axis image
   
   xlabel('Easting (m)');
   ylabel('Northing (m)');
   colorbar
   
   hold on
   plot(x,y,'k.');
   
end

