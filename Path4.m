function x=Path4(i,j,next)
 
    if dist(i,j)==inf
      x=NaN; 
      
     end

     intermediate = next(i,j);
      if intermediate ==0
         x=[]; % the direct edge from i to j gives the shortest path
      else
       
        x=[Path4(i,intermediate,next),intermediate,Path4(intermediate,j,next)];
 
      end
