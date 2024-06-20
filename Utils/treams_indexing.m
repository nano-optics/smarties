function [u] = treams_indexing(q, qmax)

pmax = qmax / 2;
s = floor(q./(pmax + 1)) + 1;
p = q - (s - 1)*pmax;
u = 2*(p-1) + s; 

end
