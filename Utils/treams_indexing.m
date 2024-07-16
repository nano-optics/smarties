function [u] = treams_indexing(q, qmax)

pmax = qmax / 2;
s = floor(q./(pmax + 1)) + 1;
p = q - (s - 1)*pmax;

% below is the TERMS -> TREAMS conversion 
% 3-s: convert magnetic/electric -> electric/magnetic
% picks odd/even electric, magnetic in alternance
u = 2*(p-1) + (3 - s); 

end
