% for i = 1:1000
%         n = randi([0,2^32-1],1,'uint32');
%         s = sprintf('%08X',n);
%         m = sscanf(s,'%X');
%         assert(m == n);
% end

n = randi([0,2^32-1],t,1,'uint32');
s = sprintf('%08X ',n);
m = sscanf(s,'%X');
assert(all(m == n));