% loop_waitbar(i,imax)
%   use the function fprintf to display a wait bar for loops.
%   
%   i     current loop iteration
%   imax  number of iterations in the loop
%

function loop_waitbar(i,imax)


if i==1
    for j=0:5
        fprintf('% 2i%%       ',20*j);
    end
    fprintf('\n ');
end
if mod(i,max(floor(imax/5),4))==1,fprintf('|'),end
if mod(i,max(floor(imax/50),2))==1,fprintf('.'),end
if i==imax,fprintf('|\n');end