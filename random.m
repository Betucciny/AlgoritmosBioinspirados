clear all
for i = 1:100
    x(i,1)=randgen(4, 50);
end


function r = randgen(min, max)
    persistent state; 
    if isempty(state)
        state = sum(100 * clock);
    end
    r = mod(1.4324672*state, 1);
    state = r;
    r = r * (max - min) + min;
end
