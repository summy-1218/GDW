function df = double_factorial(n)
if n <= 0
    df = 1;
    return;
end
if mod(n,2)==0
    df = prod(2:2:n);
else
    df = prod(1:2:n);
end
end

