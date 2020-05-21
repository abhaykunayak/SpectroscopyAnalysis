function result = dft(x)
N = length(x);
n = 1:N;
k = reshape(n,N,1);

M = exp(-2j.*pi.*k.*n./N);

result = sum(M.*x,2);
end