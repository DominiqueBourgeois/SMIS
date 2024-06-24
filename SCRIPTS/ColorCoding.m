function ColorVector = ColorCoding(N)

ColorVector = zeros(N,3);

N_1_half = round(N/2);
N_2_half = N - N_1_half;

for n = 1 : N_1_half
    ColorVector(n,1) = 0;
    ColorVector(n,2) = (n-1)/N_1_half;
    ColorVector(n,3) = 1 - (n-1)/N_1_half;
end

for n = 1 : N_2_half
    ColorVector(n+N_1_half,1) = (n-1)/N_2_half;
    ColorVector(n+N_1_half,2) = 1 - (n-1)/N_2_half;
    ColorVector(n+N_1_half,3) = 0;
end