s = tf('s');
A = 1/(0.5*s+1);
figure(1);
bode(A)
hold on