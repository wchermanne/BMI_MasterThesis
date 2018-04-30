clear all;
close all;
S1 =  1:10;
lS1 = length(S1)
S2 = -(-10:-1);
level = 2;
dwtmode('sym')
[LO_D,HI_D,LO_R,HI_R] = wfilters('db3')
A = (conv([S2 S1 S2],LO_D))
B = (conv([S2 S1 S2],HI_D)) % dyaddown
Aprime = dyaddown(A(lS1+1:end-lS1))
Bprime = dyaddown(B(lS1+1:end-lS1))

lS2 = length(Aprime);
A = (conv([flip(Aprime) Aprime flip(Aprime)],LO_D));
B = (conv([flip(Aprime) Aprime flip(Aprime)],HI_D));
Asec = dyaddown(A(lS2+1:end-lS2))
Bsec = dyaddown(B(lS2+1:end-lS2))

[C_1,L1] = wavedec(S1,level,LO_D,HI_D)

%C_1 - [Aprime Bprime]

