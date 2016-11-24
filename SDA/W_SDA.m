clc
clear variables
cvx_setup
 I = imread('boat.bmp');
 I = double(I);
Size = size(I);
k=1;
Framesize=5;
% Create k-1 lost blocks in the image.In order to assure that the random
% blocks won't exceed matrix dimensions, we demand the random numbers to be from one 
% until  the end f each line substracted by the length of the random block minus 1.
%(Therefore we set the initial parameters from 1 to Size(1)-7. We also substract 5 since the biggest
% frame size is equal to five. That requirement is set for the same reason
% as before.
% Definitions. 1) M->Mask Matrix,2) D1->DCT matrix,3) p,r the randomly
% chosen indices of the lost block.4)q the size of the framesize around the
% lost block.5)Framesize-> the highest framesize we are testing.maxPSNR is
% value which helps us in finding the highest PSNR for each concealed block 
% regarding the frame-size so we can take the right concealed block.I_cor
% is the corrupted image and I_con is the concealed image.
I_cor=I;
I_con=I;
while k<=10
p = randi([1+Framesize,Size(1)-(7+Framesize)],1);
r=  randi([1+Framesize,Size(2)-(7+Framesize)],1);
maxPSNR=0;
if I_cor(p:p+7,r:r+7)~=0
 Block=I(p:p+7,r:r+7);
 I_cor(p:p+7,r:r+7)=0;
 Block_cor=I_cor(p:p+7,r:r+7);
 error= Block - Block_cor;
 MSE(k) = (sum(sum(error .* error)))./(8*8);
 PSNR(k)= 10*log(255^2/MSE(k)) / log(10);
 o=1;
% Different frame sizes. Each frame size is equal to q.
for N=10 
M=eye(N*N);
q=(N-8)/2;
for i=q*N+1:N*N-q*N
 for j=q+1:N-q
 if mod(i,N)==j
  M(i,i)=0;
 end
 end
end
% W is the weighting matrix
vector=ones(1,N^2);
for i=1:2*q*N
 vector(i)=0.37;
end
W=diag(vector);
% FullBLock is the initial 8X8 lost block including the frame size. 
FullBlock = I_cor(p-q:p+7+q,r-q:r+7+q);
DCT1=dctmtx(N);
D2 = kron(DCT1,DCT1);
DCT_tovector= reshape(D2,(N^2)*(N^2),1);
u= reshape(FullBlock,N*N,1);
w=M*u;
d=0;
   cvx_begin 
     variable v(N*N)
     minimize norm(W*D2*v,1);
     subject to 
     norm(M*v-w,2) <= d;
    cvx_end 
 Reshape_v= reshape(v,N,N);
 Block_concealed=Reshape_v;
 error1 = Block - Block_concealed(q+1:N-q,q+1:N-q);
 MSE1(k,o) = (sum(sum(error1.*error1)))./(8*8);
%  Calculate PSNR for k-1 lost blocks for q different frame sizes. PSNR1 is a
%  matrix which in each row contains each lost block with different frame
%  sizes along the columns.
PSNR1(k,o) = 10*log(255*255/MSE1(k,o)) / log(10);
if  PSNR1(k,o)>maxPSNR 
 I_con(p:p+7,r:r+7)=Block_concealed(q+1:N-q,q+1:N-q);
 maxPSNR=PSNR1(k,o);
end
o=o+1;
end
k=k+1;
end
end
% Find the highest PSNR for each block.We take out  a matrix(PSNR1 in our
% case). Matrix in from the following command contains the indices of the
% highest PSNR in each case.
% imshow(I,[])
[B in]= max(PSNR1,[],2);
sum1= zeros(q,1);
for i=1:q
 for j=1:length(B)
  if in(j)==i
   sum1(i)=sum1(i)+1;
  end
 end
end
%  Calculate the average performance of PSNR1 for each framesize..
%  Average= mean(PSNR1);
% % Calculate the average perfrmance of the maximum PSNR1 for each block.
%  Average_max= mean(B);
% % Times for each PSNR having the best performance.
%    subplot(2,3,4),stem(sum1),axis([-1 o 0 k]),xlabel('Frame Size');
% %  Average PSNR-performance for each framesize.
%    subplot(2,3,5),stem(Average),axis([-1 o 0 50]),xlabel('Average PSNR-performance')
% % plot the original image, the corrupted image and the concealed image.
subplot(1,3,1), imshow(I,[]);
title('Original Image')
subplot(1,3,2), imshow(I_cor,[]);
title('Corrupted Image')
subplot(1,3,3),imshow(I_con,[]);
title('Concealed Image')
