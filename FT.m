clear
clc
load brain512;
data64=imresize(data,[64,64]); %��ͼ���С��Ϊ64*64
image=data64;
[M,N]=size(image);
Fuv1=zeros(M,N);%%����Ҷ�任
for u=0:M-1
    for v=0:N-1
        temp=0;
        for x=0:M-1
           for y=0:N-1
               temp=image(x+1,y+1)*exp(-1i*2*pi*(u*x/M+v*y/N))+temp;
           end
        end
        Fuv1(u+1,v+1)=temp;
    end
end

fxy1=zeros(M,N);%%����Ҷ��任
for x=0:M-1
    for y=0:N-1
        temp=0;
        for u=0:M-1
           for v=0:N-1
               temp=Fuv1(u+1,v+1)*exp(1i*2*pi*(u*x/M+v*y/N))+temp;
           end
        end
        fxy1(x+1,y+1)=temp/M/N;
    end
end

w=exp(-1i*2*pi/N);
Fmatrix=zeros(M,N);%%���ɸ���Ҷ�任����
for m=0:M-1
    for n=0:N-1
        Fmatrix(m+1,n+1)=w^(m*n);
    end
end


Fuv2=(Fmatrix*(Fmatrix*image)')';%%����Ҷ�任
fxy2=(Fuv2'/Fmatrix)'/Fmatrix;%����Ҷ��任

C=fft2(image);

subplot(3,2,1); imshow(image);title('origin image');
subplot(3,2,2); imshow(C);title('fft2 image');
subplot(3,2,3); imshow(Fuv1) ;title('for loop image');
subplot(3,2,4);imshow(fxy1);
subplot(3,2,5); imshow(Fuv2);title('matrix image');
subplot(3,2,6);imshow(fxy2);
