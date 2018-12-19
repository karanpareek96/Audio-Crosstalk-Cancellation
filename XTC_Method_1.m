%{
    Assignment 4 - 3D Audio
    Crosstalk Cancellation method using regularization
    Karan Pareek, Chengze Zuo, Aneesh Athrey
%}

% Read the binaural file
[binauralSignal,fs] = audioread('test.wav');

% read the impulse responses 
[output_L,~] = audioread('LeftChannel.wav');
[output_R,~] = audioread('RightChannel.wav');

% C(z) - response of loudspeakers at the ears
% H(z) - filters for cross-talk minimization

b=0.05; % regularization factor

HRIRs(1,1,:)=output_L(:,1)';
HRIRs(1,2,:)=output_R(:,1)';
HRIRs(2,1,:)=output_L(:,2)';
HRIRs(2,2,:)=output_R(:,2)';
len=length(HRIRs(1,1,:));

% Transfering into the frequency domain
for i = 1:2
    for j = 1:2
        C_f(i,j,:)=fft(HRIRs(i,j,:),len);
    end
end

% Regularized inversion of matrix C. For the simplicity of the z transform,
% we just compute the regularization on the frequency domain
H_f=zeros(2,2,len);

for k = 1:len
    H_f(:,:,k)=inv((C_f(:,:,k)'*C_f(:,:,k)+eye(2)*b))*C_f(:,:,k)';
end

% Back to time domain
for k = 1:2
    for m = 1:2
        H_n(k,m,:)=real(ifft(H_f(k,m,:)));
        H_n(k,m,:)=fftshift(H_n(k,m,:));
    end
end

loudspsig=[conv(reshape(H_n(1,1,:),len,1),binauralSignal(:,1)) + ...
    conv(reshape(H_n(1,2,:),len,1),binauralSignal(:,2)) ...
    conv(reshape(H_n(2,1,:),len,1),binauralSignal(:,1)) + ...
    conv(reshape(H_n(2,2,:),len,1),binauralSignal(:,2))];

%soundsc(loudspsig,fs)
audiowrite('cross_talk_signal.wav',loudspsig,fs); 