 % 507 Final Project
% Gwen Montague & Olivia Killingsworth

% load audio file
[ss, Fs] = audioread("Soft Piano Music_16000hz.wav" );

% call provided bands_cutoff function & set input values
fmin=300; %in Hz
fmax=6000; %in Hz
N=20; %N is number of bandpass filters
fco=bands_cutoff(fmin, fmax, N);

% initialize zero vector w/ length of sound signal's
z=zeros(size(ss));

for i=1:N
    % bandpass filter design
    Wn=[fco(i) fco(i+1)]/(Fs/2);
    [bb,aa]=butter(3,Wn,'bandpass'); %3rd order BFP at 300 to 6000 Hz
    % apply BPF
    BPFilter_signal=filter(bb,aa,ss);
    
    % pass band-passed signal through envelope
    % rectification
    FWRectified_signal=abs(BPFilter_signal);
    % create & apply LPF to rectified signal 
    [b,a]=butter(2,400/(Fs/2),'low'); %2nd order LPF at 400 Hz
    LPFilter_signal=filter(b,a,FWRectified_signal); %envelope of bandpassed signal
    
    Fc = (fco(i) + fco(i+1))/2; %center frequency of ith band


    % sine carrier
    t=(0:length(ss)-1)/Fs; %time vector w/ length same as signal
    T=t'; %transposed time vector
    carrier=sin(2*pi*Fc*T);
    % apply modulation via multiplication
    CS=carrier.*LPFilter_signal;

    % band sum
    [bbb,aaa]=butter(3,Wn,'bandpass');
    BPFcs=filter(bbb,aaa,CS);
    z=z+BPFcs;
    
end 

% output vocoded audio signal
Signal_out = z/max(abs(z)); %normalize signal in range of [-1,+1]
sound(Signal_out,Fs); pause(3);
sound(ss,Fs);



%----------------design of the logarithmic filter bank---------------------

% fmin: the lowest frequency to cover 
% fmax: the highest frequency to cover 
% N: number of bands or number of channels 
% The output fco contains N+1 frequency values for N bandpass filters. For 
% example, [fco(i) fco(i+1)] is for the ith bandpass filter. 
function fco=bands_cutoff(fmin, fmax, N)
xmin = log10(fmin/165.4+1)/2.1;
xmax = log10(fmax/165.4+1)/2.1; %relative value
dx = (xmax-xmin)/N;
x = xmin:dx:xmax;
fco=zeros(1,N+1);
for i=1:N+1
 fco(i)=165.4*(10^(x(i)*2.1)-1);
end
end

