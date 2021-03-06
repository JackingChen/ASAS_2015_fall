function out=sound_array_generator(fs,duration,t_note,f_flank1,f_note,f_flank2)
%%

duration = 0.2;
fs = 16000;
current_t = normrnd(0,Dt+dt);
current_f = normrnd(f1,Df+df);

t_note=current_t;
f_flank1=f1;
f_note=current_f;
f_flank2=f2;

x=(1:fs*duration);
Nn=fs*duration;
% gua=gaussmf(x,[Nn f_flank1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% YOUR CODE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%
note_length = round(fs*t_note);
out = zeros(fs*3,1);

gua=gaussmf(x,[Nn f_flank1]);
out((fs*1+1):(fs*(1+duration))) = gua.*(0.5*(sin(2*pi*f_flank1*(1:fs*duration)/fs)));

gua=gaussmf(x,[Nn f_note]);
out(fs*2+note_length+1:fs*(2+duration)+note_length) = 0.5*(sin(2*pi*f_note*(1:fs*duration)/fs)).*gua;

gua=gaussmf(x,[Nn f_flank2]);
out(fs*2+1:fs*(2+duration)) = out(fs*2+1:fs*(2+duration))+ 0.5*sin(2*pi*f_flank2*(1:fs*duration)/fs)';
%keyboard



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% YOUR CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%%%%
%%
end