%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EE6641 HW1                        %
% sound_generator.m                 %
% created by Jeffrey Huang, 10/2015 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out ucert] = sound_generator(Dt,Df,dt,df,f1,f2,opt)
%ft and Df stand for the original meanings in the paper, while dt and df
%are the detla_t and detla_f. (Here, delta is the lowercase of Greek alphabet)
duration = 0.2;% duration of every notes
fs = 16000;
correctness.time = zeros(1000,1);
correctness.frequency = zeros(1000,1);
correctness.resolution = zeros(2,1);
correctness.stage_idx = 1;

u_count=1;
uncertainty=zeros(15,1);

Dt_alt=0.0001;
Df_alt=0.0001;
if opt==1
    fprintf('Welcome to the example mode!\n');
    while 1
        current_t = normrnd(0,Dt+dt);
        current_f = normrnd(f1,Df+df);
        sound_array = sound_array_generator(fs,duration,current_t,f1,current_f,f2);
        %keyboard
        sound(sound_array,fs);
        fprintf('The note comes before or after the highest note?\n');
        query = input('0 is before, 1 is after, or press other digits to hear the sound again.\n');
        
        while query ~= 0&&query ~= 1
            sound(sound_array,fs);
            fprintf('The note comes before or after the highest note?\n');
            query = input('0 is before, 1 is after, or press other digits to hear the sound again.\n');
        end
        
        fprintf('The pitch of note is higher or lower the first one?\n');
        query2 = input('0 is lower, 1 is higher, or press other digits to hear the sound again.\n');
        
        if query2 ~= 0&&query2 ~= 1
            sound(sound_array,fs);
            fprintf('The pitch of note is higher or lower the first one?\n');
            query2 = input('0 is lower, 1 is higher, or press other digits to hear the sound again.\n');
        end
        
        if (query-0.5)*current_t >= 0
            fprintf('Your answer to the note appearance is correct! (%.3f s)\n',Dt);
            Dt=Dt-Dt_alt;
        else
            fprintf('Your answer to the note appearance is incorrect! (%.3f s)\n',Dt);
            Dt=Dt+Dt_alt;
        end
        
        if (query2-0.5)*(current_f-f1) >= 0
            fprintf('Your answer to the note pitch is correct! (%.3f hz)\n',Df);
            Df=Df-Df_alt;
        else
            fprintf('Your answer to the note pitch is incorrect! (%.3f hz)\n',Df);
            Df=Df+Dt_alt;
        end
        
        fprintf('The multiplication of the uncertainty rule is %f. \n',Dt*Df*4*pi);
        
        break_example_query = input('Exit example mode (Y/N) ?','s');
        
        if break_example_query=='Y'||break_example_query=='y'
            break;
        end
    end
    out = -1;
else
    fprintf('Welcome to the test mode!\n')
    turns=10;
    while turns>0
        fprintf('stage %u?\n',20-turns);
        turns=turns-1;
        current_t = normrnd(0,Dt+dt);
        current_f = normrnd(f1,Df+df);
        sound_array = sound_array_generator(fs,duration,current_t,f1,current_f,f2);
        %keyboard
        sound(sound_array,fs);
        fprintf('The note comes before or after the highest note?\n');
        query = input('0 is before, 1 is after, or press other digits to hear the sound again.\n');
        chance=0;
        while query ~= 0&&query ~= 1
            if chance~=0
                sound(sound_array,fs);
            else
                fprintf('You cannot repeat again, just pick an answer?\n');
            end
            sound(sound_array,fs);
            fprintf('The note comes before or after the highest note?\n');
            query = input('0 is before, 1 is after, or press other digits to hear the sound again (you have no chance to repeat).\n');
            chance=0;
        end
        
        fprintf('The pitch of note is higher or lower the first one?\n');
        query2 = input('0 is lower, 1 is higher, or press other digits to hear the sound again.\n');
        chance=0;
        while query2 ~= 0&&query2~= 1
            if chance~=0
                sound(sound_array,fs);
            else
                fprintf('You cannot repeat again, just pick an answer?\n');
            end
            
            fprintf('The pitch of note is higher or lower the first one?\n');
            query2 = input('0 is lower, 1 is higher, or press other digits to hear the sound again (you have no chance to repeat).\n');
            chance=0;
        end
        
        if (query-0.5)*current_t>=0
            fprintf('Your answer to the note appearance is correct! (%.3f s)\n',Dt);
            correctness.time(correctness.stage_idx) = abs(Dt);
            Dt = Dt*0.7; %ft = dt*0.8367;
        else
            fprintf('Your answer to the note appearance is incorrect! (%.3f s)\n',Dt);
            correctness.time(correctness.stage_idx) = -abs(Dt);
            Dt = Dt*1.4286; %ft = dt*1.1952;
        end
        
        if (query2-0.5)*(current_f-f1) >= 0
            fprintf('Your answer to the note pitch is correct! (%.3f hz)\n',Df);
            correctness.frequency(correctness.stage_idx) = abs(Df);
            Df = Df*0.9091; %ff = df*0.8367;
        else
            fprintf('Your answer to the note pitch is incorrect! (%.3f hz)\n',Df);
            
            correctness.frequency(correctness.stage_idx) = -abs(Df);
            Df = Df*1.1; %ff = df*1.1952;
        end
        
        
        
%         break_test_query = input('Exit test (Y/N) ?','s');
%         
%         if break_test_query=='Y'||break_test_query=='y'
%             correctness.resolution(1) = Dt;
%             correctness.resolution(2) = Df;
%             break;
%         end
          correctness.resolution(1) = Dt;
          correctness.resolution(2) = Df;
          if (query-0.5)*current_t>=0 && (query2-0.5)*(current_f-f1) >= 0
              uncertainty(u_count)=4*pi*correctness.resolution(1)*correctness.resolution(2);
              u_count=u_count+1;
          end
        correctness.stage_idx = correctness.stage_idx+1;
    end
    out=correctness;
    ucert=uncertainty;
end

end

