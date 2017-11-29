%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Monterrey Nuevo Leon Mexico 11/28/2017 Biomedical Engineering Department    
%Advisor: PHD Luz Maria Alonso Velardi                                         
%María Del Carmen Guzmán Martínez
%Francisco Darío Sánchez De la Cruz  
% Instructions to use this code:
% 1.Before using these codes EEGlab should be installed on Matlab, and run it  
% and close it to preload the biosig plug in
% 2. The Preprocessing lines require a GDF file on EEG structure format to run
% 3. This particular approach uses the standard 10/20 electrode configuration, 
% with 19 channels, a window will pop up assuming this acquisition setup, 
% conversely EEGlab allows the user to rename the channels.               
% 4. If the user selects to manually reject EMG artifacts, the user can       
% select multiple little areas.
% 5. When the "ICA" pop up window appears just hit OK. And wait.
% 6. 19 ICA components will be displayed just click OK.
% 7. Then on the next popping up window the program asks the user to pass the  
% ICA components not rejected from one set to another, on the display look on % the last for cells to input sets 
% click on the set selection button(left from the input cell) and select set   % 3 
% 8. Once again the 19 components will be displayed, allowing the user to      % reject the desired ones.    
                        %End of the preprocessing lines. 
% 9. On the EPOCH prompt window write "-.5 1.5" if seconds seconds are         % requested, or "-500 1500" for milliseconds.
% 10. Compute the map of each EPOCH takes 7 minutes each one, each iterations outputs 19 maps on frequency, time and power, (one for each channel).
% 11. The processing time takes up to 7 minutes X #Of EPOCHS 
% 12. The final output will be a cell structure with dimension equal to #of epochs X 1 where each epochs includes the 19 channels.  
% 13. PLVFunct  requires the user to give the cell structure and a vector indicating the EPOCHS that will be used to calculate the values. 
% 14. The MappingEpoch Function also requires the cell structure, a vector specifying the Epochs to be mapped, and finally a vector of the form [x y] where x is the lower frequency limit and  y the high frequency limit.
% Special consideration the "dcblock.m" file should be on the Matlab folder, to make the filtering function work.    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DataPREProcessing
clear
EEG=pop_biosig(); %load gdf file 
EEG.data=rmbase(EEG.data(:,:));%RemoveBaseline
EEG=pop_chanedit(EEG); %set channels
ALLEEG=[];
ALLCOM=[];% spaces used to work with structures on eeglab based functions
CURRENTSET=[];
xx=input('Quieres eliminar EMG? Y/N:','s');
if xx=='Y'
pop_eegplot( EEG, 1, 1, 1);%manual EMG component removal dont forget to hit have  dataset
pause();
EEG=EEGTMP;
clear EEGTMP;
end
EEG=BandPass(EEG,1,100); %filter 
clear xx
ALLEEG=EEG;
a=EEG;
a=BandPass(a,1,100); %filter 
CURRENTSET=2;
ALLEEG(double(CURRENTSET))=a; %store structure
a=pop_runica(a); %run ICA to eliminate unwanted components   
%unwanted components elimination
CURRENTSET=[3];
ALLEEG(double(CURRENTSET))=a;
pause();
EEGf=pop_editset(EEG);
pause();
EEGf=pop_selectcomps(EEGf);
pause();
EEGf=pop_subcomp(EEGf);
pause();
pop_eegplot(EEGf,1,1,1);
pause();
CURRENTSET=4;
ALLEEG(double(CURRENTSET))=EEGf;
%END of data preprocessing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EEGfep=pop_epoch(EEGf);%Epoch split
EEGfep=pop_rmbase(EEGfep,[-500 0]);%Remove baseline
x=EEGfep.data(); %get data outside str

  n=1000;
    ii=length(x)-n;
    hanW=hann(n);%hanning window
    gg=[1:34]; %selected epochs
    ii=length(x)-n;
    CellFtf=cell(length(gg),1);

    for jj=1:length(gg) 
        X=[]; %%clear wft vector
            for i=1:ii
                wft=x(:,1+(i-1):n+(i-1),gg(jj)).*hanW'; %multiply the truncated signal by the window on a determined epoch
                WFT=fft(wft,[],2); %wft of one epoch
                X=[X,WFT]; %concatenated the 19 channels  on time
            end
            for iii=1:19
                X22=[];
                for i=1:length(X)/n
                    X22=[X22;X(iii,1+(i-1)*n:n+(i-1)*n)];%rearrange them by time of a single channel
                end
                X22=(X22(:,1:end/2));% select only half of the spectram component
                Ftf(:,:,iii)=X22;
            end
        CellFtf{jj,1}=Ftf;
        pause()
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [CellData]=PLVFunct(a,st)
    CellData=cell(length(a),1);
    for act=1:length(a)
        Ftf=st{a(act),1};
        CellData1{1,1}=abs(Ftf);
        CellData1{1,2}='Ftf';
        CellData1{2,1}=abs(mean(Ftf./abs(Ftf),3));
        CellData1{2,2}='PLV';
        CellData1{3,1}=abs(Ftf).^2;
        CellData1{3,2}='Ptf';
        CellData1{4,1}=mean(Ptf(100:400,:,:));
        CellData1{4,2}='Rf';
        CellData1{5,1}=CellData1{3,1}-CellData1{4,1};
        CellData1{5,2}='Pstf';
        CellData{act,1}=CellData1;
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Function to map desired Epochs.
function MappingEpochs(a,st,lim)
    labl=['Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','F7','F8','T3','T4','T5','T6','Fz','Cz','Pz'];
    verg=[2,4,7,9,12,14,17,19,22,24,6,10,11,15,16,20,8,13,18];
        for act=1:length(a) 
            figure;    
            X=abs(st{a(act),1});
            for iii=1:19
                X22=X(:,:,iii);
                h=subplot(5,5,verg(iii));
                image(squeeze(X22), 'CDataMapping', 'scaled', 'XData', .5:500, 'YData', [-500 1000]);
                caxis([min(min(abs(X22))) max(max(abs(X22)))])
                title(labl(iii))
                colormap(jet)
                grid on
                xlim(lim)
                colorbar()
            end
        end  

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%EEG structure filtering function 
function [EEGout] = BandPass(EEGin,fH,fL)
% --- Variable Declaration ---
EEGout = EEGin; fs = EEGin.srate;
% high-pass coeffs
if fH ~= 0
    if fH < 2
        a = dcblock(fH, fs);
        Bh = [1, -1];
        Ah = [1, -a];
    else
        Wn = fH/(fs/2);
        [Bh, Ah] = butter(6, Wn, 'high');            
    end
end
% low-pass coeffs
if fL ~= 0
    Wn = fL/(fs/2);
    [Bl, Al] = butter(6, Wn, 'low');
end
% --- Filtering ---
data = EEGout.data(:,:)';
data = double(data);
if fH ~= 0
    data = filter(Bh, Ah, data);
end
if fL ~= 0
    data = filter(Bl, Al, data);
end
EEGout.data(:,:) = single(data)';
end
          