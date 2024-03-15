%%% Processing Skript Rodrigues: Made by Dr. rer. nat. Johannes Rodrigues, Dipl. Psych. Julius-Maximilians University of Würzburg. johannes.rodrigues@uni-wuerzburg.de; Started 2012, Latest update: 2021_05
%%% Important Inputs were given by: Prof. John J. B. Allen, PhD, Prof. Johannes Hewig
%%% Lots ot parts or programming logic were taken from scripts by Prof. John J. B. Allen (I am not able to clearly disentangle where his ideas and input ended and where my part beginns...) Thanks a lot !!! 
%%% Important steady input over the years was also given by the Wintersymposium Montafon.
%%% Other important input was given by Janir Nuno Ramos da Cruz and Makoto Miyakoshi as well as Nathan Fox in the final stages of the script
%%% IMPORTATANT NOTE: THERE IS NO WARRENTY INCLUDED ! -> GNU 
%%% PLEASE SCROLL DOWN AND ADJUST THE PATHS AND THE DATA IMPORT FUNCTION ACCCORDING TO YOUR EEG FILE FORMAT!!!
%%% PLEASE ALSO SCROLL DOWN AND ADJUST SOME LINES ACCORDING TO YOUR EEG MONTAGE !!!
%%% THERE ARE MANY THINGS THAT NEED TO BE ADJUSTED TO YOUR DATA !
%%% FIND ALL LINES THAT NEED TO BE ADJUSTED BY THE SEARCH TERM: "THE GREAT CHANGE"
%%% PLEASE ALSO KEEP IN MIND, THAT DIFFERENT MATLAB VERSIONS MIGHT HAVE SOME CRITICAL CHANGES IN THEM THAT MAY ALTER YOUR RESULTS !!! One example is the differences in the round function that changed the Baseline EEGLAB function on "older" MATLAB Version. 
%%% Therefore some steps that are implemented in EEGlab are done "by hand" in this script.
%%% PLEASE DON´T USE THIS SCRIPT WITHOUT CONTROLLING YOUR RESULTS ! CHECK FOR PLAUSIBILITY OF THE SIGNAL AND TOPOGRAPHY
%%% IF SOMETHING WENT WRONG, CHECK ALSO THE PREPROCESSING AND SEGMENTATION ! VERY ! CAREFULLY !


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Processing chain Rodrigues: Adaptation and automatization using Matlab 2011b / 2015b
%%% EEGlab (Delorme & Makeig, 2004) and code provided by Cohen, 2014
%%%%
%%%%%%  PLEASE REMEMBER: GARBAGE IN -> GARBAGE OUT ! MAKE A GOOD AND CLEAN EEG DATA RECORDING ! TAKE YOUR TIME TO GET THE ELECTRODES RIGHT ! PLACE THE ELECTRODE CAP RESPONSIBLY AND WITH GREAT CARE !
%%%%%	TAKE
%%%%% 	YOUR
%%%%%	TIME 
%%%%%	!!!!
%%
%% Processing steps:
%%
%% 9. Segment the data for analysis and create a 4 Dimensional matrix for signal and each frequency separately and a 5 Dimensional matrix for the respective single trials (of course one can also create a 5 D matrix or 6 D matrix...) 
%%		This approach has the problem of sometimes running into memory problems as the matrices are rather large. 
%%		If this happens, then please either only compute one matrix at a time or start resampling the data in the preprocessing (not only in the processing script). If you do so, then remember that the resampling should be done without unnecessary extrapolation (meaning that you should resample in a divisor of the current sampling rate)
%%		In order to go into single trial analysis of frequencies later, an adjustment was made to the frequency extraction functions provided by Cohen (2014). Otherwise, the functions are used as shown in Cohen (2014)
%%   
%%
%%      9a. Segment the data according to relevant markers
%%		    The markers need to be selected and a segmentation-file should be provided. For this chain, an example segmentation-script will be provided (in another script).
%%          Put the segmentation-script in the EEGlab directory. In the segmentation-script, the Cellarray "CASEARRAY" contains all the conditions of the experiment.
%%          Keep in mind, that the data matrix that is used for automatic peak detection and visualization assumes that all conditions are equally often present. If this is not the case, use the weighted averages examples.     
%%          Step 5. revisited: In some processing chains (mostly ERP related) there is a second bad segment detection step at this point (see Step 5 pre-processing script). If you want to perform this step with different parameters, be sure to mention them or just perform it with the same parameters than before.
%%		
%%      9b. Define the Baseline automatically without using EEGlab function 
%%		    This is done to avoid problems with the round function in differnt Matlab version that may or may not work correctly with the respective EEGlab versions and therefore may lead to a correct EEGlab baselinefunction or not.
%%          The recommendation is to chose an appropriate baseline dependent on the event one is regarding: On pure feedback a baseline from - 500 or -200 to 0 might be a good idea, while a motor reaction normally does not have a respective baseline to 0 but to ~ -200 ms due to premotor activations. Mind also, that frequency analysis might lead to some jittering in time.
%%
%%      9c. Choose the frequencies of interest
%%		    In this chain, I think you know what frequencies and ERP signals you are interested in and only want to look at them, because we first want to gather evidence concerning our hypotheses. (Of course you can look on other frequency bands for exploratory purposes)
%%          
%%	    9d. Choose a filter (if wanted for ERPs). Depending on the ERP you are interested in, you are of course familiar with the filters (if any) you need. 
%%
%%	    9e. Choose whether you want to look on and analyze single trials. I recommend using multilevel modelling if single trial analysis is wanted. Also I recommend using frequencies instead of raw (erp-related) data, because of the higher reliablity (see Rodrigues et al. 2020)
%%
%%
%% 10. (optional): Loose cases that are not present but in your segmentation file. This is relevant for free choice paradigms.
%%
%%
%% 11. Automatically detect the peak in a given time window in EEG signal (here example FRN)
%%		The peak is searched in a time window of interest on an electrode of interest for the mean signal. The respective parameters depend on the ERP. Please consult the literature but also criticize the literature (e.g. still looking for FRN on Fz might seem not appropriate in many cases if FCz is available due to more than 32 electrodes)
%%		Depending on your task, it might be worth to look at the total mean of all cases, or just to look for the peak in several cases that are different from others. Examples are provided, but have to be adjusted to your tasks.
%%
%%
%% 12. Create ERPs
%%      In this step, ERPs (figures) are created. 
%%
%%      12a: "Normal" ERP: This ERP is just a line as in the(ancient) manuscripts... follow the tradition (?)
%%
%%      12b: In addition to the "normal" ERP, also ERPs with shaded errorlines are created in order to give an impression of the distribution. In this step, only between errorlines are drawn to the figure.
%%
%%      12c: In addition to the "normal" ERP, also ERPs with shaded errorlines are created in order to give an impression of the distribution. In this step, mean within errorlines are drawn to the figure (mean within SE).
%%          Note that you need to display meaningful conditions and maybe also compute meaningful bundles of conditions (not shown in this example, but simply use nanmean)
%%
%% 13. Create topographical maps (Topoplots)
%%      
%%      13a: In this step, a topographical map of the time window of interest (peak-window) is made. 
%%           
%%      13b: Also, a GIF is created that shows different timesteps in order to see the dynamic changes in activation in the topography and verify the time-window of interest for the electrodes.
%%
%%
%% 14. Automatically detect peak in a given time window in frequency responses (here example midfrontal theta)
%%		The peak is searched in the time window of interest on an electrode of interest the frequency response. The respective parameters are the same as specified for the ERP of interest. Please consult the literature but also criticize the literature (e.g. still looking for peak midfrontal theta on Fz might seem not appropriate in many cases if FCz is available due to more than 32 electrodes)
%%		Depending on your task, it might be worth to look at the total mean of all cases, or just to look for the peak in several cases that are different from others. Examples are provided, but have to be adjusted to your tasks.
%%
%%
%% 15. Create topographical maps for frequency response(Topoplots)
%%      
%%      15a: In this step, a topographical map of the time window of interest (peak-window of frequency response) is made. 
%%           
%%      15b: Also, a GIF is created that shows different timesteps in order to see the dynamic changes in activation in the topography and verify the time-window of interest for the electrodes.
%%
%%
%% 16. Create a time-frequency plot for a specific electrode in a broad frequency window
%%      For this time frequency plot, once again a function that is based on the code of Mike X Cohen (2014) is used. It was edited by John J.B. Allen and me. This provides a log transformed power output
%%      This function gets you a time frequncy plot that displays the frequency reaction not only limited to your functional frequency, but it is meant to be used as a display of larger frequency windows. (Suggestion: 1-30 Hz, as I am filled with scepticism concerning gamma, since microsaccades have been discovered)
%%      Note that a dB to baseline output is recommended here in order to correct for the power law and therefore provide an adequate visualization of the frequencies in their spectrum related to the chosen baseline. The baseline must be part of the chosen display window.
%%      Note also, that this transformation is not done to the data per se, but only to your visualization. But if you are only interested in the frequency responses of your respective bands and do not compare power values of different frequency bands, the baseline standardization is not necessary.  
%%
%%
%% 17. Export the data to statistical software: As xlsx or matlab file. 
%%
%%		17a: Long format export for the mean signal/frequencies. This is for example used by many R packages to calculate anovas... Note that additional information is exported in other files (.txt and .mat)(like the channel names)
%%		
%%      17b: Wide format export for the mean signal/frequencies. This is for example used Jamovi or PSPP/SPSS to calculate anovas...Note that additional information is exported in other files (.txt and .mat)(like the channel names)
%%		
%%      17c: Long format export for single trial signal/frequencies. This is for example used by many R packages and SPSS to calculate multilevel mixed models. Note that additional information is exported in other files (.txt and .mat)(like the channel names)
%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% END OF PROCESSING RODRIGUES "THEORY"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING CHAIN RODRIGUES: STEP 0: PREPARATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
eeglab
close all

%"THE GREAT CHANGE":If you are not in the EEGlab directory and have not initialized it, then please do so by typing "eeglab" in matlab

%"THE GREAT CHANGE": These two Matlab Addons need to be used in order to have nice pictures of the effects and an easy way to print them. Modify to your own path:
addpath(genpath('C:\MATLAB\Plugins\export_fig'))
addpath(genpath('C:\Software\MATLAB\Plugins\boundedline-pkg'))

%read in Data (still to modify)
%"THE GREAT CHANGE"
read_dir = 'C:\data\Results\removed_automatically\ICA\automatically_rejected\CSD\';   %There is your EEG preprocessed data with the respective reference you are interested in. 

%specify where to save the .txt files of the necessary information:
printpath = 'C:\MATLAB\Plugins\eeglab2023.1';



%"THE GREAT CHANGE":Please load one eegfile or specify the sampling rate.
load('montage_for_topoplot_CSD_CoScience.mat') % montage file that should be safed during STEP 8 in preprocessing to ica file
EEG.srate = 250;
srate = EEG.srate 

%"THE GREAT CHANGE":Please load the previous information from preprocessing
load('Evaluation.mat') % evaluation file that should be safed at the End of preprocessing file




%"THE GREAT CHANGE": Insert here your desired Baseline parameters. Please consider that the start of the segment is the start of the final segmentation and not the preprocessing segmentation
%Manual Baseline: Dependent on the sampling rate and the segmentation
baselinestart = ((1000-500)/(1000/EEG.srate))+1; %from -500, as the segments start here from -1000, sampling rate set above
baselinestop =	((1000-0)/(1000/EEG.srate))+1; % to 0, as the segments start here from -1000, sampling rate set above

%"THE GREAT CHANGE": Change the values of the final segmentation (given in seconds). Here it is -1s to 2 seconds
Segementation_time_start = -1;
Segementation_time_end = 2;

%"THE GREAT CHANGE" : Adjust the number of maximal trials in a condition.
%SINGLETRIAL: Max trials in a conditions:
MAXTRIALS = 26; % Important to set this rather a bit higher than to low. if the matrix is extended later on it adds zeros and this will messup the NaN mean,( also i assume it will take more processing power, but in the end idk) This variable stores the maximum of trial that can be seen in any condition. The variable is needed in case one wants to go for single trial analysis as the matrix should be precomputed in order to save memory.

%"THE GREAT CHANGE": Change the channel number of your electrodes of interest.
%Elektrode of interest: Please select the electrode you are interested in (for ERPs). The examples here are provided for the electrodes Fz FCz and Cz. Of course you can add more electrodes of interest...
electrode_of_interest1 = 2; %here example Fz in Marias montage



%"THE GREAT CHANGE": Adjust the target search window of your ERP dependent on the relevant literature. The search window depends on your segmentation time_start. 
searchwindowlength_in_ms = 200; % 100ms (see Yeung & Sanfey, 2004, FRN)
searchwindowstart_from_zero_in_ms = 250; % 200ms (see Yeung & Sanfey, 2004, FRN)
%searchwindow
searchwindowlength = round(searchwindowlength_in_ms/(1000/EEG.srate)); % 100ms dependent on the sampling rate
searchwindowstart = round(searchwindowstart_from_zero_in_ms/(1000/EEG.srate)-(Segementation_time_start*EEG.srate)); % 200ms dependent on the starting point of segement and sampling rate


%"THE GREAT CHANGE": Adjust the target time window of your ERP dependent on the relevant literature. The time window will be centered around the automatically detected peak. Keep in mind to check if you got different tasks that are compared whether the individual task is still in the window. Else adjust window length and also peak location.
%windowlength
windowlength_in_ms = 100; %40 ms (see e.g. Rodrigues et al. 2020): 
%windowlength in points:
windowlength = windowlength_in_ms/(1000/EEG.srate); % Your extraction time window for the FRN / ERP you are interested in. This is a peak centered time window. It is dependening on the sampling rate

%"THE GREAT CHANGE": Adjust the time to be displayed in the ERPs (STEP 12) and in the broad time frequency plot (STEP 16)
%Display parameters for ERP: Set paramter
display_time_start_from_zero_in_ms = -200; %Note that one can also go into - : this means it is left from 0, being a display of the baseline. Note also that you cannot display data that is not in your segmentation.
display_time_end_from_zero_in_ms = 800; %Note also that you cannot display data that is not in your segmentation.
%display parameter for ERP:
display_time_start = display_time_start_from_zero_in_ms/(1000/EEG.srate)-(Segementation_time_start*EEG.srate); % dependent on the starting point of segement and sampling rate
display_time_end = display_time_end_from_zero_in_ms/(1000/EEG.srate)-(Segementation_time_start*EEG.srate); % dependent on the starting point of segement and sampling rate

%"THE GREAT CHANGE": Adjust the timesteps for the topographical map GIF
%Choose timesteps for topographical maps gif
timesteps_in_ms = 40;
timesteps = timesteps_in_ms/(1000/EEG.srate);

%"THE GREAT CHANGE": Adjust the time to be displayed in the topographical map GIF
%Display parameters for topoplot GIF: Set paramter (note that this is for ERP as well as frequency responses: STEP 13 and STEP 15)
display_time_start_from_zero_in_ms_GIF_topo = 0; %Note that one can also go into - : this means it is left from 0, being a display of the baseline. Note also that you cannot display data that is not in your segmentation.
display_time_end_from_zero_in_ms_GIF_topo = 600; %Note also that you cannot display data that is not in your segmentation.
%display parameter for ERP: 
display_time_start_GIF_topo = display_time_start_from_zero_in_ms_GIF_topo/(1000/EEG.srate)-(Segementation_time_start*EEG.srate); % dependent on the starting point of segement and sampling rate
display_time_end_GIF_topo = display_time_end_from_zero_in_ms_GIF_topo/(1000/EEG.srate)-(Segementation_time_start*EEG.srate); % dependent on the starting point of segement and sampling rate

%"THE GREAT CHANGE": Adjust the frequencies displayed in the broad time frequency plot. (in Hz)
%Set the frequencies for the broad frequency plot (STEP 16)
min_freq = 1;
max_freq = 30;


%select the files that are in the relevant directory:
files = dir([read_dir '*.set']);
filenames = {files.name}';
%How many participants / files are there ?
COUNTPARTICIPANT = size(filenames,1);

%"THE GREAT CHANGE" : Here you need to insert your script for the final segmentation. For inexperienced users, I recommend using EEGlab with the "eegh" commmand:
%(Edit -> extract epoches -> select your relevant markers in "time locking event types", think about baseline and buffer times in "epoche limits" in s)
%CAREFUL: SOME MATLAB VERSIONS ARE NOT COMPATIBLE WITH THE EEGLAB BASELINE CORRECTION (pop_rmbase(...)) ! IN THIS CASE, INSERT THE BASELINE CORRECTION MANUALLY AS DONE HERE !	
%Please also keep in mind that it may be the case that some conditions are not met by all participants.
%Therefore I created a "demo" script of a final segmentation: Named script_for_final_segmentation
%in this script, the Variable "CASEARRAY" defines the final segmentation markers
%load wthe final segmentation file in order to get the number of relevant cases that one is interested in:
script_for_final_segmentation % this file needs to be put into the EEGlab Folder

%"THE GREAT CHANGE": Adjust the matrices that you need. Select only the matrices that you need because of working memory (see above)
%Creating all needed arrays: Here only mean_signal array and mean_theta array, but other arrays are possible (e.g. signle trial array shown for signal here)
% Total_mean_signal_array(:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts, 'single');	%adjust number of points	    %4D array: VP,CASES,ELECTRODES,TIMES
% Total_mean_theta_array(:,:,:,:) =  nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts, 'single');			%4D array: VP,CASES,ELECTRODES,TIMES
%Total_mean_alpha_array(:,:,:,:) =  nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts, 'single');			%4D array: VP,CASES,ELECTRODES,TIMES
%Total_mean_beta_array(:,:,:,:) =   nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts, 'single');			%4D array: VP,CASES,ELECTRODES,TIMES
%Total_mean_delta_array(:,:,:,:) =  nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts, 'single');			%4D array: VP,CASES,ELECTRODES,TIMES
%Total_mean_gamma_array(:,:,:,:) =  zeros(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts, 'single');		%4D array: VP,CASES,ELECTRODES,TIMES
% Total_signal_array(:,:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts,MAXTRIALS, 'single');						%5D array: VP,CASES,ELECTRODES,TIMES,TRIALS
% Total_theta_log_array(:,:,:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,9,EEG.pnts,MAXTRIALS, 'single');						%6D array: VP,CASES,ELECTRODES, FREQS,TIMES,TRIALS  %% Attention! number of frequencies for the wavlet is hhardcoded - prbly needs to be adjsuted to the band and the settings
% Total_theta_lin_array(:,:,:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,9,EEG.pnts,MAXTRIALS, 'single');						%6D array: VP,CASES,ELECTRODES,TIMES,TRIALS

%Create an array only to visually quickly check whether a case / condition is given in a participant and if so how many times.
Casecheck (:,:) = zeros(COUNTPARTICIPANT,size(CASEARRAY,2), 'single');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING CHAIN RODRIGUES: STEP 9: Segmentation of the data for analysis and create a 4 Dimensional matrix for signal and each frequency separately and a 5 Dimensional matrix for the respective single trials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrays_initialized = false;
file_info = nan(COUNTPARTICIPANT,2); %2D array info about VP and Session number
for VP = 1:COUNTPARTICIPANT
    %get the info
    vp_exp = '[0-9]{3}(?=_s[0-9]{2}_r[0-9]{2}_TMaze.set)';
    matchStr_vp = regexp(filenames{VP},vp_exp,'match');
    vp_num = str2num(matchStr_vp{1});
    sess_exp = '(?<=[0-9]{3}_s)[0-9]{2}(?=_r[0-9]{2}_TMaze.set)';
    matchStr_sess = regexp(filenames{VP},sess_exp,'match');
    sess_num = str2num(matchStr_sess{1});
    
    file_info(VP,:)= [vp_num,sess_num];

	EEG = pop_loadset('filename',filenames{VP},'filepath',read_dir);									%load set (first time)	-> Reason for this here: in case of error, the file is not loaded every case but only if something is done correctly
	
	for CASES = 1:size(CASEARRAY,2) 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 9a. Segment the data according to relevant markers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%Making the Event segementation
		script_for_final_segmentation % script_for making Casearray -> put in eeglab directory
        
		try % see whether the relevant segements are there... else do the next iteration
            EEG = pop_epoch( EEG, {CASEARRAY{CASES}}, [Segementation_time_start Segementation_time_end ], 'newname', strcat(filenames{VP},CASEARRAY{CASES}), 'epochinfo', 'yes'); %selection here: -1 to 2 seconds 
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
            EEG = eeg_checkset( EEG );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % END OF STEP 9a.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 9b. Define the Baseline automatically
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Here the baseline is taken from the definded baseline start to the end of the baseline (parameters definded above)
            %Note that there is no separate baseline taken for the single trials vs. the mean trial, but the same baseline is taken for all approaches.
            %automatical Baseline total file:
            for i = 1:size(EEG.data,1)
                for j = 1:size(EEG.data,3)
                        EEG.data (i,:,j) = EEG.data (i,:,j) - nanmean(EEG.data(i,baselinestart:baselinestop,j),2);
                end
            end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % END OF STEP 9b.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 9c. Choose your frequency bands
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%"THE GREAT CHANGE": Change the desired frequencies and choose what frequency you are interested in. Also adjust the relevant cycles for the wavelets. (For beginners I recommend reading Cohen, 2014)
            %%%Note: The files "wavelet_power" and "wavelet_power2" need to be in the eeglab directory. They are from Cohen, 2014, commented bei John J.B. Allen or a slightly modified version (power2). Also note their output commented below.
            %%%Another possibility for frequency extraction is of course the FFT-function. Problem of FFT: Accurate in frequency, not in time -> whole time window; If parameters set right, wavelet = FFT
            %Frequency extraction:
            %Alpha = wavelet_power_2(EEG,'lowfreq', 8, 'highfreq', 13, 'log_spacing', 1, 'fixed_cycles', 3.5); % 3 dim array = Channel x Time x Trials
            % log spaced freq bins? -> log besser für Funktion, linear besser für Interpretation, besser zu Rechtfertigen 
            [Theta_log, freq_idx_log] = wavelet_power_3(EEG,'lowfreq', 4, 'highfreq', 8, 'log_spacing', 1, 'fixed_cycles', 3.5); % 4 dim array = Channel x Freqs x Time x Trials
            % linear freq bins
            [Theta_lin, freq_idx_lin] = wavelet_power_3(EEG,'lowfreq', 4, 'highfreq', 8, 'log_spacing', 0, 'fixed_cycles', 3.5); % 4 dim array = Channel x Freqs x Time x Trials



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 9d. Choose your filters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Filter according to your need...
            %EEG = pop_eegfiltnew(EEG, 1, 20, [], 0, [], 0); %%comment: Use filters if you need them, if you don´t -> don´t ! 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % END OF STEP 9d. (strangely befor 1c ends)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%"THE GREAT CHANGE": Choose your desired mean signals and frequencies
            % Total_mean_signal_array(VP,CASES,1:size(EEG.data,1),1:size(EEG.data,2))   = single(nanmean(EEG.data,3));	    %4D array: VP,CASES,ELECTRODES,TIMES
            % Total_mean_theta_array(VP,CASES,1:size(EEG.data,1),1:size(EEG.data,2))    = single(nanmean(Theta,3));			%4D array: VP,CASES,ELECTRODES,TIMES
            %Total_mean_alpha_array(VP,CASES,size(1:EEG.data,1),1:size(EEG.data,2))   = single(nanmean(Alpha,3));		    %4D array: VP,CASES,ELECTRODES,TIMES
            %Total_mean_beta_array(VP,CASES,size(1:EEG.data,1),1:size(EEG.data,2))    = single(nanmean(Beta,3));			%4D array: VP,CASES,ELECTRODES,TIMES
            %Total_mean_delta_array(VP,CASES,size(1:EEG.data,1),1:size(EEG.data,2))   = single(nanmean(Delta,3));		    %4D array: VP,CASES,ELECTRODES,TIMES
            %Total_mean_gamma_array(VP,CASES,size(1:EEG.data,1),1:size(EEG.data,2))   = single(nanmean(Gamma,3));		    %4D array: VP,CASES,ELECTRODES,TIMES
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % END OF STEP 9c. (strangely after 1d ended)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 9e Choose whether you want to look at single trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % initialize the arrays once in the correct size:  
            if ~ arrays_initialized
                Total_signal_array(:,:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,EEG.pnts,MAXTRIALS, 'single');						%5D array: VP,CASES,ELECTRODES,TIMES,TRIALS
                Total_theta_log_array(:,:,:,:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,size(freq_idx_log,2),EEG.pnts,MAXTRIALS, 'single');						%6D array: VP,CASES,ELECTRODES, FREQS,TIMES,TRIALS  %% Attention! number of frequencies for the wavlet is hhardcoded - prbly needs to be adjsuted to the band and the settings
                Total_theta_log_array_db(:,:,:,:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,size(freq_idx_log,2),EEG.pnts,MAXTRIALS, 'single');						%6D array: VP,CASES,ELECTRODES, FREQS,TIMES,TRIALS  %% Attention! number of frequencies for the wavlet is hhardcoded - prbly needs to be adjsuted to the band and the settings
                Total_theta_lin_array(:,:,:,:,:,:,:) = nan(COUNTPARTICIPANT,size(CASEARRAY,2),EEG.nbchan,size(freq_idx_lin,2),EEG.pnts,MAXTRIALS, 'single');
        		arrays_initialized = true;

            end

            %%%"THE GREAT CHANGE": Choose in case of single trial the appropriate frequencies:
            %IN CASE OF SINGLE TRIALS:     
            %check if only one trial:
            if size(EEG.data,3) == 1
            	Total_signal_array(VP,CASES,1:size(EEG.data,1),1:size(EEG.data,2),1)  = single(EEG.data);                                                                %5D array: VP,CASES,ELECTRODES,TIMES,TRIALS
            	Total_theta_log_array(VP,CASES,1:size(EEG.data,1),1:size(freq_idx_log,2),1:size(EEG.data,2),1,1:size(freq_idx_log))   = single(Theta_log);               %6D array: VP,CASES,ELECTRODES,FREQS,TIMES,TRIALS
            	Total_theta_log_array_db(VP,CASES,1:size(EEG.data,1),1:size(freq_idx_log,2),1:size(EEG.data,2),1,1:size(freq_idx_log))   = single(theta_db_transform1);   
            	Total_theta_lin_array(VP,CASES,1:size(EEG.data,1),1:size(freq_idx_lin,2),1:size(EEG.data,2),1,1:size(freq_idx_lin))   = single(Theta_lin);               %6D array: VP,CASES,ELECTRODES,FREQS,TIMES,TRIALS
            else
	            Total_signal_array(VP,CASES,1:size(EEG.data,1),1:size(EEG.data,2),1:size(EEG.data,3))  = single(EEG.data);	                                                        %5D array: VP,CASES,ELECTRODES,TIMES,TRIALS, FREQ
            	Total_theta_log_array(VP,CASES,1:size(EEG.data,1),1:size(freq_idx_log,2),1:size(EEG.data,2),1:size(EEG.data,3))   = single(Theta_log);               %6D array: VP,CASES,ELECTRODES,FREQS,TIMES,TRIALS
            	Total_theta_lin_array(VP,CASES,1:size(EEG.data,1),1:size(freq_idx_lin,2),1:size(EEG.data,2),1:size(EEG.data,3))   = single(Theta_lin);         %6D array: VP,CASES,ELECTRODES,FREQS,TIMES,TRIALS
            end

            

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % END OF STEP 9e
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Count the trials in the conditions and create a casecheck array
            try 
                if size(EEG.data) == [EEG.nbchan,EEG.pnts]
                    Casecheck (VP,CASES) =  single(1);                                  %2D array: VP,CASES
                end
            end
            try
                if size(EEG.data,3)>1
                    Casecheck (VP,CASES) = single(size(EEG.data,3));                    %2D array: VP,CASES
                end				
            end		
        
            %clear all temporary variables that could be used in this condition
            clear Alpha Alpha2 Theta Theta2 Beta Beta2 Gamma Gamma Delta Delta2
		
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];	    %clear the EEG sets
            EEG = pop_loadset('filename',filenames{VP},'filepath',read_dir);       %reload data here if something was cut from it
		end %try end: If this condition can not be found, then simply start from here -> next condition
	end
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];		    %clear the EEG sets
end

%As this takes some time: Save the results in the end. 
save('Backup_workspace', '-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF STEP 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING CHAIN RODRIGUES: STEP 10: Loose unnecessary cases: THIS IS OPTIONAL ! "THE GREAT CHANGE": comment this step out if not needed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%"THE GREAT CHANGE": Choose your desired matrices
%control the size of your matices
size(Total_signal_array)
size(Total_theta_log_array)
size(Total_theta_lin_array)

%now we got massive 4d arrays: Loose now the unnecessary cases that are not present in none of the participants.
%find the stuff to delete:
delete_array = [];
j =1;
for i = 1:size(Casecheck,2)
	if mean(Casecheck(:,i)) == 0
		delete_array(j) = i;
		j=j+1;
	end
end

NEWCASEARRAY = CASEARRAY(:,:,:,:) ;

%delete stuff...
NEWCASEARRAY(:,delete_array,:,:) = [];

%%%"THE GREAT CHANGE": Choose your desired matrices
%delete unnecessary cases:

Total_signal_array(:,delete_array,:,:,:) = [];
Total_theta_log_array(:,delete_array,:,:,:) = [];
Total_theta_lin_array(:,delete_array,:,:,:) = [];

Casecheck(:,delete_array) = []; 


%%%"THE GREAT CHANGE": Choose your desired matrices
%control the size of your matices again !
size(Total_signal_array)
size(Total_theta_log_array)
size(Total_theta_lin_array)
size(Casecheck)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF STEP 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING CHAIN RODRIGUES: STEP 14: Detect frequency peaks in relevant frequency bands on relevant electrode positions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%"THE GREAT CHANGE": Comment in / replicate with other electrodes if you want to have more electrodes... this is true for the entire step and following steps. This may not be mentioned again...
%%%"THE GREAT CHANGE": Comment in / replicate with other frequency bands if you want to have different frequency bands... this is true for the entire step and following steps. This may not be mentioned again...
%%%"THE GREAT CHANGE": Note that the timewindow for frequency responses might not be the same as for ERPs also due to time jittering
%%"THE GREAT CHANGE": If your conditions are not equally often in there, the mean has to be weighted on count of trials for an adequate peak detection
%for i = 1:size(Total_mean_theta_array,1)
%   for j = 1:size(Total_mean_theta_array,2)
%       wTotal_mean_theta_array(i,j,:,:)= Total_mean_theta_array(i,j,:,:).*Casecheck(i,j);
%   end
%end
%
%%%"THE GREAT CHANGE": If your conditions are not equally often in there, the mean has to be weighted on count of trials:
%lookforPeak1Freq = double(squeeze(nanmean(nanmean(wTotal_mean_theta_array(:,:,electrode_of_interest1,:),1),2))/mean(nonzeros(Casecheck'))); %mean over all participants and all coditions on electrode position 1

%%%"THE GREAT CHANGE": If you need to select some conditions, use this and fill in the conitions_of_interests. Note that this is only the example for the default, where the trials are not equally often
%example if only some conditions are ment to be seen together:
%lookforPeak1Freq = double(squeeze(nanmean(nanmean(wTotal_mean_theta_array(:,[conitions_of_interest1 conditions_of_interest2],electrode_of_interest1,:),1),2))/mean(nonzeros(Casecheck(:,[conitions_of_interest1 conditions_of_interest2])'))); %mean over all participants in some oditions on electrode position 1


 % calculate db bsl change

        
%%%%%%% Baseline Transform:
%%Total_theta_array: VP,CASES,ELECTRODES,FREQS,TIMES,TRIALS
% THETA = 4 dim array = Channel x Freqs x Time x Trials

baseline_mean_log = mean(Total_theta_log_array(:,:,:,:,baselinestart:baselinestop,:),5);
theta_db_transform_log = 10*log10(Total_theta_log_array);
theta_db_bsl_transform_log = 10*log10(Total_theta_log_array./baseline_mean_log);
 
baseline_mean_lin = mean(Total_theta_lin_array(:,:,:,:,baselinestart:baselinestop,:),5);
theta_db_transform_lin = 10*log10(Total_theta_lin_array);
theta_db_bsl_transform_lin = 10*log10(Total_theta_lin_array./baseline_mean_lin);

% create empty arrray which will be filled only for those trials where a
% peak was detected

%% create one 2D array to export, keep 1 6D- Matrix ? 
%Column Names      1    2         3        4    5                6                                          
% ExportDatframe: 1 VP | 2Session | 3CASE | 4trial| 5Freq_log | (Peak_time | Peak_amplitude ) for regular (6|7), in decible(8|9), in decible change to baseline(10|11) | 12 Freq_lin | (Peak_time | Peak_amplitude ) for regular (13|14), in decible(15|16), in decible change to baseline(17|18)

for peak_detection =1:4 %1: Center-of-gravity, Broaband Timing, 2: CoG, Frequencywise timed peak, 3: Peak_Window, Broadband Timining, 4: Peak_Window;  Frequencywise timed peak)
    
    data_pnts = COUNTPARTICIPANT*size(CASEARRAY,2)*MAXTRIALS*size(freq_idx_lin,2);
    export_df = nan(data_pnts,18);
    data_pnts_per_file = data_pnts/COUNTPARTICIPANT;
    data_pnts_per_case = data_pnts_per_file /size(Casecheck,2);
    data_pnts_per_trial = data_pnts_per_case / MAXTRIALS;
    
    
    
    peak_array= nan(COUNTPARTICIPANT,size(CASEARRAY,2),MAXTRIALS,size(freq_idx_log,2),12); %dimensions: VP, CASES,Trial, freq, peakinfo(2*2: timelocation + power for log an lin)
    idx = 1;
    sliding_window_length = round(50/1000*250); %50ms window
    sliding_window_step = round(10/1000*250); %10ms window step
    %iterate over 
    
    
    
    for VP = 1: COUNTPARTICIPANT %participants
        for CASES =1:size(Casecheck,2) %conditions
            for trial =1:Casecheck(VP,CASES) %trials
                freq_bins = size(freq_idx_log,2);
                export_df(idx:idx+freq_bins,1) =file_info(VP,1);
                export_df(idx:idx+freq_bins,2) =file_info(VP,2);
                export_df(idx:idx+freq_bins,3) =CASES;
                export_df(idx:idx+freq_bins,4) =trial;
                export_df(idx:idx+freq_bins-1,5) = freq_idx_log;
                export_df(idx:idx+freq_bins-1,12) = freq_idx_lin;
    
                for freq_space = ["log","lin", "log_db","lin_db","log_db_bsl","lin_db_bs"] % detect the peak in all conditons: regular, in decible, in decible change to baseline
                    if freq_space == "log"
                        theta_array = Total_theta_log_array;
                        freq_idx = freq_idx_log;
                        export_peak_idx = [6:7];
                        peak_idx = [1:2];                        
                    elseif freq_space == "log_db"
                        theta_array = theta_db_transform_log;
                        freq_idx = freq_idx_log;
                        export_peak_idx = [8:9];
                        peak_idx = [3:4];
                    elseif freq_space == "log_db_bsl"
                        theta_array = theta_db_bsl_transform_log;
                        freq_idx = freq_idx_log;
                        export_peak_idx = [10:11];
                        peak_idx = [5:6];
                    elseif freq_space == "lin"
                        theta_array = Total_theta_lin_array;
                        freq_idx = freq_idx_lin;
                        export_peak_idx = [13:14];
                        peak_idx = [7:8];
                    elseif freq_space == "lin_db"
                        theta_array = theta_db_transform_lin;
                        freq_idx = freq_idx_lin;
                        export_peak_idx = [15:16];
                        peak_idx = [9:10];
                    elseif freq_space == "lin_db_bs"
                        theta_array = theta_db_bsl_transform_lin;
                        freq_idx = freq_idx_lin;
                        export_peak_idx = [17:18];
                        peak_idx = [11:12];
                    end


                    if peak_detection ==1
                        %center of gravity
        
                        % find peak time in theta (4-8Hz); disassemble the peak in the
                        % contribution of the individual frequencies
        
                        single_trial = mean(squeeze(theta_array(VP,CASES,electrode_of_interest1,:,:,trial)),1);
                        max_amp = -Inf;
                        peak_start = 0;
                        for  start_idx = searchwindowstart:sliding_window_step:searchwindowstart+searchwindowlength-sliding_window_length
                            if mean(single_trial(1,start_idx:start_idx+sliding_window_length)) > max_amp
                                max_amp = mean(single_trial(1,start_idx:start_idx+sliding_window_length)); % not used later on
                                peak_start = start_idx;
                            end            
                        end
    
                        peak_amp = mean(squeeze(theta_array(VP,CASES,electrode_of_interest1,:,peak_start:peak_start+sliding_window_length,trial)),2);
                        peak_time = Segementation_time_start+ (peak_start+(25/(1000/srate)))/srate; 
    
        
                        export_df(idx:idx+freq_bins-1,export_peak_idx(2)) = peak_amp;
                        export_df(idx:idx+freq_bins-1,export_peak_idx(1)) = peak_time;
                    elseif peak_detection ==2
                        %center of gravity
        
                        % find peak time for each of the individual frequencies
    
                        single_trial = squeeze(theta_array(VP,CASES,electrode_of_interest1,:,:,trial));
                        max_amp(1:size(single_trial,1)) = -Inf;
                        peak_start(1:size(single_trial,1)) = 0;
                        for freq_bin = 1:size(single_trial,1)
                            for  start_idx = searchwindowstart:sliding_window_step:searchwindowstart+searchwindowlength-sliding_window_length
                                if mean(single_trial(freq_bin,start_idx:start_idx+sliding_window_length)) > max_amp(freq_bin)
                                    max_amp(freq_bin) = mean(single_trial(freq_bin,start_idx:start_idx+sliding_window_length)); % not used later on
                                    peak_start(freq_bin) = start_idx;
                                end            
                            end
                        end   

                        peak_time = Segementation_time_start+ (peak_start+(25/(1000/srate)))/srate; 
                        peak_amp = max_amp


                        export_df(idx:idx+freq_bins-1,export_peak_idx(2)) = peak_amp;
                        export_df(idx:idx+freq_bins-1,export_peak_idx(1)) = peak_time;
        
                    elseif peak_detection ==3
                        %Peak Window
        
                        % find peak time in theta (4-8Hz); disassemble the peak in the
                        % contribution of the individual frequencies
        
        
                        single_trial = mean(squeeze(theta_array(VP,CASES,electrode_of_interest1,:,:,trial)),1);
                        [peak_amp_mid, peak_mid] =  max(single_trial(:,searchwindowstart:searchwindowstart+searchwindowlength),[],2);
                        peak_start = round(peak_mid - (25/(1000/srate)) +searchwindowstart); % 50ms window around the peak
                        peak_end = round(peak_mid + (25/(1000/srate)) +searchwindowstart);
        
                        peak_amp = mean(squeeze(theta_array(VP,CASES,electrode_of_interest1,:,peak_start:peak_start+sliding_window_length,trial)),2);
                        peak_time = Segementation_time_start+ (peak_start+(25/(1000/srate)))/srate; 
                        
                        export_df(idx:idx+freq_bins-1,export_peak_idx(2)) = peak_amp;
                        export_df(idx:idx+freq_bins-1,export_peak_idx(1)) = peak_time;
        
                    elseif peak_detection ==4
                        %Peak Window
        
                        % find peak time for each of the individual frequencies
        
    
                        single_trial = squeeze(theta_array(VP,CASES,electrode_of_interest1,:,:,trial));
                        [peak_amp_mid, peak_mid] =  max(single_trial(:,searchwindowstart:searchwindowstart+searchwindowlength),[],2);
                        peak_start = round(peak_mid - (25/(1000/srate)) +searchwindowstart); % 50ms window around the peak
                        peak_end = round(peak_mid + (25/(1000/srate)) +searchwindowstart);
                        for i =1: length(peak_start)
                            peak_amp(i) =  mean(single_trial(i,peak_start(i):peak_end(i)),2);
                        end
                        peak_time = (searchwindowstart + peak_mid) /(srate)+ Segementation_time_start; 
        
        
                        export_df(idx:idx+freq_bins-1,export_peak_idx(2)) = peak_amp;
                        export_df(idx:idx+freq_bins-1,export_peak_idx(1)) = peak_time;
        
                    end
                            
    
                    % 
                    % for freq = 1: size(freq_idx_log,2) %freqs
                    % 
                    % lookforPeak1Freq_log = squeeze(theta_array_log(VP,CASES,electrode_of_interest1,freq,:,trial));
                    % max_peak_log = max(findpeaks(lookforPeak1Freq_log(searchwindowstart:searchwindowstart+searchwindowlength,1))); % find the peak in the window of interest
                    % %find peak as maximum point
                    % if ~ isempty(max_peak_log)
                    %     [r1freq_log,c1freq_log]=find(lookforPeak1Freq_log==max_peak_log);% get the time 
                    %     peak_array(VP,CASES,trial, freq,1:2) = [r1freq_log(1),max_peak_log(1)]; % fill the matrix
                    %     export_df(idx,6:7) = [r1freq_log(1),max_peak_log(1)];
                    % end
                    % 
                    % lookforPeak1Freq_lin = squeeze(theta_array_lin(VP,CASES,electrode_of_interest1,freq,:,trial));
                    % max_peak_lin = max(findpeaks(lookforPeak1Freq_lin(searchwindowstart:searchwindowstart+searchwindowlength,1))); % find the peak in the window of interest
                    % 
                    % if ~ isempty(max_peak_lin)
                    %     [r1freq_lin,c1freq_lin]=find(lookforPeak1Freq_lin==max_peak_lin);% get the time 
                    %     peak_array(VP,CASES,trial, freq,3:4) = [r1freq_lin(1),max_peak_lin(1)]; % fill the matrix
                    %     export_df(idx,9:10) = [r1freq_lin(1),max_peak_lin(1)];
                    % end


                 
                end
                idx =idx+freq_bins ;
            end
        end
    end
    
    % ExportDatframe: 1 VP | 2Session | 3CASE | 4trial| 5Freq_log | (Peak_time | Peak_amplitude ) for regular (6|7), in decible(8|9), in decible change to baseline(10|11) | 12 Freq_lin | (Peak_time | Peak_amplitude ) for regular (13|14), in decible(15|16), in decible change to baseline(17|18)
    header = {"vp","session","trial_type","trial","freq_log","peak_time_log","peak_amplitude_log","peak_time_log_db","peak_amplitude_log_db","peak_time_log_db_bsl","peak_amplitude_log_db_bsl","freq_lin","peak_time_lin","peak_amplitude_lin","peak_time_lin_db","peak_amplitude_lin_db","peak_time_lin_db_bsl","peak_amplitude_lin_db_bsl"};
    
    export_df = [header; num2cell(export_df)];
    
    if peak_detection ==1
        writecell(export_df,'CoG_broad_peak.csv') 
    elseif peak_detection ==2
        writecell(export_df,'CoG_ind_peak.csv') 
    elseif peak_detection ==3
        writecell(export_df,'Peak_Window_broad_peak.csv') 
    elseif peak_detection ==4
        writecell(export_df,'Peak_Window_ind_peak.csv')
    end

end 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF STEP 14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING CHAIN RODRIGUES: STEP 16: Create a time frequency plot in a broad frequency window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%"THE GREAT CHANGE": See whether you need to load the montage (again) because it contains necessary information for the next step like sampling rate 
%%load a montage to plot in. This montage should be saved during STEP 16 in preprocessing to ICA
%load('montage_for_topoplot.mat') % montage file

%%%"THE GREAT CHANGE": Choose the array according to your data structure
%Create an array for the display. Also create the weighted array, if necessary
%Total_freq_array(:,:,:,:) = nan(COUNTPARTICIPANT,size(NEWCASEARRAY,2),2*(max_freq-min_freq)+1,(display_time_end_from_zero_in_ms-display_time_start_from_zero_in_ms)/(1000/EEG.srate), 'single');	
%wTotal_freq_array(:,:,:,:) = nan(COUNTPARTICIPANT,size(NEWCASEARRAY,2),2*(max_freq-min_freq)+1,(display_time_end_from_zero_in_ms-display_time_start_from_zero_in_ms)/(1000/EEG.srate), 'single');	

%%%"THE GREAT CHANGE": Check whether all parameters are valid. Note that this step assumes that step 2 was performed (in order to avoid unnecessary conditions). If not, then just comment in the other line and comment out the "NEWCASEARRAY" lines
%This reloads the data for the relevant time window and performs the "big" frequency analysis. This is necessary if no single trial data is already present. If single trial data is present, you mal also start from the single trial array
for VP = 1:COUNTPARTICIPANT

	for CASES = 1:size(NEWCASEARRAY,2) %we already know what is there in the dataset
        EEG = pop_loadset('filename',filenames{VP},'filepath',read_dir);									%load set (first time)	-> Reason for this here: in case of error, the file is not loaded every case but only if something is done correctly
        script_for_final_segmentation
        %Segment the data according to relevant markers
		try % see whether the relevant segements are there... else do the next iteration
            EEG = pop_epoch( EEG, {NEWCASEARRAY{CASES}}, [Segementation_time_start Segementation_time_end ], 'newname', strcat(filenames{VP},NEWCASEARRAY{CASES}), 'epochinfo', 'yes'); %selection here: trial
            %EEG = pop_epoch( EEG, {CASEARRAY{CASES}}, [Segementation_time_start Segementation_time_end ], 'newname', strcat(filenames{VP},CASEARRAY{CASES}), 'epochinfo', 'yes'); %selection here: trial
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
            EEG = eeg_checkset( EEG );

            %Define the Baseline manually
            %Manual Baseline total file:
            for i = 1:size(EEG.data,1)
                for j = 1:size(EEG.data,3)
                        EEG.data (i,:,j) = EEG.data (i,:,j) - nanmean(EEG.data(i,baselinestart:baselinestop,j),2);
                end
            end




            EEG = pop_select( EEG,'time',[display_time_start_from_zero_in_ms/1000 display_time_end_from_zero_in_ms/1000] );
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 

            [GraphicsPlot frequ] = wavelet_power_plot(EEG,'lowfreq', min_freq, 'highfreq', max_freq, 'log_spacing', 1, 'fixed_cycles', 3.5); % 4 dim array = Channelx Freq x Time x Trials (Epochs)

            baseline_end = display_time_start_from_zero_in_ms/1000 *-250;
            baseline_mean = mean(GraphicsPlot(2,:,1:baseline_end,:),3); %average over timepoints of baseline for each VP,case,electrode, freq and trial
            db_transform= 10*log10(GraphicsPlot(2,:,:,:)./baseline_mean); %elementwise division of each datapoint by the corresponding average baseline
            case_average_db = squeeze(mean(db_transform,4));
            %%%"THE GREAT CHANGE": Select how many steps you want
            %Steps
            frequencystepparamter=30;
            %plot
            figure
            %%%"THE GREAT CHANGE": Change your stile of Plot ?
            contourf(EEG.times,frequ,case_average_db,frequencystepparamter,'linecolor','none')
            %imagesc(EEG.times,frequ,eegpower) %different stile
            %set(gca,'YDir','normal')          %different stile line 2: This is needed too, because of the y axis flip of imagesc


            %%%"THE GREAT CHANGE": Choose your scaling: Adjust it to the data: If you see all green, blue or all red, it is not appropriate. As a hint where to start from the mean of the frequency response is suggested. Also choose which frequency is chosen and whether you want raw power or log power
            %Choose the power limit. 
            %powerlimit = nanmean(nanmean(eegpower,1),2)*1.96; %inside-joke... but it is really an arbitraty suggestion... Note that with this formular, if the mean is negative (which could be) there will be an error !
            %own powerlimits
            powerlimit = 3;
            %%%"THE GREAT CHANGE": Choose linear or log spacing of the frequencies. Also adjust your display limit. In db to baseline, negative and positive is ok. Otherwise (if not to baseline), negative values will not appear
            %linear spacing plot
            set(gca,'clim',[-powerlimit powerlimit],'xlim',[display_time_start_from_zero_in_ms display_time_end_from_zero_in_ms],'xtick',display_time_start_from_zero_in_ms:(display_time_end_from_zero_in_ms-display_time_start_from_zero_in_ms)/5:display_time_end_from_zero_in_ms,'yscale','linear','ytick',min_freq:5:max_freq,'yticklabel',min_freq:5:max_freq)
            %%log spacing plot
            %%set(gca,'clim',[-powerlimit powerlimit],'xlim',[display_time_start_from_zero_in_ms display_time_end_from_zero_in_ms],'xtick',display_time_start_from_zero_in_ms:(display_time_end_from_zero_in_ms-display_time_start_from_zero_in_ms)/5:display_time_end_from_zero_in_ms,'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
            ylabel('Frequency (Hz)')
            xlabel('Time (ms)')
            %%%"THE GREAT CHANGE": Give it a a good title
            title_str = sprintf('TFR at Fz for VP%i Session%i %s trials', file_info(VP,1), file_info(VP,2),string(NEWCASEARRAY(CASES)))
            title(title_str )
            h = colorbar;
            %%%"THE GREAT CHANGE": Check your units: Is it a CSD or normal reference ? Remember that you made a frequency power transformation. Also remember whether you displayed log power or raw power.
            %ylabel(h, 'log power (µV²/m²)', 'Rotation',90)
            ylabel(h, 'dB change from baseline', 'Rotation',90)
            %%%"THE GREAT CHANGE": Change size of figure ?
            set(gcf, 'Units', 'points', 'OuterPosition', [0, 0, 400, 400]);
            set(gcf, 'color', [1 1 1])
            %%%"THE GREAT CHANGE": Name the file an change format or resolution
            filename_str =fullfile('.\data\\Results\\Plots\\', sprintf('TFR_Fz_VP%i_Sess%i_%s', file_info(VP,1), file_info(VP,2),string(NEWCASEARRAY(CASES))));
            export_fig(filename_str, '-tif', '-r400');
            %Total_freq_array(VP,CASES,:,:) =  squeeze(nanmean(GraphicsPlot(electrode_of_interest1,:,:,:),4)); % 4 D array: VP, Cases, Frequencies, times (on electrode of interest): be careful: this may be not weighted due to differnt count of trials in conditions
            %wTotal_freq_array(VP,CASES,:,:) =  squeeze(nanmean(GraphicsPlot(electrode_of_interest1,:,:,:),4).*Casecheck(VP,CASES));
            clear GraphicsPlot
		    close all
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];	    %clear the EEG sets
		end %try end: If this condition can not be found, then simply start from here -> next condition
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];		    %clear the EEG sets
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF STEP 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF PROCESSING CHAIN RODRIGUES: That´s all folks: Still to come in another script (postprocessing/premiumprocessing): Cross frequency coupling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Altman, Y (2020). export_fig (https://www.github.com/altmany/export_fig), GitHub. Retrieved March 9, 2020.
%%%%%% Cohen, M. X. (2014). Analyzing Neural Time Series Data Theory and Practice. Cambridge, Massachusetts, London, England.
%%%%%% Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of Neuroscience Methods, 134(1), 9–21. https://doi.org/10.1016/j.jneumeth.2003.10.009
%%%%%% Kearney, K. (2020). boundedline.m (https://www.github.com/kakearney/boundedline-pkg), GitHub. Retrieved March 9, 2020. 
%%%%%% Rodrigues, J., Liesner, M., Reutter, M., Mussel, P., & Hewig, J. (2020). It’s costly punishment, not altruistic: Low midfrontal theta and state anger predict punishment. Psychophysiology. https://doi.org/10.1111/PSYP.13557
%%%%%% Yeung, N., & Sanfey, A. G. (2004). Independent Coding of Reward Magnitude and Valence in the Human Brain. Journal of Neuroscience, 24(28), 6258–6264. https://doi.org/10.1523/JNEUROSCI.4537-03.2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
