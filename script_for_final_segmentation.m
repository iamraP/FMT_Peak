%%% Skript Rodrigues: Made by Dr. rer. nat. Johannes Rodrigues, Dipl. Psych. Julius-Maximilians University of W�rzburg. johannes.rodrigues@uni-wuerzburg.de; Started 2012, Latest update: 2021_05
%%% IMPORTATANT NOTE: THERE IS NO WARRENTY INCLUDED ! -> GNU 
%%% THERE ARE MANY THINGS THAT NEED TO BE ADJUSTED TO YOUR DATA !
%%% PLEASE ALSO KEEP IN MIND, THAT DIFFERENT MATLAB VERSIONS MIGHT HAVE SOME CRITICAL CHANGES IN THEM THAT MAY ALTER YOUR RESULTS !!! One example is the differences in the round function that changed the Baseline EEGLAB function on "older" MATLAB Version. 
%%% PLEASE DON�T USE THIS SCRIPT WITHOUT CONTROLLING YOUR RESULTS ! CHECK FOR PLAUSIBILITY OF THE SIGNAL AND TOPOGRAPHY

%script_for_final_segmentation
%This script is a simple example for setting the "CASEARRAY" that is used for the final segmentation in the processing script.
%Note that the CASEARRAY contains all relevant conditions. These are not necessarily the same as were used for the first segmentation.
%In this very basic example, we are only interested in Target markers, that are after event 1 or after event 2

 %%%"THE GREAT CHANGE": Adjust the marker names
%now tell the script the relevant marker names. 
%Relevant_Markers = {'S_T3_ambiguous_0' 'S_T3_ambiguous_1' 'S_T3_ambiguous_2' 'S_T3_ambiguous_3' 'S_T3_approach_0' 'S_T3_approach_1' 'S_T3_avoidance_0' 'S_T3_avoidance_1' 'S_T3_conflict_0' 'S_T3_conflict_1'};
for i=1:size(EEG.event,2)
	if contains(EEG.event(1,i).type,'ambiguous')
        EEG.event(1,i).type = 'ambiguous';
    elseif contains(EEG.event(1,i).type,'approach')
        EEG.event(1,i).type = 'approach';
    elseif contains(EEG.event(1,i).type,'avoidance')
        EEG.event(1,i).type = 'avoidance';
    elseif contains(EEG.event(1,i).type,'conflict')
        EEG.event(1,i).type = 'conflict';
    end
end 


Relevant_Markers ={'ambiguous','approach','avoidance','conflict'};
%see whether there are some double or triple mentioned markers because one might get carried away if there are many conditions...
Relevant_Markers = unique(Relevant_Markers);

%Create the Casearray:
CASEARRAY = Relevant_Markers;