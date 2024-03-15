library(tidyverse)
library(readr)
library(readxl)
library(SimplyAgree)
library(ggplot2)

# set the working directory
############### enter working directory where the files are located
setwd("C:/XXXXXX")


#Option A)
# 1.	Peak-Frequenz für jeden Trial
# 2.	MITTELN der Amplitude der Peaks für jedes Frequenzband
# 3.	Frequenz mit dem höchsten Peak als ITF

# Option B)
# 1.	Peak-Frequenz für jeden Trial
# 2.	SUMMIEREN der Amplitude der Peaks für jedes Frequenzband
# 3.	Frequenz mit dem höchsten Peak als ITF

# Option C) (meines Erachtens nicht wirklich sinnig..)
# 1.	Peak-Frequenz für jeden Trial
# 2.	Frequenz mit den meisten Trials mit höchstem Peak als ITF


# Option D)
# 1.	Mitteln der Amplitude der Peaks über die Trials für jedes Frequenzband
# 2.	Frequenz mit dem höchsten Peak als ITF

freq_space = c("log","lin", "log_db","lin_db","log_db_bsl","lin_db_bsl") # c("log","lin", "log_db","lin_db","log_db_bsl","lin_db_bsl")
#  Trial Types: 1: ambigous 2: approach 3: avoidance, 4: conflict
trial_types = c( "ambi", "approach", "avoidance","conflict", "conf_ambi")

peak_detection_algo = c('CoG_broad_peak','CoG_ind_peak','Peak_Window_broad_peak','Peak_Window_ind_peak') 

for (detect in 1:length(peak_detection_algo)){
  data_file = paste(peak_detection_algo[detect],".csv",sep="")

  for (i in 1:length(freq_space)){
    data <- read.csv(data_file,sep =",",fileEncoding="UTF-8-BOM", header = TRUE)
    # data$session[data$session== 1] <- "Initial Session (Day 0) "
    # data$session[data$session==8] <- "after 7 NF sessions (Day 10-14)"
    # data$session[data$session==9] <- "1 week later (Day 17-21 )"
    data$peak_amplitude_log_db <- data$peak_amplitude_log_db*-1
    data$peak_amplitude_lin_db <- data$peak_amplitude_lin_db*-1
    
    data$session <- factor(data$session)
    #data$session <- factor(data$session, levels=c("Initial Session (Day 0) ","after 7 NF sessions (Day 10-14)","1 week later (Day 17-21 )"))
    
    #initialize peak_df
    if(!exists("peak_df")){
      peak_df <- unique(data[,c('vp','session')])
      peak_df <- drop_na(peak_df)
    }
    
    
    data$vp <- as.factor(data$vp)
    
  
    if (i==1|i==2){ #log
      unit = "[µV/m²]"
      plot_lim_a = c(0,0.0025)
      plot_lim_b = c(0,0.7)
      plot_lim_d = plot_lim_a
    } else if(i==3|i==4){ #log_db
      unit = "[dB]*-1"
      plot_lim_a = c(0,50)
      plot_lim_b = c(400,1400)
      plot_lim_d = plot_lim_a
    } else if(i==5|i==6){ #lin_db_bsl
      unit = "[dB change to bsl]"
      plot_lim_a = c(0,20)
      plot_lim_b = c(0,150)
      plot_lim_d = plot_lim_a
    }
    # scalings lin / log : 0-20
    data$vp <- as.factor(data$vp)
    for (k in 1:length(trial_types)){
    
      if (k==5){
        data$trial_type[data$trial_type==4 | data$trial_type==1] = 5
      }     
      
      if (grepl("log",freq_space[i])){
        freq_bin <- "freq_log"
        bin_space = "log"
      } else if (grepl("lin",freq_space[i])){
        freq_bin <- "freq_lin"
        bin_space = "linear"
      }
      
      
      ################  Option A)  ###########################
      
      # 1.	Peak-Frequenz für jeden Trial
      # 2.	MITTELN der Amplitude der Peaks für jedes Frequenzband
      # 3.	Frequenz mit dem höchsten Peak als ITF
      
  
      peak_time <- paste("peak_time",freq_space[i],sep="_")
      peak_amp <- paste("peak_amplitude",freq_space[i],sep="_")
      
      current_df <- data[c("vp","session","trial_type","trial",freq_bin,peak_time,peak_amp)]
      colnames(current_df)<- c("vp","session","trial_type","trial","freq_bin","peak_time","peak_amp")
      
  
      
      max_avg_sess <-current_df %>% 
        drop_na(peak_amp) %>%  # get rid of trials where no peak was detected
        group_by(vp,session,trial_type,freq_bin)%>% 
        mutate(summed_amp=sum(peak_amp)) # get the sum amplitude for each freq in each trial_type of each session 
      
      num_trials <- count(max_avg_sess)  
      max_avg_sess<- max_avg_sess[c("vp","session","trial_type","freq_bin","summed_amp")]
      max_avg_sess  <-   max_avg_sess[!duplicated(max_avg_sess), ]
      max_avg_sess['n'] <- num_trials$n
      max_avg_sess['avg_amp'] <- max_avg_sess$summed_amp/max_avg_sess$n
      max_avg_sess <- max_avg_sess%>%  
        group_by(vp,session,trial_type)%>% 
        mutate(max_amp=max(avg_amp))%>% 
        subset(avg_amp== max_amp)
      
  
      conflict =  subset(max_avg_sess,max_avg_sess$trial_type==k)
      
      legend_title_size = paste("Average Amplitude\n", unit, sep=" ")
      
      ggplot(data = conflict,aes(x = session, y = freq_bin, group=vp, color= vp)) +
        ylim(4,8)+
        geom_point(aes(size=max_amp),alpha = 0.5)+
        geom_line()+
        labs(y= paste("Frequency [Hz]      ",bin_space, "- spcaced bins ", sep =" "))+
        #ggtitle(paste("A: Avg amp of trials where freq. was peak-freq.       Trials:", trial_types[k])) + 
        guides(size=guide_legend(title=legend_title_size))+
        scale_size_binned(
          limits = plot_lim_a,
          range = c(1, 12)
        )+ theme_classic()+     theme(
          # Set the plot title and axis titles styles
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text( size = 25),
          legend.title = element_text(size = 20),
          
          # Set the axis text styles
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20))
      
      ggsave(filename = paste("Results/Plots/Freq_Stability/Current/", trial_types[k],"_", freq_space [i], "_", peak_detection_algo[detect],"_A.png", sep=""),
                              plot = last_plot())
      
      if(k==5){
        peak_df[paste(freq_space[i],peak_detection_algo[detect],"A",sep="_")] = conflict["freq_bin"]
      }
      
      ################  Option B)  ###########################
      # 1.	Peak-"Sieger"-Frequenz für jeden Trial
      # 2.	SUMMIEREN der Amplitude der Sieger-Peaks für jedes Frequenzband
      # 3.	Frequenz mit dem höchsten Peak als ITF
      
      sum_max_log <-current_df %>% 
        drop_na(peak_amp) %>%  # get rid of trials where no peak was detected
        group_by(vp,session,trial_type,freq_bin)%>% 
        mutate(summed_amp=sum(peak_amp))%>%  
        group_by(vp,session,trial_type)%>% 
        mutate(max_amp=max(summed_amp))%>% 
        group_by(vp,session,trial_type)%>% 
        subset(summed_amp== max_amp)
      sum_max_log<- sum_max_log[c("vp","session","trial_type","freq_bin","summed_amp")]
      sum_max_log  <-   sum_max_log[!duplicated(sum_max_log), ]
      
     legend_title_size = paste("Summed Amplitude\n", unit, sep=" ")
      
  
      conflict =  subset(sum_max_log,sum_max_log$trial_type==k)
      
      
      ggplot(data = conflict,aes(x = session, y = freq_bin, group=vp, color= vp))+
        ylim(4,8)+ 
        geom_point(aes(size=summed_amp),alpha = 0.5)+
        geom_line()+
        labs(y= paste("Frequency [Hz]      ",bin_space, "- spcaced bins ", sep =" "))+
        #ggtitle(paste("B: Summed amp of trials where freq. was peak-freq.       Trials:", trial_types[k])) + 
        guides(size=guide_legend(title=legend_title_size))+
        scale_size_binned(
          limits = plot_lim_b,
          range = c(1, 12)
        )+theme_classic()+    theme(
          # Set the plot title and axis titles styles
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text( size = 25),
          legend.title = element_text(size = 20),
          
          # Set the axis text styles
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20))
      
      
      if(k==5){
        peak_df[paste(freq_space[i],peak_detection_algo[detect],"B",sep="_")] = conflict["freq_bin"]
      }
      
      ggsave(filename = paste("Results/Plots/Freq_Stability/Current/", trial_types[k],"_", freq_space [i], "_", peak_detection_algo[detect],"_B.png", sep=""),
                              plot = last_plot())
      
      ################  Option C)  ###########################
      # 1.	Peak-Frequenz für jeden Trial
      # 2.	Frequenz mit den meisten Trials mit höchstem Peak als ITF
      
      max_trial_count <-current_df %>% 
        drop_na(peak_amp) %>%  # get rid of trials where no peak was detected
        group_by(vp,session,trial_type,trial)%>% 
        mutate(max_freq=max(peak_amp)) %>% # get the max amplitude for each trial 
        subset(peak_amp == max_freq)  
      
      max_trial_count <-  max_trial_count%>% # only keep the freq for each trial where the amp was max
        group_by(vp,session,trial_type,freq_bin)%>% 
        count() %>%
        group_by(vp,session,trial_type)%>% 
        mutate(most_occ_freq=max(n), sum_amp = sum(freq_bin)) %>% 
        subset(n== most_occ_freq)
      
  
      legend_title_size = paste("Number of Trials")
      conflict =  subset(max_trial_count,max_trial_count$trial_type==k)
      {if (k==5) plot_lim_c = c(0,52) else plot_lim_c = c(0,26)}
    
      
      ggplot(data = conflict,aes(x = session, y = freq_bin, group=vp, color= vp)) +
        ylim(4,8)+
        geom_point(aes(size=most_occ_freq),alpha = 0.5)+
        geom_line()+
        labs(y= paste("Frequency [Hz]      ",bin_space, "- spcaced bins ", sep =" "))+
        #ggtitle(paste("C:Number of trials where freq. was peak-freq.       Trials:", trial_types[k])) + 
        guides(size=guide_legend(title=legend_title_size))+
        scale_size_binned(
          limits = plot_lim_c,
          range = c(1, 12)
        )+theme_classic()+   theme(
          # Set the plot title and axis titles styles
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text( size = 25),
          legend.title = element_text(size = 20),
          
          # Set the axis text styles
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20))
      
      
      
      
      if(k==5){
        peak_df[paste(freq_space[i],peak_detection_algo[detect],"C",sep="_")] = conflict["freq_bin"]
      }
      
      ggsave(filename = paste("Results/Plots/Freq_Stability/Current/", trial_types[k],"_", freq_space [i], "_", peak_detection_algo[detect],"_C.png", sep=""),
                              width = 7.47, height = 6.37, dpi = 150, units = "in",
                         plot = last_plot())
  
      
      ################  Option D)  ###########################
      # 1.	Mitteln der Amplitude der Peaks über die Trials für jedes Frequenzband
      # 2.	Frequenz mit dem höchsten Peak als ITF
      
      avg_sess_max <-current_df %>% 
        drop_na(peak_amp) %>%  # get rid of trials where no peak was detected
        group_by(vp,session,trial_type,freq_bin)%>% 
        mutate(avg_peak=mean(peak_amp)) # get the avg5 amplitude for each freq in each trial_type of each session 
      avg_sess_max<- avg_sess_max[c("vp","session","trial_type","freq_bin","avg_peak")]
      avg_sess_max  <-   avg_sess_max[!duplicated(avg_sess_max), ]
      avg_sess_max <- avg_sess_max%>%  
        group_by(vp,session,trial_type)%>% 
        mutate(max_peak=max(avg_peak))%>% 
        subset(avg_peak== max_peak)
      
      legend_title_size = paste("Average Amplitude\n", unit, sep=" ")
  
      conflict =  subset(avg_sess_max,avg_sess_max$trial_type==k)
      
     # %>% group_by(vp)
      ggplot(data = conflict,aes(x = session, y = freq_bin, group=vp, color= vp)) +
        ylim(4,8)+
        geom_point(aes(size=max_peak),alpha = 0.5)+
        labs(y= paste("Frequency [Hz]      ",bin_space, "- spcaced bins ", sep =" "))+
        geom_line()+
        #ggtitle(paste("D: Avg amp in peak freq across all trials.       Trials:", trial_types[k])) + 
        guides(size=guide_legend(title=legend_title_size))+
        scale_size_binned(
          limits = plot_lim_d,
          range = c(1, 12)
        )+theme_classic()+  theme(
          # Set the plot title and axis titles styles
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text( size = 25),
          legend.title = element_text(size = 20),
          
          # Set the axis text styles
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20))
      
          
      
      if(k==5){
        peak_df[paste(freq_space[i],peak_detection_algo[detect],"D",sep="_")] = conflict["freq_bin"]
      }
      
      ggsave(filename = paste("Results/Plots/Freq_Stability/Current/", trial_types[k],"_", freq_space [i], "_", peak_detection_algo[detect],"_D.png", sep=""),
        plot = last_plot())
    }
  }
}  
  
  

  for (i in 1:length(freq_space)){
    data <- read.csv(data_file,sep =",",fileEncoding="UTF-8-BOM", header = TRUE)
    data <- na.omit(data)
    data$freq_log <- as.factor(data$freq_log)
    
    
    
    if (grepl("log",freq_space[i])){
      freq_bin <- "freq_log"
      bin_space = "log"
    } else if (grepl("lin",freq_space[i])){
      freq_bin <- "freq_lin"
      bin_space = "linear"
    }
    peak_time <- paste("peak_time",freq_space[i],sep="_")
    
    current_df <- data[c("vp","session","trial_type","trial",freq_bin,peak_time)]
    colnames(current_df)<- c("vp","session","trial_type","trial","freq","time")
    
    for (k in 1:length(trial_types)){
      if (k==5){
        current_df$trial_type[current_df$trial_type==4 | current_df$trial_type==1] = 5
      }     
  
      for (j in 1:length(unique(data$vp))){
      sub_data <- subset(current_df,vp==j & trial_type==k)
      sub_data$freq <- as.factor(sub_data$freq)
      ggplot(data = sub_data ,aes(x = freq, y = time, fill= freq))+
        geom_violin()+
        ggtitle(paste("VP", unique(data$vp)[j], "Trials:", trial_types[k]))
        
        ggsave(filename = paste("Results/Plots/Freq_Stability/",peak_detection_algo[detect],"/Timing/VP_", unique(data$vp)[j],"_",trial_types[k],"_", freq_space [i],".png", sep=""),
               plot = last_plot())
      }
    }
  }





write.csv(peak_df, "peak_df.csv", row.names=FALSE)


## peak timings 


for (i in 1:length(freq_space)){
  data <- read.csv(data_file,sep =",",fileEncoding="UTF-8-BOM", header = TRUE)
  data <- na.omit(data)
  data$freq_log <- as.factor(data$freq_log)



  if (grepl("log",freq_space[i])){
    freq_bin <- "freq_log"
    bin_space = "log"
  } else if (grepl("lin",freq_space[i])){
    freq_bin <- "freq_lin"
    bin_space = "linear"
  }
  peak_time <- paste("peak_time",freq_space[i],sep="_")

  current_df <- data[c("vp","session","trial_type","trial",freq_bin,peak_time)]
  colnames(current_df)<- c("vp","session","trial_type","trial","freq","time")

  current_df <- current_df %>% group_by(vp,session,trial_type,trial) %>%
    mutate(std =sd(time))

  for (k in 1:length(trial_types)){
    if (k==5){
      current_df$trial_type[current_df$trial_type==4 | current_df$trial_type==1] = 5
    }

    for (j in 1:length(unique(data$vp))){
      sub_data <- subset(current_df,vp==j & trial_type==k)
      sub_data$freq <- as.factor(sub_data$freq)
      current_df$vp <- as.factor(current_df$vp)

      ggplot(data = current_df ,aes(x= vp, y = std, fill= vp))+
        geom_violin()+
        #geom_histogram()
        ggtitle(paste("Standard Deviation of peak timing for different frequencyies  Trials:", trial_types[k]))

      ggsave(filename = paste("Results/Plots/Freq_Stability/",peak_detection_algo[detect],"/Timing/VP_", unique(data$vp)[j],"_",trial_types[k],"_", freq_space [i],".png", sep=""),
             plot = last_plot())

    }
  }
}



# Intraclasscorrelation 


peak_df= read.csv("peak_df.csv",sep =",",fileEncoding="UTF-8-BOM", header = TRUE)


# for each pipeline
res_df = data.frame(pipeline = character() ,icc = numeric(), lower_ci = numeric(), upper_ci = numeric())
for(i in 3:length(colnames(peak_df))){
  #test auf ungleichheit: Deshalb 90% KI
  print( colnames(peak_df)[i])
  test <- reli_stats(measure=colnames(peak_df)[i], item="session",id="vp",data=peak_df,wide = FALSE, conf.level = .90, cv_calc = "MSE")
  res_df[i-2, ] <- c(colnames(peak_df)[i], test$icc$icc[2], test$icc$lower.ci[2], test$icc$upper.ci[2])
}

sum(is.na(res_df$icc))


res_df$icc = as.numeric(res_df$icc)  

ggplot(data=res_df, aes(x=icc))+
  geom_histogram(color = "white", fill = "grey30")


# ICC across the pipelines
peak_df_long <- gather(peak_df, pipeline, peak, 3:98, factor_key=TRUE)
peak_df_long <- mutate(peak_df_long, vp_sess= paste(vp, session))


test <- reli_stats(measure="peak", item="pipeline" ,id="vp_sess",data=peak_df_long,wide = FALSE, conf.level = .90, cv_calc = "MSE")




ind_cols <- grep("_ind_peak", names(peak_df), perl = TRUE)
ind_df = peak_df[ind_cols]
perc_ind_edge = 1- (sum(ind_df<9) - sum(ind_df>7.5) - sum(ind_df<4.5))/sum(ind_df<9)

broad_cols <- grep("_broad_peak", names(peak_df), perl = TRUE)
broad_df = peak_df[broad_cols]
perc_broad_edge = 1- (sum(broad_df<9) - sum(broad_df>7.5) - sum(broad_df<4.5))/sum(broad_df<9)



