 clear all
 clc

 Subject = {'Bryan', 'Thomas', 'Hamid', 'Kayla', 'Sevack', 'Jim', 'Tarun', 'Anastasia', 'Shivam', 'Brian', 'Angela', 'Darya', 'Dave', 'Dorisha', 'Gary', 'Greg', 'Larry R', 'Larry D', 'Marcus', 'Alycia'};


 for sbjct = 1:20     
 
     MIK_S_load = load(['MIK_S', char(Subject(sbjct))]);
     M_i_k_S = MIK_S_load.M_i_k_S;

     MIK_E_load = load(['MIK_E', char(Subject(sbjct))]);  
     M_i_k_E = MIK_E_load.M_i_k_E;
          
     for i = 1:length(M_i_k_S)
         k = 1;
         for j = [1 2 4 6 8]
            M_i_k_S_mean{1,sbjct}(i,k) = mean(M_i_k_S{i,j});
            M_i_k_S_peak{1,sbjct}(i,k) = max(M_i_k_S{i,j});
            k = k + 1;
         end
     end
     
     MikS_mean(sbjct,:) = mean(M_i_k_S_mean{1,sbjct});
     MikS_peak(sbjct,:) = max(M_i_k_S_peak{1,sbjct});
     
     for i = 1:length(M_i_k_E)
         k = 1;
         for j = [1 2 4 6 8]
            M_i_k_E_mean{1,sbjct}(i,k) = mean(M_i_k_E{i,j});
            M_i_k_E_peak{1,sbjct}(i,k) = max(M_i_k_E{i,j});
            k = k + 1;
         end
     end
     
     MikE_mean(sbjct,:) = mean(M_i_k_E_mean{1,sbjct});
     MikE_peak(sbjct,:) = max(M_i_k_E_peak{1,sbjct});
     
 end     
 