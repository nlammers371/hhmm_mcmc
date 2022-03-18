function initializePool(mcmcInfo)

    % determine max number of workers available
    myCluster = parcluster('local');
    maxWorkers = myCluster.NumWorkers;
    mcmcInfo.NumWorkers = min([mcmcInfo.NumWorkers maxWorkers]);
    
    % initialize pool if necessary
    pool = gcp('nocreate');
    if isempty(pool)
      parpool(mcmcInfo.NumWorkers);  
    elseif  pool.NumWorkers ~= mcmcInfo.NumWorkers     
      delete(pool)
      parpool(mcmcInfo.NumWorkers);  
    end  
   