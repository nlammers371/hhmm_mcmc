function initializePool(mcmcInfo)

    pool = gcp('nocreate');
    if isempty(pool)
      parpool(mcmcInfo.NumWorkers);  
    elseif  pool.NumWorkers ~= mcmcInfo.NumWorkers     
      delete(pool)
      parpool(mcmcInfo.NumWorkers);  
    end  
   