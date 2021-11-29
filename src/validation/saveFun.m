function saveFun(mcmcInfo, outPath, saveString)

    save([outPath 'mcmcInfo_' saveString '.mat'], 'mcmcInfo')