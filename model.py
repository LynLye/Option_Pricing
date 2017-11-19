"""
Created on Sat Apr  8 10:18:27 2017
@author: Lyn Ye
"""

import numpy as np
import scipy.stats as stats








def PriceByBSFormula(T, K, S0, r, sigma, OptType, q = 0) :
    d1 = ((np.log(S0 / K) + (r - q) * T) / (sigma * np.sqrt(T))+ 0.5 * sigma * np.sqrt(T))
    d2 = d1 - (sigma * np.sqrt(T))
    if OptType == "CALL":
        return (S0 * np.exp(- q * T) * stats.norm.cdf(d1)- K * np.exp(-r * T) * stats.norm.cdf(d2))
    else:
        return (-S0 * np.exp(- q * T) *stats.norm.cdf(-d1) + K * np.exp(-r * T) * stats.norm.cdf(-d2))
    
def implied_volatility(T, K, S0, r, price, OptType, q):
    x0 = np.sqrt(2 * abs((np.log(S0 / K) + (r - q) * T) / T))
    error= 1
    i = 0            
    while ( error > 1e-10 or i < 700):
        d1 = ((np.log(S0 / K) + (r - q) * T) / (x0 * np.sqrt(T))+ 0.5 * x0 * np.sqrt(T))
        vega = 1 / np.sqrt(2 * np.pi) * S0 *np.exp(-q * T)* np.exp(-0.5 * d1 ** 2)
        error_ = ( PriceByBSFormula(T, K, S0, r, x0, OptType, q ) - price)/vega
        x0 -= error_
        error = abs(error_)
        i += 1
    return x0





   
    

def PriceBySnell(T, K, S0, r, sigma, OptType, N):
    val_num_steps = int(N)
    val_cache = {}

    dt = float(T) / val_num_steps
    up_fact= np.exp(sigma * (dt ** 0.5))
    down_fact = 1.0 / up_fact
    prob_up = (np.exp((r ) * dt) - down_fact) / (up_fact - down_fact)    
    def spot_price(num_ups, num_downs):
        return (S0) * (up_fact ** (num_ups - num_downs))
    def node_value(n, num_ups, num_downs):
        value_cache_key = (n, num_ups, num_downs)
        if value_cache_key not in val_cache:
            spot = spot_price(num_ups, num_downs)
            if  OptType == "CALL":
                exer_profit = max(0, spot - S0)
            else:  
                exer_profit = max(0, S0 - spot)

            if n >= val_num_steps:
                val = exer_profit
            else:
                fv = prob_up * node_value(n + 1, num_ups + 1, num_downs) + \
                    (1 - prob_up) * node_value(n + 1, num_ups, num_downs + 1)
                pv = np.exp(r * dt) * fv
                val = max(pv, exer_profit)
            val_cache[value_cache_key] = val
        return val_cache[value_cache_key]

    val = node_value(0, 0, 0)
    return val


def PriceByClosedFormulaForGmtrAsian(T, K, N,S0, r, sigma, OptType):
    a = S0 * np.exp(-r * T) * np.exp(((N + 1) * T) *(1/ (2 * N)) *(r + (sigma * sigma / 2) * ((2 * N + 1) / (3 * N) - 1)))
    b = sigma * np.sqrt((N + 1) * (2 * N + 1) / 6 / (N * N))
    return PriceByBSFormula(T, K, a, r, b, OptType)
    
    
    
def PriceByClosedFormulaForGmtBasket(T, K, S1, S2, r, sigma1, sigma2, rho, OptType):
    sigma_b =( np.sqrt((sigma1 * sigma1 +sigma2 * sigma2 +2 *rho * sigma1 * sigma2)) / 2)
    niu_b =( r - (sigma1 * sigma1 + sigma2 * sigma2) / 4 + sigma_b * sigma_b / 2)
    S0_b = np.sqrt(S1 * S2) * np.exp(-r * T) * np.exp(niu_b *T)
    return PriceByBSFormula(T, K, S0_b, r, sigma_b, OptType)


def PriceByMCSimulationForArthmAsian(T, K, N,S0, r, sigma, OptType, NumPath, isGmtrCon = 0):
    dT = T/N
    sqrtdT = np.sqrt(dT)
    drift = np.exp((r - 0.5 * sigma * sigma) * dT)
    SPath = np.zeros((NumPath,N), dtype = float)
    df = np.exp( -r * T)
    np.random.seed()
    for i in range(NumPath):
        deltaZ = np.random.standard_normal()
        growthFactor = drift * np.exp(sigma * sqrtdT * deltaZ)
        SPath[i][0] = S0 * growthFactor
        for j in range(1,N):
            deltaZ = np.random.standard_normal()
            growthFactor = drift * np.exp(sigma * sqrtdT * deltaZ)
            SPath[i][j] = SPath[i][j-1] * growthFactor
    SPathArithMean = SPath.mean(1)
    if OptType == "CALL":
        OptSample = np.maximum(SPathArithMean - K, 0) * df
    else:
        OptSample = np.maximum(K - SPathArithMean , 0) * df
    if isGmtrCon == 0:        
        meanOptSample = OptSample.mean()
        stdOptSample = OptSample.std()
        OptVal =  meanOptSample
        LCI = meanOptSample - 1.96 * stdOptSample/np.sqrt(NumPath)
        RCI =meanOptSample + 1.96 * stdOptSample/np.sqrt(NumPath)
    else: 
        GmtrConVal = PriceByClosedFormulaForGmtrAsian(T, K, N,S0, r, sigma, OptType)   
        SPathGmtrMean = SPath.prod(1) ** (1/N)        
        if OptType == "CALL":
            OptSampleGmtr= np.maximum(SPathGmtrMean - K , 0)*df
             
        else:
            OptSampleGmtr = np.maximum(K - SPathGmtrMean, 0)*df
        theta = ((OptSample*OptSampleGmtr).mean() -OptSample.mean()*OptSampleGmtr.mean())/(OptSampleGmtr.var())
        OptSampleErrRed = OptSample + (GmtrConVal-OptSampleGmtr) * theta
        mean_OptSampleErrRed = OptSampleErrRed.mean()
        std_OptSampleErrRed = OptSampleErrRed.std()
        
        OptVal =  mean_OptSampleErrRed
        LCI = mean_OptSampleErrRed - 1.96 * std_OptSampleErrRed/np.sqrt(NumPath)
        RCI = mean_OptSampleErrRed + 1.96 * std_OptSampleErrRed/np.sqrt(NumPath)

    return OptVal, [LCI, RCI]




def PriceByMCSimulationForArthmBasket(T, K, S1, S2, r, sigma1, sigma2, rho, OptType,NumPath, isGmtrCon=0):
    df= np.exp( -r *T)
    sqrtT = np.sqrt(T)
    drift1 = np.exp((r - 0.5 * sigma1 * sigma1) * T)
    drift2 = np.exp((r - 0.5 * sigma2 * sigma2) * T)
    np.random.seed()
    sample_1 = np.random.normal(size=(NumPath,1))#(3,1)
    sample_ = np.random.normal(size=(NumPath,1))
    sample_2 = sample_1 * rho + np.sqrt(1 - rho * rho) * sample_
    sample_1_ = np.exp(sample_1 * sqrtT * sigma1) * drift1 * S1
    sample_2_ = np.exp(sample_2 * sqrtT * sigma2) * drift2 * S2
    sample = (sample_1_ +sample_2_)/2 
    if OptType == "CALL":        
        sample = np.maximum( sample - K, 0)*df   
    else:
        sample = np.maximum( K - sample, 0)*df     
    if isGmtrCon == 0: 
        mean_sample = sample.mean()
        std_sample = sample.std()
        OptVal =  mean_sample
        LCI = mean_sample - 1.96 * std_sample/np.sqrt(NumPath)
        RCI = mean_sample + 1.96 * std_sample/np.sqrt(NumPath)
    else: 
        GmtrConVal = PriceByClosedFormulaForGmtBasket(T, K, S1, S2, r, sigma1, sigma2, rho, OptType)
        sampleGmtr_ = np.sqrt(sample_1_ * sample_2_ )
        if OptType == "CALL": 
            sampleGmtr = np.maximum( sampleGmtr_ - K, 0)*df
        else:
            sampleGmtr = np.maximum( K - sampleGmtr_ , 0)*df
        theta = ((sample * sampleGmtr).mean() -sample.mean() * sampleGmtr.mean())/(sampleGmtr.var())
        sampleErrRed = sample + (GmtrConVal - sampleGmtr)*theta 
        mean_sampleErrRed = sampleErrRed.mean()
        std_sampleErrRed = sampleErrRed.std()
        
        
        OptVal =  mean_sampleErrRed
        LCI = mean_sampleErrRed - 1.96 * std_sampleErrRed/np.sqrt(NumPath)
        RCI = mean_sampleErrRed + 1.96 * std_sampleErrRed/np.sqrt(NumPath)
    return OptVal, [LCI, RCI]

if __name__=="__main__":
    r = 0.05
    T = 3
    S0 = 100 
    print("============================European=====================")
    print(PriceByBSFormula(T, 100, S0, r, 0.3, 'PUT',0) )
    print(PriceByBSFormula(T, 100, S0, r, 0.3, 'CALL',0) )
    print("============================AMERICAN=====================")
    print(PriceBySnell(T, 100, S0, r, 0.3, 'PUT', 100) )
    print(PriceBySnell(T, 100, S0, r, 0.3, 'CALL', 100) )
    print("============================GmtrAsian=====================")
    
    print (PriceByClosedFormulaForGmtrAsian(T, 100 , 50, S0, r, 0.3, 'PUT'))
    print (PriceByClosedFormulaForGmtrAsian(T, 100 , 100, S0, r, 0.3, 'PUT'))
    print (PriceByClosedFormulaForGmtrAsian(T, 100 , 50, S0, r, 0.4, 'PUT'))
    print("Kay")
    print (PriceByClosedFormulaForGmtrAsian(T, 100 , 50, S0, r, 0.3, 'CALL'))
    print (PriceByClosedFormulaForGmtrAsian(T, 100 , 100, S0, r, 0.3, 'CALL'))
    print (PriceByClosedFormulaForGmtrAsian(T, 100 , 50, S0, r, 0.4, 'CALL'))
    print("===========================ArthmAsian NoCont=====================")
    print("Kay")
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.3, 'PUT', 100000))
    print(PriceByMCSimulationForArthmAsian(T,100,100,S0, r, 0.3, 'PUT', 30))
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.4, 'PUT', 30))
    print("Kay")
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.3, 'CALL', 100000))
    print(PriceByMCSimulationForArthmAsian(T,100,100,S0, r, 0.3, 'CALL', 30))
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.4, 'CALL', 30))
    print("===========================ArthmAsian Cont=====================")
    print("Kay")
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.3, 'PUT', 30,1))
    print(PriceByMCSimulationForArthmAsian(T,100,100,S0, r, 0.3, 'PUT', 30,1))
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.4, 'PUT', 30,1))
    print("Kay")
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.3, 'CALL', 30,1))
    print(PriceByMCSimulationForArthmAsian(T,100,100,S0, r, 0.3, 'CALL', 30,1))
    print(PriceByMCSimulationForArthmAsian(T,100,50,S0, r, 0.4, 'CALL', 30,1))
    print("============================GmtrBasket=====================")
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.5, 'PUT'))
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.9, 'PUT'))
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.1, 0.3, 0.5, 'PUT'))
    print(PriceByClosedFormulaForGmtBasket(T, 80, 100, 100, r, 0.3, 0.3, 0.5, 'PUT'))
    print(PriceByClosedFormulaForGmtBasket(T, 120, 100, 100, r, 0.3, 0.3, 0.5, 'PUT'))
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.5, 0.5, 0.5, 'PUT'))
    print("Kay")
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.5, 'CALL'))
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.9, 'CALL'))
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.1, 0.3, 0.5, 'CALL'))
    print(PriceByClosedFormulaForGmtBasket(T, 80, 100, 100, r, 0.3, 0.3, 0.5, 'CALL'))
    print(PriceByClosedFormulaForGmtBasket(T, 120, 100, 100, r, 0.3, 0.3, 0.5, 'CALL'))
    print(PriceByClosedFormulaForGmtBasket(T, 100, 100, 100, r, 0.5, 0.5, 0.5, 'CALLl'))
    print("============================ArthmBasket NOCont=====================")
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.5, 'PUT',100000))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.9, 'PUT',100000))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.1, 0.3, 0.5, 'PUT',100000))
    print(PriceByMCSimulationForArthmBasket(T, 80, 100, 100, r, 0.3, 0.3, 0.5, 'PUT',100000))
    print(PriceByMCSimulationForArthmBasket(T, 120, 100, 100, r, 0.3, 0.3, 0.5, 'PUT',100000))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.5, 0.5, 0.5, 'PUT',100000))
    print("Kay")
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.5, 'CALL',100000))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.9, 'CALL',100000))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.1, 0.3, 0.5, 'CALL',100))
    print(PriceByMCSimulationForArthmBasket(T, 80, 100, 100, r, 0.3, 0.3, 0.5, 'CALL',100))
    print(PriceByMCSimulationForArthmBasket(T, 120, 100, 100, r, 0.3, 0.3, 0.5, 'CALL',100))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.5, 0.5, 0.5, 'CALLl',100))    
    print("============================ArthmBasket  Cont=====================")
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.5, 'PUT',100000,1))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.9, 'PUT',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.1, 0.3, 0.5, 'PUT',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 80, 100, 100, r, 0.3, 0.3, 0.5, 'PUT',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 120, 100, 100, r, 0.3, 0.3, 0.5, 'PUT',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.5, 0.5, 0.5, 'PUT',100,1))
    print("Kay")
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.5, 'CALL',100000,1))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.3, 0.3, 0.9, 'CALL',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.1, 0.3, 0.5, 'CALL',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 80, 100, 100, r, 0.3, 0.3, 0.5, 'CALL',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 120, 100, 100, r, 0.3, 0.3, 0.5, 'CALL',100,1))
    print(PriceByMCSimulationForArthmBasket(T, 100, 100, 100, r, 0.5, 0.5, 0.5, 'CALLl',100,1))

    
    
       
   
