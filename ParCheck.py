# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:35:33 2017

@author: Lyn Ye
"""



from model import *

def implied_volatility_inter(T, K, S0, r, price, OptType, q):
    if (r < 0) or (q < 0) or (S0 <= 0) or (T <= 0) or (price <= 0):
        return "Invalid input !"
    return implied_volatility(T, K, S0, r, price, OptType, q)

def PriceBySnell_inter(T, K, S0, r, sigma, OptType, N):
    if (r < 0) or (S0 <= 0) or (K <= 0) or ( T <= 0) or (sigma <= 0) or (N <= 0):
        return "Invalid input !"
    if N >= 2500:
        return "The maximum number of steps is 2500 !"
    return PriceBySnell(T, K, S0, r, sigma, OptType, N)

def  PriceByBSFormula_inter(T, K, S0, r, sigma, OptType, q) :
    if (r < 0) or (S0 <= 0) or (K <= 0) or (T <= 0) or (sigma <= 0) or (q < 0):
        return "Invalid input !"
    return  PriceByBSFormula(T, K, S0, r, sigma, OptType, q) 

def PriceByClosedFormulaForGmtrAsian_inter(T, K, N,S0, r, sigma, OptType):
    if (r < 0) or (S0 <= 0) or (K <= 0) or (T <= 0) or (sigma <= 0) or (N <= 0):
        return "Invalid input !"
    return PriceByClosedFormulaForGmtrAsian(T, K, N,S0, r, sigma, OptType)

def PriceByClosedFormulaForGmtBasket_inter(T, K, S1, S2, r, sigma1, sigma2, rho, OptType):
    if (r < 0) or (S1 <= 0) or (S2 <= 0) or (K <= 0) or (T <= 0) or (sigma1 <= 0) or (sigma2 <= 0) or (rho < -1) or (rho > 1):
        return "Invalid input !"
    return PriceByClosedFormulaForGmtBasket(T, K, S1, S2, r, sigma1, sigma2, rho, OptType)

def PriceByMCSimulationForArthmAsian_inter(T, K, N,S0, r, sigma, OptType, NumPath, isGmtrCon):
    if (r < 0) or (S0 <= 0) or (K <= 0) or (T <= 0) or (sigma <= 0) or (N <= 0) or (NumPath <= 0):
        return "Invalid input !"
    return PriceByMCSimulationForArthmAsian(T, K, N,S0, r, sigma, OptType, NumPath, isGmtrCon = 0)

def  PriceByMCSimulationForArthmBasket_inter(T, K, S1, S2, r, sigma1, sigma2, rho, OptType,NumPath, isGmtrCon=0):
    if (r < 0) or (S1 <= 0) or (S2 <= 0) or (K <= 0) or (T <= 0) or (sigma1 <= 0) or (sigma2 <= 0) or (rho < -1) or (rho > 1) or (NumPath <= 0):
        return "Invalid input !"
    return PriceByMCSimulationForArthmBasket(T, K, S1, S2, r, sigma1, sigma2, rho, OptType,NumPath, isGmtrCon)

