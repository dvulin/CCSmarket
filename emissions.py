# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:58:12 2022

@author: dvulin
"""

import numpy as np
class Emissions(object):
    """
    some variables shared at class level - for all instances (objects)
    """
    allowance_EUA_price = 80.0              # fixed for each time_step *** DA LI STAVITI DA NEKI KUPUJU CO2 U STARTU ISPOD OVE CIJENE?
    domestic_CO2_price = 80.0               # domestic CO2 price should decrease as sellers appear
    trading_volume = []                     # tCO2 = trading volume
    time_steps = None
    trading_history = None                  # all history of trades ([timestep, volume, price])
    token_number = 0
    
    def __init__(self, ts = np.arange(0,30), cbr = 0.1, ip = 0, mb = 0):
        """
        
        """
        self.core_business_r = cbr          # core business interest rate
        self.investment_potential = ip      # money that can be available for new investments
        self._released_CO2 = None
        self._reduced_CO2 = None
        self._core_cash_flow = None
        self._money_balance = None    # money balance (used for more specific profitability analysis)
        self._allowance_wallet = 0    # array, CO2 for trading (1t CO2 in wallet = 1 token)
              
        self.setTimeSteps(ts)
      
    @classmethod
    def setTimeSteps(cls, values):
        if type(values) == int:
            cls.time_steps = np.arange(values)
        else:
            cls.time_steps = values
    
    @classmethod
    def setAllowanceEUAPrice(cls, values):
        cls.allowance_EUA_price = values
    
    @classmethod
    def setDomesticCO2Price(cls, value):
        cls.domestic_CO2_price = value
    
    @classmethod
    def setTradingVolume(cls, value):
        cls.trading_volume = value

    @classmethod
    def updateTradingHistory(cls, values):
        cls.trading_history.append(values)
        
    @classmethod
    def updateTokenNumber(cls, values):
        cls.token_number += values
        
    """
    @property
    def time_steps(self):
        return self._time_steps
    
    @time_steps.setter
    def time_steps(self, value):
        if type(value) == int:
            self._time_steps = np.arange(value)
        else:
            self._time_steps = value
    """
    
    @property
    def released_CO2(self):
        return self._released_CO2
    
    @released_CO2.setter
    def released_CO2(self, values):
        self.updateTokenNumber(values)          # tokens are minted when more CO2 is emitted
        self._released_CO2 = values
        self.setDomesticCO2Price(self.domestic_CO2_price - 
                                 self.domestic_CO2_price*(values/self.token_number.cumsum()))
        
    @property
    def reduced_CO2(self):
        return self._reduced_CO2
    
    @reduced_CO2.setter
    def reduced_CO2(self, values):
        self._allowance_wallet += values
        self._reduced_CO2 = values
        self.setDomesticCO2Price(self.domestic_CO2_price + 
                                self.domestic_CO2_price*(values/self.token_number.cumsum()))
           
    @property
    def free_allowances(self):
      return self._free_alowances    # allowances are available for trading
    
    @free_allowances.setter
    def free_allowances(self, tonnes_CO2):
        self.updateTokenNumber(-tonnes_CO2)         # burns tokens and rises token price
        self._free_alowances = tonnes_CO2
        
    @property
    def allowance_wallet(self):
        return self._allowance_wallet
    
    @allowance_wallet.setter
    def allowance_wallet(self, value):
        self._allowance_wallet = value
    """
    @property
    def allowance_EUA_price(self, price):
        return self._allowance_EUA_price
    
    @allowance_EUA_price.setter
    def allowance_EUA_price(self, price):
        self._allowance_EUA_price = price

    @property
    def domestic_CO2_price(self):
        return self._domestic_CO2_price
    
    @domestic_CO2_price.setter
    def domestic_CO2_Price(self, price):
        self._domestic_CO2_price = price
    """
        
    @property
    def money_balance(self):
        return self._money_balance
    
    @money_balance.setter
    def money_balance(self, value):
        self._money_balance = value 
 
    def sellAllowances(self, time_step, tonnes_CO2, unit_price):
        """
        Parameters
        ----------
        time_step : int
            time step when CO2 is sold
        tonnes_CO2 : float
            ammount of CO2
        unit_price : float
            price of 1t of CO2 at market
        Returns
        -------
        (1) checks if allowances are available for trading
        (2) adds money to the balance at given time_step
        """
        if self._allowance_wallet.cumsum()[time_step] < tonnes_CO2:
            return 'not enough allowances in wallet. transaction unsuccessful'
        else:
            self._allowance_wallet[time_step] -= tonnes_CO2
            self._money_balance[time_step] += tonnes_CO2*unit_price
            ### TODO: prilagoditi cijenu domestic_CO2_price ovisno o prodanoj kolicini
            self.setDomesticCO2Price()
            self.updateTradingHistory([time_step, tonnes_CO2, unit_price])
            
        return ("""sold {tCO2} allowances for {up} per 1t of CO2.
                    balance increased for {total} at timestep {ts}.
                """.format(tCO2 = tonnes_CO2, up = unit_price, total = tonnes_CO2*unit_price, ts = time_step))
        
    def buyAllowances(self, time_step, tonnes_CO2, unit_price):
        if unit_price*tonnes_CO2 > self._money_balance[time_step]:
            return 'not enough money. transaction unsuccessful'
        else:
            self._free_allowances[time_step] += tonnes_CO2
            self._money_balance[time_step] -= tonnes_CO2*unit_price
        return ("""bought {tCO2} allowances for {up} per 1t of CO2.
                    balance reduced for {total} at timestep {ts}.
                """.format(tCO2 = tonnes_CO2, up = unit_price, total = tonnes_CO2*unit_price, ts = time_step))    
    
    @property
    def core_cash_flow(self):
        return (self._core_cash_flow)
    
    @core_cash_flow.setter
    def core_cash_flow(self, value):
        self._core_cash_flow = value
              
    def calculateCoreCashFlow(self, CAPEX, OPEX, income, r):
        """
        Parameters
        ----------
        CAPEX : float
            capital expenditure, can be an array (if reinvestments are considered).
        OPEX : float
            DESCRIPTION.
        income : float
            DESCRIPTION.
        r : TYPE
            DESCRIPTION.
        Returns
        -------
        cash flow time-series for production of core production
        """
        import numpy_financial as npf
        if type(CAPEX) == int:
            CAPEX = np.ones(len(self.time_steps))      # if it is not array, CAPEX appears only at first time_step
            CAPEX[0] = CAPEX                           # meaning there are no reinvestments
        if type(OPEX) == int:
            OPEX = np.ones(len(self.time_steps))*OPEX  # if it is not array, OPEX is constant through the time
        if type(income) == int:
            income = np.ones(len(self.time_steps))*income
        cf = income - (CAPEX + OPEX)
        self._core_cash_flow = npf.pv(r, np.arange(len(self.time_steps)), 0, -cf )
        
    def geometricBrownianMotion(self, time_horizon= 30, n_time_steps = 30, initial_value = 0, sigma = 0.8, mu = 0):
        dt = time_horizon/n_time_steps
        vty = np.exp((mu - sigma ** 2 / 2) * dt 
                     + sigma * np.random.normal(0, np.sqrt(dt), 
                    size=n_time_steps))
        vty[0] = 1
        GBM = initial_value*vty.cumprod()
        return (GBM)
    
    def fitGeometricBrownianMotion(self, values, dt=1/12, log_returns = False):
        """
        Parameters
        ----------
        values : np.array(float)
            historical prices
        dt : float, optional
            time period. The default is 1/12 = 1 month.
        log_returns : True/False, optional
            if log-returns are used. The default is False.
        
        Returns
        -------
        mu : float
            GBM mean = drift
        sigma : float
            standard deviation
        """           
        if log_returns:
            ri = np.log(values[1:]/values[:-1])         # log returns
        else:
            ri = (values[1:]-values[:-1])/values[:-1]    # returns
        m = np.mean(ri)                                 # mean
        v = np.sum((ri-m)**2)/(len(ri)-1)*dt            # variance
        sigma = v**0.5                                  # standard deviation
        mu = m+v/2       
        return (mu, sigma)
        
    def trendModel(self, time_series, initial, change, model = 'linear'):
        """
        Parameters
        ----------
        time_series : array
            DESCRIPTION.
        initial : TYPE
            initial value.
        percent_change : float
            percent change is defined as part of 1, eg. 2% reduction each timestep = -0.02

        Returns
        -------
        numpy array
            trend values.
        """
        p = [] 
        if model == 'percent':
            for i, t in enumerate(time_series):
                p.append(initial*(1+change)**i)
            return (np.array(p))
        if model == 'linear':
            for i, t in enumerate(time_series):
                p.append(initial+change*i)
            return (np.array(p))       
            
class Industry(Emissions):
    def __init__(self):
        super().__init__()
      
class BECCS(Emissions):
    def __init__(self):
        pass

class CarbonFarm(Emissions):
    def __init__(self):
        pass
    
class StorageOperator(Emissions):
    def __init__(self):
        pass
    
class Transport(Emissions):
    def __init__(self):
        pass

# initialize parent object - it will store general class-level variables:
    # allowance_EUA_price = 85.0              
    # domestic_CO2_price = 85.0
    # trading_volume = 0
    # time_steps = None
    # trading_history = []
    
e = Emissions()
e.setAllowanceEUAPrice(e.trendModel(time_series = e.time_steps, 
                              initial = 85., 
                              change = 0.05, 
                              model = 'percent'))


"----------------------------------------------------------------------------------------------------------"
# (1) initialize stakeholder
ina = Industry()

# (2) initialize timeseries of:
#   - released_CO2 (emissions)
#   - reduced_CO2 (stored, sequestered)
#   - free_allowances (to calculate penalties charged by EUA_price)

ina.released_CO2 = e.trendModel(time_series = ina.time_steps, 
                              initial = 1.2e6, 
                              change = -0.1, 
                              model = 'percent')
ina.reduced_CO2 = e.trendModel(time_series = ina.time_steps, 
                              initial = 0.15e6, 
                              change = 0.1, 
                              model = 'percent')
ina.free_allowances = e.trendModel(time_series = ina.time_steps, 
                              initial = 300000, 
                              change = -0.1, 
                              model = 'percent')


# iz https://www.poslovni.hr/kompanije/neto-dobit-ine-u-prvom-kvartalu-znatno-veca-nego-lani-4334651
ina_neto_prihod = 4*6.24E9/7.54                 #eur, u 2022., na tamelju kvartala
ina.CAPEX = 4*846e6/7.54           
ina_neto_dobit = 4*586e6/7.54                   # na temelju kvartalne dobiti
ina.OPEX = (ina_neto_prihod-ina_neto_dobit-ina.CAPEX)
ina.core_business_r = 0.09
ina.calculateCoreCashFlow(ina.CAPEX, ina.OPEX, ina_neto_prihod, ina.core_business_r)

# print (ina.allowance_wallet)        # wallet is initialized for all steps when CO2_reduced is set
# run simulation
for i, ts in enumerate(ina.time_steps):
    if i>0: ina.allowance_wallet[i] += ina.released_CO2[i-1]-ina.released_CO2[i]   #reduction is already included, this are the emissions
# print (ina.allowance_wallet)

"----------------------------------------------------------------------------------------------------------"
NEXE = Industry()
NEXE.released_CO2 = e.trendModel(time_series = NEXE.time_steps, 
                              initial = 0.7e6, 
                              change =  0.005, 
                              model = 'percent')
NEXE.reduced_CO2 = e.trendModel(time_series = NEXE.time_steps, 
                              initial = 0.15e6, 
                              change = 0.0, 
                              model = 'percent')
NEXE.free_allowances = e.trendModel(time_series = NEXE.time_steps, 
                              initial = 250000, 
                              change = -0.1, 
                              model = 'percent')

# iz https://www.poslovni.hr/kompanije/neto-dobit-ine-u-prvom-kvartalu-znatno-veca-nego-lani-4334651
NEXE_neto_prihod = 1*6.24E9/7.54                 #eur, u 2022., na tamelju kvartala
NEXE.CAPEX = 0.1*846e6/7.54           
NEXE_neto_dobit = 0.1*586e6/7.54                   # na temelju kvartalne dobiti
NEXE.OPEX = (NEXE_neto_prihod-NEXE_neto_dobit-ina.CAPEX)
NEXE.core_business_r = 0.08
NEXE.calculateCoreCashFlow(NEXE.CAPEX, NEXE.OPEX, NEXE_neto_prihod, NEXE.core_business_r)

# print (NEXE.allowance_wallet)        # wallet is initialized for all steps when CO2_reduced is set
# run simulation
for i, ts in enumerate(NEXE.time_steps):
    if i>0: NEXE.allowance_wallet[i] += NEXE.released_CO2[i-1]-NEXE.released_CO2[i]   #reduction is already included, this are the emissions
# print (NEXE.allowance_wallet)


"""
example of GBM fitting
- number of values should be checked in historical data (h_CO2_price)) in order to correctly define dt
- drift might vary a lot, depending on time-window used ("1d". "1mo")

for CO2, we took daily closing price and got the results:
mu, sigma =  (0.001117098100358904, 0.00227731379179647)

optimal time window for similar analyses has been prposed in:
    Vulin, D., Arnaut, M., & Karasalihović Sedlar, D. (2020). 
    Forecast of long‐term EUA price probability using momentum strategy and GBM simulation. 
    Greenhouse Gases: Science and Technology, 10(1), 230-248.
"""
# import yfinance as yf
# ticker = yf.Ticker("CO2.L")
# h_CO2_price = ticker.history(start = "2018-12-31", end = "2022-10-13", interval = '1d')['Close']
# mu, sigma = e.fitGeometricBrownianMotion(dt = 1/248, values = h_CO2_price.to_numpy(), log_returns = False)
# paths = []
# for i in range(50000):
#     paths.append(e.geometricBrownianMotion(time_horizon = 30, 
#                                            n_time_steps = 30, 
#                                            initial_value = 70., 
#                                            sigma = sigma, mu = mu))
# paths = np.array(paths).T
# ppath = np.percentile(paths, q = [5, 25, 50, 75, 95], axis=1)
# import matplotlib.pyplot as plt
# plt.plot(ppath.T)
# plt.legend(['p5', 'p25', 'p50', 'p75', 'p95'])

import pandas as pd
input = pd.read_excel('input.xlsx', sheet_name=None)
stakeholders = []
gi_columns = ['variable', 'value', 'comment']
ts_columns = ['year','released_CO2','reduced_CO2','free_allowances','core_cash_flow']
for name, sheet in input.items():
    stakeholders.append({'sheet_name':name, 
                        'gi' : input[name][gi_columns].dropna(axis='rows').set_index('variable'),
                        'ts' : input[name][ts_columns].dropna(axis='columns')
                    })

# after all inputs are loaded they can be accessed as list elements:
stakeholders[0]['gi'].loc['r']['value']

"""
gi is general information table, and r value is called
ts is time series value
"""

