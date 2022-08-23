import xspec
from preflowscp import *
from lmodel_preflow import *
xspec.AllModels.addPyMod(preflowscp, lmodel_preflowscp_info, 'add')
xspec.AllData.dummyrsp(0.01, 50., 1000, scaleType='log')
xspec.Model('preflowscp')
xspec.Plot.device='/xw'
xspec.Plot('emo')
