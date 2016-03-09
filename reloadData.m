function [yield, tf,succes]=reloadData(filename)

load(filename)
yield=evaluateData(time,ri*1E-6,XX,param,true);
tf=param.tf;