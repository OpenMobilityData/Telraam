
clear all
close all
clc

config = setupAnalysisConfig();
config.plots.showAdjustedCounts = false;
config.plots.showWeather = false;
config.plots.plotTypes = {'daily'};
runTelraamAnalysis();