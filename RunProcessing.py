from email.mime import base
from socket import gaierror
from traitlets import Instance
from ROOT import TF1, TH1F, gPad, TCanvas, TFile, TGraph
from tracemalloc import start
import numpy as np
import array as arr
import codecs
import math
import scipy 
import json
import os 
import natsort



configfile = 'config.json'
with open(configfile) as datafile:
    config = json.load(datafile)

def getInfo(files): #Scoate informatiile de care ai nevoie din fisierel e tip info pe care le scoate Janus
    infoFiles = [os.path.relpath(x) for x in os.scandir(files)]
    infoFiles = natsort.natsorted(infoFiles)
    holdDelay = np.array([])
    for i in range(len(infoFiles)):
        with codecs.open(infoFiles[i], 'r', 'utf-8', errors='ignore') as file:
            for i, line in enumerate(file):
                line = line.strip()
                info = line.split()
                if(len(info) == 0):
                 i = i+1
                else:
                    if info[0] == 'HoldDelay':
                        holdDelay = np.append(holdDelay, info[1])
    return holdDelay

def getPeak(hist, start, end): #Scoate niste parametrii preliminari pentru a fita gausiene
    mean = 0
    height = 0
    fhwm = 0
    sigma = 0
    baseline = np.array([])
    par = np.array([]), 
    print ('start = ', start, "end = ", end)
    for i in range(start, start + 10):
        baseline = np.append(baseline, hist.GetBinContent(i))
    bslValue = np.sum(baseline) / len(baseline)
    offset = 50
    hist.GetXaxis().SetRange(start, end)
    treshold = float(hist.GetMaximum())*0.8
    print('treshold= ', treshold)
    for i in range(start, end):
        dif1 = offset + hist.GetBinContent(i) - bslValue
        dif2 = offset + hist.GetBinContent(i+1) - bslValue +  0.05*(offset + hist.GetBinContent(i+1) - bslValue)
        print('dif1 = ', dif1, "dif2 = ", dif2)
        if dif2 < dif1 and dif1 > treshold and dif2 > treshold:
            height = hist.GetBinContent(i)
            mean = hist.GetBinCenter(i)
            par = np.append(par, height)
            par = np.append(par, mean)
            break
    print("mean = ", mean, 'height = ', height)
    for i in range(hist.FindBin(mean), end):
        h = hist.GetBinContent(i)
        c = hist.GetBinCenter(i)
        uplim = height/2 + 0.15 * (height/2)
        downlim = height/2 - 0.15 * (height/2)
        if downlim < h and h < uplim:
            fhwm = (c - mean) * 2
            sigma = fhwm/2.355 
            par = np.append(par, sigma)
            par = np.append(par, fhwm) 
            break
    print('fhwm = ', fhwm, 'sigma = ', sigma)
    return par

def fit(hist, par, start, end): #fiteaza gausiana cu un vector de parametrii dati
    fit = TF1("fitfunc", "gaus(0) + pol1(3)", start, end)
    fit.SetParameters(par[0], par[1], par[2])
    fit.SetParLimits(0, par[0] - par[0]*0.3, par[0] + par[0]*0.3)
    fit.SetParLimits(1, par[1] - par[1]*0.3, par[1] + par[1]*0.3)
    fit.SetParLimits(2, par[2] - par[2]*0.5, par[2] + par[2]*0.5)
    hist.Fit(fit, 'RNQ')
    return fit

def integrate(hist, start, end): #O mica functie de integrare pentru binii histogramei
    integration = 0
    for i in range(hist.FindBin(start), hist.FindBin(end)):
        integration += hist.GetBinContent(i)
    return integration

def gaussianArea(func): #Aria gausienei
    mu = func.GetParameter(1)
    sigma = func.GetParameter(2)
    return mu * sigma / 0.3989

def plot(files): #Ploeaza histogramele generate din fisierele de tip list pe care le scoate Janus
    histsCh0 = []
    histsCh32 =[]
    dataFiles = [os.path.relpath(x) for x in os.scandir(files)]
    dataFiles = natsort.natsorted(dataFiles)
    Ch0 = TCanvas('c1', "Channel 0")
    Ch32 = TCanvas('c2', 'Channel 32')
    for  i in range(len(dataFiles)):
        print(i)
        with codecs.open(dataFiles[i], 'r', 'utf-8', errors='ignore') as file:
            histCh0 = TH1F(f"h0_{i}", f"Histogram Channel 0_{i}", 2**13, 0, 2**13-1)
            histCh32 = TH1F(f"h32_{i}", f"Histogram Channel 32_{i}", 2**13, 0, 2**13-1)
            for i, line in enumerate(file):
                if (i>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        histCh32.Fill(int(event[-2]))
            histsCh0.append(histCh0)
            histsCh32.append(histCh32)
            # print(f'Max hist {i} = ', histCh0.GetMaximum(), f'Max hist {i} = ', histCh32.GetMaximum())
    print(len(histCh32), " ", len(histCh0))
    

    file = TFile('hists.root', 'RECREATE')

    for hist in histsCh0:
        if hist.Integral():
            hist.Write()
    for hist in histsCh32:
        if hist.Integral():
            hist.Write()

    for i in range(len(histCh0)):   
        Ch0.cd()
        histsCh0[i].Draw('hist')
        Ch32.cd()
        histsCh32[i].SetLineColor(2)
        histsCh32[i].Draw('hist')
        gPad.WaitPrimitive('ggg')

def analyseFiles(files): #Scoate parametrii doriti din fisierele de tip list de la Janus
    parCh0 = np.array([])
    parCh32 = np.array([])
    histsCh0 = []
    histsCh32 =[]
    dataFiles = [os.path.relpath(x) for x in os.scandir(files)]
    dataFiles = natsort.natsorted(dataFiles)
    for  i in range(len(dataFiles)):
        print(i)
        with codecs.open(dataFiles[i], 'r', 'utf-8', errors='ignore') as file:
            pCh0 = np.array([])
            pCh32 = np.array([])
            histCh0 = TH1F(f"h0_{i}", f"Histogram Channel 0_{i}", 2**13, 0, 2**13-1)
            histCh32 = TH1F(f"h32_{i}", f"Histogram Channel 32_{i}", 2**13, 0, 2**13-1)
            for i, line in enumerate(file):
                if (i>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        histCh32.Fill(int(event[-2]))
            histCh0Co1 = histCh0.Clone()
            histCh0Co2 = histCh0.Clone()
            histCh32Co1 = histCh32.Clone()
            histCh32Co2 = histCh32.Clone()

            par0 = getPeak(histCh0, histCh0.FindBin(config['start']), histCh0.FindBin(config['end']))
            par32 = getPeak(histCh32, histCh32.FindBin(config['start']), histCh32.FindBin(config['end']))
            parCh0Co1 = getPeak(histCh0Co1, histCh0.FindBin(config['startCo1']), histCh0.FindBin(config['endCo1']))
            parCh0Co2 = getPeak(histCh0Co2, histCh0.FindBin(config['startCo2']), histCh0.FindBin(config['endCo2']))
            parCh32Co1 = getPeak(histCh32Co1, histCh0.FindBin(config['startCo1']), histCh0.FindBin(config['endCo1']))
            parCh32Co2 = getPeak(histCh32Co2, histCh0.FindBin(config['startCo2']), histCh0.FindBin(config['endCo2']))
 
            print(par0)
            print(par32)
            print(parCh0Co1)
            print(parCh0Co2)
            print(parCh32Co1)
            print(parCh32Co2)

            fitFunc0 = fit(histCh0, par0, config['start'], config['end'])
            fitFunc32 = fit(histCh32, par32, config['start'], config['end'])
            fitfuncCh0Co1 = fit(histCh0Co1, parCh0Co1, config['startCo1'], config['endCo1'])
            fitfuncCh0Co2 = fit(histCh0Co2, parCh0Co2, config['startCo2'], config['endCo2'])
            fitfuncCh32Co1 = fit(histCh32Co1, parCh32Co1, config['startCo1'], config['endCo1'])
            fitfuncCh32Co2 = fit(histCh32Co2, parCh32Co2, config['startCo2'], config['endCo2'])

            totalCounts0 = integrate(histCh0, config['specStart'], config['specEnd'])
            totalCounts32 = integrate(histCh0, config['specStart'], config['specEnd'])

            Ch0 = TCanvas('c1', "Channel 0")
            Ch32 = TCanvas('c2', 'Channel 32')

            Ch0.cd()
            histCh0.Draw('hist')
            histCh0.Fit(fitFunc0, "R")
            fitFunc0.Draw("same")

            histCh0Co1.Draw('same')
            histCh0Co1.Fit(fitfuncCh0Co1, "R")
            fitfuncCh0Co1.Draw('same')

            histCh0Co2.Draw('same')
            histCh0Co2.Fit(fitfuncCh0Co2, "R")
            fitfuncCh0Co2.Draw('same')

            Ch32.cd()
            histCh32.Draw('hist')
            histCh32.Fit(fitFunc32, "R")
            fitFunc32.Draw('same')

            histCh32Co1.Draw('same')
            histCh32Co1.Fit(fitfuncCh32Co1, "R")
            fitfuncCh32Co1.Draw('same')

            histCh32Co2.Draw('same')
            histCh32Co2.Fit(fitfuncCh32Co2, "R")
            fitfuncCh32Co2.Draw('same')    
           
            gPad.WaitPrimitive('ggg')
            
            rez0 = fitFunc0.GetParameter(2) 
            rez32 = fitFunc32.GetParameter(2) 

            CsArea0 = gaussianArea(fitFunc0)
            CsArea32 = gaussianArea(fitFunc32)

            eff0 = totalCounts0/CsArea0
            eff32 = totalCounts32/CsArea32 

            pCh0 = np.append(pCh0, rez0)
            # pCh0 = np.append(pCh0, eff0)
            print(pCh0)

            pCh32 = np.append(pCh32, rez32)
            # pCh32 = np.append(pCh32, eff32)
            print(pCh32)
            
            parCh0 = np.append(parCh0, pCh0)
            parCh32 = np.append(parCh32, pCh32)

            histsCh0.append(histCh0)
            histsCh32.append(histCh32)
        
    print(parCh32)
    print(parCh0)


    return parCh0, parCh32 

def plotGraphs(parms1, parms2, parms3): #Ploteaza mai multe grafice 
    c1 = TCanvas('graphs0', 'Hold Delay vs Sigma graphs Ch0')
    c2 = TCanvas('graphs32', 'Hold Delay vs Sigma graphs Ch32')
    graph1 = TGraph()
    graph1.SetName('Hold Delay vs Sigma Ch0')
    graph1.SetMarkerSize(1.5)
    graph1.SetMarkerStyle(20)
    graph1.GetXaxis().SetTitle("Hold Delay")
    graph1.GetYaxis().SetTitle("Sigma")

    graph2 = TGraph()
    graph2.SetName('Hold Delay vs Sigma Ch32')
    graph2.SetMarkerSize(1.5)
    graph2.SetMarkerStyle(20)
    graph2.GetXaxis().SetTitle("Hold Delay")
    graph2.GetYaxis().SetTitle("Sigma")
    graph2.SetMarkerColor(2)
    for i in range(len(parms1)):
        graph1.AddPoint(float(parms1[i]), float(parms2[i]))
        graph2.AddPoint(float(parms1[i]), float(parms3[i]))
    file = TFile('graphs.root', 'RECREATE')
    c1.cd()
    graph1.Draw('ap')
    c2.cd()
    graph2.Draw('ap')
    c1.Write()
    c2.Write()
    gPad.WaitPrimitive('ggg')

# plot(config['dataFiles'])
info = getInfo(config['infoFiles'])
parameters = analyseFiles(config['dataFiles'])
# plotGraphs(info, parameters[0], parameters[1])




 



     

