#!/usr/bin/env python
#------------------------------------------------------------------
# Description: complete example of an unbinned fit to CMS H . gg
#              discovery data - HBP
#$Revision:$
#------------------------------------------------------------------
import os,sys
from histutil import setStyle
from time import sleep
from ROOT import *
#------------------------------------------------------------------		
def main():

	masspoints = [200,250,300,320,340,360,380,400,420]

	for masspoint in masspoints:
		# Suppress all messages except those that matter
		RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
	
		file = open('FitParameters.txt', 'a')

		# Make the workspace
		wspace = RooWorkspace('ThreePhoton')

		# define fit range
		massmin=masspoint - masspoint/30.0
		massmax=masspoint + masspoint/30.0
		massrange = massmax-massmin

		#define ridiculously large weight range
		wmin=-5000.0
		wmax=5000.0
	
		# Read data
		wspace.factory('x[%f, %f]' % (massmin, massmax))
		wspace.factory('w[%f, %f]' % (wmin, wmax))
		w = wspace.var('w')
		w.SetTitle('weights')
		x = wspace.var('x')
		x.SetTitle('m_{#gamma#gamma} (GeV)')
		# 7 TeV
		data = RooDataSet.read('/Users/diamond/2012Tuples/OUTPUTHIST/PhoIDHistsMC%s_InvMassUnbinned.txt' %(masspoint), RooArgList(x,w))
		# troubleshooting to see contents of data variable
		#data.write('dataout')
		#print

		#instruct wdata to interpret w as event weight rather than observable
		wdata = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"",w.GetName())
	
		getattr(wspace,'import')(wdata)

		ndata = wdata.numEntries()
		print "number of entries: %d" % ndata

		# create background parameters
		wspace.factory('background[%d, 0.0, %d]' % (ndata, 3*ndata))	
		wspace.factory('A[1., %d, %d]' % (-2*ndata,2*ndata))
		wspace.factory('B[0., %d, %d]' % (-2*ndata,2*ndata))
	
		# Background model
		#wspace.factory('Chebychev::bmodel(x, {A, B})')
		#wspace.factory('Gaussian::bmodel(x, A, B)')
		wspace.factory('GenericPdf::bmodel("A*x*x + B*x + 1",{x, A, B})')
		bmodel  = wspace.pdf('bmodel')
	
	
		# Signal model:CB Parameters
		wspace.factory('mass[%f, %f, %f]' % (massmin, massmin, massmax)) 
		wspace.factory('width[2.0,0.5,5.0]')
		wspace.factory('a[1.0,0.5,3.0]')
		wspace.factory('n[1.0,0.5,500.0]')

		wspace.factory('signal[%f, 0, %f]' % (ndata/2, 2*ndata))
	#	wspace.factory('mass[250.0, %f, %f]' % (massmin, massmax))
	#	wspace.factory('width[2.0,0.0,10.0]')
	# 	wspace.factory('a[0.5,0,2.0]')
	# 	wspace.factory('n[0.5,0,10.0]')
		wspace.factory('CBShape::smodel(x, mass, width, a, n)'); 
	#	wspace.factory('Gaussian::smodel(x, mass, width)')
		smodel   =  wspace.pdf('smodel')
	
	

		print 'signal title %s' % wspace.var('signal').GetTitle()
		print 'signal min = %f, max = %f' % (wspace.var('signal').getMin(),
							 wspace.var('signal').getMax())
		print ndata
		# According to Wouter, the following automatically
		# produces an extended pdf
		# p(x|...) = exp(-(background+signal)) * PROD p(x_i|...)
		# B + S = avg number of events
		wspace.factory('SUM::model(background*bmodel, signal*smodel)')
		model  = wspace.pdf('model')
	
		print "="*80
		wspace.Print()
		print "="*80

		# ----------------------------
		# Fit
		# ----------------------------
		print "fit extended model to data"
		print "="*80

		swatch = TStopwatch()
		swatch.Start()
		results = model.fitTo(wdata, RooFit.Save())
	# -2 ln l(theta)/l(max) ~ X^2
		print "real time: %10.3f s" % swatch.RealTime()
	
		vbkg = wspace.var('background')
		vsig = wspace.var('signal')
		vmass= wspace.var('mass')
		vwidth=wspace.var('width')

		# Get results and compute a simple measure of signal significance
		bkg  = vbkg.getVal()
		ebkg = vbkg.getError()
		sig  = vsig.getVal()
		esig = vsig.getError()
		mass  = vmass.getVal()
		emass = vmass.getError()
		width = vwidth.getVal()
		ewidth= vwidth.getError()
		zvalue= sig / esig
	
		print "="*80
		print "background: %10.1f +\-%-5.1f GeV" % (bkg, ebkg)
		print "signal:     %10.1f +\-%-5.1f" % (sig, esig)
		print "mass:       %10.1f +\-%-4.1f GeV" % (mass, emass)
		print "width:      %10.1f +\-%-4.1f GeV" % (width, ewidth)
		print "sig/esig:   %10.1f" % zvalue
		print
	
		# now plot results of fit

		setStyle()
		xframe = x.frame(RooFit.Bins(50))
		xframe.GetXaxis().SetNdivisions(505)
	
		wdata.plotOn(xframe)
		model.plotOn(xframe)
	
	
		# Hack to prevent display of background parameters
		wspace.var('A').setConstant()
		wspace.var('B').setConstant()
		wspace.var('n').setConstant()
		wspace.var('a').setConstant()
		model.paramOn(xframe)

		chi2 = xframe.chiSquare() 
		print "chiSq/NDF :      %f" % chi2
		print
		file.write('{0:10} {1:10.2f} {2:10.2f} {3:10.2f} {4:10.2f} {5:10.2f}\n'.format(masspoint, mass, emass, width, ewidth, chi2))

		model.plotOn(xframe,RooFit.Components('smodel'), RooFit.LineColor(kRed+1))
		model.plotOn(xframe,RooFit.Components('bmodel'), RooFit.LineColor(kMagenta+1))
	
		c1 = TCanvas('fig_InvMass_fit_' + str(masspoint), 'fit', 10, 10, 800, 800)
		ymax = xframe.GetMaximum()
		xframe.SetMaximum(ymax*1.5)
	#	.Integral()
		xframe.Draw()
		c1.SaveAs('.pdf')
		c1.SaveAs('.gif')


		sleep(5)
#------------------------------------------------------------------
try:
	main()
except KeyboardInterrupt:
	print
	print "bye cruel world!"
	print
	


