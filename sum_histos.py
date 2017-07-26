import ROOT as R
R.gROOT.SetBatch(True)                        # go batch!


histos = ['n_mu', 'n_mu6ab', 'n_pairs_mu6_0dr15', 'n_pairs_0dr15', 'n_pairs_mu6ab', 'n_cand_pairs', 'dr_min', 'dr_0dr15', 'dr_mu6', 'dr_mu6_min', 'dr_mu6_0dr15', 'dr_any', 'Phi_mu6', 'Eta_mu6', 'pt_0dr15', 'pt_any']
outcomes = ['pass_em_pass_hw', 'pass_em_fail_hw', 'fail_em_pass_hw', 'fail_em_fail_hw']
histo_names = [name+'_'+outcome for outcome in outcomes for name in histos]

c = R.TCanvas('c')
for name in histos:
	i=1
	c.Clear()
	c.Divide(2,2)
	for outcome in outcomes:
		#get first histogram
		canvasfile = R.TFile('./Histograms/file0/'+name+'.root')
		C = R.TCanvas('C')
		canvasfile.GetObject('c', C)
		pad = C.GetPad(i)
		sum_hist = pad.GetPrimitive(name+'_'+outcome)

		for j in range(1,6): #sum over all histograms
			canvasfile = R.TFile('./Histograms/file'+str(j)+'/'+name+'.root')
		 	C = R.TCanvas('C')
			canvasfile.GetObject('c', C)
			pad = C.GetPad(i)
			h = pad.GetPrimitive(name+'_'+outcome)
			sum_hist.Add(h)
		#print result
		c.cd(i)
		sum_hist.Draw('h text')
		c.Update()
		i+=1
		sum_hist.Print("all")
		print("- Mean     = %.3f +- %.3f"%(sum_hist.GetMean(),sum_hist.GetMeanError()))
		print("- Std Dev  = %.3f +- %.3f"%(sum_hist.GetStdDev(), sum_hist.GetStdDevError()))
		print("- Skewness = %.3f"%(sum_hist.GetSkewness()))
		print("- Kurtosis = %.3f"%(sum_hist.GetKurtosis()))
	
	print('\n\n')
	c.SaveAs(name+'.png')
	c.SaveAs(name+'.root')

	i=1
	c.Clear()
	c.Divide(2,2)
	for outcome in outcomes:
		#get first histogram
		canvasfile = R.TFile('./Histograms/file0/PhiEta_mu6.root')
		C = R.TCanvas('C')
		canvasfile.GetObject('c', C)
		pad = C.GetPad(i)
		sum_hist = pad.GetPrimitive('PhiEta_mu6_'+outcome)

		for j in range(1,6): #sum over all histograms
			canvasfile = R.TFile('./Histograms/file'+str(j)+'/PhiEta_mu6.root')
			C = R.TCanvas('C')
			canvasfile.GetObject('c', C)
			pad = C.GetPad(i)
			h = pad.GetPrimitive('PhiEta_mu6_'+outcome)
			sum_hist.Add(h)
		c.cd(i)
		sum_hist.Draw('Colz')
		c.Update()
		i+=1
		sum_hist.Print("all")
	
	print('\n')

	c.SaveAs('PhiEta_mu6.png')
	c.SaveAs('PhiEta_mu6.root')