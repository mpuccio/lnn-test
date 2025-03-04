import ROOT

pt_range = [4,8]

file = ROOT.TFile.Open("bw/H3L_BWFit.root")
bw_0_10 = file.Get("BlastWave_H3L_0_10")
bw_10_30 = file.Get("BlastWave_H3L_10_30")
bw_30_50 = file.Get("BlastWave_H3L_30_50")

bw_integral_010 = bw_0_10.Integral(pt_range[0],pt_range[1])
bw_integral_10_30 = bw_10_30.Integral(pt_range[0],pt_range[1])
bw_integral_30_50 = bw_30_50.Integral(pt_range[0],pt_range[1])

print("0-10%: ", bw_integral_010)
print("10-30%: ", bw_integral_10_30)
print("30-50%: ", bw_integral_30_50)


bw_integral_050 = (bw_integral_010 + 2 * bw_integral_10_30 + 2 * bw_integral_30_50) / 5
print("0-50%: ", bw_integral_050)
print('dN/dydpt: ',bw_integral_050 / (pt_range[1]-pt_range[0]))
