import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("/home/cczhu/GitHubTemp/czerny_wd/")
import StarModSteve as sms
import pandas as pd
import cPickle as cP
import copy
import rhoTcontours as rtc
sys.path.append("/home/cczhu/Runaway/")
from runaway_prod_anal_func import *
import gc
import h5py
import matplotlib.lines as mlines

def schwab_vals():

	# rho_c read off of table 3, rho_max and temp_max read off of their respective figures
	# from figures, temp_c is always tiny.  Mtot_new is the sum of Mc and Mtp in table 3

	schwab = {"Md": np.array([0.2, 0.3, 0.5, 0.3, 0.6, 0.2, 0.3, 0.9]),
			"Ma": np.array([0.8, 1.1, 1.2, 0.6, 0.9, 0.3, 0.4, 1.2]),
			"Tc_after": 1e8*np.ones(8),
			"rho(Tc)": np.array([8.8e6, 4.7e7, 9.5e7, 3.8e6, 2.8e7, 6.4e5, 1.5e6, 3.3e8]),
			"Tmax_after": np.array([0., 7.5e8, 0., 3e8, 7.25e8, 1e8, 2e8, 0.]),
			"rho(Tmax)": np.array([0., 9e4, 0., 1.2e5, 7e5, 5e2, 1.5e3, 0.]),
			"Mtot_new": np.array([0.71+0.1, 0.98+0.12, 1.05+0.16, 0.53+0.13, 0.84+0.2, 0.28+0.08, 0.38+0.12,1.11+0.24])}

 	schwab = pd.DataFrame(schwab)
	schwab["Mtot"] = schwab["Ma"] + schwab["Md"]

	return schwab

def getscalarmap(vmin, vmax, cmap_name, logv=False, fullreturn=False):
	"""obtains scalarMap instance for using colormap on lines and points"""
	my_cmap = plt.get_cmap(cmap_name)
	if logv:
		cNorm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
	else:
		cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	if fullreturn:
		return [matplotlib.cm.ScalarMappable(norm=cNorm, cmap=my_cmap), cNorm, my_cmap]
	else:
		return matplotlib.cm.ScalarMappable(norm=cNorm, cmap=my_cmap)


def get_runaway_line(od, ax, td, od_wwk04_i=False, wwkmfc="k", color="red", wwkmarker="o", wwksize=5, label="", lw=1, ls='-', apls = '-.'):
	od_ign = find_igniter(od, td=td)
	if not od_wwk04_i:
		[od_wwk04_i, od_dyn_i, od_sound_i] = find_exploder(od, td=td)
	ax.plot(od["dens_c"][:od_ign + 1], od["temp_c"][:od_ign + 1], ls=apls, color=color, lw=lw, label=label)
	if od_wwk04_i == -1:
		ax.plot(od["dens_c"][od_ign:], od["temp_c"][od_ign:], ls=ls, color=color, lw=lw, label=label)
	else:
		ax.plot(od["dens_c"][od_ign:od_wwk04_i + 1], od["temp_c"][od_ign:od_wwk04_i + 1], ls=ls, color=color, lw=lw, label=label)
		ax.plot(od["dens_c"][od_wwk04_i + 1:], od["temp_c"][od_wwk04_i + 1:], ls=apls, color=color, lw=lw, label=label)
		ax.plot(od["dens_c"][od_wwk04_i], od["temp_c"][od_wwk04_i], 'ro', marker=wwkmarker, markersize=wwksize, markerfacecolor=wwkmfc, markeredgecolor="k")


def getzhu13pme_sub(files, filedir="/home/cczhu/SunnyvaleOutput/zhu16/"):

	outp = {"Mtot": np.zeros(len(files)), "Mce": np.zeros(len(files)), "Mdisk": np.zeros(len(files)), 
			"Mtot_new": np.zeros(len(files)), "Txy": np.zeros([len(files),2]), "Tz": np.zeros([len(files),2]),
			"Tc_before": np.zeros([len(files),2]), "Tc_after": np.zeros([len(files),2]),
			"Txy_before": np.zeros([len(files),2])}

	for i in range(len(files)):
		data = cP.load(open(filedir + files[i], 'r'))
		outp["Mtot"][i] = data["Mtot"]
		outp["Mce"][i] = data["Mce"]
		outp["Mdisk"][i] = data["Mdisk"]
		outp["Mtot_new"][i] = data["Mtot_new"]
		outp["Txy_before"][i,0] = data["Tmax"][0]; outp["Txy_before"][i,1] = data["Tmax"][1]
		outp["Txy"][i,0] = data["Tmax_new"][0]; outp["Txy"][i,1] = data["Tmax_new"][1]
		outp["Tz"][i,0] = data["zTmax_new"][0]; outp["Tz"][i,1] = data["zTmax_new"][1]
		outp["Tc_before"][i,0] = data["dens_before"][0]; outp["Tc_before"][i,1] = data["temp_before"][0]
		outp["Tc_after"][i,0] = data["dens_after"][0]; outp["Tc_after"][i,1] = data["temp_after"][0]

	return outp


def load_ji():

	ji13r = h5py.File('../../PaperRunaway/ji_radial_profile.hdf5', 'r')
	ji13 = {"r": np.array(ji13r["r"]),
			"omega": np.array(ji13r["omega"]),
			"temp": np.array(ji13r["temp"]),
			"dens": np.array(ji13r["dens"])}

	mystar = sms.mhs_steve(0., 1e8, mintemp=3e7, stop_mindenserr=1e0, dontintegrate=True)
	
	ji13["eint"] = np.zeros(len(ji13["r"]))
	ji13["edeg"] = np.zeros(len(ji13["r"]))
	ji13["press"] = np.zeros(len(ji13["r"]))
	ji13["ent"] = np.zeros(len(ji13["r"]))
	for i in range(len(ji13["r"])):
		ji13["press"][i], ji13["ent"][i] = mystar.getpress_rhoT(ji13["dens"][i], ji13["temp"][i])
		ji13["eint"][i], ji13["edeg"][i], dummy = mystar.gethelmeos_energies(ji13["dens"][i], ji13["temp"][i])

	ji13["m"] = 4.*np.pi*scipyinteg.cumtrapz(ji13["r"]**2*ji13["dens"], x=ji13["r"], initial=0.)

	return ji13


data = pd.read_csv("../../PaperRunaway/zhu13ichart.csv", sep="\t", skiprows=2)

files = ["pt4pt4o.p","pt4pt5o.p","pt5pt5o.p","pt4pt55o.p","pt5pt55o.p","pt55pt55o.p","pt4pt6o.p","pt5pt6o.p","pt55pt6o.p","pt6pt6o.p",
		"pt4pt65o.p","pt5pt65o.p","pt55pt65o.p","pt575pt65o.p","pt6pt65o.p","pt625pt65o.p","pt64pt65o.p","pt65pt65o.p",
		"pt4pt7o.p","pt5pt7o.p","pt55pt7o.p","pt6pt7o.p","pt65pt7o.p","pt7pt7o.p","pt4pt8o.p","pt5pt8o.p","pt55pt8o.p",
		"pt6pt8o.p","pt65pt8o.p","pt7pt8o.p","pt8pt8o.p","pt4pt9o.p","pt5pt9o.p","pt55pt9o.p","pt6pt9o.p","pt65pt9o.p",
		"pt7pt9o.p","pt8pt9o.p","pt9pt9o.p","pt4pt10o.p","pt5pt10o.p","pt55pt10o.p","pt6pt10o.p","pt65pt10o.p","pt7pt10o.p",
		"pt8pt10o.p","pt9pt10o.p","pt10pt10o.p"]

outp = getzhu13pme_sub(files)
outp["Mtot_new"] /= 1.9891e33

for item in ["Mtot", "Mce", "Mdisk", "Mtot_new"]:
	data[item] = pd.Series(outp[item])
data["mrho_xy"] = pd.Series(outp["Txy"][:,0])
data["mT_xy"] = pd.Series(outp["Txy"][:,1])
data["mrho_xy_before"] = pd.Series(outp["Txy_before"][:,0])
data["mT_xy_before"] = pd.Series(outp["Txy_before"][:,1])
data["mrho_z"] = pd.Series(outp["Tz"][:,0])
data["mT_z"] = pd.Series(outp["Tz"][:,1])
data["crho_xy_before"] = pd.Series(outp["Tc_before"][:,0])
data["cT_xy_before"] = pd.Series(outp["Tc_before"][:,1])
data["crho_xy"] = pd.Series(outp["Tc_after"][:,0])
data["cT_xy"] = pd.Series(outp["Tc_after"][:,1])

colorlist = ['red','orange','lime','cyan','blue','magenta','purple','brown','black']
labellist = ["0.4", "0.5", "0.55", "0.6", "0.65", "0.7", "0.8", "0.9", "1.0"]
Macc = np.array([ 0.40231300000000003,  0.502649,  0.553285,  0.603529,  0.653632,  0.703745, 0.80457,  0.904897,  1.005611])

p_ar = cP.load(open("/home/cczhu/SunnyvaleOutput/zhu16/pt625pt65AREPO.p", 'r'))
ji13 = load_ji()
schw12 = schwab_vals()

width=6.; height=5.; axes=[10**5.5,10**8.5,10**7.8,10**9.75]

pylab.rc('font', family="serif")
pylab.rc('font', size=14)
fig = pylab.figure(figsize=(width,height))
ax = fig.add_axes([0.8/width, 0.75/height, (width - 1.2)/width, (height - 1.)/height])

zhu13_options = {"marker": "o", "ls": "None"}
arepo_options = {"mfc": "b", "mec": "r", "marker": "*", "ls": "None", "mew": 2}
schwab_options = {"mec": "m", "marker": "x", "ls": "None", "mew": 3}
ji_options = {"mec": "c", "marker": "x", "ls": "None", "mew": 3}
fitline_options = {"color": "k", "ls": ":", "lw": "2"}

for i in range(len(Macc)):
	args = abs(data["M2 (Msun) "] - Macc[i])/Macc[i] < 1e-6
	ax.plot(data.loc[args, "Mtot"], data.loc[args, "Mtot_new"], mfc=colorlist[i], **zhu13_options)

ax.plot(p_ar["Mtot"], p_ar["Mtot_new"]/1.9891e33, markersize=20, **arepo_options)
ax.plot(schw12["Mtot"], schw12["Mtot_new"], markersize=10, **schwab_options)
ax.plot(np.array([1.2]), np.array(ji13["m"][max(np.where(ji13["r"] <= 2e9)[0])]/1.9891e33), markersize=10, **ji_options)
fit_out = np.polyfit(data["Mtot"], data["Mtot_new"], 1)
print "Fit y = {0:.3e}x + {1:.3e}".format(fit_out[0],fit_out[1])
fitx = np.arange(0.25,2.2,0.01)
fitfunc = lambda x: fit_out[1] + fit_out[0]*x
fity = fitfunc(fitx)
ax.plot(fitx, fity, **fitline_options)
line_artists = (mlines.Line2D([],[], mfc="k", ms=10, **zhu13_options),
				mlines.Line2D([],[], **fitline_options),
				mlines.Line2D([],[], ms=20, **arepo_options),
				mlines.Line2D([],[], ms=10, **schwab_options),
				mlines.Line2D([],[], ms=10, **ji_options))
leg = ax.legend(line_artists, ("Zhu+13 Estimate", "Zhu+13 Fit", "Zhu+15 Arepo", "Schwab+12", "Ji+13"), loc=2, fontsize=13, numpoints = 1, markerscale=1)
leg.draw_frame(False)
#ax.set_xlim(0.1,2.2); ax.set_ylim(0.0, 1.6)
ax.set_xlabel(r"$M_\mathrm{tot}$ ($M_\odot$)")
ax.tick_params(axis='x', which='major', pad=7)
ax.set_ylabel(r"$M_\mathrm{c, pv}$ ($M_\odot$)")

print "mean(Mtot_new - Mtot)_schw12 = {0:.3e}".format(np.mean(abs(schw12["Mtot_new"] - fitfunc(schw12["Mtot"]))))

pylab.savefig("figures/c2a_mcpv.pdf", dpi=150)
