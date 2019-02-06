
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

class FreysoldtPlotter(object):
    """
    This class plots Freysoldt planar-averaged potential alignment
    """
    def __init__(self, x, Vr, dft_diff, final_shift, check):
        self.x = x
        self.Vr = Vr
        self.dft_diff = dft_diff
        self.final_shift = final_shift
        self.check = check

    def plot(self, title='freysoldt'):
        plt.figure()
        plt.clf()
        plt.plot(self.x, self.Vr, c="green", zorder=1,
                 label="long range from model")
        plt.plot(self.x, self.dft_diff, c="red", label="DFT locpot diff")
        plt.plot(self.x, self.final_shift, c="blue",
                 label="short range (aligned)")
        tmpx = [self.x[i] for i in range(self.check[0], self.check[1])]
        plt.fill_between(tmpx, -100, 100, facecolor='red', alpha=0.15,
                         label='sampling region')

        plt.xlim(round(self.x[0]), round(self.x[-1]))
        ymin = min(min(self.Vr), min(self.dft_diff), min(self.final_shift))
        ymax = max(max(self.Vr), max(self.dft_diff), max(self.final_shift))
        plt.ylim(-0.2+ymin, 0.2+ymax)
        plt.xlabel('distance along axis ($\AA$)', fontsize=20)
        plt.ylabel('Potential (V)', fontsize=20)
        plt.legend(loc=9)
        plt.axhline(y=0, linewidth=0.2, color='black')
        plt.title(str(title) + ' defect potential')
        plt.xlim(0, max(self.x))

        return plt

    def as_dict(self):
        d = {"x": self.x, "Vr": self.Vr, "dft_diff": self.dft_diff,
             "final_shift": self.final_shift, "check": self.check
             }
        return d

    @staticmethod
    def from_dict( d):
        return FreysoldtPlotter( d["x"], d["Vr"], d["dft_diff"],
                                 d["final_shift"], d["check"])


class KumagaiPlotter(object):
    """
    This class plots Freysoldt planar-averaged potential alignment

    takes site_dict, which is a dictionary with keys as index of site
    and values:
        'Vpc' = point charge potential
        'Vqb' = difference in DFT electrostatic potential between bulk and defect
        'dist_to_defect' = cartesian distance to defect
    """

    def __init__(self, site_dict, sampling_radius, potalign):
        self.site_dict = site_dict
        self.sampling_radius = sampling_radius
        self.potalign = potalign

    def plot(self, title='kumagai'):
        plt.figure()
        plt.clf()

        distances, sample_region  =  [], []
        Vqb_list, Vpc_list, diff_list = [], [], []
        for site_ind, site_dict in self.site_dict.items():
            dist = site_dict["dist_to_defect"]
            distances.append( dist)

            Vqb = site_dict["Vqb"]
            Vpc = site_dict["Vpc"]

            Vqb_list.append( Vqb)
            Vpc_list.append( Vpc)
            diff_list.append( Vqb - Vpc)

            if dist > self.sampling_radius:
                sample_region.append( Vqb - Vpc)

        plt.plot( distances, Vqb_list,
                  color='r', marker='^', linestyle='None',
                  label='$V_{q/b}$')

        plt.plot( distances, Vpc_list,
                  color='g', marker='o', linestyle='None',
                  label='$V_{pc}$')

        plt.plot( distances, diff_list, color='b', marker='x', linestyle='None',
                 label='$V_{q/b}$ - $V_{pc}$')

        x = np.arange(self.sampling_radius, max(distances) * 1.05, 0.01)
        y_max = max( max( Vqb_list), max( Vpc_list), max( diff_list)) + .1
        y_min = min( min( Vqb_list), min( Vpc_list), min( diff_list)) - .1
        plt.fill_between(x, y_min, y_max, facecolor='red',
                         alpha=0.15, label='sampling region')
        plt.axhline(y=self.potalign, linewidth=0.5, color='red',
                    label='pot. align. / -q')

        plt.legend(loc=0)
        plt.axhline(y=0, linewidth=0.2, color='black')

        plt.ylim([y_min, y_max])
        plt.xlim([0, max(distances)*1.1 ])

        plt.xlabel('Distance from defect ($\AA$)',fontsize=20)
        plt.ylabel('Potential (V)',fontsize=20)
        plt.title('%s atomic site potential plot' % title, fontsize=20)

        return plt

    def as_dict(self):
        d = {"site_dict": self.site_dict, "sampling_radius": self.sampling_radius,
             "potalign": self.potalign
             }
        return d

    @staticmethod
    def from_dict(d):
        return KumagaiPlotter(d["site_dict"], d["sampling_radius"], d["potalign"])

