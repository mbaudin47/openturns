# -*- coding: utf-8 -*-
"""
Created on Thu Jun 01 11:48:42 2017

@author: c61372, dumas

Analyse la précision des estimateurs des intervalles de confiance 
pour la fonction G-Sobol et Ishigami
"""

#! /usr/bin/env python

from __future__ import print_function
import matplotlib as mpl
mpl.use('agg')
import openturns as ot
import numpy as np
from openturns.viewer import View
from gsobollib import (gsobolSAExact, 
                       gsobolDistribution, gsobol)
from ishigamilib import (ishigamiSAExact, ishigamiDistribution, 
                         ishigamiAB, ishigamiGSymbolic)
from numpy import zeros, sqrt, array, linspace
from pylab import plot, show, xlabel, ylabel, xscale, yscale, legend, title, savefig, subplot, hist
import pylab as pl
from multiprocessing import Pool
import pathlib
from sobol_variance_estimators import (SaltelliSensitivityAlgorithm,
                                       MartinezSensitivityAlgorithm,
                                       JansenSensitivityAlgorithm,
                                       MauntzKucherenkoSensitivityAlgorithm)


class SensitivityConfidenceTest():
    """
    Class pour valider la variance des estimateurs de sobols par la méthode Janon.
    La distribution asymptotique est comparé à la distribution empirique obtenue
    par répétitions du calcul des estimateurs de Sobol.
    """
    def __init__(self, model, distribution, sobol_estimator, FOexact=None,
                 TOexact=None, sampleSize=5000, seed=154681,
                 nrepetitions=1000, alpha=0.95, savefig=False, plot_figure=False):

        # Le modèle test : g-sobol ou ishigami
        self.model = model
        # la distribution des paramètres du modèle test
        self.distribution = distribution
        self.dim = distribution.getDimension()
        # classe améliorée pour le calcul des indices (par ex. SaltelliSensitivityAlgorithm) 
        self.sobol_estimator = sobol_estimator
        # Taille du plan d'expérience de base pour estimer S et ST
        self.sampleSize = sampleSize
        # Nombre de répétition de l'expérience
        self.nrepetitions = nrepetitions
        # bootstrap confidence level
        self.alpha = alpha
        # flag pour sauvegarder les figures
        self.savefig = savefig
        # first order exact indices if known
        self.FOexact = FOexact
        # total order exact indices if known
        self.TOexact = TOexact
        # seed
        self.seed = seed

        sampleFirst, sampleTotal, foInterval, toInterval, distFirstCol, distTotalCol = self.compute_sample_indices()
        if plot_figure:
            self.compare_last_repetition(sampleFirst, sampleTotal, foInterval, toInterval, distFirstCol, distTotalCol)
            self.plot_indices_histogram(sampleFirst, sampleTotal, distFirstCol, distTotalCol)
            self.plot_indices_histogram(sampleFirst, sampleTotal, distFirstCol, distTotalCol, True)

    def compute_sample_indices(self):

        # Estimations des indices du premier ordre
        sampleFirst = zeros((self.nrepetitions,self.dim))

        # Estimations des indices totaux
        sampleTotal = zeros((self.nrepetitions,self.dim))

        # loi asymptotique
        distFirstCol = [object] * self.nrepetitions
        distTotalCol = [object] * self.nrepetitions

        # set seed of the random generator
        ot.RandomGenerator.SetSeed(self.seed)
        for i in range(self.nrepetitions):
            sobolexperiment = ot.SobolIndicesExperiment(self.distribution,
                                                        int(self.sampleSize),
                                                        False)
            inputDesign = sobolexperiment.generate()
            outputDesign = self.model(inputDesign)
            self.sensitivity_algorithm = self.sobol_estimator(
                inputDesign, outputDesign, int(self.sampleSize))
            # self.sensitivity_algorithm = self.sobol_estimator(self.distribution,
            #                             int(self.sampleSize), self.model)
            self.sensitivity_algorithm.setBootstrapConfidenceLevel(self.alpha)
            fo = self.sensitivity_algorithm.getAggregatedFirstOrderIndices()
            to = self.sensitivity_algorithm.getAggregatedTotalOrderIndices()
            # Récupère les distributions asymptotiques
            distFirstCol[i] = self.sensitivity_algorithm.getFirstOrderAsymptoticDistribution()
            distTotalCol[i] = self.sensitivity_algorithm.getTotalOrderAsymptoticDistribution()
            for j in range(self.dim):
                sampleFirst[i, j] = fo[j]
            for j in range(self.dim):
                sampleTotal[i, j] = to[j]

        # Récupère l'intervalle de confiance bootstrap pour le dernier échantillon
        foInterval = self.sensitivity_algorithm.getFirstOrderIndicesInterval()
        toInterval = self.sensitivity_algorithm.getTotalOrderIndicesInterval()

        # compute empirical variance
        self.std_first_empirical = ot.Sample(sampleFirst).computeStandardDeviation()
        self.std_total_empirical = ot.Sample(sampleTotal).computeStandardDeviation()
        return sampleFirst, sampleTotal, foInterval, toInterval, distFirstCol, distTotalCol

    def compare_last_repetition(self, sampleFirst, sampleTotal, foInterval, toInterval, distFirstCol, distTotalCol):

        # récupère les valeurs des min et max des intervalles
        foIntervalMin = foInterval.getLowerBound()
        foIntervalMax = foInterval.getUpperBound()
        toIntervalMin = toInterval.getLowerBound()
        toIntervalMax = toInterval.getUpperBound()
        # Compare les intervalles bootstrap pour le dernier échantillon 
        # et les quantiles issus des répétitions
        for j in range(self.dim):
            # Calcule les quantiles empiriques
            sampleFirstPerDim = ot.Sample(sampleFirst[:,j],1)
            foMinj = sampleFirstPerDim.computeQuantile((1-self.alpha)/2)[0]
            foMaxj = sampleFirstPerDim.computeQuantile(1-(1-self.alpha)/2)[0]
            sampleTotalPerDim = ot.Sample(sampleTotal[:,j],1)
            toMinj = sampleTotalPerDim.computeQuantile((1-self.alpha)/2)[0]
            toMaxj = sampleTotalPerDim.computeQuantile(1-(1-self.alpha)/2)[0]
            foAsympMinj = np.mean([distFirstCol[i][j].computeQuantile((1-self.alpha)/2)[0] 
                                   for i in range(self.nrepetitions)])
            foAsympMaxj = np.mean([distFirstCol[i][j].computeQuantile(1-(1-self.alpha)/2)[0] 
                                   for i in range(self.nrepetitions)])
            toAsympMinj = np.mean([distTotalCol[i][j].computeQuantile((1-self.alpha)/2)[0] 
                                   for i in range(self.nrepetitions)])
            toAsympMaxj = np.mean([distTotalCol[i][j].computeQuantile(1-(1-self.alpha)/2)[0] 
                                   for i in range(self.nrepetitions)])
            meanAsympFO = np.mean([distFirstCol[i][j].getStandardDeviation()[0] 
                                   for i in range(self.nrepetitions)])
            meanAsympTO = np.mean([distTotalCol[i][j].getStandardDeviation()[0]
                                   for i in range(self.nrepetitions)])

            print("X%d" % (j))
            print("   First standard deviation, Sample=%.5f, Asymptotic=%.5f" % \
                    (sampleFirstPerDim.computeStandardDeviation()[0, 0], meanAsympFO))
            print("   Total standard deviation, Sample=%.5f, Asymptotic=%.5f" % \
                    (sampleTotalPerDim.computeStandardDeviation()[0, 0], meanAsympTO))
            print("   First, Bootstrap=[%.4f,%.4f], Sample=[%.4f,%.4f], Asymptotic=[%.4f,%.4f]" % \
                    (foIntervalMin[j],foIntervalMax[j],foMinj,foMaxj, foAsympMinj, foAsympMaxj))
            print("   Total, Bootstrap=[%.4f,%.4f], Sample=[%.4f,%.4f], Asymptotic=[%.4f,%.4f]" % \
                    (toIntervalMin[j],toIntervalMax[j],toMinj,toMaxj, toAsympMinj, toAsympMaxj))
            print("")


    def plot_indices_histogram(self, sampleFirst, sampleTotal, distFirstCol,
                              distTotalCol, mean_distribution=False):
        fig, ax = pl.subplots(2,self.dim, figsize=(4*self.dim, 8))
        fig.suptitle("%s - N=%d - Repetitions = %d" % \
        (self.sobol_estimator.__name__, self.sampleSize,self.nrepetitions))

        # Pour chaque estimateur, compare la répartition empirique et la loi exacte
        for j in range(self.dim):
            # Indice du premier ordre
            sampleJ = sampleFirst[:,j]
            ax[0, j].hist(sampleJ,histtype="step",normed=True, label='empirique - std=%.2e' % np.std(sampleJ))
            if mean_distribution:
                # moyenne et ecart-type = moyenne des moyenne et écart-types des lois asymptotiques
                mu = np.mean([distFirstCol[i][j].getMean()[0] 
                                       for i in range(self.nrepetitions)]) # valeur moyenne
                sigma = np.mean([distFirstCol[i][j].getStandardDeviation()[0] 
                                       for i in range(self.nrepetitions)]) # écart-type
                label = 'asymptotique - std(mean)=%.2e' % sigma
            else:
                mu = distFirstCol[-1][j].getMean()[0]
                sigma = distFirstCol[-1][j].getStandardDeviation()[0]
                label = 'asymptotique - std=%.2e' % sigma
            loiFo = ot.Normal(mu,sigma)
            View(loiFo.drawPDF(), axes=[ax[0, j]], plot_kwargs={"label":label})
            if self.FOexact is not None:
                ax[0, j].vlines(self.FOexact[j], 0, ax[0,j].get_ylim()[1], label="exact - %0.3f" % self.FOexact[j])
            ax[0, j].set_xlabel("S%d" % (j))
            ax[0, j].set_ylabel("Density")
            ax[0, j].legend(loc='lower right')
            # Indice du total
            sampleJ = sampleTotal[:,j]
            ax[1, j].hist(sampleJ,histtype="step",normed=True, label='empirique - std=%.2e' % np.std(sampleJ))
            if mean_distribution:
                # moyenne et ecart-type = moyenne des moyenne et écart-types des lois asymptotiques
                mu = np.mean([distTotalCol[i][j].getMean()[0] 
                                       for i in range(self.nrepetitions)]) # valeur moyenne
                sigma = np.mean([distTotalCol[i][j].getStandardDeviation()[0] 
                                       for i in range(self.nrepetitions)]) # écart-type
                label = 'asymptotique- std(mean)=%.2e' % sigma
            else:
                mu = distTotalCol[-1][j].getMean()[0]
                sigma = distTotalCol[-1][j].getStandardDeviation()[0]
                label = 'asymptotique- std=%.2e' % sigma
            loiFo = ot.Normal(mu,sigma)
            View(loiFo.drawPDF(), axes=[ax[1, j]], plot_kwargs={"label":label})
            if self.TOexact is not None:
                ax[1, j].vlines(self.TOexact[j], 0, ax[1,j].get_ylim()[1], label="exact - %0.3f" % self.TOexact[j])
            ax[1, j].set_xlabel("ST%d" % (j))
            ax[1, j].set_ylabel("Density")
            ax[1, j].legend(loc='lower right')
        # fig.show()
        if mean_distribution:
            directory = 'graphe_validation_mean_distribution/'
        else:
            directory = 'graphe_validation/'
        # create directory if not exist
        pathlib.Path(directory).mkdir(parents=True, exist_ok=True) 
        if self.savefig:
            fig.savefig(directory+"%s-%s.png" % (self.model.getName(),
                                                   self.sobol_estimator.__name__),
                        transparent=True, bbox_inches="tight")


if __name__ == '__main__':

    ################################################################################
    #################                  GSOBOL                  #####################
    ################################################################################

    a = array([0,9,99])
    # Distribution uniforme associée au cas-test GSobol
    distribution_gsobol = gsobolDistribution(len(a))
    model_gsobol = ot.PythonFunction(len(a), 1, func_sample=lambda X: gsobol(X,a))
    model_gsobol.setName("G-Sobol")
    # Indices de sensibilité exacts
    [muexact,vexact,sexact,stexact] = gsobolSAExact(a)

    sTest_gsobol_saltelli = SensitivityConfidenceTest(model_gsobol, distribution_gsobol,
                                    SaltelliSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)

    sTest_gsobol_jansen = SensitivityConfidenceTest(model_gsobol, distribution_gsobol,
                                    JansenSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)

    sTest_gsobol_mauntz = SensitivityConfidenceTest(model_gsobol, distribution_gsobol,
                                    MauntzKucherenkoSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)

    sTest_gsobol_martinez = SensitivityConfidenceTest(model_gsobol, distribution_gsobol,
                                    MartinezSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)

    ################################################################################
    ################                  ISHIGAMI                  ####################
    ################################################################################

    model_ishigami = ishigamiGSymbolic()
    model_ishigami.setName("Ishigami")
    distribution_ishigami = ishigamiDistribution()
    # Indices de sensibilité exacts
    a, b = ishigamiAB()
    meanY,varY,S1,S2,S3,ST1,ST2,ST3 = ishigamiSAExact(a,b)
    sexact = ot.Point([S1, S2, S3])
    stexact = ot.Point([ST1, ST2, ST3])

    sTest_ishigami_saltelli = SensitivityConfidenceTest(model_ishigami, distribution_ishigami,
                                    SaltelliSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)

    sTest_ishigami_jansen = SensitivityConfidenceTest(model_ishigami, distribution_ishigami,
                                    JansenSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)

    sTest_ishigami_mauntz = SensitivityConfidenceTest(model_ishigami, distribution_ishigami,
                                    MauntzKucherenkoSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)

    sTest_ishigami_martinez = SensitivityConfidenceTest(model_ishigami, distribution_ishigami,
                                    MartinezSensitivityAlgorithm,
                                    FOexact=sexact, TOexact=stexact,
                                    savefig=True, plot_figure=True)


    ################################################################################
    #################                  POUTRE                  #####################
    ################################################################################


    model_poutre = ot.SymbolicFunction(['L', 'b', 'h', 'E', 'F'],
                                       ['F * L^3 / (48 * E * b * h^3 / 12)'])
    model_poutre.setName("poutre")
    L = ot.LogNormal()
    L.setParameter(ot.LogNormalMuSigmaOverMu()([5., .02, 0.]))
    b = ot.LogNormal()
    b.setParameter(ot.LogNormalMuSigmaOverMu()([.2, .05, 0.]))
    h = ot.LogNormal()
    h.setParameter(ot.LogNormalMuSigmaOverMu()([.4, .05, 0.]))
    E = ot.LogNormal()
    E.setParameter(ot.LogNormalMuSigmaOverMu()([3e4, .12, 0.]))
    F = ot.LogNormal()
    F.setParameter(ot.LogNormalMuSigmaOverMu()([.1, .20, 0.]))
    distribution_poutre = ot.ComposedDistribution([L, b, h, E, F])

    sTest_poutre_saltelli = SensitivityConfidenceTest(model_poutre, distribution_poutre,
                                    SaltelliSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)

    sTest_poutre_jansen = SensitivityConfidenceTest(model_poutre, distribution_poutre,
                                    JansenSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)

    sTest_poutre_mauntz = SensitivityConfidenceTest(model_poutre, distribution_poutre,
                                    MauntzKucherenkoSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)

    sTest_poutre_martinez = SensitivityConfidenceTest(model_poutre, distribution_poutre,
                                    MartinezSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)



    ################################################################################
    ###################            Aggregated Sobol            #####################
    ################################################################################


    model_aggregated = ot.SymbolicFunction(['X1', 'X2', 'X3'],
                                       ['2*X1 + X2 - 3*X3 + 0.3*X1*X2', 
                                        '-5*X1 + 4*X2 - 0.8*X2*X3 + 2*X3'])
    model_aggregated.setName("AggregatedSobol")
    distribution_aggregated = ot.ComposedDistribution([ot.Uniform()]*3)

    sTest_poutre_saltelli = SensitivityConfidenceTest(model_aggregated, distribution_aggregated,
                                    SaltelliSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)

    sTest_poutre_jansen = SensitivityConfidenceTest(model_aggregated, distribution_aggregated,
                                    JansenSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)

    sTest_poutre_mauntz = SensitivityConfidenceTest(model_aggregated, distribution_aggregated,
                                    MauntzKucherenkoSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)

    sTest_poutre_martinez = SensitivityConfidenceTest(model_aggregated, distribution_aggregated,
                                    MartinezSensitivityAlgorithm,
                                    savefig=True, plot_figure=True)
