import sys
import time

import numpy as np
import scipy.stats as stats
import sklearn.mixture as mixture


class MembershipModel:
    # inputdata is going to be NxM
    def calculateMembership(inputdata):
        return np.zeros((inputdata.shape[0]), np.bool)


class GMMMembershipModel(MembershipModel):
    # model id should be an iteratable (list/tuple/array) of ids
    def __init__(self, gmm, model_id):
        self.gmm = gmm
        self.model_id = model_id

    def calculateMembership(self, spikeset):
        inputdata = spikeset.featureByName('PCA').data
        labels = self.gmm.predict(inputdata)
        retval = np.zeros((inputdata.shape[0],), dtype=np.bool)
        for cid in self.model_id:
            retval = retval | (labels == cid)
        return retval

class PrecalculatedLabelsMembershipModel(MembershipModel):
    def __init__(self, labels, model_id):
        self.labels = labels
        self.model_id = model_id

    def calculateMembership(self, spikeset):
        retval = np.zeros((spikeset.N,), dtype=np.bool)
        if spikeset.spikes.shape[0] != self.labels.size:
            print 'size mismatch', spikeset.spikes.shape[0], self.labels.size
        else:
            for cid in self.model_id:
                retval = retval | (self.labels == cid)
        return retval


def fitGMMMembershipModel(spikes, gui_object):
    inputdata = spikes.featureByName('PCA').data
    # Fit a bunch of models
    bic = np.inf
    gmms = dict()

    Nbest = 0
    N = 0

    t1 = time.clock()
    print "Fitting GMMs"
    print 'Best model so far:',
    while True:
        N = N+1
        gmms = mixture.GMM(n_components=N, covariance_type='full')
        gmms.fit(inputdata)
        temp = gmms.bic(inputdata)
        if temp < bic:
            bic = temp
            Nbest = N
            gmm = gmms
            print N,
        else: # if bic has been increasing for 3 models then weve probably
              # passed the minimum
            print '('+str(N)+')',
            if Nbest < N-3:
                print ''
                break
        sys.stdout.flush()
    print "Settled on model with", Nbest, "components"
    N = Nbest

#    for n in xrange(num_comps.size):
#        print num_comps[n],
#        sys.stdout.flush()
#        gmms[n] = mixture.GMM(n_components=num_comps[n], covariance_type='full')
#        gmms[n].fit(inputdata)
#        bic[n] = gmms[n].bic(inputdata)
    t2 = time.clock()
    print "Took", (t2-t1), "seconds."

    # Pick the model with best BIC
#    idx = np.argmin(bic)
#    print 'BIC:', np.log10(bic)
#    print 'Best model:', num_comps[idx]

    # Now check if some clusters are so overlapping that they should
    # be merged - a future version of this procedure might involve looking
    # for joing refractory periods, etc. for now we use a crude mahal
    # threshold
#    N = num_comps[idx]
    p = spikes.featureByName('Peak').data
#    labels = gmms[idx].predict(inputdata)
    labels = gmm.predict(inputdata)

    # generate pairwise mahal distances
    dist = np.zeros((N,N))
    for i in xrange(N):
        for j in xrange(N):
            if i == j: continue

            gen_set = p[labels==i,:]
            cv = np.cov(gen_set.T)
            mu = np.mean(gen_set, axis=0)

            test_set = p[labels==j,:]

            x = test_set - mu
            temp = np.linalg.solve(cv.T, x.T) # equivalent to temp = x * inv(C)
            mhl = np.sum(x.T * temp, axis=0)
            dist[i,j] = np.mean(mhl)
    dist = 0.5 * (dist + dist.T)
    dist = stats.chi2.cdf(dist, p.shape[1])

    # now merge labels, this code works, but its yucky, should do it better
    threshold = 0.95
    cluster_labels = [[i,] for i in range(N)]
    for i in xrange(N):
        for j in xrange(i+1,N):
            if np.isnan(dist[i,j]): continue

            if dist[i,j] < threshold:
                print 'Merging components', (i,j)
                dist[j,:] = np.nan
                dist[:,j] = np.nan
                cluster_labels.remove([j,])
                for cluster in cluster_labels:
                    if i in cluster:
                        cluster.append(j)

    print "Creating", len(cluster_labels), "cluters from mixture model with",
    print N, "components"

    for label in cluster_labels:
        update = True
        cluster = gui_object.add_cluster()
        cluster.membership_model.append(GMMMembershipModel(gmm, label))
        cluster.calculateMembership(spikes)
        if cluster.stats['refr_frac'] > 0.01:
            print 'Deleting auto-cluster with > 1% refractory'
            gui_object.delete_cluster(cluster)
            update = False

    if update:
        gui_object.update_active_cluster()
    gui_object.updateFeaturePlot()




