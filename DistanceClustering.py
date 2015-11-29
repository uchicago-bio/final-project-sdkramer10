import numpy as np

class DistanceClustering:
    
    def __init__(self, num_clusters, distance_function, distance_minimization_function):
        self.num_clusters = num_clusters
        self.distance_function = distance_function
        self.distance_minimization_function = distance_minimization_function
    
    def fit(self, data):

        clusters = np.zeros(len(data))
        centers = np.copy(data[0:self.num_clusters])

        nextPossibleStartingCenter = 1
        centersChosen = 1
        while centersChosen < self.num_clusters:
            centerAlreadyIncluded = False

            for i in range(centersChosen):
                centerAlreadyIncluded |= np.array_equal(centers[i], data[nextPossibleStartingCenter])

            if not centerAlreadyIncluded:
                centers[centersChosen] = data[nextPossibleStartingCenter]
                centersChosen += 1

            nextPossibleStartingCenter += 1

        while True:
            distances = self.distance_function(data, centers)
            newClusters = np.argmin(distances, axis=1)

            if np.array_equal(newClusters, clusters):
                return clusters, centers

            clusters = newClusters
            centers = np.array([self.distance_minimization_function(data[np.where(clusters == i)], axis=0) for i in range(self.num_clusters)])
