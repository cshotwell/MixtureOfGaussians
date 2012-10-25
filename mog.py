# Carl Shotwell, Mixture of Gaussians EM
# To run, type ~ python mog.py "infile" in the command line 

import math
import random
import sys 


class MOG():
    """Mixture of Gaussians learning class."""
    def __init__(self, mixes, data_in):
    """Constructor accepts number of gaussians, and 1d list of data."""
        self.data = data_in
        self.gaussians = []
        self.gaussians.append(self.__initMeanStdDev())
        #stores stores the remaining gaussians with random changes.
        #[mean, stddev]
        for i in range(1,mixes):
            self.gaussians.append([self.gaussians[0][0]+random.uniform(-1,1),self.gaussians[0][1]+random.uniform(-1,1)])
               
    #cleans up some of the __init__ function
    def __initMeanStdDev(self):
        #calc mean
        mean = 0.0
        for i in self.data:
            mean += float(i)
        mean = mean/len(self.data)
        #calc stddev
        variance = 0
        for i in self.data:
            variance += (float(i) - mean)**2
        stddev = math.sqrt(variance/float(len(self.data)))
        return [mean, stddev]

    def pprint(self):
    """Pretty print the mean and std_dev found for each gaussian"""
        for i in range(len(self.gaussians)):
            print "Class " + str(i+1)
            print "Parameters learned"
            print "Mean:   " + self.truncate(self.gaussians[i][0], 4)
            print "Stddev: " + self.truncate(self.gaussians[i][1], 4)
            print " "

    
    def truncate(self,n,d):
        temp = str(n).split('.')
        return ".".join([temp[0], temp[1][:d]])    

    def run(self,iterations):
    """Run the algorithm for the specified number of iterations."""
        for i in range(iterations):
            self.step()

    def step(self):
    """Perform one complete EM step."""
        expectations = self.E_Step()
        self.M_Step(expectations)

    def E_Step(self):
    """Perform the E step, return expectations for each distrubution."""
        # nested list representing expectations for each distribution.
        expectations = []
        for x in self.data:
            temp = []
            for parameters in self.gaussians:
                temp.append(self.Gaussian(x, parameters))
            expectations.append(self.normalize(temp))   
        return expectations
    
    def M_Step(self, expectations):
    """Perform the M step on a list of expectations."""
        count = 0
        new_values = []
        #loop for gaussian
        for i in range(len(self.gaussians)):
            n_mean = 0
            n_stddev = 0
            Ni = 0
            for x in range(len(self.data)):
                n_mean += expectations[x][count] * self.data[x]
                n_stddev += expectations[x][count] * (self.data[x]**2) 
                Ni += expectations[x][count]
            n_stddev = n_stddev/Ni
            n_mean = n_mean/Ni

            n_stddev = math.sqrt(math.fabs(n_stddev - n_mean**2))
            self.gaussians[count] = [n_mean, n_stddev]
            count += 1 
        
    def normalize(self, vector):
    """Normalize a vector."""
        scalar = 1.0 / sum(vector)
        return [scalar * i for i in vector]


    #compute the gaussian
    def Gaussian(self,x, parameters):
    """Compute the value of a gaussian PDF given x, parameters."""
        (mean, stddev) = parameters 
        return (1/(stddev * math.sqrt(2*math.pi))) * math.e**(-.5 * ((x-mean)/stddev)**2)


if __name__ == "__main__":
    #data collection step
    data_in = [float(i.strip()) for i in open(sys.argv[1])]    
    #Constructor takes number of gaussians and the data set
    N = 3
    iterations = 100

    EM = MOG(N,data_in)
    EM.run(iterations)
    EM.pprint()
