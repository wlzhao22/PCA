# PCA
Project of Principle Content Analysis in C++. Given an input matrix 'mat' in 'nRow' rows and 'nCol' columns,  if nRow >= nCol the PCA mapping matrix is trained with covariance matrix otherwise when nRow < nCol, the PCA mapping matrix is trained by Gram method.

# Required library 
### lapack
### blas
Under Linux, one can install them as follows

```
sudo apt-get install libblas-dev liblapack-dev
```

# Compile
``` 
cd PCA/
make release
```

# Command
```
pcatrn -pc mat -d dstfn
```

A sample input matrix 'sift4k.txt' is found from 'data/'. Since we have to calculate the covariance matrix during the training, one should not provide a matrix in many rows. The training could be run out of memory. The number rows should be around 20 times of vector dimension. 

# Author
Wan-Lei Zhao
