/*
Copyright (c) 2012 Daniel Forsberg
dannne.forsberg@outlook.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "mex.h"
#include "math.h"

void normalize_distribution(double *distribution, int numberOfElements);

void compute_intensity_change(double *conditionalDistribution, double pixelValue,
        double min, double max, int numberOfChannels,
        double *intensityChange, int *certainty);

double entropy_diff(double *channelDistribution, double pixelValue,
        int numberOfChannels, double min, double max);

void entropy_channelize(double x, double min, double max,
        int numberOfChannels, int *firstChannel,
        double *channelWeight1,
        double *channelWeight2,
        double *channelWeight3, int order);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Input arguments
    double *fixedImage;
    double *movingImage;
    int numberOfChannels;
    double *mask = NULL;
    
    // Output arguments
    double *deltaC;
    double *deltaCMask;
    double MI;
    
    int numberOfPixels;
    int numberOfDimensionsFixedImage;
    int numberOfDimensionsMovingImage;
    const int *dimensionsFixedImage;
    const int *dimensionsMovingImage;
    
    //****************************************************
    // Check the input and output arguments.
    //****************************************************
    
    // Require at least three input arguments
    if (nrhs < 3) {
        mexErrMsgTxt("Too few arguments.");
    }
    
    // First we expect an array containing fixedImage image
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
            mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0])) {
        mexErrMsgTxt("fixedImage is expected to be an array.");
    }
    
    fixedImage = mxGetPr(prhs[0]);
    numberOfDimensionsFixedImage = mxGetNumberOfDimensions(prhs[0]);
    dimensionsFixedImage = mxGetDimensions(prhs[0]);
    numberOfPixels = 1;
    for (int k = 0; k < numberOfDimensionsFixedImage; k++) {
        numberOfPixels *= dimensionsFixedImage[k];
    }
    
    // Second we expect an array containing movingImage image
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
            mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1])) {
        mexErrMsgTxt("movingImage is expected to be an array.");
    }
    
    movingImage = mxGetPr(prhs[1]);
    numberOfDimensionsMovingImage = mxGetNumberOfDimensions(prhs[1]);
    if (numberOfDimensionsMovingImage != numberOfDimensionsFixedImage) {
        mexErrMsgTxt("fixedImage and movingImage must have the same size.");
    }
    
    dimensionsMovingImage = mxGetDimensions(prhs[1]);
    for (int k = 0; k < numberOfDimensionsFixedImage; k++) {
        if (dimensionsMovingImage[k] != dimensionsFixedImage[k]) {
            mexErrMsgTxt("fixedImage and movingImage must have the same size.");
        }
    }
    
    // Third we expect a scalar containing the number of channels for estimating MI
    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) ||
            mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]) ||
            mxGetNumberOfDimensions(prhs[2]) > 2 ||
            mxGetDimensions(prhs[2])[0] != 1 ||
            mxGetDimensions(prhs[2])[1] != 1) {
        mexErrMsgTxt("N is expected to be a scalar.");
    }
    
    numberOfChannels = (int) mxGetScalar(prhs[2]);
    if ((double) numberOfChannels != mxGetScalar(prhs[2])) {
        mexErrMsgTxt("N is expected to be an integer.");
    }
    
    if (numberOfChannels < 3) {
        mexErrMsgTxt("N must be at least 3.");
    }
    
    // Fourth we expect an array containing a mask
    if (nrhs > 3) {
        if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
                mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3])) {
            mexErrMsgTxt("mask is expected to be an array.");
        }
        
        mask = mxGetPr(prhs[3]);
        numberOfDimensionsMovingImage = mxGetNumberOfDimensions(prhs[3]);
        if (numberOfDimensionsMovingImage != numberOfDimensionsFixedImage) {
            mexErrMsgTxt("fixedImage and mask must have the same size.");
        }
        
        dimensionsMovingImage = mxGetDimensions(prhs[3]);
        for (int k = 0; k < numberOfDimensionsFixedImage; k++) {
            if (dimensionsMovingImage[k] != dimensionsFixedImage[k]) {
                mexErrMsgTxt("fixedImage and mask must have the same size.");
            }
        }
    }
    
    // Require at least one output argument
    if (nlhs < 1) {
        mexErrMsgTxt("Too few output arguments.");
    }
    
    plhs[0] = mxCreateNumericArray(numberOfDimensionsFixedImage,
            mxGetDimensions(prhs[0]),
            mxDOUBLE_CLASS, mxREAL);
    deltaC = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(numberOfDimensionsFixedImage,
            mxGetDimensions(prhs[0]),
            mxDOUBLE_CLASS, mxREAL);
    deltaCMask = mxGetPr(plhs[1]);
    
    // Get min and max of movingImage and fixedImage images
    double minFixedImage = fixedImage[0];
    double maxFixedImage = fixedImage[0];
    double minMovingImage = movingImage[0];
    double maxMovingImage = movingImage[0];
    
    for (int k = 1; k < numberOfPixels; k++) {
        if (fixedImage[k] < minFixedImage) {
            minFixedImage = fixedImage[k];
        }
        if (fixedImage[k] > maxFixedImage) {
            maxFixedImage = fixedImage[k];
        }
        if (movingImage[k] < minMovingImage) {
            minMovingImage = movingImage[k];
        }
        if (movingImage[k] > maxMovingImage) {
            maxMovingImage = movingImage[k];
        }
    }
    
    double *jointDistribution = (double *)mxCalloc(numberOfChannels * numberOfChannels, sizeof(double));
    double *distributionFixedImage = (double *)mxCalloc(numberOfChannels, sizeof(double));
    double *distributionMovingImage = (double *)mxCalloc(numberOfChannels, sizeof(double));
    
    // Update distributions of fixedImage and deformed images and also update joint distribution
    for (int k = 1; k < numberOfPixels; k++) {
        double fixedValue = fixedImage[k];
        double movingValue =  movingImage[k];
        
        if (mask && !mask[k]) {
            continue;
        }
        
        int firstChannelFixed, firstChannelMoving;
        double channelWeight1Fixed, channelWeight2Fixed, channelWeight3Fixed;
        double channelWeight1Moving, channelWeight2Moving, channelWeight3Moving;
        
        // Perform channel coding of deformed and fixedImage image values
        entropy_channelize(fixedValue,
                minFixedImage, maxFixedImage,
                numberOfChannels, &firstChannelFixed,
                &channelWeight1Fixed, &channelWeight2Fixed, &channelWeight3Fixed, 0);
        
        entropy_channelize(movingValue,
                minMovingImage, maxMovingImage,
                numberOfChannels, &firstChannelMoving,
                &channelWeight1Moving, &channelWeight2Moving, &channelWeight3Moving, 0);
        
        // Update distributions
        jointDistribution[ firstChannelFixed      * numberOfChannels + firstChannelMoving    ] += channelWeight1Fixed * channelWeight1Moving;
        jointDistribution[(firstChannelFixed + 1) * numberOfChannels + firstChannelMoving    ] += channelWeight2Fixed * channelWeight1Moving;
        jointDistribution[(firstChannelFixed + 2) * numberOfChannels + firstChannelMoving    ] += channelWeight3Fixed * channelWeight1Moving;
        jointDistribution[ firstChannelFixed      * numberOfChannels + firstChannelMoving + 1] += channelWeight1Fixed * channelWeight2Moving;
        jointDistribution[(firstChannelFixed + 1) * numberOfChannels + firstChannelMoving + 1] += channelWeight2Fixed * channelWeight2Moving;
        jointDistribution[(firstChannelFixed + 2) * numberOfChannels + firstChannelMoving + 1] += channelWeight3Fixed * channelWeight2Moving;
        jointDistribution[ firstChannelFixed      * numberOfChannels + firstChannelMoving + 2] += channelWeight1Fixed * channelWeight3Moving;
        jointDistribution[(firstChannelFixed + 1) * numberOfChannels + firstChannelMoving + 2] += channelWeight2Fixed * channelWeight3Moving;
        jointDistribution[(firstChannelFixed + 2) * numberOfChannels + firstChannelMoving + 2] += channelWeight3Fixed * channelWeight3Moving;
        distributionFixedImage[firstChannelFixed    ] += channelWeight1Fixed;
        distributionFixedImage[firstChannelFixed + 1] += channelWeight2Fixed;
        distributionFixedImage[firstChannelFixed + 2] += channelWeight3Fixed;
        distributionMovingImage[firstChannelMoving    ] += channelWeight1Moving;
        distributionMovingImage[firstChannelMoving + 1] += channelWeight2Moving;
        distributionMovingImage[firstChannelMoving + 2] += channelWeight3Moving;
    }
    
    // Normalize the distributions
    normalize_distribution(jointDistribution, numberOfChannels * numberOfChannels);
    normalize_distribution(distributionFixedImage, numberOfChannels);
    normalize_distribution(distributionMovingImage, numberOfChannels);
    
    double entropyFixedImage = 0.0;
    double entropyMovingImage = 0.0;
    double jointEntropy = 0.0;
    
    // Estimate entropy
    for (int i = 0; i < numberOfChannels; i++) {
        if (distributionFixedImage[i] > 0.0)
            entropyFixedImage -= distributionFixedImage[i] * log(distributionFixedImage[i]);
        if (distributionMovingImage[i] > 0.0)
            entropyMovingImage -= distributionMovingImage[i] * log(distributionMovingImage[i]);
        for (int j = 0; j < numberOfChannels; j++) {
            double x = jointDistribution[i * numberOfChannels + j];
            if (x > 0)
                jointEntropy -= x * log(x);
        }
    }
    
    MI = entropyFixedImage + entropyMovingImage - jointEntropy;
    
    double *conditionalDistribution = (double *)mxCalloc(numberOfChannels, sizeof(double));
    
    // Estimate delta C
    for (int k = 1; k < numberOfPixels; k++) {
        deltaC[k] = 0;
        deltaCMask[k] = 0;
        
        if (mask && !mask[k]) {
            continue;
        }
        
        double movingValue = movingImage[k];
        double fixedValue = fixedImage[k];
        
        int firstChannel;
        double channelWeight1, channelWeight2, channelWeight3;
        
        // Perform channel coding of deformed pixel value
        entropy_channelize(movingValue, minMovingImage, maxMovingImage,
                numberOfChannels, &firstChannel, &channelWeight1, &channelWeight2, &channelWeight3, 0);
        
        // Estimate the conditional distribution for the fixedImage image given the deformed pixel
        for (int channel = 0; channel < numberOfChannels; channel++) {
            conditionalDistribution[channel] =
                    channelWeight1 * jointDistribution[channel * numberOfChannels + firstChannel] +
                    channelWeight2 * jointDistribution[channel * numberOfChannels + firstChannel + 1] +
                    channelWeight3 * jointDistribution[channel * numberOfChannels + firstChannel + 2];
        }
        
        // Normalized the marginal distribution
        normalize_distribution(conditionalDistribution, numberOfChannels);
        
        double intensityChangeFixed;
        int certaintyFixed;
        
        // Estimate the intensity change
        compute_intensity_change(conditionalDistribution, fixedValue,
                minFixedImage, maxFixedImage,
                numberOfChannels, &intensityChangeFixed, &certaintyFixed);
        
        if (certaintyFixed) {
            deltaC[k] = -intensityChangeFixed;
            deltaCMask[k] = 1;
        }
    }
    
    // Return MI
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[2])[0] = MI;
    }
}

// Call this function to normalize the distribution with number of elements
void normalize_distribution(double *distribution, int numberOfElements) {
    double sum = 0.0;
    
    for (int k = 0; k < numberOfElements; k++) {
        sum += distribution[k];
    }
    
    if (sum > 0) {
        for (int k = 0; k < numberOfElements; k++) {
            distribution[k] /= sum;
        }
    }
}

// Call this function to compute the intensity change between the pixel value
// of the deformed image, with the given the conditional distribution, and
// the given pixel value of the fixed image
// conditionalDistribution of the deformed pixel value
// pixelValue corresponding fixed pixel value
// min of fixed image
// max of fixed image
// numberOfChannels
// intensityChange
// certainty
void compute_intensity_change(double *conditionalDistribution, double pixelValue,
        double min, double max, int numberOfChannels,
        double *intensityChange, int *certainty) {
    
    double entDiff = entropy_diff(conditionalDistribution, pixelValue, numberOfChannels, min, max);
    int channel = floor((pixelValue - min) / (max - min) * (numberOfChannels - 2));
    double delta = 0.0;
    double channelDistance = (max - min) / (numberOfChannels - 2);
    double lowLowChanVal, lowChanVal, upChanVal, upUpChanVal;
    double entDiffLowLowChanVal, entDiffLowChanVal, entDiffUpChanVal, entDiffUpUpChanVal;
    
    // Update channel if too close to the borders
    if (channel == numberOfChannels - 2)
        channel--;
    if (channel == -1)
        channel = 0;
    
    lowChanVal = min + channel * channelDistance;
    upChanVal = lowChanVal + channelDistance;
    entDiffLowChanVal = entropy_diff(conditionalDistribution, lowChanVal,
            numberOfChannels,
            min, max);
    entDiffUpChanVal = entropy_diff(conditionalDistribution, upChanVal,
            numberOfChannels,
            min, max);
    *certainty = 1;
    
    if (entDiffLowChanVal < 0 && entDiffUpChanVal > 0) {
        if (entDiff < 0)
            delta = (upChanVal - pixelValue) * (-entDiff) /
                    (entDiffUpChanVal - entDiff);
        else
            delta = -(pixelValue - lowChanVal) * entDiff /
                    (entDiff - entDiffLowChanVal);
    } else if (entDiffUpChanVal > entDiffLowChanVal) {
        if (entDiff < 0) {
            if (channel == numberOfChannels - 3) {
                delta = upChanVal - pixelValue;
            } else {
                upUpChanVal = pixelValue + channelDistance;
                entDiffUpUpChanVal = entropy_diff(conditionalDistribution, upUpChanVal, numberOfChannels, min, max);
                if (entDiffUpUpChanVal <= 0)
                    delta = channelDistance;
                else
                    delta = upChanVal - pixelValue + (upUpChanVal - upChanVal) *
                            (-entDiffUpChanVal) / (entDiffUpUpChanVal - entDiffUpChanVal);
            }
        } else {
            if (channel == 0) {
                delta = lowChanVal - pixelValue;
            } else {
                lowLowChanVal = pixelValue - channelDistance;
                entDiffLowLowChanVal = entropy_diff(conditionalDistribution, lowLowChanVal,
                        numberOfChannels, min, max);
                if (entDiffLowLowChanVal >= 0)
                    delta = -channelDistance;
                else
                    delta = lowChanVal - pixelValue + (lowLowChanVal - lowChanVal) *
                            entDiffLowChanVal / (entDiffLowChanVal - entDiffLowLowChanVal);
            }
        }
    } else {
        *certainty = 0;
    }
    
    *intensityChange = delta;
}

// Call this function to estimate the entropy difference
// channelDistribution
// pixelValue
// numberOfChannels
// min
// max
double entropy_diff(double *channelDistribution, double pixelValue,
        int numberOfChannels, double min, double max) {
    int firstChannel;
    double channelWeight1, channelWeight2, channelWeight3;
    
    // Estimate channel code for pixel value
    entropy_channelize(pixelValue, min, max, numberOfChannels, &firstChannel,
            &channelWeight1, &channelWeight2, &channelWeight3, 1);
    
    // Return entropy difference
    return -(channelWeight1 * (0 + log(channelDistribution[firstChannel] + 0.0001)) +
            channelWeight2 * (0 + log(channelDistribution[firstChannel + 1] + 0.0001)) +
            channelWeight3 * (0 + log(channelDistribution[firstChannel + 2] + 0.0001)));
}

// Call this function to perform channel coding on pixelValue
// pixelValue
// min
// max
// numberOfChannels
// firstChannel
// channelWeight1
// channelWeight2
// channelWeight3
// order of spline for channel weight approximation
void entropy_channelize(double pixelValue, double min, double max,
        int numberOfChannels, int *firstChannel,
        double *channelWeight1, double *channelWeight2, double *channelWeight3, int order) {
    double channel = (pixelValue - min) / (max - min) * (numberOfChannels - 2);
    double lowChannel = floor(channel);
    double channelDifference = channel - lowChannel;
    
    // Set first channel number
    *firstChannel = (int)lowChannel;
    
    // Update first channel and channel difference if to close to last channel
    if (*firstChannel == numberOfChannels - 2) {
        (*firstChannel)--;
        channelDifference += 1.0;
    }
    
    // Estimate channel weights
    if (order == 0) {
        // Approximate weights using a spline of order 2
        *channelWeight1 = 0.5 * (1 - channelDifference) * (1 - channelDifference);
        *channelWeight2 = 0.75 - (channelDifference - 0.5) * (channelDifference - 0.5);
        *channelWeight3 = 0.5 * channelDifference * channelDifference;
    }
    else {
        // Approximate weights using a spline of order 1
        *channelWeight1 = channelDifference - 1;
        *channelWeight2 = 1 - 2 * channelDifference;
        *channelWeight3 = channelDifference;
    }
}
