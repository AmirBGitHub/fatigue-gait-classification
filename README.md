# fatigue-gait-classification

# Guideline

download the package and run Read_and_Classify_CrossValSubj.m file.

This repository contains the data and code for reproducing the results from:

Amir Baghdadi, Fadel M. Megahed, Ehsan T. Esfahani & Lora A. Cavuoto (2018) A machine learning approach to detect changes in gait parameters following a fatiguing occupational task, Ergonomics, 61:8, 1116-1129, DOI: 10.1080/00140139.2018.1442936

The purpose of this study is to provide a method for classifying non-fatigued vs. fatigued states following manual material handling. A method of template matching pattern recognition for feature extraction ($1 Recognizer) along with the support vector machine model for classification were applied on the kinematics of gait cycles segmented by our stepwise search-based segmentation algorithm. A single inertial measurement unit on the ankle was used, providing a minimally intrusive and inexpensive tool for monitoring. The classifier distinguished between states using distance based scores from the recogniser and the step duration. The results of fatigue detection showed an accuracy of 90% across data from 20 recruited subjects. This method utilises the minimum amount of data and features from only one low-cost sensor to reliably classify the state of fatigue induced by a realistic manufacturing task using a simple machine learning algorithm that can be extended to real-time fatigue monitoring as a future technology to be employed in the manufacturing facilities.
