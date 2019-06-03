# MLEwt
We applied Gaussian function weights to commonly used maximum likelihood estimation (MLE) based single-molecule localization algorithm for both 2D and astigmatism 3D imaging (abbreviated as MLEwt).

We found that the MLEwt is able to suppress the signal contamination from neighboring molecules, thus improve localization precision and improve tolerable activation density compared with MLE. 


# Results of three ROIs (cropped from experimental data set)

![](https://github.com/SRMLabHUST/QC-STORM/blob/master/MLEwt%20MATLAB/image%20for%20display/figure%201s%20results%20on%20three%20ROI%202D.png)


![](https://github.com/SRMLabHUST/QC-STORM/blob/master/MLEwt%20MATLAB/image%20for%20display/figure%201s%20results%20on%20three%20ROI%203D.png)




# How to use
Download the codes and open MATLAB software (the author run the codes by MATLAB R2017a).

Run Single Molecule Fitting 2D/BFGS_2D_main.m or Single Molecule Fitting AS3D/BFGS_3D_main.m

The codes will load a test data and give the fitting results.

You can switch between MLEwt and MLE to watch the difference, and can replace the test data by your own data.



# Declaration
This program is free software: you can redistribute it and/or modify it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU LESSER GENERAL PUBLIC LICENSE for more details.

You should have received a copy of the GNU LESSER GENERAL PUBLIC LICENSE along with this program.  If not, see <https://www.gnu.org/licenses/>.

For more questions, please contact Prof. Zhenli-Huang at "leo@mail.hust.edu.cn" or author Luchang Li at "luchangli1993@163.com".
