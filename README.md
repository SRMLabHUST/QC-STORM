# QC-STORM
QC-STORM is an online processing plugin (Micro-Manager & ImageJ) for quality-controlled super-resolution localization imaging.

# Key features

Ultrahigh efficient MLE localization without sacrificing localization precition: online processing for sCMOS camera based localization microscopy at full FOV and fastest frame rate (localization speed > 5x10^6 for 7x7 pixels ROI based on NVidia Titan xp graphics card). 

Most processing parts of QC-STORM are GPU accelerated, and their performance increase almost linearly with GPU’s single floating performance (GFLOPS). Note the speed may be reduced by CPU memory bandwidth and low speed image reading from hard disk.

The localization performance evaluation are available at previous work of Sage, Daniel, et al. "Super-resolution fight club: A broad assessment of 2D & 3D single-molecule localization microscopy software." bioRxiv (2018): 362517.

MLE Localization type: 2D, Astigmatism 3D

Online super-resolution image rendering and statistical information analyzing (photon number, background intensity, SNR, PSF width, localization density, on time)

Online localization precision (CRLB) calculation

Online Nyquist resolution and convolved spatial resolution calculation

Online spatial resolution guided acquisition number

Drift correction by cross-correlation (post-processing)

Merging molecules emitting in consecutive frames (localization precision weighted average)

Online feedback control for localization density stabilizing and axial drift correction (hardware dependent)

Sequential Multi-FOV acquisition (hardware dependent)


# System requirements and installation
1, Windows 7 sp1 or newer, x64.

2, NVidia CUDA enabled GPU with compute capability no less than 3.5.

3, ImageJ/Fiji, Micro-Manager 2.0 (beta 3).

4, Download and install Microsoft Visual C++ 2015 Redistributable Update 3 (x64) (https://www.microsoft.com/en-us/download/details.aspx?id=53587).

5, Download QC-STORM release version from https://github.com/SRMLabHUST/QC-STORM/releases.

6, Copy the downloaded .dll files into installation directory of ImageJ or Micro-Manager, and .jar files into "plugins" or "mmplugins" folder for ImageJ and Micro-Manager respectively.


# Recompile the source codes
1, The Java GUI for ImageJ (v1.52g) and Micro-Manager (MMSetup_64bit_2.0.0-beta3 from nightly builds of Micro-Manager 2.0) are develped by NetBeans IDE 7.3.1 with Java jdk1.6.0_45.

2, The C++ core image processing algorithms are developed by Visual Studio 2015 Update 3, and accelerated by NVidia CUDA 9.0.

3, Download the codes, open the projects by corresponding software and recompile the codes.


# Declaration
This program is free software: you can redistribute it and/or modify it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU LESSER GENERAL PUBLIC LICENSE for more details.

You should have received a copy of the GNU LESSER GENERAL PUBLIC LICENSE along with this program.  If not, see <https://www.gnu.org/licenses/>.

For more questions, pleast contact Prof. Zhenli-Huang at "leo@mail.hust.edu.cn" or author Luchang Li at "luchangli1993@163.com".
