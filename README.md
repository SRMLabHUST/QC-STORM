# QC-STORM
Online processing plugin (ImageJ & Micro-Manager) for quality-controlled super-resolution localization imaging.

# Features
MLE Localization type: 2D, Astigmatism 3D

Online processing for sCMOS camera based imaging at full FOV and fastest frame rate (localization speed > 5x10^6 for 7x7 pixels ROI based on NVidia Titan xp graphics card)

Merging molecules emitting consecutively in adjacent frames (by localization precision weighted average)

Online super-resolution image rendering and statistical information analyzing (photon number, background intensity, SNR, PSF width, localization density, localization precision (CRLB for 2D), on time)

Online Nyquist resolution and convolved spatial resolution estimation

Background noise filtering

Drift correction by cross-correlation

Online feedback control

Spatial resolution controlled acquisition

Sequential Multi-FOV acquisition



# System requirements and installation
1, Windows 7 sp1 or newer, x64.

2, NVidia CUDA enabled GPU with compute capability no less than 3.5.

3, ImageJ/Fiji, Micro-Manager 2.0 (beta 3).

4, Download and install Microsoft Visual C++ 2015 Redistributable Update 3.

5, Download QC-STORM release version from https://github.com/SRMLabHUST/QC-STORM/releases.

6, Copy .dll files into installation directory of ImageJ or Micro-Manager, and .jar files into folder plugins or mmplugins for ImageJ and Micro-Manager respectively.

# Recompile the source codes
1, The Java GUI for ImageJ and Micro-Manager are develped by NetBeans IDE 7.3.1 with Java jdk1.6.0_45.

2, The DLL wrapped image processing are developed by Visual Studio 2015 end accelerated by CUDA 9.0.



# Declaration
This program is free software: you can redistribute it and/or modify it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU LESSER GENERAL PUBLIC LICENSE for more details.

You should have received a copy of the GNU LESSER GENERAL PUBLIC LICENSE along with this program.  If not, see <https://www.gnu.org/licenses/>.
