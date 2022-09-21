#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
function [DS_GRACE_fields,RMSdiff] = statistical_downscaling_grids_TWS(nGRACE, nPrec, nET, nRunoff, model_TWS, k, time_vec_M, basinMASK, tvec_ds,nm)
%UNTITLED3 Summary of this function goes here

% This function uses MA-PLR to downscale the GRACE catchment scale products
%% Input:
% nGRACE: number of GRACE products
% nPrec: number of precipitation products
% nET: number of Evapotranspiration products
% nRunnoff: number of Runoff products
% model_TWS: modelled estimate of global TWS change from a model (mm), such as WGHM. It should be in a cell format and the model should be at a resoluion of 0.5 degree grid resolution.
% k: number of time shifts
% time_vec_M: time vector of model_TWS
% basinMASK: mask of the catchment we are interested in (360 X 720 matrix) or a cell structure with 1 column and each row with a 360 X 720 field (1 inside and 0 outside the catchment).
% tvec_ds: time period of the output, it should be less than the length of GRACE time series and the length of model time series. See paper.
% nm : number of modes to use for reconstruction
%% Output
% DS_GRACE: a cell fromat with downscaled GRACE product for each time point in tvec_assm. It is a cell fromat with each cell containing a downscaled field.
% RMSdiff: RMS of the difference between GRACE time series and the downscaled catchment average
%% Required functions:
% naninterp
% plr
% inputDATA
%% Author
% Dr. Bramha Dutt Vishwakarma
% School of Geographical sciences, University of Bristol, BS8 1SS, Bristol, UK.
% Date: 25/02/2020
%% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.

%    see <http://www.gnu.org/licenses/>
% -------------------------------------------------------------------------


%% here are the details of step by step processing
% Prepare datasets (pre-processing)
"""

#GRID AREA CALCULATION
x = np.linspace(0, 359.5, 720);
y = np.linspace(0, 179.5, 360);
x1 = np.linspace(0.5, 360, 720);
y1 = np.linspace(0.5, 180, 360);
lambdd,theta = np.meshgrid(x,y)  
lambdd1,theta1 = np.meshgrid(x1,y1)  
a = np.sin((90-theta)/180*np.pi)-np.sin((90-theta1)/180*np.pi)
b = lambdd1 - lambdd
area = (6378.137**2)*pow(10,6)*np.pi/180*(np.multiply(a,b)) 