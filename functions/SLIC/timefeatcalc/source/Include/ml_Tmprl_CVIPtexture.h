/***************************************************************************
* ======================================================================
* Computer Vision/Image Processing Tool Project - Dr. Scott Umbaugh SIUE
* ======================================================================
*
*             File Name: CVIPtexture.h
*           Description: contains function prototypes, type names, constants,
*			 etc. related to libdataserv (Data Services Toolkit.)
*         Related Files: Imakefile, cvip_pgmtexture.c
*   Initial Coding Date: 6/19/96
*           Portability: Standard (ANSI) C
*             Credit(s): Steve Costello
*                        Southern Illinois University @ Edwardsville
*
** Copyright (C) 1993 SIUE - by Gregory Hance.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
**
****************************************************************************/
/* Yanhua Hu: modify the struct and decrease elements to 1 in each array  June 25, 03*/

typedef struct  {
/* [0] -> 0 degree, [1] -> 45 degree, [2] -> 90 degree, [3] -> 135 degree,
   [4] -> average, [5] -> range (max - min)  yh: no degree*/
	float ASM[1];		/*  (1) Angular Second Moment */
	float contrast[1];	/*  (2) Contrast */
	float correlation[1];	/*  (3) Correlation */
	float variance[1];	/*  (4) Variance */
	float IDM[1];		/*  (5) Inverse Diffenence Moment */
	float sum_avg[1];	/*  (6) Sum Average */
	float sum_var[1];	/*  (7) Sum Variance */
	float sum_entropy[1];	/*  (8) Sum Entropy */
	float entropy[1];	/*  (9) Entropy */
	float diff_var[1];	/* (10) Difference Variance */
	float diff_entropy[1];	/* (11) Diffenence Entropy */
	float meas_corr1[1];	/* (12) Measure of Correlation 1 */
	float meas_corr2[1];	/* (13) Measure of Correlation 2 */
	float max_corr_coef[1]; /* (14) Maximal Correlation Coefficient */
	} TEXTURE;


typedef struct {
/* Allows the user to choose which features to extract, a zero will cause 
   the feature to be ignored, the returned feature value will be 0.0 */
	int ASM;		/*  (1) Angular Second Moment */
	int contrast;		/*  (2) Contrast */
	int correlation;	/*  (3) Correlation */
	int variance;		/*  (4) Variance */
	int IDM;		/*  (5) Inverse Diffenence Moment */
	int sum_avg;		/*  (6) Sum Average */
	int sum_var;		/*  (7) Sum Variance */
	int sum_entropy;	/*  (8) Sum Entropy */
	int entropy;		/*  (9) Entropy */
	int diff_var;		/* (10) Difference Variance */
	int diff_entropy;	/* (11) Diffenence Entropy */
	int meas_corr1;		/* (12) Measure of Correlation 1 */
	int meas_corr2;		/* (13) Measure of Correlation 2 */
	int max_corr_coef; 	/* (14) Maximal Correlation Coefficient */
	} TEXTURE_FEATURE_MAP;

