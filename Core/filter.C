// provide function to perform filtering, all in public name space
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "ifm_common.h"
#include "engine.h"

/* 
 generate a 1D gaussian filter with given sigma
 */
int64_t gen_gaussian1d(int64_t sz, real_t sigma, real_t *space)
{
    int64_t jj, cntr;
    real_t sig2, sum, tmp;

    if (sz <= 0)
        return -1;
    sig2 = 2.0*sigma*sigma;
    sum = 0.0;
    cntr = 0;
    for (jj=-sz; jj<=sz; jj++) {
        tmp = exp( - ( (real_t)(jj*jj) )/sig2);
        sum += tmp;
        space[cntr++] = tmp;
    }
    for (jj=0; jj<cntr; jj++)
    {  space[jj] /= sum;
    }
     
    return cntr;
}

// convolution 1D in X direction
void convolutionx1d(real_t *target, real_t *source, real_t *filter, real_t nnd, int64_t filter_sz, uint64_t num_row, uint64_t num_col)
{
	int64_t offset, jj, ii, kk;
	real_t tmp, sum;

	offset = (filter_sz-1)/2;  /* filter size has to be odd! */
	for (jj=0; jj<num_row; jj++) {
		/* first few items , partially patched with zero */
		for (kk=0; kk<offset; kk++) {
			tmp = 0.0;
			sum = 0.0;
			for (ii=0; ii<filter_sz; ii++) {
				if ( kk + ii - offset >= 0 ) {
					if ( !(FLOAT_EQ( source[kk-offset+ii+num_col*jj], nnd))) {
						tmp += filter[ii] * source[ kk - offset + ii + num_col*jj];
						sum += filter[ii];
					}
				}
			}
			if ( !(FLOAT_EQ( source[kk+num_col*jj], nnd)) ) {
				if (sum>0)  target[ kk + num_col*jj] = tmp/sum;
			}
		}

		/* the middle trunk, */
		for (kk=offset; kk<num_col-offset; kk++) {
			tmp = 0.0;
			sum = 0.0;
			for (ii=0; ii<filter_sz; ii++) {

				if ( !(FLOAT_EQ( source[kk-offset+ii+num_col*jj],nnd))) {
					tmp += filter[ii] * source[ kk - offset + ii + num_col*jj];
					sum += filter[ii];
				}
			}
			if ( !(FLOAT_EQ( source[kk + num_col*jj],nnd)) ) {
				if (sum>0) target[ kk + num_col*jj] = tmp/sum;
			}
		}

		/* the last few of them, again padded with zero  */
		for (kk=num_col-offset; kk<num_col; kk++) {
			tmp = 0.0;
			sum = 0.0;
			for (ii=0; ii<filter_sz; ii++) {
				if ( kk + ii - offset < num_col ) { 
					if ( !(FLOAT_EQ(source[kk-offset+ii+num_col*jj], nnd))) {
						tmp += filter[ii] * source[ kk - offset + ii + num_col*jj];
						sum += filter[ii];
					}
				}
			}
			if ( !(FLOAT_EQ( source[ kk+num_col*jj], nnd))) {
				if (sum>0)  target[ kk + num_col*jj] = tmp/sum;
			}
		}


	}

}

/* Ditto, but in Y-direction */
void convolutiony1d(real_t *target, real_t *source, real_t *filter, real_t nnd, int64_t filter_sz, int64_t num_row, int64_t num_col)
{
	int64_t offset, jj, ii, kk;
	real_t tmp, sum;

	offset = (filter_sz-1)/2;  /* filter size has to be odd! */

	for (jj=0; jj<num_col; jj++) {
		/* first few items , partially zero */
		for (kk=0; kk<offset; kk++) {
			tmp = 0.0;
			sum = 0.0;
			for (ii=0; ii<filter_sz; ii++) {
				if ( kk + ii - offset >= 0 ) { 
					if ( !(FLOAT_EQ(source[ (kk-offset + ii)*num_col+jj], nnd)) ) {
						tmp += filter[ii] * source[ (kk - offset + ii)*num_col + jj];
						sum += filter[ii];
					}
				}
			}
			if ( !(FLOAT_EQ( source[kk*num_col + jj], nnd))) {
				if ( sum>0) target[ kk*num_col + jj] = tmp/sum;
			}
		}

		/* the middle */
		for (kk=offset; kk<num_row-offset; kk++) {
			tmp = 0.0;
			sum = 0.0;
			for (ii=0; ii<filter_sz; ii++) {
				if ( !(FLOAT_EQ(source[ (kk-offset+ii)*num_col + jj], nnd))) {
					tmp += filter[ii] * source[ (kk - offset + ii)*num_col + jj];
					sum += filter[ii];
				}
			}
			if ( !(FLOAT_EQ( source[kk*num_col + jj], nnd))) {
				if ( sum>0) target[ kk*num_col + jj] = tmp/sum;
			}
		}

		/* the last few of them  */

		for (kk=num_row-offset; kk<num_row; kk++) {
			tmp = 0.0;
			sum = 0.0;
			for (ii=0; ii<filter_sz; ii++) {
				if ( kk + ii - offset < num_row ) { 
					if ( !(FLOAT_EQ(source[ (kk-offset+ii)*num_col+jj], nnd))) {
						tmp += filter[ii] * source[ (kk - offset + ii)*num_col + jj];
						sum += filter[ii];
					}
				}
			}
			if ( !(FLOAT_EQ( source[kk*num_col + jj], nnd))) {
				if (sum>0) target[ kk*num_col + jj] = tmp/sum;
			}
		}
	}
}

// perform 2D Gaussian convolution, overwrite the original
int gaussian_filter(Grid *grid, int64_t sz, real_t sigma) {
	real_t *tmp1, *tmp2, *filter;
	int64_t fsz;

	assert(sz>0);   // at least 1
	assert(grid);

	real_t *dem = grid->data;
	real_t nnd = grid->nodata;
	uint64_t row = grid->rows;
	uint64_t col = grid->cols;

	fsz = 2*sz+1;  // got to be odd 

	filter = (real_t*)malloc(sizeof(real_t)*fsz);
	tmp1 = (real_t*)malloc(sizeof(real_t)*col*row);
	tmp2 = (real_t*)malloc(sizeof(real_t)*col*row);

	memset(tmp1, 0, sizeof(real_t)*row*col);
	memcpy(tmp2, dem, sizeof(real_t)*row*col); // we have the original

	if (gen_gaussian1d(sz, sigma, filter) < 0) {
		printf("Bummer: error in pre-processing.\n");
		return -1;
	}

	convolutionx1d(tmp1, dem,  filter, nnd, fsz, row, col);

        convolutiony1d(tmp2, tmp1, filter, nnd, fsz, row, col);
        
	memcpy(dem, tmp2,sizeof(real_t)*row*col);  // we overwrite the original 

	if (tmp2) free (tmp2);
	if (tmp1) free (tmp1);
	if (filter) free (filter);

	return 0;
}

// fill the NODATA pixels with fixed values D meters below the lowest four neighbors
// (default M=10m)
// modifies dem and mask files
//
int fix_boundaries(real_t *dem, short *mask, uint64_t row, uint64_t col, real_t nodata_tag, real_t offset)
{
	real_t *tmp_dem;
	real_t buffer[4];
	real_t junk;
	uint64_t cntr;
	uint64_t cur, lft, rgt, bot, top;
	uint64_t ii, jj, kk;
	uint64_t nsink=0;

	assert(offset<0);

	tmp_dem = (real_t*)malloc(sizeof(real_t)*row*col);
	memcpy(tmp_dem, dem, sizeof(real_t)*row*col);

	///* take care of the four corners */
	//if (FLOAT_EQ( dem[0], nodata_tag)) {
	//	cntr = 0;
	//	if ( !(FLOAT_EQ(dem[1],nodata_tag)) ) buffer[cntr++] = dem[1];
	//	if ( !(FLOAT_EQ(dem[col],nodata_tag)) ) buffer[cntr++] = dem[col];
	//	if ( cntr>0) {
	//		junk = buffer[0];
	//		for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//		tmp_dem[0] = junk + offset;
	//		mask[0] = -1;
	//		nsink++;
	//	}
	//}

	//if (FLOAT_EQ( dem[col-1], nodata_tag)) {
	//	cntr = 0;
	//	if ( !(FLOAT_EQ(dem[col-2],nodata_tag)) ) buffer[cntr++] = dem[col-2];
	//	if ( !(FLOAT_EQ(dem[2*col-1],nodata_tag)) ) buffer[cntr++] = dem[2*col-1];
	//	if ( cntr>0) {
	//		junk = buffer[0];
	//		for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//		tmp_dem[col-1] = junk + offset;
	//		mask[col-1] = -1;
	//		nsink++;

	//	}
	//}

	//if (FLOAT_EQ( dem[(row-1)*col], nodata_tag) ) {
	//	cntr = 0;
	//	if ( !(FLOAT_EQ(dem[(row-1)*col+1],nodata_tag)) ) buffer[cntr++] = dem[(row-1)*col+1];
	//	if ( !(FLOAT_EQ(dem[(row-2)*col],nodata_tag)) ) buffer[cntr++] = dem[(row-2)*col];
	//	if ( cntr>0) {
	//		junk = buffer[0];
	//		for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//		tmp_dem[(row-1)*col] = junk + offset;
	//		mask[(row-1)*col] = -1;
	//		nsink++;

	//	}
	//}

	//if (FLOAT_EQ( dem[ row*col-1], nodata_tag)) {
	//	cntr = 0;
	//	if ( !(FLOAT_EQ(dem[(row-1)*col-1],nodata_tag)) ) buffer[cntr++] = dem[(row-1)*col-1];
	//	if ( !(FLOAT_EQ(dem[row*col-2],nodata_tag)) ) buffer[cntr++] = dem[row*col-2];
	//	if ( cntr>0) {
	//		junk = buffer[0];
	//		for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//		tmp_dem[row*col-1] = junk + offset;
	//		mask[row*col-1] = -1;
	//		nsink++;

	//	}
	//}

	///* the first row */
	//for (ii=1; ii<col-1; ii++) {
	//	if ( FLOAT_EQ( dem[ii], nodata_tag) ) {
	//		cntr=0;
	//		if ( !(FLOAT_EQ(dem[ii-1], nodata_tag))   ) buffer[cntr++]=dem[ii-1];
	//		if ( !(FLOAT_EQ(dem[ii+1], nodata_tag))   ) buffer[cntr++]=dem[ii+1];
	//		if ( !(FLOAT_EQ(dem[col+ii], nodata_tag)) ) buffer[cntr++]=dem[row+ii];

	//		if (cntr>0) {
	//			junk = buffer[0];
	//			for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//			tmp_dem[ii] = junk + offset;
	//			mask[ii] = -1;
	//			nsink++;

	//		}	
	//	}
	//}

	///* the last row */
	//for (ii=1; ii<col-1; ii++) {
	//	if ( FLOAT_EQ( dem[ (row-1)*col + ii], nodata_tag) ) {
	//		cntr = 0;
	//		if ( !(FLOAT_EQ(dem[(row-1)*col+ii-1], nodata_tag)) ) buffer[cntr++]=dem[(row-1)*col+ii-1];
	//		if ( !(FLOAT_EQ(dem[(row-1)*col+ii+1], nodata_tag)) ) buffer[cntr++]=dem[(row-1)*col+ii+1];
	//		if ( !(FLOAT_EQ(dem[(row-2)*col+ii],   nodata_tag)) ) buffer[cntr++]=dem[(row-2)*col+ii];

	//		if (cntr>0) {
	//			junk = buffer[0];
	//			for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//			tmp_dem[(row-1)*col+ii] = junk + offset;
	//			mask[(row-1)*col+ii] = -1;
	//			nsink++;

	//		}	
	//	}
	//}

	///* the first column */
	//for (jj=1; jj<row-1; jj++) {
	//	if ( FLOAT_EQ( dem[jj*col], nodata_tag) ) {
	//		cntr = 0;
	//		if ( !(FLOAT_EQ(dem[(jj-1)*col], nodata_tag)) ) buffer[cntr++]=dem[(jj-1)*col];
	//		if ( !(FLOAT_EQ(dem[(jj+1)*col], nodata_tag)) ) buffer[cntr++]=dem[(jj+1)*col];
	//		if ( !(FLOAT_EQ(dem[jj*col+1],  nodata_tag)) ) buffer[cntr++]=dem[jj*col+1];
	//		if (cntr>0) {
	//			junk = buffer[0];
	//			for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//			tmp_dem[jj*col] = junk + offset;
	//			mask[jj*col] = -1;
	//			nsink++;

	//		}	
	//	}
	//}

	///* the last column */
	//for (jj=1; jj<row-1; jj++) {
	//	if ( FLOAT_EQ( dem[(jj+1)*col-1], nodata_tag) ) {
	//		cntr = 0;
	//		if ( !(FLOAT_EQ(dem[(jj)*col-1],   nodata_tag)) ) buffer[cntr++]=dem[(jj-1)*col-1];
	//		if ( !(FLOAT_EQ(dem[(jj+2)*col-1], nodata_tag)) ) buffer[cntr++]=dem[(jj+2)*col-1];
	//		if ( !(FLOAT_EQ(dem[(jj+1)*col-2], nodata_tag)) ) buffer[cntr++]=dem[(jj+1)*col-2];
	//		if (cntr>0) {
	//			junk = buffer[0];
	//			for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//			tmp_dem[(jj+1)*col-1] = junk + offset;
	//			mask[(jj+1)*col-1] = -1;
	//			nsink++;

	//		}	
	//	}
	//}

	///* the rest */
	//for (jj=1; jj<row-1; jj++) {
	//	for (ii=1; ii<col-1; ii++) {
	//		if ( FLOAT_EQ( dem[ jj*col+ii ], nodata_tag) ) {
	//			cntr = 0;
	//			if ( !(FLOAT_EQ(dem[(jj-1)*col+ii], nodata_tag)) ) buffer[cntr++]=dem[(jj-1)*col+ii];
	//			if ( !(FLOAT_EQ(dem[(jj+1)*col+ii], nodata_tag)) ) buffer[cntr++]=dem[(jj+1)*col+ii];
	//			if ( !(FLOAT_EQ(dem[jj*col+ii+1],   nodata_tag)) ) buffer[cntr++]=dem[jj*col+ii+1];
	//			if ( !(FLOAT_EQ(dem[jj*col+ii-1],   nodata_tag)) ) buffer[cntr++]=dem[jj*col+ii-1];
	//			if (cntr>0) {
	//				junk = buffer[0];
	//				for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
	//				tmp_dem[jj*col+ii] = junk + offset;
	//				mask[jj*col+ii] = -1;
	//				nsink++;
	//			}
	//		}
	//	}
	//}


	//optimize
    for (jj=0; jj<row; jj++) {
		for (ii=0; ii<col; ii++) {
			cur = jj*col+ii;
            if ( FLOAT_EQ( dem[cur], nodata_tag) ) {
				cntr = 0;
				if(ii>0) {
					lft = cur-1;
					if (!(FLOAT_EQ(dem[lft], nodata_tag)) ) buffer[cntr++]=dem[lft];
				}

				if(ii<col-1) {
					rgt = cur+1;
					if (!(FLOAT_EQ(dem[rgt], nodata_tag)) ) buffer[cntr++]=dem[rgt];
				}

				if(jj>0) {
					bot = cur-col;
					if (!(FLOAT_EQ(dem[bot], nodata_tag)) ) buffer[cntr++]=dem[bot];
				}

				if(jj<row-1) {
					top = cur+col;
					if (!(FLOAT_EQ(dem[top], nodata_tag)) ) buffer[cntr++]=dem[top];
				}


				if (cntr>0) {
				  junk = buffer[0];
				  for (kk=1; kk<cntr; kk++) junk = MIN2(junk, buffer[kk]);
				  tmp_dem[cur] = junk + offset;
				  mask[cur] = -1;
				  nsink++;
				}
			}
		}
	}

	memcpy(dem, tmp_dem, sizeof(real_t)*row*col);
	if (tmp_dem) free (tmp_dem);

	printf("Found %ld water surface points.\n", nsink);
	return 0;
}

/* end */
