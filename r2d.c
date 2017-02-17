/*
 *
 *		r2d.c
 *		
 *		Devon Powell
 *		31 August 2015
 *		
 *		See r2d.h for usage.
 *
 *		This program was prepared by Los Alamos National Security, LLC at Los Alamos National
 *		Laboratory (LANL) under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). 
 *		All rights in the program are reserved by the DOE and Los Alamos National Security, LLC.  
 *		Permission is granted to the public to copy and use this software without charge, provided that 
 *		this Notice and any statement of authorship are reproduced on all copies.  Neither the U.S. 
 *		Government nor LANS makes any warranty, express or implied, or assumes any liability 
 *		or responsibility for the use of this software.
 *
 */

#include "r2d.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// useful macros
#define ONE_THIRD 0.333333333333333333333333333333333333333333333333333333
#define ONE_SIXTH 0.16666666666666666666666666666666666666666666666666666667
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)
#define wav(va, wa, vb, wb, vr) {			\
	vr.x = (wa*va.x + wb*vb.x)/(wa + wb);	\
	vr.y = (wa*va.y + wb*vb.y)/(wa + wb);	\
}
#define norm(v) {					\
	r2d_real tmplen = sqrt(dot(v, v));	\
	v.x /= (tmplen + 1.0e-299);		\
	v.y /= (tmplen + 1.0e-299);		\
}


void r2d_clip(r2d_poly* poly, r2d_plane* planes, r2d_int nplanes) {

	// variable declarations
	r2d_int v, p, np, onv, vstart, vcur, vnext, numunclipped; 

        // Number of starting verts.
	r2d_int* nverts = &poly->nverts; 
	if(*nverts <= 0) return;

	// signed distances to the clipping plane
	r2d_real* sdists = (r2d_real*) malloc((*nverts)*sizeof(r2d_real));
	r2d_real smin, smax;

	// for marking clipped vertices
	r2d_int* clipped = (r2d_int*) malloc((*nverts)*sizeof(r2d_int));

        // Be prepared for reallocating with new pointers.
	r2d_vertex* verts_new;
        r2d_real* sdists_new;
        r2d_int* clipped_new;

	// loop over each clip plane
	for(p = 0; p < nplanes; ++p) {
	
		// calculate signed distances to the clip plane
		onv = *nverts;
		smin = 1.0e30;
		smax = -1.0e30;
		memset(clipped, 0, onv*sizeof(r2d_int));
		for(v = 0; v < onv; ++v) {
			sdists[v] = planes[p].d + dot(poly->verts[v].pos, planes[p].n);
			if(sdists[v] < smin) smin = sdists[v];
			if(sdists[v] > smax) smax = sdists[v];
			if(sdists[v] < 0.0) clipped[v] = 1;
		}

		// skip this face if the poly lies entirely on one side of it 
		if(smin >= 0.0) continue;
		if(smax <= 0.0) {
			*nverts = 0;
                        free(sdists);
                        free(clipped);
			return;
		}

		// check all edges and insert new vertices on the bisected edges 
		for(vcur = 0; vcur < onv; ++vcur) {
			if(clipped[vcur]) continue;
			for(np = 0; np < 2; ++np) {
				vnext = poly->verts[vcur].pnbrs[np];
				if(!clipped[vnext]) continue;
                                verts_new = realloc(poly->verts, (*nverts + 1)*sizeof(r2d_vertex));
                                sdists_new = realloc(sdists, (*nverts + 1)*sizeof(r2d_real));
                                clipped_new = realloc(clipped, (*nverts + 1)*sizeof(r2d_int));
                                if (verts_new == NULL || sdists_new == NULL || clipped_new == NULL) {
                                  printf("r2d_clip ERROR: unable to increase vertex buffer size.\n");
                                  abort();
                                }
                                poly->verts = verts_new;
                                sdists = sdists_new;
                                clipped = clipped_new;
				poly->verts[*nverts].pnbrs[1-np] = vcur;
				poly->verts[*nverts].pnbrs[np] = -1;
				poly->verts[vcur].pnbrs[np] = *nverts;
                                clipped[*nverts] = 0;
				wav(poly->verts[vcur].pos, -sdists[vnext],
					poly->verts[vnext].pos, sdists[vcur],
					poly->verts[*nverts].pos);
				(*nverts)++;
			}
		}

		// for each new vert, search around the poly for its new neighbors
		// and doubly-link everything
		for(vstart = onv; vstart < *nverts; ++vstart) {
			if(poly->verts[vstart].pnbrs[1] >= 0) continue;
			vcur = poly->verts[vstart].pnbrs[0];
			do {
				vcur = poly->verts[vcur].pnbrs[0]; 
			} while(vcur < onv);
			poly->verts[vstart].pnbrs[1] = vcur;
			poly->verts[vcur].pnbrs[0] = vstart;
		}

		// go through and compress the vertex list, removing clipped verts
		// and re-indexing accordingly (reusing `clipped` to re-index everything)
		numunclipped = 0;
		for(v = 0; v < *nverts; ++v) {
			if(!clipped[v]) {
				poly->verts[numunclipped] = poly->verts[v];
				clipped[v] = numunclipped++;
			}
		}
		*nverts = numunclipped;
		for(v = 0; v < *nverts; ++v) {
			poly->verts[v].pnbrs[0] = clipped[poly->verts[v].pnbrs[0]];
			poly->verts[v].pnbrs[1] = clipped[poly->verts[v].pnbrs[1]];
		}	
	}

        // Clean up.
        free(sdists);
        free(clipped);
}

void r2d_reduce(r2d_poly* poly, r2d_real* moments, r2d_int polyorder) {

	// var declarations
	r2d_int vcur, vnext, m, i, j, corder;
	r2d_real twoa;
	r2d_rvec2 v0, v1; 

	// direct access to vertex buffer
	r2d_vertex* vertbuffer = poly->verts; 
	r2d_int* nverts = &poly->nverts; 

	// zero the moments
	for(m = 0; m < R2D_NUM_MOMENTS(polyorder); ++m)
		moments[m] = 0.0;

	if(*nverts <= 0) return;

	// Storage for coefficients
	// keep two layers of the triangle of coefficients
	r2d_int prevlayer = 0;
	r2d_int curlayer = 1;
	r2d_real D[polyorder+1][2];
	r2d_real C[polyorder+1][2];

	// iterate over edges and compute a sum over simplices 
	for(vcur = 0; vcur < *nverts; ++vcur) {

		vnext = vertbuffer[vcur].pnbrs[0];
		v0 = vertbuffer[vcur].pos;
		v1 = vertbuffer[vnext].pos;
		twoa = (v0.x*v1.y - v0.y*v1.x); 

		// calculate the moments
		// using the fast recursive method of Koehl (2012)
		// essentially building a set of Pascal's triangles, one layer at a time

		// base case
		D[0][prevlayer] = 1.0;
		C[0][prevlayer] = 1.0;
		moments[0] += 0.5*twoa;

		// build up successive polynomial orders
		for(corder = 1, m = 1; corder <= polyorder; ++corder) {
			for(i = corder; i >= 0; --i, ++m) {
				j = corder - i;
				C[i][curlayer] = 0; 
				D[i][curlayer] = 0;  
				if(i > 0) {
					C[i][curlayer] += v1.x*C[i-1][prevlayer];
					D[i][curlayer] += v0.x*D[i-1][prevlayer]; 
				}
				if(j > 0) {
					C[i][curlayer] += v1.y*C[i][prevlayer];
					D[i][curlayer] += v0.y*D[i][prevlayer]; 
				}
				D[i][curlayer] += C[i][curlayer]; 
				moments[m] += twoa*D[i][curlayer];
			}
			curlayer = 1 - curlayer;
			prevlayer = 1 - prevlayer;
		}
	}

	// reuse C to recursively compute the leading multinomial coefficients
	C[0][prevlayer] = 1.0;
	for(corder = 1, m = 1; corder <= polyorder; ++corder) {
		for(i = corder; i >= 0; --i, ++m) {
			j = corder - i;
			C[i][curlayer] = 0.0; 
			if(i > 0) C[i][curlayer] += C[i-1][prevlayer];
			if(j > 0) C[i][curlayer] += C[i][prevlayer];
			moments[m] /= C[i][curlayer]*(corder+1)*(corder+2);
		}
		curlayer = 1 - curlayer;
		prevlayer = 1 - prevlayer;
	}
}

r2d_int r2d_is_good(r2d_poly* poly) {

	r2d_int v;
	r2d_int* nverts = &poly->nverts; 
	r2d_int vct[*nverts];

	// direct access to vertex buffer
	r2d_vertex* vertbuffer = poly->verts; 

	// consistency check
	memset(&vct, 0, sizeof(vct));
	for(v = 0; v < *nverts; ++v) {

		// return false if vertices share an edge with themselves 
		// or if any edges are obviously invalid
		if(vertbuffer[v].pnbrs[0] == vertbuffer[v].pnbrs[1]) return 0;
		if(vertbuffer[v].pnbrs[0] >= *nverts) return 0;
		if(vertbuffer[v].pnbrs[1] >= *nverts) return 0;

		vct[vertbuffer[v].pnbrs[0]]++;
		vct[vertbuffer[v].pnbrs[1]]++;
	}
	
	// return false if any vertices are pointed to 
	// by more or fewer than two other vertices
	for(v = 0; v < *nverts; ++v) if(vct[v] != 2) return 0;

	return 1;
}

void r2d_rotate(r2d_poly* poly, r2d_real theta) {
	r2d_int v;
	r2d_rvec2 tmp;
	r2d_real sine = sin(theta);
	r2d_real cosine = cos(theta);
	for(v = 0; v < poly->nverts; ++v) {
		tmp = poly->verts[v].pos;
		poly->verts[v].pos.x = cosine*tmp.x - sine*tmp.y; 
		poly->verts[v].pos.x = sine*tmp.x + cosine*tmp.y; 
	}
}

void r2d_translate(r2d_poly* poly, r2d_rvec2 shift) {
	r2d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		poly->verts[v].pos.x += shift.x;
		poly->verts[v].pos.y += shift.y;
	}
}

void r2d_scale(r2d_poly* poly, r2d_real scale) {
	r2d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		poly->verts[v].pos.x *= scale;
		poly->verts[v].pos.y *= scale;
	}
}

void r2d_shear(r2d_poly* poly, r2d_real shear, r2d_int axb, r2d_int axs) {
	r2d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		poly->verts[v].pos.xy[axb] += shear*poly->verts[v].pos.xy[axs];
	}
}

void r2d_affine(r2d_poly* poly, r2d_real mat[3][3]) {
	r2d_int v;
	r2d_rvec2 tmp;
	r2d_real w;
	for(v = 0; v < poly->nverts; ++v) {
		tmp = poly->verts[v].pos;

		// affine transformation
		poly->verts[v].pos.x = tmp.x*mat[0][0] + tmp.y*mat[0][1] + mat[0][2];
		poly->verts[v].pos.y = tmp.x*mat[1][0] + tmp.y*mat[1][1] + mat[1][2];
		w = tmp.x*mat[2][0] + tmp.y*mat[2][1] + mat[2][2];
	
		// homogeneous divide if w != 1, i.e. in a perspective projection
		poly->verts[v].pos.x /= w;
		poly->verts[v].pos.y /= w;
	}
}


r2d_poly r2d_init_empty_poly() {
	r2d_poly result = {NULL, 0};
	return result;
}


void r2d_init_box(r2d_poly* poly, r2d_rvec2 rbounds[2]) {

	// Allocate memory
	r2d_int* nverts = &poly->nverts; 
	*nverts = 4;
        if (poly->verts != NULL) free(poly->verts);
	poly->verts = malloc((*nverts)*sizeof(r2d_vertex));
	if (poly->verts == NULL) {
		printf("r2d_init_box ERROR: unable to allocate verts\n");
                abort();
	}

	// direct access to vertex buffer
	r2d_vertex* vertbuffer = poly->verts; 
	
	vertbuffer[0].pnbrs[0] = 1;	
	vertbuffer[0].pnbrs[1] = 3;	
	vertbuffer[1].pnbrs[0] = 2;	
	vertbuffer[1].pnbrs[1] = 0;	
	vertbuffer[2].pnbrs[0] = 3;	
	vertbuffer[2].pnbrs[1] = 1;	
	vertbuffer[3].pnbrs[0] = 0;	
	vertbuffer[3].pnbrs[1] = 2;	
	vertbuffer[0].pos.x = rbounds[0].x; 
	vertbuffer[0].pos.y = rbounds[0].y; 
	vertbuffer[1].pos.x = rbounds[1].x; 
	vertbuffer[1].pos.y = rbounds[0].y; 
	vertbuffer[2].pos.x = rbounds[1].x; 
	vertbuffer[2].pos.y = rbounds[1].y; 
	vertbuffer[3].pos.x = rbounds[0].x; 
	vertbuffer[3].pos.y = rbounds[1].y; 

}


void r2d_init_poly(r2d_poly* poly, r2d_rvec2* vertices, r2d_int numverts) {

	// Allocate memory
	r2d_int* nverts = &poly->nverts; 
	*nverts = numverts;
        if (poly->verts != NULL) free(poly->verts);
	poly->verts = malloc((*nverts)*sizeof(r2d_vertex));
	if (poly->verts == NULL) {
		printf("r2d_init_poly ERROR: unable to allocate verts\n");
                abort();
	}

	// direct access to vertex buffer
	r2d_vertex* vertbuffer = poly->verts; 

	// init the poly
	r2d_int v;
	for(v = 0; v < *nverts; ++v) {
		vertbuffer[v].pos = vertices[v];
		vertbuffer[v].pnbrs[0] = (v+1)%(*nverts);
		vertbuffer[v].pnbrs[1] = (*nverts+v-1)%(*nverts);
	}
}


void r2d_free_poly(r2d_poly* poly) {
	if (poly->verts != NULL) {
		free(poly->verts);
	}
}

void r2d_copy_poly(r2d_poly* topoly, r2d_poly* frompoly) {
	topoly->nverts = frompoly->nverts;
	if (topoly->verts != NULL) {
		free(topoly->verts);
	}
	topoly->verts = malloc((topoly->nverts)*sizeof(r2d_vertex));
	if (topoly->verts == NULL) {
		printf("r2d_copy_poly ERROR: unable to allocate target vertex buffer.\n");
                abort();
	}
	memcpy(topoly->verts, frompoly->verts, (topoly->nverts)*sizeof(r2d_vertex));
}

void r2d_box_faces_from_verts(r2d_plane* faces, r2d_rvec2* rbounds) {
	faces[0].n.x = 0.0; faces[0].n.y = 1.0; faces[0].d = rbounds[0].y; 
	faces[1].n.x = 1.0; faces[1].n.y = 0.0; faces[1].d = rbounds[0].x; 
	faces[2].n.x = 0.0; faces[2].n.y = -1.0; faces[2].d = rbounds[1].y; 
	faces[3].n.x = -1.0; faces[3].n.y = 0.0; faces[3].d = rbounds[1].x; 
}

void r2d_poly_faces_from_verts(r2d_plane* faces, r2d_rvec2* vertices, r2d_int numverts) {

	// dummy vars
	r2d_int f;
	r2d_rvec2 p0, p1;

	// calculate a centroid and a unit normal for each face 
	for(f = 0; f < numverts; ++f) {

		p0 = vertices[f];
		p1 = vertices[(f+1)%numverts];

		// normal of the edge
		faces[f].n.x = p0.y - p1.y;
		faces[f].n.y = p1.x - p0.x;

		// normalize the normals and set the signed distance to origin
		norm(faces[f].n);
		faces[f].d = -dot(faces[f].n, p0);

	}
}

r2d_real r2d_orient(r2d_rvec2 pa, r2d_rvec2 pb, r2d_rvec2 pc) {
	return 0.5*((pa.x - pc.x)*(pb.y - pc.y) - (pb.x - pc.x)*(pa.y - pc.y)); 
}

void r2d_print(r2d_poly* poly) {
	r2d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		printf("  vertex %d: pos = ( %.10e , %.10e ), nbrs = %d %d\n", 
				v, poly->verts[v].pos.x, poly->verts[v].pos.y, poly->verts[v].pnbrs[0], poly->verts[v].pnbrs[1]);
	}
}

